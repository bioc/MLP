#' Graphical Representation of GO Based MLP Results
#' @param object object of class MLP (as produced by the MLP function)
#' @param nRow number of GO IDs for which to produce the plot
#' @param main main title of the graph; if NULL (default) the main title is set to 'GO graph'
#' @param nCutDescPath number of characters at which the pathway description should be cut 
#' (inserted in a new line), 30 by default
#' @return GO graph is plotted to the current device
#' @examples if (require(GO.db) && require(Rgraphviz)){
#'   pathExampleMLPResult <- system.file("exampleFiles", "exampleMLPResult.rda", package = "MLP")
#'   load(pathExampleMLPResult)
#'   plotGOgraph(exampleMLPResult, main = "GO Graph")
#' }
#' @importFrom graphics legend
#' @importFrom utils getFromNamespace
#' @export
plotGOgraph <- function (object, nRow = 5, main = NULL, 
	nCutDescPath = 30) {

  if (!inherits(object, "MLP"))
    stop("The 'object' argument should be an object of class 'MLP' as produced by the MLP function")
  if (is.data.frame(attributes(object)$geneSetSource))
    stop("Plotting a GO graph is only possible for MLP results based om geneSetSource 'GOBP', 'GOMF', or 'GOCC'")

  main <- if (is.null(main)) "Go graph" else main

  requireNamespace("GO.db")
  requireNamespace("Rgraphviz")
  requireNamespace("GOstats")
  requireNamespace("annotate")

  goids <- rownames(object)[seq.int(nRow)]
  ontology <- sub("GO", "", attributes(object)$geneSetSource)
  graphEnv <- getFromNamespace(x = paste("GO", ontology, "PARENTS", sep = ""), ns = "GO.db")
  basicGraph <- GOstats::GOGraph(
		x = goids, 
		dataenv = graphEnv
	)
  basicGraph <- removeNode("all", basicGraph)
  basicGraph <- removeNode(setdiff(nodes(basicGraph), rownames(object)),
      basicGraph)
  basicGraph <- Rgraphviz::layoutGraph(basicGraph)
  pvalues <- object[nodes(basicGraph), "geneSetPValue"]
  names(pvalues) <- nodes(basicGraph)
  pvalues <- pvalues[!is.na(pvalues)]
  pvalues[pvalues == 0] <- min(pvalues[pvalues != 0])/10
  scores <- -log10(pvalues)
  scores[scores <= 0.1] <- 0.1
  nColors <- round(max(scores) * 10)
  gocolors <- gplots::colorpanel(nColors, low = "lightyellow", high = "olivedrab")
  nodeFillColor <- rep("white", length(nodes(basicGraph)))
  names(nodeFillColor) <- nodes(basicGraph)
  nodeFillColor[names(scores)] <- gocolors[trunc(scores * 10)]
  nodeRenderInfo(basicGraph) <- list(fill = nodeFillColor)
  inbin <- object$totalGeneSetSize
  names(inbin) <- rownames(object)
  onchip <- object$testedGeneSetSize
  names(onchip) <- rownames(object)
  allcounts <- matrix(ncol = 2, nrow = length(inbin), dimnames = list(c(names(inbin)),
          c("onchip", "inbin")))
  allcounts[names(onchip), "onchip"] <- onchip
  allcounts[names(inbin), "inbin"] <- inbin
  counts <- apply(allcounts, 1, function(x) {
        paste(x[1], x[2], sep = " - ")
      })
  terms <- annotate::getGOTerm(nodes(basicGraph))
  goTerm <- terms[[1]]
  goTerm <- sapply(goTerm, function(x) {
       paste(
		substr(x, start = 1, stop = nCutDescPath), 
		substr(x, start =  nCutDescPath+1, stop = nCutDescPath*2), 
		sep = "\n")
      })
  nodeLabel <- paste(nodes(basicGraph), goTerm[nodes(basicGraph)],
      counts[nodes(basicGraph)], sep = "\n")
  names(nodeLabel) <- nodes(basicGraph)
  nodeRenderInfo(basicGraph) <- list(label = nodeLabel)
  edgeRenderInfo(basicGraph)$arrowhead <- "none"
  nodeRenderInfo(basicGraph) <- list(shape = "ellipse")
  nodeRenderInfo(basicGraph) <- list(cex = 0.5)
  nodeRenderInfo(basicGraph) <- list(lWidth = 60)
  nodeRenderInfo(basicGraph) <- list(labelJust = "c")
  graphRenderInfo(basicGraph)$main <- main
  
  # if graph without edge, renderGraph returns error:
  #Error in text.default(labelX, labelY, label, col = textCol, cex = cex *  :  no coordinates were supplied
  if(length(unlist(edges(basicGraph))) > 0){
	  
	  Rgraphviz::renderGraph(basicGraph)
	  legend("right", "bottom", legend = paste(c(" least", 
	              " medium", " most"), " (scores ", round(max(scores)) *
	              c(2, 5, 8)/10, ")", sep = ""), fill = gocolors[round(max(scores)) *
	              c(2, 5, 8)], cex = 0.7)

  }else{
	stop("The function cannot deal currently with a set of pathways without any links.")
  }
}
