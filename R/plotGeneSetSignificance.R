#' Plot the Significance for the Genes of a Given Gene Set 
#' @param geneSet object of class 'geneSetMLP' as produced by function getGeneSets
#' @param geneSetIdentifier identifier of the gene set for which a significance plot should be produced;
#'   character of length one  
#' @param geneStatistic named vector of gene statistics (e.g. p values); the names of the vector
#'   are Entrez Gene identifiers 
#' @param annotationPackage name of the annotation package to be used (without .db extension);
#'   character of length one
#' @param barColors named vector of colors to use for the bars of the barplot; the 
#' names of the vector are Entrez Gene identifiers and the vector should be of length equal
#' to the length of the geneStatistic vector 
#' defaults to NULL in which case 'grey50' is used 
#' @return no return value 
#' @examples
#' pathExamplePValues <- system.file("exampleFiles", "examplePValues.rda", package = "MLP")
#' pathExampleGeneSet <- system.file("exampleFiles", "exampleGeneSet.rda", package = "MLP")
#' pathExampleMLPResult <- system.file("exampleFiles", "exampleMLPResult.rda", package = "MLP")
#' load(pathExampleGeneSet)
#' load(pathExamplePValues)
#' load(pathExampleMLPResult) 
#' # annotationPackage <- if (require(mouse4302mmentrezg.db)) "mouse4302mmentrezg" else "mouse4302"
#' annotationPackage <- "mouse4302"
#' geneSetID <- rownames(exampleMLPResult)[1]
#' dev.new(width = 10, height = 10)
#' op <- par(mar = c(25, 10, 6, 2))
#' plotGeneSetSignificance(
#'     geneSet = exampleGeneSet, 
#'     geneSetIdentifier = geneSetID, 
#'     geneStatistic = examplePValues, 
#'     annotationPackage = annotationPackage
#' )
#' par(op)
#' @export
plotGeneSetSignificance <- function(geneSet, geneSetIdentifier, geneStatistic, annotationPackage, barColors = NULL, descriptionInMainTitle = TRUE){
  
  if (!inherits(geneSet, "geneSetMLP"))
    stop("geneSet should be an object of class 'geneSetMLP' as produced by 'getGeneSets'")
  
  if (!geneSetIdentifier %in% names(geneSet))
    stop("Please provide as 'geneSetName' a gene set name belonging to 'geneSets', i.e. the group of gene sets specified.")
  
  if (!is.null(barColors) && (length(barColors) != length(geneStatistic)))
    stop("'colorBars' should be a vector of the same length as 'geneStatistic'")
  
  if (!is.null(barColors) && is.null(names(barColors)))
    stop("'colorBars' should be a named vector of the same length as 'geneStatistic'")
  
  
  require(annotate)
  require(paste(annotationPackage, ".db", sep = ""), character.only = TRUE)
  
  descriptions <- attr(geneSet, "descriptions")
  
  entrezids <- geneSet[[geneSetIdentifier]]
  entrezids <- entrezids[entrezids %in% names(geneStatistic)]
  genePValues <- geneStatistic[entrezids]
  genePValues <- sort(genePValues)
  psids <- if (length(grep("_at", names(genePValues))) == 0){
    paste(names(genePValues), "_at", sep = "")
  } else {
    names(genePValues)
  }
  barColors <- if (is.null(barColors)) "grey50" else barColors[names(genePValues)]
  
  names(genePValues) <- paste(
      unlist(lookUp(psids, annotationPackage, "SYMBOL")),
      unlist(lookUp(psids, annotationPackage, "GENENAME")),
      sep = ":")
  names(genePValues) <- substr(names(genePValues), 1, 60)
  
  title = paste("Significance of tested genes involved in gene set", geneSetIdentifier)
  if(descriptionInMainTitle){
    title = paste("Significance of tested genes involved in gene set", geneSetIdentifier, descriptions[geneSetIdentifier])
  }
  
  barplot(-log10(genePValues), xlab = "", 
      main = title, 
      border = "white", col = barColors,
      las = 3, ylab = "Significance")
}
