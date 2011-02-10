#' Draw a Barplot for MLP Results 
#' @param object object of class MLP 
#' @param nRow number of rows of the MLP data frame to depict in the barplot; defaults to 20.
#' @param barColors vector of colors to use for the bars of the barplot; defaults to NULL; 
#'   if NULL, three gray shades are used reflecting the proportion of tested genes of a gene set
#'   versus the total number of genes in a geneset. If the proportion exceeds 75\%, the darkest
#'   shade is used; between 50 and 75\% a moderately dark shade is used; below 50\% a lighter gray
#'   shade is used.
#' @param main main title; if NULL (default) "Effect of the treatment on <geneSetSource> gene sets"
#' will be used
#' @return the midpoints of all the bars are returned invisibly (using the conventions of barplot); 
#'   an MLP-specific barplot is drawn to the current device;  
#' @seealso barplot
#' @examples pathExampleMLPResult <- system.file("exampleFiles", "exampleMLPResult.rda", package = "MLP")
#' load(pathExampleMLPResult)
#' dev.new(width = 10, height = 10)
#' op <- par(mar = c(30, 10, 6, 2))
#' mlpBarplot(exampleMLPResult)
#' par(op)
#' @export
mlpBarplot <- function (object, nRow = 20, barColors = NULL, main = NULL) {
  
  if (!inherits(object, "MLP")) 
    stop("'object' should be an object of class 'MLP' as produced by the MLP function")
  
  geneSetSource <- attr(object, "geneSetSource")
  
  if (is.null(object$geneSetDescription)) {
    object <- addGeneSetDescription(object, geneSetSource = geneSetSource)
  }
  mlpResults <- head(object, nRow)
  dat <- -log10(mlpResults$geneSetPValue)
  #Fix Inf values
  dat[is.infinite(dat)] <- max(dat[is.finite(dat)]) + 50
  names(dat) <- rownames(mlpResults)
  
  if (is.null(barColors)){
    barColors <- rep("grey", length(dat))
    percentTested <- 100*(object$testedGeneSetSize/object$totalGeneSetSize)
    percentTested <- percentTested[1:length(dat)]
    
    barColors[percentTested >= 75] <- "grey50"
    barColors[percentTested < 75 & percentTested > 50] <- "grey60"
    barColors[percentTested <= 50] <- "grey70"
  } else {
    barColors <- barColors
  }
  
  descr <- mlpResults$geneSetDescription
  descriptionLength <- 60
  descr <- substr(descr, 1, descriptionLength)
  names(dat) <- paste(descr, " (", mlpResults$testedGeneSetSize, 
      "-", mlpResults$totalGeneSetSize, ")", sep = "")
  bottomMar <- 30
  op <- par(mar = c(bottomMar, 10, 6, 2))
  mp <- barplot(dat, xlab = "", main = "", border = "white", 
      las = 3, ylab = "", col = barColors)
  if (is.null(main)) {
    if (is.data.frame(geneSetSource)){ # external data
      mainTitle <- "Effect of the treatment on the gene sets"
    } else { # geneSetSource %in% c("GOBP", "GOMF", "GOCC", "KEGG")
      mainTitle <- paste("Effect of the treatment on", geneSetSource, "gene sets")
    }
  } else {
    mainTitle <- main
  }
  mtext(mainTitle, side = 2, line = 5)
  par(op)
  invisible(mp)
}
