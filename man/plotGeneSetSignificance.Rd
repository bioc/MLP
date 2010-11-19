\name{plotGeneSetSignificance}
\alias{plotGeneSetSignificance}
\title{Plot the Significance for the Genes of a Given Gene Set...}
\usage{plotGeneSetSignificance(geneSet, geneSetIdentifier, geneStatistic, annotationPackage, barColors)
}
\description{Plot the Significance for the Genes of a Given Gene Set}
\value{no return value}
\arguments{\item{geneSet}{object of class 'geneSetMLP' as produced by function getGeneSets}
\item{geneSetIdentifier}{identifier of the gene set for which a significance plot should be produced;
character of length one}
\item{geneStatistic}{named vector of gene statistics (e.g. p values); the names of the vector
are Entrez Gene identifiers}
\item{annotationPackage}{name of the annotation package to be used (without .db extension);
character of length one}
\item{barColors}{named vector of colors to use for the bars of the barplot; the 
names of the vector are Entrez Gene identifiers and the vector should be of length equal
to the length of the geneStatistic vector 
defaults to NULL in which case 'grey50' is used}
}
\examples{pathExamplePValues <- system.file("exampleFiles", "examplePValues.rda", package = "MLP")
pathExampleGeneSet <- system.file("exampleFiles", "exampleGeneSet.rda", package = "MLP")
pathExampleMLPResult <- system.file("exampleFiles", "exampleMLPResult.rda", package = "MLP")
load(pathExampleGeneSet)
load(pathExamplePValues)
load(pathExampleMLPResult) 
annotationPackage <- if (require(mouse4302mmentrezg.db)) "mouse4302mmentrezg" else "mouse4302"
geneSetID <- rownames(exampleMLPResult)[1]
dev.new(width = 10, height = 10)
op <- par(mar = c(25, 10, 6, 2))
plotGeneSetSignificance(
geneSet = exampleGeneSet, 
geneSetIdentifier = geneSetID, 
geneStatistic = examplePValues, 
annotationPackage = annotationPackage
)
par(op)}

