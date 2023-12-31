% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addGeneSetDescription.R
\name{addGeneSetDescription}
\alias{addGeneSetDescription}
\title{Utility function which adds the biological description of the gene sets as
a column to the return value of the MLP function (data frame)}
\usage{
addGeneSetDescription(object, geneSetSource = NULL)
}
\arguments{
\item{object}{object of class 'MLP' as produced by the 'MLP' function}

\item{geneSetSource}{source to be used to construct the list of pathway categories; 
for public data sources, the user can specify a string (one of 'GOBP', 'GOMF', 'GOCC', 'KEGG' or 'REACTOME')
and BioC packages will be used to construct the list of pathway categories; 
for non-public data sources, the user can pass the pathway data as a dataframe with (at least) 
the following four columns: PATHWAYID, TAXID, PATHWAYNAME and GENEID. It is assumed all columns
are of type character. The 'geneSetSource' argument should be the same as the argument
provided to the getGeneSets function; defaults to NULL}
}
\value{
the data frame as returned by MLP enriched with an additional column geneSetDescription, providing
a concise description of the gene set
}
\description{
Utility function which adds the biological description of the gene sets as
a column to the return value of the MLP function (data frame)
}
\examples{
if (require(GO.db)){
  pathExamplePValues <- system.file("exampleFiles", "examplePValues.rda", package = "MLP")
  load(pathExamplePValues)
	 geneSet <- getGeneSets(species = "Mouse", geneSetSource = "GOBP", entrezIdentifiers = names(examplePValues))
  mlpResult <- MLP(geneSet = geneSet, geneStatistic = examplePValues, addGeneSetDescription = FALSE)
  head(mlpResult)
  mlpResultsWithGSDescr <- addGeneSetDescription(object = mlpResult, geneSetSource = "GOBP")
  head(mlpResultsWithGSDescr)
}
}
\seealso{
\link{MLP}
}
