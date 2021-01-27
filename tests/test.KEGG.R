
pathPValues <- system.file("exampleFiles", "examplePValues.rda", package = "MLP")
load(pathPValues)

pvalues <- examplePValues[seq.int(1000)]

system.time(geneSet <- getGeneSets(
	species = "Mouse", 
	geneSetSource = "KEGG", 
	entrezIdentifiers = names(pvalues)
))

set.seed(111)
mlpOut <- MLP(
	geneSet = geneSet, 
	geneStatistic = pvalues
) 	

mlpOutWithGeneSetDescr <- addGeneSetDescription(object = mlpOut, geneSetSource = "KEGG")
	