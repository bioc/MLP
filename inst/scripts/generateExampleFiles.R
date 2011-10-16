# Script to generate example files of the MLP package
# 
# Author: Tobias Verbeke
###############################################################################

library(MLP)

# in inst/exampleFiles
load("examplePValues.rda")
exampleGeneSet <- getGeneSets(species = "Mouse", geneSetSource = "GOBP", 
    entrezIdentifiers = names(examplePValues))
save(exampleGeneSet, file = "exampleGeneSet.rda")
