# Script to generate example files of the MLP package
# 
# Author: Tobias Verbeke
###############################################################################

library(MLP)

# set working directory to inst/exampleFiles
load("examplePValues.rda")
exampleGeneSet <- getGeneSets(species = "Mouse", geneSetSource = "GOBP", 
    entrezIdentifiers = names(examplePValues))
save(exampleGeneSet, file = "exampleGeneSet.rda")

