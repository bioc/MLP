# Script to generate example files of the MLP package
# 
# Author: Tobias Verbeke
###############################################################################

library(MLP)

# make sure the R version is the one that corresponds to BioC-devel
# see also http://www.bioconductor.org/developers/useDevel/

# use the BioC-devel version of GO.db 

# set working directory to inst/exampleFiles
load("examplePValues.rda")
exampleGeneSet <- getGeneSets(species = "Mouse", geneSetSource = "GOBP", 
    entrezIdentifiers = names(examplePValues))
save(exampleGeneSet, file = "exampleGeneSet.rda")

