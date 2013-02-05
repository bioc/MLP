#' Prepare Pathway Data for the MLP Function
#' 
#' The return value of the getGeneSets function has as primary use
#' to serve as geneSet argument for the MLP function
#' @param species character vector of length one indicating the species, one of
#' 'Mouse', 'Human' or 'Rat'; defaults to 'Mouse'. 
#' @param geneSetSource source to be used to construct the list of pathway categories; 
#' for public data sources, the user can specify a string (one of 'GOBP', 'GOMF', 'GOCC', 'KEGG' or 'REACTOME')
#' and BioC packages will be used to construct the list of pathway categories; 
#' for non-public data sources, the user can pass the pathway data as a dataframe with (at least) 
#' the following four columns: PATHWAYID, TAXID, PATHWAYNAME and GENEID. It is assumed all columns
#' are of type character.
#' @param entrezIdentifiers Entrez Gene identifiers used to subset the relevant gene set
#' @return object of class geneSetMLP which is essentially a named list of pathway categories. 
#' Each list component contains a vector of Entrez Gene identifiers related to that particular pathway
#' @import AnnotationDbi
#' @examples if (require(GO.db) && require(org.Mm.eg.db)){
#'   pathExampleData <- system.file("exampleFiles", "expressionSetGcrma.rda", package = "MLP")
#'   pathExamplePValues <- system.file("exampleFiles", "examplePValues.rda", package = "MLP")
#'   load(pathExampleData)
#'   load(pathExamplePValues)
#'   geneSet <- getGeneSets(species = "Mouse", geneSetSource = "GOBP", entrezIdentifiers = names(examplePValues)[1:2000])
#'   head(geneSet)
#' }
#' @export
getGeneSets <- function (species = "Mouse", geneSetSource = NULL, entrezIdentifiers){
  descriptions = "blah"
  
  ### checks
  if (!species %in% c("Mouse", "Human", "Rat", "Dog")) 
    stop("The 'species' argument should be one of 'Mouse', 'Human', 'Rat' or 'Dog')")

  if (is.null(geneSetSource)) 
    stop("Please provide a source of gene sets. For more info, see the help page.")
  
  if (!is.data.frame(geneSetSource) && length(geneSetSource) != 1 && !(geneSetSource %in% 
        c("GOBP", "GOMF", "GOCC", "KEGG", "REACTOME"))) 
    stop("The 'geneSetSource' argument should be one of 'GOBP', 'GOMF', 'GOCC', 'KEGG', 'REACTOME' or a data.frame.  More info, see help.")
  
  if (is.data.frame(geneSetSource) && any(!(c("PATHWAYID", "TAXID", "PATHWAYNAME", "GENEID") %in% colnames(geneSetSource))))
      stop("The geneSetSource as data.frame should have at least the 4 columns 'PATHWAYID', 'TAXID', 'PATHWAYNAME' and 'GENEID'. More info on their content, see help.")
  
  ### character
  if (is.character(geneSetSource)) {
    if (geneSetSource %in% c("GOBP", "GOMF", "GOCC")) {
      require(GO.db)
      ontology <- sub("GO", "", geneSetSource)
      switch(species, Mouse = {
            require(org.Mm.eg.db)
            geneSetToEntrez <- as.list(org.Mm.egGO2ALLEGS)
          }, Human = {
            require(org.Hs.eg.db)
            geneSetToEntrez <- as.list(org.Hs.egGO2ALLEGS)
          }, Rat = {
            require(org.Rn.eg.db)
            geneSetToEntrez <- as.list(org.Rn.egGO2ALLEGS)
          }, Dog = {
            require(org.Cf.eg.db)
            geneSetToEntrez <- as.list(org.Cf.egGO2ALLEGS)
          })
      switch(ontology, BP = {
            GOs <- names(as.list(GOBPANCESTOR))
          }, MF = {
            GOs <- names(as.list(GOMFANCESTOR))
          }, CC = {
            GOs <- names(as.list(GOCCANCESTOR))
          })
      go <- geneSetToEntrez[names(geneSetToEntrez) %in% 
              GOs]
      geneSets <- lapply(go, unique)
      anyGenesInGeneSet <- ifelse(unlist(lapply(geneSets, 
                  length)) > 0, TRUE, FALSE)
      geneSets <- geneSets[anyGenesInGeneSet]
      goterms <- Term(GOTERM)
      descriptions <- goterms[names(geneSets)]
    }
    if (geneSetSource == "KEGG") {
      require(KEGG.db)
      geneSetToEntrez <- as.list(KEGGPATHID2EXTID)
      switch(species, 
          Mouse = {
            prefix <- "mmu"
          },
          Human = {
            prefix <- "hsa"
          }, 
          Rat = {
            prefix <- "rno"
          }, 
          Dog = {
            prefix <- "cfa"
          })
      geneSets <- geneSetToEntrez[grep(prefix, names(geneSetToEntrez))]
      descriptions <- unlist(mget(sub("\\D+?(\\d+)", "\\1", names(geneSetToEntrez)), KEGGPATHID2NAME))
      names(descriptions) <- paste(prefix, names(descriptions), sep = "")
    }
    if (geneSetSource == "REACTOME") {
      require(reactome.db)
      pathways <- toTable(reactomePATHNAME2ID)
      allReactomeIDs <- ls(reactomePATHID2EXTID)
      switch(species, Mouse = {
            pathwaysSelectedSpecies <- pathways[grep("Mus musculus: ", iconv(pathways$path_name)), ]
            #Neem eerst de reactomeids met een mapping beperk deze vervolgens met de geselecteerde species.
            geneSets <- mget(pathwaysSelectedSpecies$reactome_id[pathwaysSelectedSpecies$reactome_id %in% allReactomeIDs], reactomePATHID2EXTID)
          }, Human = {
            pathwaysSelectedSpecies <- pathways[grep("Homo sapiens: ", iconv(pathways$path_name)), ]
            geneSets <- mget(pathwaysSelectedSpecies$reactome_id[pathwaysSelectedSpecies$reactome_id %in% allReactomeIDs], reactomePATHID2EXTID)
          }, Rat = {
            pathwaysSelectedSpecies <- pathways[grep("Rattus norvegicus: ", iconv(pathways$path_name)), ]
            geneSets <- mget(pathwaysSelectedSpecies$reactome_id[pathwaysSelectedSpecies$reactome_id %in% allReactomeIDs], reactomePATHID2EXTID)
          }, Dog = {
            pathwaysSelectedSpecies <- pathways[grep("Canis familiaris: ", iconv(pathways$path_name)), ]
            geneSets <- mget(pathwaysSelectedSpecies$reactome_id[pathwaysSelectedSpecies$reactome_id %in% allReactomeIDs], reactomePATHID2EXTID)
          })
      descriptions <- pathwaysSelectedSpecies$path_name
      names(descriptions) <- pathwaysSelectedSpecies$reactome_id
    }
  } else { ### data frame
    switch(species, Mouse = {
          TAXID <- "10090"
        }, Human = {
          TAXID <- "9606"
        }, Rat = {
          TAXID <- "10116"
        }, Dog = {
          TAXID <- "9615"
        })
    if (!(TAXID %in% geneSetSource$TAXID)) 
      stop(paste("Please check the available TAXIDs in your geneSetSource. The geneSetSource object does not contain gene sets for ", 
              species, sep = ""))
    geneSetSource <- geneSetSource[geneSetSource$TAXID == 
            TAXID, ]
    geneSets <- by(geneSetSource$GENEID, INDICES = geneSetSource$PATHWAYID, 
        FUN = list)
    geneSets <- lapply(geneSets, as.character)
    # Get the descriptions
    tempDescriptions <- geneSetSource[, c("PATHWAYID", "PATHWAYNAME")]
    tempDescriptions$PATHWAYID <- as.character(tempDescriptions$PATHWAYID)
    tempDescriptions$PATHWAYNAME <- as.character(tempDescriptions$PATHWAYNAME)
    tempDescriptions <- unique(tempDescriptions)
    descriptions <- tempDescriptions[, "PATHWAYNAME"]
    names(descriptions) <- tempDescriptions$PATHWAYID
  }
  
  # entrezIds <- sub("_at", "", featureNames(eset))
  tfidx <- sapply(geneSets, function(geneSet) {
        sum(geneSet %in% entrezIdentifiers) > 0
      })
  geneSets <- geneSets[tfidx]
  attr(geneSets, "species") <- species
  attr(geneSets, "geneSetSource") <- geneSetSource
  attr(geneSets, "descriptions") <- descriptions

  class(geneSets) <- c("geneSetMLP", class(geneSets))
  return(geneSets)
}
