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
#' @return object of class geneSetMLP which is essentially a named 
#' list of pathway categories. \cr
#' Each list component is named with the pathway ID and 
#' contains a vector of Entrez Gene identifiers 
#' related to that particular pathway.\cr
#' The object contains additionally the attributes:
#' \itemize{
#' \item{'species' and 'geneSetSource': }{\code{species} and \code{geneSetSource}
#' (as provided as input)}
#' \item{'descriptions': }{named character vector with pathway descriptions.
#' The vector is named with the pathway ID.}
#' }
#' @examples if (require(GO.db) && require(org.Mm.eg.db)){
#'   pathExamplePValues <- system.file("exampleFiles", "examplePValues.rda", package = "MLP")
#'   load(pathExamplePValues)
#'   geneSet <- getGeneSets(species = "Mouse", geneSetSource = "GOBP", entrezIdentifiers = names(examplePValues)[1:2000])
#' 	 geneSet <- getGeneSets(species = "Mouse", geneSetSource = "KEGG", entrezIdentifiers = names(examplePValues)[1:2000])
#' }
#' @importFrom utils head
#' @importFrom AnnotationDbi Term toTable
#' @export
getGeneSets <- function (species = "Mouse", geneSetSource = NULL, entrezIdentifiers){
  
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
		
		if(!requireNamespace("GO.db")){
			stop("Package 'GO.db' should be available ",
				"for the extraction of pathways from: ", 
				geneSetSource, ".")
		}
		
      switch(species, Mouse = {
			requireNamespace("org.Mm.eg.db")
            geneSetToEntrez <- as.list(org.Mm.eg.db::org.Mm.egGO2ALLEGS)
          }, Human = {
			requireNamespace("org.Hs.eg.db")
            geneSetToEntrez <- as.list(org.Hs.eg.db::org.Hs.egGO2ALLEGS)
          }, Rat = {
			requireNamespace("org.Rn.eg.db")
            geneSetToEntrez <- as.list(org.Rn.eg.db::org.Rn.egGO2ALLEGS)
          }, Dog = {
			requireNamespace("org.Cf.eg.db")
            geneSetToEntrez <- as.list(org.Cf.eg.db::org.Cf.egGO2ALLEGS)
          })
  
 	  ontology <- sub("GO", "", geneSetSource)
      switch(ontology, BP = {
            GOs <- names(as.list(GO.db::GOBPANCESTOR))
          }, MF = {
            GOs <- names(as.list(GO.db::GOMFANCESTOR))
          }, CC = {
            GOs <- names(as.list(GO.db::GOCCANCESTOR))
          })
  
      go <- geneSetToEntrez[names(geneSetToEntrez) %in% 
              GOs]
      geneSets <- lapply(go, unique)
      anyGenesInGeneSet <- ifelse(unlist(lapply(geneSets, 
                  length)) > 0, TRUE, FALSE)
      geneSets <- geneSets[anyGenesInGeneSet]
      goterms <- Term(GO.db::GOTERM)
      descriptions <- goterms[names(geneSets)]
    }else	if (geneSetSource == "KEGG") {
		
		if(!requireNamespace("KEGGREST")){
			stop("Package 'KEGGREST' should be available ",
				"for the extraction of pathways from: ", 
				geneSetSource, ".")
		}
		
	  prefix <- switch(species, 
          Mouse = "mmu",
		  Human = "hsa", 
          Rat = "rno", 
          Dog = "cfa"
      )
	  
	  # extract pathway IDs
  	  keggPathways <- KEGGREST::keggList(database = "pathway", organism = prefix)
	  
	  # extract all info for these pathways
	  # is there a more efficient way to only query the gene IDs?
	  # keggGet only gets info for 10 pathways maximum
	  idxPathways <- c(seq(from = 1, to = length(keggPathways), by = 10), length(keggPathways)+1)
	  keggDbList <- lapply(head(seq_along(idxPathways), -1), function(i){
		 idxPathSel <- seq(from = idxPathways[i], to = idxPathways[i+1]-1)	  
		 keggDbSel <- KEGGREST::keggGet(dbentries = names(keggPathways)[idxPathSel])
		 lapply(keggDbSel, "[", c("ENTRY", "NAME", "GENE"))
	  }) 
	  keggDb <- do.call(c, keggDbList)
  
	  # extract gene IDs for each pathway
	  # 'GENE' includes gene ID and description -> retain IDs only
      geneSets <- lapply(keggDb, function(x) grep("^\\d+$", x$GENE, value = TRUE))
	  
	  # extract pathway name
	  descriptions <- sapply(keggDb, "[[", "NAME")
	  # remove specie name in description:
	  org <- KEGGREST::keggList("organism")
	  idxOrg <- which(org[, which(colnames(org) == "organism")] == prefix)
	  keggSpecie <- org[idxOrg, which(colnames(org) == "species")]
	  descriptions <- sub(paste(" -", keggSpecie), "", descriptions, fixed = TRUE)
	  
	  # extract pathway IDs
	  names(geneSets) <- names(descriptions) <- sapply(keggDb, "[[", "ENTRY")
	  
    }else	if (geneSetSource == "REACTOME") {
		
		if(!requireNamespace("reactome.db")){
			stop("Package 'reactome.db' should be available ",
					"for the extraction of pathways from: ", 
					geneSetSource, ".")
		}
      pathways <- toTable(reactome.db::reactomePATHNAME2ID)
      allReactomeIDs <- ls(reactome.db::reactomePATHID2EXTID)
      switch(species, Mouse = {
            pathwaysSelectedSpecies <- pathways[grep("Mus musculus: ", iconv(pathways$path_name)), ]
            #Neem eerst de reactomeids met een mapping beperk deze vervolgens met de geselecteerde species.
            geneSets <- mget(pathwaysSelectedSpecies$reactome_id[pathwaysSelectedSpecies$reactome_id %in% allReactomeIDs], reactome.db::reactomePATHID2EXTID)
          }, Human = {
            pathwaysSelectedSpecies <- pathways[grep("Homo sapiens: ", iconv(pathways$path_name)), ]
            geneSets <- mget(pathwaysSelectedSpecies$reactome_id[pathwaysSelectedSpecies$reactome_id %in% allReactomeIDs], reactome.db::reactomePATHID2EXTID)
          }, Rat = {
            pathwaysSelectedSpecies <- pathways[grep("Rattus norvegicus: ", iconv(pathways$path_name)), ]
            geneSets <- mget(pathwaysSelectedSpecies$reactome_id[pathwaysSelectedSpecies$reactome_id %in% allReactomeIDs], reactome.db::reactomePATHID2EXTID)
          }, Dog = {
            pathwaysSelectedSpecies <- pathways[grep("Canis familiaris: ", iconv(pathways$path_name)), ]
            geneSets <- mget(pathwaysSelectedSpecies$reactome_id[pathwaysSelectedSpecies$reactome_id %in% allReactomeIDs], reactome.db::reactomePATHID2EXTID)
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
  descriptions <- descriptions[names(geneSets)]
  
  attr(geneSets, "species") <- species
  attr(geneSets, "geneSetSource") <- geneSetSource
  attr(geneSets, "descriptions") <- descriptions

  class(geneSets) <- c("geneSetMLP", class(geneSets))
  return(geneSets)
}
