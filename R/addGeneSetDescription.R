#' Utility function which adds the biological description of the gene sets as
#' a column to the return value of the MLP function (data frame)
#' @param object object of class 'MLP' as produced by the 'MLP' function 
#' @param geneSetSource source to be used to construct the list of pathway categories; 
#' for public data sources, the user can specify a string (one of 'GOBP', 'GOMF', 'GOCC', 'KEGG' or 'REACTOME')
#' and BioC packages will be used to construct the list of pathway categories; 
#' for non-public data sources, the user can pass the pathway data as a dataframe with (at least) 
#' the following four columns: PATHWAYID, TAXID, PATHWAYNAME and GENEID. It is assumed all columns
#' are of type character. The 'geneSetSource' argument should be the same as the argument
#' provided to the getGeneSets function; defaults to NULL 
#' @return the data frame as returned by MLP enriched with an additional column geneSetDescription, providing
#' a concise description of the gene set
#' @seealso \link{MLP} 
#' @export
addGeneSetDescription <- function (object, geneSetSource = NULL){
    
    ### checks
    if (!inherits(object, "MLP")) 
        stop("'object' should be an object of class 'MLP' as produced by the MLP function")
    
    if (any(names(object) == "geneSetDescription"))
        warning("The MLP object already contains a column 'geneSetDescription'")
    
    if (is.null(geneSetSource)) 
        stop("Please provide the same source of gene sets as provided to the getGeneSets function. More info, see help.")
    
    if (!is.data.frame(geneSetSource) && (length(geneSetSource) != 1) && !(geneSetSource %in% c("GOBP", "GOMF", "GOCC", "KEGG", "REACTOME")))
        stop("Please provide the same source of gene sets as provided to the getGeneSets function. More info, see help.")  
    
    if (is.data.frame(geneSetSource) && any(!(c("PATHWAYID", "TAXID", 
                                "PATHWAYNAME", "GENEID") %in% colnames(geneSetSource)))) 
        stop("Please provide the same source of gene sets as provided to the getGeneSets function. More info, see help.")
    
    species <- attr(object, "species")
    
    ### deal with character vectors
    if (!is.data.frame(geneSetSource)){
        
        if (geneSetSource %in% c("GOBP", "GOMF", "GOCC")) {
            allGOTerms <- as.list(Term(GOTERM))
            geneSetNames <- rownames(object)
            if (!all(geneSetNames %in% names(allGOTerms))) 
                stop("Check the geneSetSource parameter and compare it to the one used in the getGeneSets function, they should be the same!")
            returnValue <- data.frame(object, geneSetDescription = unlist(allGOTerms[geneSetNames]), stringsAsFactors = FALSE)
        }
        
        if (geneSetSource == "KEGG") {
            allKEGGterms <- as.list(KEGGPATHID2NAME)
            geneSetNames <- gsub("^[[:alpha:]]{3}", "", rownames(object))
            if (!all(geneSetNames %in% names(allKEGGterms))) 
                stop("Check the geneSetSource parameter and compare it to the one used in the getGeneSets function, they should be the same!")
            returnValue <- data.frame(object, geneSetDescription = unlist(allKEGGterms[geneSetNames]), stringsAsFactors = FALSE)
        }
        
        if (geneSetSource == "REACTOME") {
            
            require(reactome.db)
            pathways <- toTable(reactomePATHNAME2ID)
            switch(species, 
                    Mouse = {pathwaysSelectedSpecies <- pathways[grep("Mus musculus: ", iconv(pathways$path_name)), ]
                        pathwaysSelectedSpecies$path_name <- gsub("Mus musculus: ", "", iconv(pathwaysSelectedSpecies$path_name))}, 
                    Human = {pathwaysSelectedSpecies <- pathways[grep("Homo sapiens: ", iconv(pathways$path_name)), ]
                        pathwaysSelectedSpecies$path_name <- gsub("Homo sapiens: ", "", iconv(pathwaysSelectedSpecies$path_name))}, 
                    Rat = {pathwaysSelectedSpecies <- pathways[grep("Rattus norvegicus: ", iconv(pathways$path_name)), ]
                        pathwaysSelectedSpecies$path_name <- gsub("Rattus norvegicus: ", "", iconv(pathwaysSelectedSpecies$path_name))}, 
                    Dog = {pathwaysSelectedSpecies <- pathways[grep("Canis familiaris: ", iconv(pathways$path_name)), ]
                        pathwaysSelectedSpecies$path_name <- gsub("Canis familiaris: ", "", iconv(pathwaysSelectedSpecies$path_name))})
            geneSetNames <- rownames(object)
            if (!all(geneSetNames %in% pathwaysSelectedSpecies$reactome_id)) 
                stop("Check the geneSetSource parameter and compare it to the one used in the getGeneSets function, they should be the same!")
            returnValue <- data.frame(object, geneSetDescription = sapply(geneSetNames, function(x) pathwaysSelectedSpecies[pathwaysSelectedSpecies$reactome_id == x, 2]), stringsAsFactors = FALSE)
        }
        ### deal with a data frame
    } else {
        if (!all(rownames(object) %in% geneSetSource$PATHWAYID)) 
            stop("Check the geneSetSource parameter and compare it to the one used in the getGeneSets function, they should be the same!")
        idx <- match(rownames(object), geneSetSource$PATHWAYID)
        returnValue <- data.frame(object, geneSetDescription = geneSetSource$PATHWAYNAME[idx], stringsAsFactors = FALSE)
    }
    
    class(returnValue) <- c("MLP", class(returnValue))
    return(returnValue)
}

