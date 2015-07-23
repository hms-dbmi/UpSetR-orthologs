library(UpSetR)
library(biomaRt)
library(data.table)
# UpSetRensembl <- function(organisms){
#   species <- read.csv("genes.csv")
#   species <- species[which(species$description %in% organisms),]
#   species <- head(as.character(species$dataset), 5)

  
#   sixSpecies <- function(species){
#     if(length(species) > 6){warning("6 species maximum")}
#     else{
#       datasets <- list()
#       for(i in seq_along(species)){
#         datasets[[i]] <- getEnsemblMarts(species[i])
#       }
#     }
#     return(datasets)
#   }
  
  # data <- sixSpecies(species)
  
  
  # namesData <- getNames(species)
  
#   grabData <- function(data, namesData, species){
#     for(i in seq_along(data)){
#       data[[i]] <- getBM(attributes = c("ensembl_gene_id",namesData[[i]]), mart = data[[i]])
#       new_colname <- paste(unlist(strsplit(species[i], split = "_"))[1], "_homolog_ensembl_gene", sep = "", collapse="")
#       colnames(data[[i]])[1] <- new_colname
#       data[[i]] <- data[[i]][ ,order(colnames(data[[i]]))]
#     }
#     alldata <- rbindlist(data)
#     alldata <- data.frame(alldata[!duplicated(alldata),])
#     return(alldata)
#   }
#   data <- grabData(data, namesData, species)
  
  
#   data <- BioMartConverter(data)
#   
#   upset(data, order.matrix="freq", empty.intersections = "on")
# }

getEnsemblMarts <- function(dSet){
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL",host = "useast.ensembl.org", dataset=dSet)
}

getNames <- function(species){
  allnames <- list()
  cnames <- c()
  for(i in seq_along(species)){
    current <- species[i]
    names <- species[which(!species %in% current)]
    for(j in seq_along(names)){
      name <- strsplit(names[j], split = "_")
      cnames[j] <- unlist(name)[1]
      cnames[j] <- paste(cnames[j], "_homolog_ensembl_gene", sep = "", collapse = "")
    }
    allnames[[i]] <- cnames
  }
  return(allnames)
}

# grabData <- function(data, namesData, species){
#   for(i in seq_along(data)){
#     data[[i]] <- getBM(attributes = c("ensembl_gene_id",namesData[[i]]), mart = data[[i]])
#     new_colname <- paste(unlist(strsplit(species[i], split = "_"))[1], "_homolog_ensembl_gene", sep = "", collapse="")
#     colnames(data[[i]])[1] <- new_colname
#     data[[i]] <- data[[i]][ ,order(colnames(data[[i]]))]
#   }
#   alldata <- rbindlist(data)
#   alldata <- data.frame(alldata[!duplicated(alldata),])
#   return(alldata)
# }
ManualBioMartConverter <- function(data){
  data <- data[which(!duplicated(data)), ]
  sets <- grep(("gene"), tolower(colnames(data)))
  setnames <- colnames(data[c(sets)]) 
  attnames <- NULL
  if(length(setnames) != length(data)){
    attnames <-colnames(data[-c(sets)])
  }
  names <- c(attnames, setnames)
  data <- data.frame(apply(data[sets], 1:2, function(x) if(isTRUE(x == "")){x <- 0} else{x=1}))
  #cbind(data[ ,!(colnames(data) %in% colnames(data[,sets]))],
  colnames(data) <- names
  return(data)
}

BioMartConverter <- function(data){
  data <- data[which(!duplicated(data)), ]
  sets <- grep("homolog", tolower(colnames(data)))
  setnames <- colnames(data[c(sets)]) 
  attnames <- NULL
  if(length(setnames) != length(data)){
    attnames <-colnames(data[-c(sets)])
  }
  names <- c(attnames, setnames)
  data <- data.frame(apply(data[sets], 1:2, function(x) if(isTRUE(x == "")){x <- 0} else{x=1}))
  #cbind(data[ ,!(colnames(data) %in% colnames(data[,sets]))],
  colnames(data) <- names
  return(data)
}

orthologs_manual <- function(datalist){
  if(length(datalist) > 1){
  names <- list()
  for(i in seq_along(datalist)){
    colname <- colnames(datalist[[i]])
    names[[i]] <- colname
  }
  
  names <- unique(unlist(names))
  names <- names[-match("Ensembl.Gene.ID", names)]
  for(i in seq_along(datalist)){
    col <- match("Ensembl.Gene.ID", colnames(datalist[[i]]))
    sets <- grep("Gene", colnames(datalist[[i]]))
    sets <- colnames(datalist[[i]])[sets]
    sets <- sets[-col]
    missing <- names[!(names %in% sets)]
    colnames(datalist[[i]])[col] <- missing
    datalist[[i]] <- datalist[[i]][,order(colnames(datalist[[i]]))]
  }
  
  datalist <- rbindlist(datalist)
  datalist <- data.frame(datalist[!duplicated(datalist),])
  }
  else{
    datalist <- datalist[!duplicated(datalist), ]
    # datalist <- ManualBioMartConverter(data.frame(datalist))
   return(data.frame(datalist)) 
  }
}