rm(list = ls())
gc()

setwd("~/Shared/Data-Science/Data-Source-Model-Repository/Monarch/scripts/")

library(XML)
library(parallel)
library(jsonlite)
library(rdflib)
library(redland)
source("../../00-Utils/df2rdf.R")

##
mc.cores <- 55
sdir <- "../sources"
ddir <- "../data"

###############################################################################@
## Source information ----
###############################################################################@

sfi <- read.table(
   file.path(sdir, "ARCHIVES/ARCHIVES.txt"),
   sep="\t",
   header=T,
   stringsAsFactors=FALSE
)
Monarch_sourceFiles <- sfi[which(sfi$inUse), c("url", "current")]


###############################################################################@
## Data from mondo.owl ----
###############################################################################@
## Convert OWL to JSON
if(!file.exists(file.path(sdir,"mondo.json"))){
  Sys.setenv(PATH = paste(Sys.getenv("PATH"),"/home/lfrancois/bin/",sep = ":"))
  system(paste("robot convert --input ",file.path(sdir,"mondo.owl"),
               " --output ",file.path(sdir,"mondo.json"), sep = ""))
}
readJson <- jsonlite::fromJSON(txt = file.path(sdir,"mondo.json"))

###########################################
## nodes (id, def, name, xref, label)
nodesJson <- do.call(rbind,
                     lapply(1:nrow(readJson$graphs$nodes[[1]]),
                            function(i){
                              ## id
                              id <- gsub("_",":",gsub(".*/","",readJson$graphs$nodes[[1]]$id[[i]]))
                              ## def
                              def <- paste(readJson$graphs$nodes[[1]]$meta$definition[i,"val"], collapse = ", ")
                              ## Cross-ref
                              Xref <- paste(readJson$graphs$nodes[[1]]$meta$xrefs[[i]]$val,collapse = ", ")
                              ## synonym
                              name <- paste(readJson$graphs$nodes[[1]]$meta$synonyms[[i]]$val,collapse = ", ")
                              ## Definition
                              lbl <- readJson$graphs$nodes[[1]]$lbl[[i]]
                              ## df
                              df <- data.frame(
                                id = id,
                                Xref = Xref,
                                def = def,
                                name = name,
                                label = lbl,
                                stringsAsFactors = FALSE)
                              return(df)
                            }
                     )
)

## edges (parents)
edgesJson <- readJson$graphs$edges[[1]]
edgesJson <- edgesJson[which(edgesJson$pred %in% c("is_a")),]
edgesJson$sub <- gsub("_",":",gsub(".*/","",edgesJson$sub))
edgesJson$obj <- gsub("_",":",gsub(".*/","",edgesJson$obj))
edgesJson$obj <- gsub("NCIT","NCIt",edgesJson$obj)
edgesJson$sub <- gsub("NCIT","NCIt",edgesJson$sub)

######################################
## crossId
crossId <- unique(nodesJson[,c("id","Xref")])
crossIdList <- strsplit(crossId$Xref, split = ", ")
names(crossIdList) <- crossId$id
crossId <- stack(crossIdList)
names(crossId) <- c("id2","id1")
crossId$id1 <- as.character(crossId$id1)
crossId$id2 <- as.character(crossId$id2)
crossId <- crossId[!grepl("#",crossId$id1),]
crossId$id2 <- gsub("MESH","MeSH",crossId$id2)
crossId$id2 <- gsub("ORDO","ORPHA",crossId$id2)
crossId$id2 <- gsub("NCIT","NCIt",crossId$id2)
crossId$id2 <- gsub("NCiT","NCIt",crossId$id2)
crossId$id2 <- gsub("SNOWMEDCT","SNOMEDCT",crossId$id2)
crossId$id2 <- gsub("SCTID","SNOMEDCT",crossId$id2)
crossId$id2 <- gsub("UMLS","MedGen",crossId$id2)
crossId$id2 <- gsub("Orphanet","ORPHA",crossId$id2)
crossId$id2 <- gsub("MEDGEN","MedGen",crossId$id2)
# crossId$id2 <- gsub("NCI Metathesaurus","NCIt",crossId$id2)
crossId$DB2 <- gsub(":.*","",crossId$id2)
crossId$DB1 <- gsub(":.*","",crossId$id1)

######################################
## CrossId and masterNode
masterNode <- data.frame(subject = paste0("MN_",1:length(unique(crossId$id1))),
                  predicate = "is_masterNode",
                  object = unique(crossId$id1),
                  datatype = "uid",
                  stringsAsFactors = F)
dgraph_masterNode <- data.frame(subject = crossId$id1,
                                predicate = "is_masterNode",
                                object = crossId$id2,
                                datatype = "uid",
                                stringsAsFactors = F)
dgraph_masterNode$subject <- masterNode$subject[match(dgraph_masterNode$subject,masterNode$object)]
dgraph_masterNode <- rbind(masterNode,dgraph_masterNode)
attr(dgraph_masterNode, "object_type") <- "uid"

######################################
## entryId
entryId <- nodesJson[c("id")]
entryId <- entryId[grep("#",entryId$id, invert = T, value = F),,drop = FALSE]
entryId$DB <- gsub(":.*","",entryId$id)
entryId <- entryId[,c("DB","id")]
entryId$definition <- nodesJson$def[match(entryId$id,nodesJson$id)]
entryId$definition <- ifelse(entryId$definition == "NA",NA,entryId$definition)
entryId$definition <- tolower(entryId$definition)
entryId$definition <- gsub("[[:punct:]]"," ",entryId$definition)
entryId$definition <- iconv(x = entryId$definition,to="ASCII//TRANSLIT")
entryId$definition <- gsub("\n"," ",entryId$definition)
# entryId <- unique(entryId)

######################################
## idNames
idNames <- unique(nodesJson[,c("id","name")])
idNames <- idNames[grep("#",idNames$id,invert = T, value = F),]
idNamesList <- strsplit(idNames$name, split = ", ")
names(idNamesList) <- idNames$id
idNames <- stack(idNamesList)
names(idNames) <- c("name","id")
## Labels
lbl <- unique(nodesJson[,c("label","id")])
lbl <- lbl[grep("#",lbl$id,invert = T, value = F),]
## 
idNames <- rbind(idNames,setNames(lbl, nm = names(idNames)))
idNames$DB <- gsub(":.*","",idNames$id)
idNames <- idNames[!is.na(idNames$name),]
idNames$canonical <- ifelse(idNames$name %in% lbl$label, TRUE, FALSE)
idNames$name <- tolower(idNames$name)
idNames$name <- gsub("[[:punct:]]"," ",idNames$name)
idNames$name <- iconv(x = idNames$name,to="ASCII//TRANSLIT")
idNames$name <- gsub("\n"," ",idNames$name)
idNames <- unique(idNames)

######################################
## parentId
parentId <- edgesJson[,c("sub","obj")]
names(parentId) <- c("id","parent")
parentId$DB <- gsub(":.*","",parentId$id)
parentId$pDB <- gsub(":.*","",parentId$parent)
parentId <- parentId[which(parentId$id %in% entryId$id & parentId$parent %in% entryId$id),]
parentId <- parentId[parentId$id %in% nodesJson$id,]
parentId <- parentId[parentId$parent %in% nodesJson$id,]

#######################################
crossId$id1 <- gsub(".*:","",crossId$id1)
crossId$id2 <- gsub(".*:","",crossId$id2)
entryId$id <- gsub(".*:","",entryId$id)
parentId$id <- gsub(".*:","",parentId$id)
parentId$parent <- gsub(".*:","",parentId$parent)
idNames$id <- gsub(".*:","",idNames$id)

############################
Monarch_idNames <- idNames[,c("DB","id","name","canonical")]
Monarch_parentId <- parentId[,c("DB","id","pDB","parent")]
Monarch_crossId <- crossId[,c("DB1","id1","DB2","id2")]
Monarch_entryId <- entryId[,c("DB","id","definition")]

############################
## Write tables
toSave <- grep("^Monarch[_]", ls(), value=T)
for(f in toSave){
  message(paste("Saving", f))
  ## Ensure unicity
  assign(f, get(f))
  if(length(names(f))==0){
    f <- unique(f)
  }
  ##
  write.table(
    get(f),
    file=file.path(ddir, paste(f, ".txt", sep="")),
    sep="|",
    row.names=FALSE, col.names=TRUE,
    quote=FALSE
  )
}
