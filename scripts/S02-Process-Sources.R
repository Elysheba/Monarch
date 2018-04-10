setwd("~/Shared/Data-Science/Data-Source-Model-Repository/Monarch/scripts/")

library(XML)
library(parallel)

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
## Data from efo.owl ----
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
crossIdList <- strsplit(crossId$Xref, split = ",")
names(crossIdList) <- crossId$id
crossId <- stack(crossIdList)
names(crossId) <- c("id2","id1")
crossId$id2 <- gsub("MESH","MeSH",crossId$id2)
crossId$id2 <- gsub("ORDO","ORPHA",crossId$id2)
crossId$id2 <- gsub("NCIT","NCIt",crossId$id2)
crossId$id2 <- gsub("NCiT","NCIt",crossId$id2)
crossId$id2 <- gsub("SNOWMEDCT","SNOMEDCT",crossId$id2)
crossId$id2 <- gsub("SCTID","SNOMEDCT",crossId$id2)
crossId$id2 <- gsub("UMLS","MedGen",crossId$id2)
crossId$id2 <- gsub("Orphanet","ORPHA",crossId$id2)
crossId$id2 <- gsub("MEDGEN","MedGen",crossId$id2)
crossId$id2 <- gsub("NCI Metathesaurus","NCIt",crossId$id2)
crossId$DB2 <- gsub(":.*","",crossId$id2)
crossId$DB1 <- gsub(":.*","",crossId$id1)

######################################
## entryId
entryId <- crossId[!duplicated(crossId$id1),c("DB1","id1")]
names(entryId) <- c("DB","id")
entryId$definition <- nodesJson$def[match(entryId$id,nodesJson$id)]

######################################
## idNames
idNames <- unique(nodesJson[,c("id","name")])
idNamesList <- strsplit(idNames$name, split = ",")
names(idNamesList) <- idNames$id
idNames <- stack(idNamesList)
names(idNames) <- c("name","id")
## Labels
lbl <- unique(nodesJson[,c("label","id")])
## 
idNames <- rbind(idNames,setNames(lbl, nm = names(idNames)))
idNames$DB <- gsub(":.*","",idNames$id)
idNames$canonical <- ifelse(idNames$name %in% lbl$label, TRUE, FALSE)

######################################
## parentId
parentId <- edgesJson[,c("sub","obj")]
names(parentId) <- c("id","parent")
parentId$DB <- gsub(":.*","",parentId$id)
parentId$pDB <- gsub(":.*","",parentId$parent)

#######################################
crossId$id1 <- gsub(".*:","",crossId$id1)
crossId$id2 <- gsub(".*:","",crossId$id2)
entryId$id <- gsub(".*:","",entryId$id)
parentId$id <- gsub(".*:","",parentId$id)
parentId$parent <- gsub(".*:","",parentId$parent)
idNames$id <- gsub(".*:","",idNames$id)

############################
Monarch_idNames <- idNames[,c("DB","id","name")]
Monarch_parentId <- parentId[,c("DB","id","pDB","parent")]
Monarch_crossId <- crossId[,c("DB1","id1","DB2","id2")]
Monarch_entryid <- entryId[,c("DB","id")]

############################
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
    sep="\t",
    row.names=FALSE, col.names=TRUE,
    quote=FALSE
  )
}




# 
# ###############################################################################@
# ## Data from mondo.owl ----
# ###############################################################################@
# # owl <- readLines(file.path(sdir,sfi$file))
# owl <- readLines(file.path("../../../Therapeutic Landscape/Data sources/mondo.owl"))
# 
# ###########################
# ## Basic information
# starts <- grep("<owl:Class rdf:about=\"http.*>",owl)
# ends <- c(starts[-1]-1, length(owl))
# ends <- apply(data.frame(starts,ends),1,
#               function(x){
#                 s <- x[1]
#                 e <- grep("<owl:Axiom>",owl[x[1]:x[2]])[1] ## Skipping Axiom class input
#                 if(is.na(e) == TRUE){ e <- x[2]
#                 }else{e <- x[1] + e -1 }
#                 # print(c(x,e))
#                 return(e)
#               })
# 
# MonarchDef <- do.call(rbind, apply(
#   data.frame(starts, ends),
#   1,
#   function(x){
#     termDesc <- trimws(owl[as.numeric(x[1]):as.numeric(x[2]-1)])
#     ## ID
#     fn <- paste("<owl:Class rdf:about=\"http://www.orpha.net/ORDO/",
#                 "<owl:Class rdf:about=\"http://purl.obolibrary.org/obo/",sep="|")
#     id <- sub(paste("\">","\\\"/>",sep="|"),"", sub(fn, "", grep(fn, termDesc, value=T)))
#     # if(length(id)==0) id <- NA
#     ## Equivalent ID
#     fn <- "^<owl:equivalentClass rdf"
#     equiId <- sub("\"/>","",sub(".*obo/", "", grep(fn, termDesc, value=T)))
#     if(length(equiId) == 0){ equiId <- NA
#     } else {
#       equiId <- paste(unique(equiId), collapse=", ")}
#     ## Alternative ID
#     fn <- "^<oboInOwl:hasDbXref rdf"
#     crossId <- sub(":","_",sub("</oboInOwl:hasDbXref>","",gsub(".*\\\">","",grep(fn,termDesc,value = T))))
#     if(length(crossId) == 0){ crossId <- NA
#     } else {
#       crossId <- paste(unique(crossId), collapse=", ")}
#     ## Labels
#     fn <- paste("^<rdfs:label rdf")
#     labels <- gsub("</rdfs:label>","",
#                    gsub(".*\\\">","",grep(fn,termDesc,value = T)))
#     if(length(labels) == 0){ labels <- NA 
#     } else{
#       labels <- paste(unique(labels), collapse=", ")}  
#     ## Exact synonym
#     fn <- "<oboInOwl:hasExactSynonym rdf"
#     esyn <- sub("</oboInOwl:hasExactSynonym>",
#                 "",
#                 sub(".*\\\">","",grep(fn,termDesc,value = T)))
#     if(length(esyn) == 0){esyn <- NA 
#     } else{
#       esyn <- paste(unique(esyn), collapse=", ")}
#     ## Related synonym
#     fn <- "^<oboInOwl:hasRelatedSynonym rdf"
#     rsyn <- sub("</oboInOwl:hasRelatedSynonym>",
#                 "",
#                 sub(".*\\\">","",grep(fn,termDesc,value = T)))
#     if(length(rsyn) == 0){rsyn <- NA 
#     } else{
#       rsyn <- paste(unique(rsyn), collapse=", ")}
#     ## Parents - matching only for MONDO IDs
#     fn <- "<rdfs:subClassOf rdf:resource=\"http://purl.obolibrary.org/obo/"
#     parents <- sub("\"/>","",sub(".*obo/","",grep(fn,termDesc,value = T)))
#     if(length(parents) == 0){ parents <- NA
#     }else{
#       parents <- paste(unique(parents),collapse = ", ")
#     }
#     ## Definition
#     fn <- c("<obo:IAO_0000115 rdf")
#     def <- sub("</obo:IAO_0000115>","",sub(".*\\\">","",grep(fn,termDesc,value = T)))
#     if(length(def) == 0){ def <- NA 
#     }else{
#       def <- paste(unique(def), collapse=", ")}
#     # return(def)
#     return(data.frame(
#       id = id,
#       equiId = equiId,
#       crossId = crossId,
#       exactSyn = esyn,
#       relSyn = rsyn,
#       def = def,
#       labels = labels, 
#       parents = parents, 
#       stringsAsFactors=F)
#     )
#   }
# ))
# 
# 
# ## There are 50.030 unique entries, however there are both MONDO entries as other ontology IDs. 
# ## The other ontologies report no or one MONDO ID as equivalent ID (one reports a GO id as equivalent id (CHEBI:8058 - GO:0005623)).
# ## In addition, they may report original DB parent information.
# ## All these other ontology entries appear as equivalent and/or alternative ID for a MONDO entry except for (CHEBI:8058 - GO:0005623)
# # a <- MonarchDef[!grepl("MONDO", MonarchDef$id),]
# # b <- strsplit(a$equiId,", ")
# # names(b) <- a$id
# # a <- stack(b)
# # table(sub("[:].*","",a$values), useNA = c("ifany"))
# # a <- a[is.na(a$values) == FALSE,]
# # 
# # mondo <- MonarchDef[grepl("MONDO", MonarchDef$id),]
# # b <- strsplit(mondo$crossId,", ")
# # names(b) <- mondo$id
# # mondo <- stack(b)
# # unique(mondo$values[grep("MONDO",mondo$values)]) ## MONDO is in alternative IDs
# # mondo[which(mondo$values == "MONDO:0005600"),] ## Some MONDO IDs refer to other MONDO ids that don't exist in the db
# # table(sub("[:].*","",mondo$values),useNA = c("ifany")) ## 110941 entries for alternative IDs
# # 
# # ## check whether all non-Mondo ID entries also occur as alternative ID for a Mondo ID
# # mondo$ind <- as.character(mondo$ind)
# # c <- a[which(a$values %in% mondo$ind),]
# # aa <- a[-which(a$values %in% c$values), ]
# # c <- c[is.na(c$values) == FALSE,]
# 
# # crossId <- unique(MonarchDef[, c("id", "crossId")])
# # MonarchDef <- MonarchDef[, setdiff(colnames(MonarchDef), "crossId")]
# 
# 
# ###############################################################################@
# ## Custom information ----
# ###############################################################################@
# 
# ############################
# ## crossId
# crossId <- unique(MonarchDef[, c("id", "crossId")])
# crossIdList <- strsplit(crossId$crossId, ", ")
# names(crossIdList) <- crossId$id
# crossId <- stack(crossIdList)
# colnames(crossId) <- c("id2", "id1")
# crossId$id2 <- as.character(crossId$id2)
# crossId$id1 <- as.character(crossId$id1)
# # crossId <- crossId[grepl(paste("UMLS","ICD10","Orphanet","DOID","MESH","OMIM","SCTID","NCIT","GARD","ICD9","EFO",
# #                                "SNOMEDCT","SNOWMEDCT","ORDO","OMIMPS","MeSH",sep = "_|"), crossId$id2, ignore.case = TRUE),]
# crossId$id2 <- gsub("MESH","MeSH",crossId$id2)
# crossId$id2 <- gsub("ORDO","ORPHA",crossId$id2)
# crossId$id2 <- gsub("NCIT","NCIt",crossId$id2)
# crossId$id2 <- gsub("NCiT","NCIt",crossId$id2)
# crossId$id2 <- gsub("SNOWMEDCT","SNOMEDCT",crossId$id2)
# crossId$id2 <- gsub("SCTID","SNOMEDCT",crossId$id2)
# crossId$id2 <- gsub("UMLS","MedGen",crossId$id2)
# crossId$id2 <- gsub("Orphanet","ORPHA",crossId$id2)
# crossId$DB1 <- gsub("_.*","",crossId$id1)
# crossId$DB2 <- gsub("_.*","",crossId$id2)
# crossId <- crossId[is.na(crossId$id2) == FALSE,c("DB1","id1","DB2","id2")]
# 
# MonarchDef <- MonarchDef[which(MonarchDef$id %in% crossId$id1),]
# MonarchDef$DB <- gsub("_.*","",MonarchDef$id)
# 
# ############################
# ## entryId
# entryId <- MonarchDef[!duplicated(MonarchDef$id),c("id","DB")]
# 
# ############################
# ## Parents
# parentId <- unique(MonarchDef[, c("id", "parents")])
# parentList <- strsplit(parentId$parents, ", ")
# names(parentList) <- parentId$id
# parentId <- stack(parentList)
# colnames(parentId) <- c("parent", "id")
# parentId$parent <- as.character(parentId$parent)
# parentId$id <- as.character(parentId$id)
# parentId$DB <- gsub("_.*","",parentId$id)
# parentId$pDB <- gsub("_.*","",parentId$parent)
# parentId <- parentId[which(parentId$id %in% crossId$id1 & is.na(parentId$parent) == FALSE),c("DB","id","pDB","parent")]
# 
# ############################
# ## Definition and name
# descrId <- unique(MonarchDef[, c("id", "labels","def")])
# names(descrId) <- c("id","name","definition")
# descrId$id <- as.character(descrId$id)
# descrId$name <- as.character(descrId$name)
# descrId$definition <- as.character(descrId$definition)
# descrId$DB <- gsub("_.*","",descrId$id)
# descrId <- descrId[which(descrId$id %in% crossId$id1),]
# entryId$definition <- descrId$def[match(entryId$id,descrId$id)]
# 
# ############################
# ## Names
# exactSynId <- unique(MonarchDef[, c("id", "exactSyn")])
# exactSynList <- strsplit(exactSynId$exactSyn, ", ")
# names(exactSynList) <- exactSynId$id
# exactSynId <- stack(exactSynList)
# colnames(exactSynId) <- c("name", "id")
# exactSynId$id <- as.character(exactSynId$id)
# exactSynId$name <- as.character(exactSynId$name)
# exactSynId$DB <- gsub("_.*","",exactSynId$id)
# 
# relSynId <- unique(MonarchDef[, c("id", "relSyn")])
# relSynList <- strsplit(relSynId$relSyn, ", ")
# names(relSynList) <- relSynId$id
# relSynId <- stack(relSynList)
# colnames(relSynId) <- c("name", "id")
# relSynId$id <- as.character(relSynId$id)
# relSynId$name <- as.character(relSynId$name)
# relSynId$DB <- gsub("_.*","",relSynId$id)
# 
# ## idNames
# idNames <- rbind(exactSynId[is.na(exactSynId$name) == FALSE,c("DB","id","name")],
#                  relSynId[is.na(relSynId$name) == FALSE,c("DB","id","name")],
#                  descrId[is.na(descrId$name) == FALSE,c("DB","id","name")])
# idNames$canonical <- ifelse(idNames$name %in% descrId$name,TRUE, FALSE)
# 
# ## Ids
# idNames <- apply(idNames,2,function(x) gsub(".*_","",x))
# parentId <- apply(parentId,2,function(x) gsub(".*_","",x))
# entryId <- apply(entryId,2,function(x) gsub(".*_","",x))
# crossId <- apply(crossId,2,function(x) gsub(".*_","",x))
# 
# ############################
# Monarch_parentId <- parentId
# Monarch_crossId <- crossId
# Monarch_entryId <- entryId 
# Monarch_idNames <- idNames
# ############################
# 
# 
# ###############################################################################@
# ## Writing tables ----
# ###############################################################################@
# message("Writing tables...")
# message(Sys.time())
# toSave <- grep("^Monarch[_]", ls(), value=T)
# for(f in toSave){
#   message(paste("Saving", f))
#   ## Ensure unicity
#   assign(f, get(f))
#   if(length(names(f))==0){
#     f <- unique(f)
#   }
#   ##
#   write.table(
#     get(f),
#     file=file.path(ddir, paste(f, ".txt", sep="")),
#     sep="\t",
#     row.names=FALSE, col.names=TRUE,
#     quote=FALSE
#   )
# }
# message(Sys.time())
# message("... Done\n")
