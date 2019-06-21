rm(list = ls())
gc()

collectOntology <- function(database = c()){
  if(!all(database %in% listDB()$db)){
    stop('database not in DOD, use listDB() to check')
  }
  if(length(database) >1){
    stop("When providing a (list of) ids only one database can be used")
  }
  q <- sprintf('{q(func:eq(type,"database")) @filter(eq(name, /%s/i)){~is_in{name}}}', database)
  # q <- sprintf('{q(func:regexp(name,/%s/i)) @filter(eq(type, "disease") OR eq(type, "phenotype")){name}}', database)
  q <- paste(q,collapse="\n")
  rq <- dodCall(dgraphRequest,postText = q)
  ids <- unlist(rq$result$data$q)
  ##
  return(buildDisNet(ids = ids, seed = ids))
}

setwd("~/Shared/Data-Science/Data-Source-Model-Repository/Monarch/scripts/")
source("../../00-Utils/writeLastUpdate.R")
library(XML)
library(parallel)
library(jsonlite)
library(tidyr)
library(tibble)
library(dplyr)
library(readr)
library(here)
# library(stringr)
# library(rdflib)
# library(redland)
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
# Sys.setenv(PATH = paste(Sys.getenv("PATH"),"/home/lfrancois/bin/",sep = ":"))
# system(paste("robot convert --input ",file.path(sdir,"mondo.owl"),
#              " --output ",file.path(sdir,"mondo.json"), sep = ""))
readJson <- jsonlite::fromJSON(txt = file.path(sdir,"mondo.json"))

checkSyn <- do.call(rbind,lapply(readJson$graphs$nodes[[1]]$meta$synonyms, function(x) x))
table(checkSyn$pred)

###########################################
## nodes (id, def, name, xref, label)
nodesJson <- lapply(1:nrow(readJson$graphs$nodes[[1]]),
                    function(i){
                      ## id
                      id <- gsub("_",":",gsub(".*/","",readJson$graphs$nodes[[1]]$id[[i]]))
                      ## def
                      def <- readJson$graphs$nodes[[1]]$meta$definition[i,"val"]
                      ## Cross-ref
                      Xref <- readJson$graphs$nodes[[1]]$meta$xrefs[[i]]$val
                      ## synonym
                      name <- readJson$graphs$nodes[[1]]$meta$synonyms[[i]]$val
                      ## Definition
                      lbl <- readJson$graphs$nodes[[1]]$lbl[[i]]
                      ## df
                      df1 <- data.frame(
                        id = id,
                        def = def,
                        label = lbl,
                        stringsAsFactors = FALSE)
                      if(length(Xref) == 0){
                        df2 <- NULL
                      }else{
                        df2 <- data.frame(
                          id = id,
                          xref = Xref,
                          stringsAsFactors = FALSE)
                      }
                      if(length(name) == 0){
                        df3 <- NULL
                      }else{
                        df3 <- data.frame(
                          id = id,
                          syn = name,
                          stringsAsFactors = FALSE)
                      }
                      return(list(id=df1,xref=df2,syn=df3))
                    }
             )
id <- do.call(rbind,lapply(nodesJson,function(x){x$id} ))
xref <- do.call(rbind,lapply(nodesJson,function(x) x$xref))
syn <- do.call(rbind,lapply(nodesJson,function(x) x$syn))

## edges (parents)
edgesJson <- readJson$graphs$edges[[1]]
table(edgesJson$pred)
edgesJson <- edgesJson[which(edgesJson$pred %in% c("is_a")),]
edgesJson$sub <- gsub("_",":",gsub(".*/","",edgesJson$sub))
edgesJson$obj <- gsub("_",":",gsub(".*/","",edgesJson$obj))
dim(edgesJson)

getDescendants <- function(sp){
  direct <- edgesJson[which(edgesJson$obj==sp),"sub"]
  descendants <- direct
  level <- 0
  dLev <- c()
  for(d in direct){
    dDesc <- getDescendants(d)
    dLev <- c(dLev, dDesc$level)
    descendants <- c(descendants, dDesc$descendants)
  }
  if(length(dLev)>0){
    level <- max(dLev)+1
  }
  return(list(descendants=unique(c(descendants,sp)), level=level))
}
## Not only disease in Orphanet, also genes, etc. Keep only child terms of "phenome" = Orphanet_C001, see https://www.ebi.ac.uk/ols/ontologies/ordo
disease <- getDescendants("MONDO:0000001")
lapply(disease,length)
dim(edgesJson)
"MONDO:0000001" %in% disease$descendants

# rm <- c("MONDO_0024571","MONDO_0021125","MONDO_0021007")

edgesJson$obj <- gsub("NCIT","NCIt",edgesJson$obj)
edgesJson$sub <- gsub("NCIT","NCIt",edgesJson$sub)
table(gsub(":.*","",edgesJson$obj))
table(gsub(":.*","",edgesJson$sub))
edgesJson <- edgesJson[!(grepl("#",edgesJson$sub) | grepl("#",edgesJson$obj)),]

######################################
## crossId
crossId <- xref[xref$id %in% disease$descendants,]
names(crossId) <- c("dbid1","dbid2")
crossId$dbid1 <- as.character(crossId$dbid1)
crossId$dbid2 <- as.character(crossId$dbid2)
table(gsub(":.*","",crossId$dbid2))
table(gsub(":.*","",crossId$dbid1))
dim(crossId)
crossId <- crossId[!(grepl("#",crossId$dbid1) | grepl("#",crossId$dbid2)),]
## Remove crossids with colon and space ": "
head(grep(": ",crossId$dbid2,value = T))
head(grep(": ",crossId$dbid1,value = T))
crossId$dbid1 <- gsub(" ", "", crossId$dbid1)
crossId$dbid2 <- gsub(" ", "", crossId$dbid2)
##
dim(crossId)

crossId$DB2 <- gsub(":.*","",crossId$dbid2)
crossId$DB1 <- gsub(":.*","",crossId$dbid1)
crossId$id2 <- gsub(".*:","",crossId$dbid2)
crossId$id1 <- gsub(".*:","",crossId$dbid1)

## Remove crossIds without a colon (e.g. definitions, ...)
head(grep(":",crossId$dbid1,invert = T,value = T))
head(grep(":",crossId$dbid2,invert = T,value = T))
crossId <- crossId[grepl(":",crossId$dbid2) & grepl(":",crossId$dbid1) ,]
dim(crossId)
## an integer is a correct disease ID
table(!is.na(as.numeric(crossId$id2)))
table(!is.na(as.numeric(crossId$id1)))
toKeep <- crossId[which(!is.na(as.numeric(crossId$id2)) &
                          !is.na(as.numeric(crossId$id1))),]
dim(toKeep)
toCheck <- crossId[-which(!is.na(as.numeric(crossId$id2)) &
                            !is.na(as.numeric(crossId$id1))),]
dim(toCheck)
## When removing prefix, an integer is a correct disease ID
table(!is.na(as.numeric(sub("^[^[:digit:]]*", "", toCheck$id2))))
table(!is.na(as.numeric(sub("^[^[:digit:]]*", "", toCheck$id1))))
toKeep <- rbind(toKeep, 
                toCheck[which(!is.na(as.numeric(sub("^[^[:digit:]]*", "", toCheck$id2))) &
                                !is.na(as.numeric(sub("^[^[:digit:]]*", "", toCheck$id1)))),])
dim(toKeep)
toCheck <- toCheck[-which(!is.na(as.numeric(sub("^[^[:digit:]]*", "", toCheck$id2))) &
                            !is.na(as.numeric(sub("^[^[:digit:]]*", "", toCheck$id1)))),]
dim(toCheck)

## Remove any DBs that are not disease DBs and DB1 can only be "EFO" or "Orphanet"
## check wrong IDs, remove weird ones still
table(toCheck$DB2)
table(toCheck$DB1)
toCheck[toCheck$DB2 == "CSP",]
toCheck[toCheck$DB2 == "https",]
toCheck[toCheck$DB2 == "ICD10",]
toCheck[toCheck$DB2 == "ICD9",]
toCheck[toCheck$DB2 == "MESH",]
toCheck[toCheck$DB2 == "UMLS",]
toCheck[toCheck$DB2 == "Wikipedia",]
toCheck[toCheck$DB2 == "ONCOTREE",]
toCheck[toCheck$DB2 == "Orphanet",]
toCheck[toCheck$DB2 == "SCTID",]

table(toKeep$DB2)
table(toKeep$DB1)
crossId <- setNames(toKeep[,c("dbid1","dbid2")],c("id1","id2"))
dim(crossId)
head(crossId)

crossId <- crossId[grep(paste("http","url","Wikidata",sep = "|"),crossId$id2,invert = T),]
table(gsub(":.*","",crossId$id2))

crossId$id2 <- gsub("MEDGEN","UMLS",crossId$id2)
crossId$id2 <- gsub("MESH","MeSH",crossId$id2)
crossId$id2 <- gsub("ORDO","ORPHA",crossId$id2)
crossId$id2 <- gsub("\\bNCI\\b","NCIt",crossId$id2)
crossId$id2 <- gsub("NCiT","NCIt",crossId$id2)
crossId$id2 <- gsub("NCIT","NCIt",crossId$id2)
crossId$id2 <- gsub("SNOWMEDCT","SNOMEDCT",crossId$id2)
crossId$id2 <- gsub("SCTID","SNOMEDCT",crossId$id2)
crossId$id2 <- gsub("UMLS_CUI","UMLS",crossId$id2)
crossId$id2 <- gsub("Orphanet","ORPHA",crossId$id2)
crossId$id2 <- gsub("MEDDRA","MedDRA",crossId$id2)
# crossId$id2 <- gsub("NCI Metathesaurus","NCIt",crossId$id2)
crossId$DB2 <- gsub(":.*","",crossId$id2)
crossId$DB1 <- gsub(":.*","",crossId$id1)
table(crossId$DB2)
table(crossId$DB1)

## Remove self references
crossId[which(crossId$id1 == crossId$id2),]
# crossId <- crossId[-which(crossId$id1 == crossId$id2),]

######################################
## entryId
entryId <- id[id$id %in% disease$descendants,]
table(gsub(":.*","",entryId$id))
grep(":",entryId$id,invert = T,value = T)
entryId <- entryId[grepl(":",entryId$id),,drop = F]
table(gsub(":.*","",entryId$id))
entryId <- entryId[grep("#",entryId$id, invert = T, value = F),,drop = FALSE]
table(gsub(":.*","",entryId$id))
entryId$DB <- gsub(":.*","",entryId$id)
entryId <- entryId[,c("DB","id","def")]
## Empty definition to NA
nc <- nchar(entryId$def)
head(table(nc), n = 20)
entryId[which(nc < 4),]
entryId[which(nc < 4),"def"] <- NA
## Check characters for \t, \n, \r and put to ASCII
entryId$def <- iconv(x = entryId$def,to="ASCII//TRANSLIT")
entryId$def <- gsub(paste("\n","\t","\r", sep = "|")," ",entryId$def)
entryId$def <- gsub("\"","'",entryId$def)
entryId$def <- gsub("\\\\","",entryId$def)
table(unlist(sapply(entryId$def, strsplit, split = "")))

## Check duplicated records
dim(entryId)
length(unique(entryId[,"id"]))

## all crossId$id1 in entryId
table(crossId$id1 %in% entryId$id)

######################################
## idNames
idNames <- syn[syn$id %in% disease$descendants,]
table(gsub(":.*","",idNames$id))
grep(":",idNames$id,invert = T,value = T)
idNames <- idNames[grepl(":",idNames$id),,drop = F]
table(gsub(":.*","",idNames$id))
grep("#",idNames$id,value = T)
idNames <- idNames[grep("#",idNames$id,invert = T, value = F),]
idNames <- idNames[!is.na(idNames$syn),]
idNames$canonical <- FALSE
## Labels
lbl <- id[id$id %in% disease$descendants,c("id","label")]
table(gsub(":.*","",lbl$id))
grep(":",lbl$id,invert = T,value = T)
lbl <- lbl[grepl(":",lbl$id),,drop = F]
table(gsub(":.*","",lbl$id))
grep("#",lbl$id,value = T)
lbl <- lbl[grep("#",lbl$id,invert = T, value = F),]
table(gsub(":.*","",lbl$id))
lbl <- lbl[!is.na(lbl$label),]
lbl$canonical <- TRUE

## 
idNames <- idNames %>%
  as_tibble() %>%
  bind_rows(lbl %>% select(id, syn = label, canonical)) %>%
  mutate(DB = gsub(":.*","", id))
## unique
dim(unique(idNames))
idNames <- idNames[order(idNames$canonical,decreasing = T),]
idNames <- unique(idNames)
dim(idNames)

## Check characters for \t, \n, \r and put to ASCII
idNames$syn <- iconv(x = idNames$syn,to="ASCII//TRANSLIT")
idNames$syn <- gsub(paste("\n","\t","\r", sep = "|")," ",idNames$syn)
idNames$syn <- gsub("\"","'",idNames$syn)
table(unlist(sapply(idNames$syn, strsplit, split = "")))

## Remove empty syn
table(is.na(idNames$syn))
## idNames <- idNames[!is.na(idNames$syn),]
dim(idNames)

## Remove empty names, ifany
nc <- nchar(idNames$syn)
table(nc)
head(idNames[which(nc == 6),])
head(idNames[which(nc < 3 & idNames$canonical == FALSE),])
## Remove names of 0 or 1 character long
idNames[which(nc == 0),]
idNames[which(nc == 1),]
# idNames <- idNames[-which(nc == 0),]

## All idnames in entryid
table(idNames$id %in% entryId$id)

######################################
## parentId
parentId <- edgesJson[which(edgesJson$obj %in% disease$descendants),c("sub","obj")]
table(gsub(":.*","",parentId$sub))
table(gsub(":.*","",parentId$obj))
grep(":",parentId$sub,invert = T,value = T)
grep(":",parentId$obj,invert = T,value = T)
grep("#",parentId$sub,value = T)
grep("#",parentId$obj,value = T)
names(parentId) <- c("id","parent")
parentId$DB <- gsub(":.*","",parentId$id)
parentId$pDB <- gsub(":.*","",parentId$parent)

## all parentId in entryId
table(parentId$id %in% entryId$id)
table(parentId$parent %in% entryId$id)
## "Phenome" not in entryId --> OK
parentId[!(parentId$parent %in% entryId$id),]

# parentId <- parentId[which(parentId$id %in% entryId$id & parentId$parent %in% entryId$id),]
# parentId <- parentId[parentId$id %in% nodesJson$id,]
# parentId <- parentId[parentId$parent %in% nodesJson$id,]

######################################
## Mondo to Phenotype
mondoHp <- read_tsv(here("sources","disease_phenotype.all.tsv"), col_names = T, col_types = cols(.default = "c")) 
mondoHp <- mutate(mondoHp,
                  DB = gsub(":.*","", subject),
                  id = gsub(".*:","", subject),
                  hp = object) %>%
  select(DB, id, hp) %>%
  filter(grepl("HP", hp))
mondoHp$DB <- gsub("Orphanet","ORPHA",mondoHp$DB)
table(mondoHp$DB)

#######################################
crossId$id1 <- gsub(".*:","",crossId$id1)
crossId$id2 <- gsub(".*:","",crossId$id2)
entryId$id <- gsub(".*:","",entryId$id)
parentId$id <- gsub(".*:","",parentId$id)
parentId$parent <- gsub(".*:","",parentId$parent)
idNames$id <- gsub(".*:","",idNames$id)
mondoHp$hp <- gsub(".*:", "", mondoHp$hp)

############################
Monarch_idNames <- idNames[,c("DB","id","syn","canonical")]
Monarch_parentId <- parentId[,c("DB","id","pDB","parent")]
Monarch_crossId <- crossId[,c("DB1","id1","DB2","id2")]
Monarch_entryId <- entryId[,c("DB","id","def")]
Monarch_hp <- mondoHp[,c("DB","id","hp")]

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
    sep="\t",
    row.names=FALSE, col.names=TRUE,
    quote=TRUE,
    qmethod = "double"
  )
}
writeLastUpdate()

##############################################################
## Check model
source("../../00-Utils/autoCheckModel.R")
