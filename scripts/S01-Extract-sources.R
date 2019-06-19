setwd("~/Shared/Data-Science/Data-Source-Model-Repository/Monarch/scripts/")

library(RJSONIO)
source("../../00-Utils/downloadSourceFiles.R")

desc <- readJSONStream("../DESCRIPTION.json")

sourceFiles <- desc$"source files"
urls <- unlist(lapply(
   sourceFiles,
   function(sf){
      toRet <- sf$"URL template"
      names(toRet) <- sf$"name"
      return(toRet)
   }
))
# urls["mondo.owl"] <- sprintf(
#    urls["mondo.owln"],
#    format(Sys.Date(), "%Y-%m")
# )
srcDir <- "../sources"

downloadSourceFiles(urls, srcDir, httpForce = T)
