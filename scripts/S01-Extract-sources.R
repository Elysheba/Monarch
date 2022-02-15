library(here)
library(RJSONIO)
source(here("../00-Utils/downloadSourceFiles.R"))

desc <- readJSONStream(here("DESCRIPTION.json"))

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
srcDir <- here("sources")

downloadSourceFiles(urls, srcDir, httpForce = T)
