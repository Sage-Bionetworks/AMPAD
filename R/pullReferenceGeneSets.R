pullReferenceGeneSets <- function(url){
  #pull gene set from Enrichr
  foo <- RCurl::getURL(url)

  #split by carriage returns
  fooSplit <- strsplit(foo,"\n")

  #unlist
  fooSplit <- unlist(fooSplit)

  #split by tab
  fooSplit <- strsplit(fooSplit,"\t")

  #name them by the first element
  catNames <- lapply(fooSplit,function(x) x[1])
  names(fooSplit) <- catNames

  #drop the first two elements
  fooSplit <- lapply(fooSplit,function(x) x[-c(1,2)])

  fxn1 <- function(x){
    library(dplyr)
    strsplit(x,',') %>%
      sapply(function(y) y[1]) %>%
      return
  }
  #remove 1.0 from elements
  fooSplit2 <- lapply(fooSplit,fxn1)
  names(fooSplit2) <- names(fooSplit)
  return(fooSplit2)
}
