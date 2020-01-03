computePairwiseMatrixOfEnrichments = function(key,value){
  #get unique keys
  unique_keys <- unique(key)
  unique_values <- unique(value)
  #create list of key value
  #df <- data.frame(keys=key,values=value,stringsAsFactors=F)

  #listify
  dfList <- lapply(unique_keys,
                   AMPAD::listify,
                   value,
                   key)
  names(dfList) <- unique_keys
  res <- list()

  library(dplyr)

  res$pval <- AMPAD::fisherWrapperPval %>%
    AMPAD::outerSapplyParallel(dfList,
                                          dfList,
                                          unique_values)

  res$or <- AMPAD::fisherWrapperOR %>%
    AMPAD::outerSapplyParallel(dfList,
                                          dfList,
                                          unique_values)
  diag(res$or) <- 0

  return(res)
}
