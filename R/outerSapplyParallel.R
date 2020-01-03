outerSapplyParallel <- function(FUN,X,Y,...){
  library(dplyr)
  library(parallel)
  no_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(no_cores)
  #pb <- progress::progress_bar$new(total=(length(X)/no_cores))
  internal <- function(X,Y,FUN,...){
    library(dplyr)
    #pb$tick()
    Y %>%
      sapply(FUN,X,...) %>%
      return
  }
  cl %>%
    parallel::parSapply(X,internal,Y,FUN,...) %>%
    return
}
