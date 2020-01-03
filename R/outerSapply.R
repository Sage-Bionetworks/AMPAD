outerSapply <- function(FUN,X,Y,...){
  require(dplyr)
  pb <- progress::progress_bar$new(total=length(X))
  internal <- function(X,Y,FUN,...){
    pb$tick()
    return(Y%>% sapply(FUN,X,...))
  }
  return(X %>% sapply(internal,Y,FUN,...))
}
