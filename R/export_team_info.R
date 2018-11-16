export_team_info <- function(){
  synapser::synLogin()
  obj <- synapser::synGet('syn12548902')
  dataBlob <- rjson::fromJSON(file=obj$path)
  fxn1 <- function(x){
    if(length(x$nominatedtarget) > 0){
      x$nominatedtarget$ensembl_gene_id <- x$ensembl_gene_id
      return(x$nominatedtarget)
    }else{
      return(NULL)
    }
  }
  res <- lapply(dataBlob,fxn1)
  len1 <- sapply(res,length)
  res <- res[len1>0]

  res2 <- unlist(res)
  cat(res2,file='~/Desktop/targetcopy.txt',sep='\n')

}
