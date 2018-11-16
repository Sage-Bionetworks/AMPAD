meta_analysis_summary_table <- function(){
  synapser::synLogin()
  foo <- synapser::synGet('syn11914808')
  bar <- data.table::fread(foo$path,data.table=F)
  bar3 <- dplyr::select(bar,ensembl_gene_id,TE.random,fdr.random)
  bar3$ad_direction <- mapply(function(x,y){
    if(y<=0.05){
      if(x>0){
        return('Up')
      }else{
        return('Down')
      }
    }else{
      return('Unchanged')
    }},bar3$TE.random,bar3$fdr.random)
  return(bar3)
}
