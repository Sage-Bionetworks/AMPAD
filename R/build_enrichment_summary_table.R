build_enrichment_summary_table <- function(moduleName){
  synapseClient::synapseLogin()
  foobar <- synapseClient::synTableQuery('select * from syn11954640')@values

  foobar <- dplyr::mutate(foobar,adj.pval = p.adjust(fisherPval,method='fdr'))
  foobar <- dplyr::filter(foobar,adj.pval <= 0.05)
  foobar <- dplyr::arrange(foobar,desc(fisherOR))
  foobar2 <- dplyr::filter(foobar,ModuleNameFull == moduleName)
  foobar3 <- dplyr::group_by(foobar2,geneSet)
  foobar4 <- dplyr::summarise(foobar3,top5 = paste0(category[1:min(5,length(category))],collapse=','))

}
