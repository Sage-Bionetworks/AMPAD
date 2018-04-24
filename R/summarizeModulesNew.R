summarizeModulesNew <- function(moduleName,gs = NULL){
  #pull enrichments from synapse
  synapseClient::synapseLogin()

  enrichmentTable <- rSynapseUtilities::loadFullTable('syn11954640')
  if(!is.null(gs)){
    enrichmentTable2 <- dplyr::filter(enrichmentTable,geneSet == gs)
  }else {
    enrichmentTable2 <- enrichmentTable
  }
  enrichmentTable3 <- dplyr::mutate(enrichmentTable2,adj.pval = p.adjust(fisherPval,method='fdr'))
  enrichmentTable3 <- dplyr::filter(enrichmentTable3,adj.pval<=0.05)
  enrichmentTable3 <- dplyr::filter(enrichmentTable3,ModuleNameFull == moduleName)
  enrichmentTable3 <- dplyr::arrange(enrichmentTable3,desc(fisherOR))
  #f
}
