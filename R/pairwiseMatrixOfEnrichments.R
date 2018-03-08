pairwiseMatrixOfEnrichments = function(synId){
  synapseClient::synapseLogin()

  ##pull aggregate modules
  aggregateModules <- rSynapseUtilities::loadFullTable(synId)

  #View(aggregateModules)

  mats <- utilityFunctions::pairwiseMatrixOfEnrichments(key=aggregateModules$Module,
                                                        value=aggregateModules$GeneID)

  #get rid of underflow issues
  mats$pval <- mats$pval + 1e-300

  pheatmap::pheatmap(-log10(mats$pval))
}
