makeDEGAD <- function(synId = 'syn10338156'){
  synapseClient::synapseLogin()
  moduleSet <- synapseClient::synTableQuery(paste0("SELECT DISTINCT ModuleNameFull, Module, method, brainRegion from ",synId))@values
  colnames(moduleSet)[c(3:4)] <- c('ModuleMethod','ModuleBrainRegion')
  degResObj <- synapseClient::synGet("syn10496554")
  load(degResObj@filePath)
  #source('enrichmentAnalysis/run_amp_ad_enrichment.R')
  keep <- grep('^Diagnosis',names(amp.ad.de.geneSets))
  keep <- intersect(keep,grep('AD-CONTROL',names(amp.ad.de.geneSets)))
  amp.ad.de.geneSets <- amp.ad.de.geneSets[keep]
  return(amp.ad.de.geneSets)
}
