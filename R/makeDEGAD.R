makeDEGAD <- function(synId = 'syn10338156'){
  synapser::synLogin()
  moduleSet <- synapser::synTableQuery(paste0("SELECT DISTINCT ModuleNameFull, Module, method, brainRegion from ",synId))$asDataFrame()
  colnames(moduleSet)[c(3:4)] <- c('ModuleMethod','ModuleBrainRegion')
  degResObj <- synapser::synGet("syn10496554")
  load(degResObj$path)
  #source('enrichmentAnalysis/run_amp_ad_enrichment.R')
  keep <- grep('^Diagnosis',names(amp.ad.de.geneSets))
  keep <- intersect(keep,grep('AD-CONTROL',names(amp.ad.de.geneSets)))
  amp.ad.de.geneSets <- amp.ad.de.geneSets[keep]
  return(amp.ad.de.geneSets)
}
