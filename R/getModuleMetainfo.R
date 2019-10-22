getModuleMetainfo <- function(synId){
  synapser::synLogin()
  moduleSet <- synapser::synTableQuery(paste0("SELECT DISTINCT ModuleNameFull, Module, method, brainRegion from ",synId))$asDataFrame()
  colnames(moduleSet)[c(3:4)] <- c('ModuleMethod','ModuleBrainRegion')
  return(moduleSet)
}
