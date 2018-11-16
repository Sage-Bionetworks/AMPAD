pullFullManifests = function(degSyn,moduleSyn){
  synapseClient::synapseLogin()
  moduleTable <- synapseClient::synTableQuery(paste0("SELECT * FROM ",moduleSyn))@values
  load(synapseClient::synGet(degSyn)@filePath)
  dataList <- list()
  dataList$IndividualModules <- moduleTable
  dataList$DEGs <- amp.ad.de.geneSets
  return(dataList)
}
