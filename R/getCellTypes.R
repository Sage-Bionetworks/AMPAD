getCellTypes <- function(synId = 'syn10338156'){
  #source('enrichmentAnalysis/run_amp_ad_enrichment.R')
  moduleSet <- synapseClient::synTableQuery(paste0("SELECT DISTINCT ModuleNameFull, Module, method, brainRegion from ",synId))@values
  colnames(moduleSet)[c(3:4)] <- c('ModuleMethod','ModuleBrainRegion')
  genesets1 <- synapseClient::synGet('syn5923958')
  load(synapseClient::getFileLocation(genesets1))
  cell_enrich <- AMPAD::run_amp_ad_enrichment(GeneSets$Cell_Markers,
                                       'cell_type',
                                       manifestId=synId)
  cell_enrich <- dplyr::select(cell_enrich,
                               ModuleNameFull,
                               category,
                               geneSet,
                               fisherPval,
                               fisherOR)
  colnames(cell_enrich) <- c('ModuleNameFull',
                             'GeneSetName',
                             'GeneSetCategoryName',
                             'GeneSetAssociationStatistic',
                             'GeneSetEffect')
  cell_enrich$GeneSetADLinked <- rep(FALSE,nrow(cell_enrich))
  cell_enrich$GeneSetBrainRegion <- rep(NA,nrow(cell_enrich))
  cell_enrich$GeneSetDirectionAD <- rep(NA,nrow(cell_enrich))
  cell_enrichSummary <- dplyr::left_join(moduleSet,
                                         cell_enrich)
  cell_enrichSummary <- AMPAD::splitByBrainRegionAdjustPvalue(cell_enrichSummary)
  #View(cell_enrichSummary)
  return(cell_enrichSummary)
  #moduleSummary <- rbind(moduleSummary,cell_enrichSummary)
  #moduleSummary$EvidenceClass <- 'cell_type'
}
