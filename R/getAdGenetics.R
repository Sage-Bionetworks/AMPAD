getAdGenetics <- function(synId='syn10338156'){
  moduleSet <- synapseClient::synTableQuery(paste0("SELECT DISTINCT ModuleNameFull, Module, method, brainRegion from ",synId))@values
  colnames(moduleSet)[c(3:4)] <- c('ModuleMethod','ModuleBrainRegion')


  #magma enrichments
  magmaResults <- synapseClient::synTableQuery("SELECT * FROM syn10380432")@values
  magmaResults <- dplyr::select(magmaResults,SET,BETA,P)
  colnames(magmaResults) <- c('ModuleNameFull',
                              'GeneSetEffect',
                              'GeneSetAssociationStatistic')
  magmaResults$GeneSetName <- rep('MAGMA',nrow(magmaResults))
  magmaResults$GeneSetCategoryName <- rep('genetics',nrow(magmaResults))
  magmaResults$GeneSetADLinked <- rep(TRUE,nrow(magmaResults))
  magmaResults$GeneSetBrainRegion <- rep(NA,nrow(magmaResults))
  magmaResults$GeneSetDirectionAD <- rep(NA,nrow(magmaResults))

  magmaResults <- dplyr::select(magmaResults,
                                ModuleNameFull,
                                GeneSetName,
                                GeneSetCategoryName,
                                GeneSetAssociationStatistic,
                                GeneSetEffect,
                                GeneSetBrainRegion,
                                GeneSetDirectionAD,
                                GeneSetADLinked)
  #####merge magma results with the module set given the schema
  moduleSummary <- dplyr::left_join(moduleSet,
                                    magmaResults)

  #add adjustd pvalue by brain region and gene set category
  moduleSummary <- AMPAD::splitByBrainRegionAdjustPvalue(moduleSummary)

  #####igap, dbgap, and genecards for ad genes enrichments
  #source('enrichmentAnalysis/run_amp_ad_enrichment.R')
  genesets1 <- synapseClient::synGet('syn5923958')
  load(synapseClient::getFileLocation(genesets1))
  dbgap <- AMPAD::pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=dbGaP")
  #kegg, wikipathways, biocarta, panther, jensen_diseases, omim_disease, omim_expanded
  kegg <- AMPAD::pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2016")
  wikipathways <- AMPAD::pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=WikiPathways_2016")
  jensen_diseases <- AMPAD::pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Jensen_DISEASES")
  biocarta <- AMPAD::pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=BioCarta_2016")
  panther <- AMPAD::pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Panther_2016")
  omim_disease <- AMPAD::pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=OMIM_Disease")
  omim_expanded <- AMPAD::pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=OMIM_Expanded")

  adList <- GeneSets$Alzheimers$`AD:GeneticLoci`
  adList <- c(adList,'HLA-DRB5','HLA-DRB1')
  adList <- adList[-which(adList=='HLA-DRB5-DRB1')]
  adList <- list(igap = adList)
  adList$dbgap <- dbgap$`Alzheimer Disease`
  adList$kegg <- kegg$`Alzheimer's disease_Homo sapiens_hsa05010`
  adList$wikipathwaysMouse <- wikipathways$`Alzheimers Disease_Mus musculus_WP2075`
  adList$wikipathwaysHuman <- wikipathways$`Alzheimers Disease_Homo sapiens_WP2059`
  adList$jensenDisease <- jensen_diseases$`Alzheimer's_disease`
  adList$biocarta <- biocarta$`Deregulation of CDK5 in Alzheimers Disease_Homo sapiens_h_p35alzheimersPathway`
  adList$pantherAmyloid <- panther$`Alzheimer disease-amyloid secretase pathway_Homo sapiens_P00003`
  adList$pantherPresenilin <- panther$`Alzheimer disease-presenilin pathway_Homo sapiens_P00004`
  adList$omim <- omim_disease$`alzheimer_disease`
  adList$omimExpanded <- omim_expanded$`alzheimer_disease`

  genecardsObj <- synapseClient::synGet('syn10507702')
  genecards <- data.table::fread(synapseClient::getFileLocation(genecardsObj),data.table=F)

  adList$genecards <- genecards$`Gene Symbol`



  adTest <- AMPAD::run_amp_ad_enrichment(adList,
                                  'genetics',
                                  manifestId=synId)

  adTest <- dplyr::select(adTest,
                          ModuleNameFull,
                          category,
                          geneSet,
                          fisherPval,
                          fisherOR)
  colnames(adTest) <- c('ModuleNameFull',
                        'GeneSetName',
                        'GeneSetCategoryName',
                        'GeneSetAssociationStatistic',
                        'GeneSetEffect')

  adTest$GeneSetADLinked <- rep(TRUE,nrow(adTest))
  adTest$GeneSetBrainRegion <- rep(NA,nrow(adTest))
  adTest$GeneSetDirectionAD <- rep(NA,nrow(adTest))
  adTestSummary <- dplyr::left_join(moduleSet,
                                    adTest)
  adTestSummary <- AMPAD::splitByBrainRegionAdjustPvalue(adTestSummary)
  #View(adTestSummary)



  moduleSummary <- rbind(moduleSummary,adTestSummary)
  moduleSummary$EvidenceClass <- 'genetics'

  return(moduleSummary)
}
