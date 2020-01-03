collateEnrichmentSets = function(){
  #synapseClient::synapseLogin()
  synapser::synLogin()
  enrichmentSets <- list()

  #1) pull ad genesets
  genesets1 <- synapser::synGet('syn5923958')
  load(genesets1$path)
  dbgap <- AMPAD::pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=dbGaP")
  #kegg, wikipathways, biocarta, panther, jensen_diseases, omim_disease, omim_expanded
  kegg <- AMPAD::pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2016")
  wikipathways <- AMPAD::pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=WikiPathways_2016")
  jensen_diseases <- AMPAD::pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Jensen_DISEASES")
  biocarta <- AMPAD::pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=BioCarta_2016")
  panther <- AMPAD::pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Panther_2016")
  omim_disease <- AMPAD::pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=OMIM_Disease")
  omim_expanded <- AMPAD::pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=OMIM_Expanded")
  go_cc <- AMPAD::pullReferenceGeneSets("http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Cellular_Component_2018")

  adList <- GeneSets$Alzheimers$`AD:GeneticLoci`
  adList <- c(adList,'HLA-DRB5','HLA-DRB1')
  adList <- adList[-which(adList=='HLA-DRB5-DRB1')]
  adList <- adList[-which(adList=='PSEN1' | adList =='PSEN2' | adList =='APP' | adList == 'P2K2B')]
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

  genecardsObj <- synapser::synGet('syn10507702')
  genecards <- data.table::fread(genecardsObj$path,data.table=F)

  adList$genecards <- genecards$`Gene Symbol`

  #load mckenzie oligo data
  mckenzieObj1 <- synapser::synGet('syn21482836')
  mckenzie <- data.table::fread(mckenzieObj1$path,data.table = F)
  adList$MSSM <- dplyr::filter(mckenzie,Module == 'Red' | Module == 'List green' | Module == 'Green')$Gene_Symbol
  enrichmentSets$MSSM2 <- AMPAD::listifyWrapper(mckenzie$Module,mckenzie$Gene_Symbol)


  sczNetworkObj <- synapser::synGet('syn7118804')
  sczNetwork <- data.table::fread(sczNetworkObj$path,data.table=F)
  enrichmentSets$scz <- AMPAD::listifyWrapper(sczNetwork$ModuleColor,sczNetwork$Ensembl)
  controlNetworkObj <- synapser::synGet('syn7118802')
  controlNetwork <- data.table::fread(controlNetworkObj$path,data.table=F)
  enrichmentSets$control <- AMPAD::listifyWrapper(controlNetwork$ModuleColor,controlNetwork$Ensembl)


  #load allen oligo data
  #allen <- data.table::fread('allen2018.csv',data.table=F)
  #adList$`AMP-AD Allen et al. 2018 AD+PSP TCX40.CS` <- dplyr::filter(allen,module=='TCX40')$gene_symbol

  #adList$`AMP-AD Allen et al. 2018 AD+PSP TCX10.CS` <- dplyr::filter(allen,module=='TCX10')$gene_symbol


  allen_simple_obj1 <- synapser::synGet('syn21482838')
  allen_simple <- data.table::fread(allen_simple_obj1$path,data.table=F)

  allen_comprehensive_obj1 <- synapser::synGet('syn21482837')
  allen_comprehensive <- data.table::fread(allen_comprehensive_obj1$path,data.table=F)

  allen_simple$module <- paste0(allen_simple$module,
                                '_simple')
  allen_comprehensive$module <- paste0(allen_comprehensive$module,
                                       '_comprehensive')
  allen_list_s <- sapply(unique(allen_simple$module),
                         AMPAD::listify,
                         allen_simple$gene_symbol,
                         allen_simple$module)
  names(allen_list_s) <- unique(allen_simple$module)
  adList$Mayo_simple <- unlist(allen_list_s)


  allen_list_c <- sapply(unique(allen_comprehensive$module),
                         AMPAD::listify,
                         allen_comprehensive$gene_symbol,
                         allen_comprehensive$module)

  names(allen_list_c) <- unique(allen_comprehensive$module)
  adList$Mayo_comprehensive <- unlist(allen_list_c)

  #load amp-ad targets
  adtargetsObj<-synapser::synGet('syn12540368',version = 17)

  adtargets <- data.table::fread(adtargetsObj$path,data.table=F)
  adList$`Nominated_targets` <- adtargets$hgnc_symbol

  #load module109
  adList$`Columbia_Broad_Rush_m109` <- unique(synapser::synTableQuery('select * from syn5321231 where speakeasyModule = 109')$asDataFrame()$hgncName)


  #load johnson 2018
  johnsonObj1 <- synapser::synGet('syn21482834')
  johnson <- data.table::fread(johnsonObj1$path,data.table=F)
  adList$Emory <- dplyr::filter(johnson,Module == WGCNA::labels2colors(17) | Module == WGCNA::labels2colors(18) | Module == WGCNA::labels2colors(10) | Module == WGCNA::labels2colors(15) | Module == WGCNA::labels2colors(29) | Module == WGCNA::labels2colors(40))$geneName
  enrichmentSets$ad <- lapply(adList,unique)

  #pull degs
  enrichmentSets$deg <- lapply(AMPAD::makeDEGAD(),unique)
  #pull cell types
  enrichmentSets$cell <- lapply(GeneSets$Cell_Markers,unique)

  df11Obj1 <- synapser::synGet('syn21482835')
  df11 <- data.table::fread(df11Obj1$path,data.table=F)
  enrichmentSets$cell2 <- AMPAD::listifyWrapper(df11$Cluster,df11$Gene)



  fob1 <- synapser::synGet('syn11914811')
  load(fob1$path)
  enrichmentSets$degMeta <- lapply(all.gs,unique)

  enrichmentSets$mito <- go_cc[grep('mitochondri',names(go_cc))]

  enrichmentSets$targetedPathways <- lapply(AMPAD::getCuratedPathways(),unique)
  return(enrichmentSets)
}
