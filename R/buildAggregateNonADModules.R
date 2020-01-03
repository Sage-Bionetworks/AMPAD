buildAggregateNonADModules <- function(){
  #syn11914811
  #replace this with the appropriate synapse call soon
  fob1 <- synapser::synGet('syn11914811')
  load(fob1$path)
  ad_names <- grep('^ad',names(all.gs))
  all.gs <- all.gs[ad_names]




  #replace with indivudal module manifest
  #indMods <- fullManifestOfDAta$IndividualModules
  indMods <- synapser::synTableQuery("select * from syn10309369")$asDataFrame()
  indMods <- indMods[,-c(1,2)]
  #re-run enrichment analysis of full module set for these lists

  indModsList <- lapply(unique(indMods$ModuleNameFull),
                        AMPAD::listify,
                        indMods$GeneID,
                        indMods$ModuleNameFull)
  names(indModsList) <- unique(indMods$ModuleNameFull)
  allgenes <- unique(indMods$GeneID)

  enrichmentMatrix <- AMPAD::outerSapply( AMPAD::fisherWrapperPval,
                                                     indModsList,
                                                     all.gs,
                                                     allgenes)
  #tidyr it
  enrichmentMatrix2 <- data.frame(enrichmentMatrix,stringsAsFactors = FALSE)
  enrichmentMatrix2$geneSet <- rownames(enrichmentMatrix2)
  enrichmentMatrix2 <- tidyr::gather(enrichmentMatrix2,key='ModuleNameFull',value='Pvalue',-geneSet)

  enrichmentMatrix2 <- dplyr::mutate(enrichmentMatrix2,adj.pval=p.adjust(Pvalue,method='bonferroni'))
  #enrichmentMatrix3 <- dplyr::filter(enrichmentMatrix2,adj.pval<=0.05)
  indModsSumm <- dplyr::group_by(indMods,ModuleNameFull,brainRegion,method)
  indModsAgg <- dplyr::summarise(indModsSumm,sizeOfMod=length(GeneID))
  enrichmentMatrix4 <- dplyr::left_join(enrichmentMatrix2,indModsAgg)
  enrSum <- dplyr::group_by(enrichmentMatrix4,brainRegion)
  enrAgg <- dplyr::summarise(enrSum,uniqueMods = length(unique(ModuleNameFull)))
  indSum <- dplyr::group_by(indModsAgg,brainRegion)
  indAgg <- dplyr::summarise(indSum,uniqueMods = length(unique(ModuleNameFull)))

  #for each brain region build module comparison graph
  #replace, using dlpfc as an example now...
  enrMatDLPFC <- dplyr::filter(enrichmentMatrix4,brainRegion == 'DLPFC')
  pairMods <- synapser::synTableQuery("SELECT * FROM syn10339153 where ModuleNameFull like '%DLPFC' and category like '%DLPFC'")$asDataFrame()
  pairMods <- pairMods[,-c(1,2)]
  library(dplyr)
  pairMods <- pairMods %>% AMPAD::removeSwappedDupKeyValueDf() %>%
    dplyr::mutate(adj=p.adjust(fisherPval,method='bonferroni')) %>%
    dplyr::filter(adj<=0.05) %>%
    dplyr::filter(from%in%enrMatDLPFC$ModuleNameFull & to %in% enrMatDLPFC$ModuleNameFull) %>%
    dplyr::mutate(weight = 1/fisherOR)

  dlpfc_mods <- AMPAD::buildModuleGraph(pairMods,indMods,'DLPFC')

  ######TCX
  enrMatTCX <- dplyr::filter(enrichmentMatrix4,brainRegion == 'TCX')
  pairMods <- synapser::synTableQuery("SELECT * FROM syn10339153 where ModuleNameFull like '%TCX' and category like '%TCX'")$asDataFrame()
  pairMods <- pairMods[,-c(1,2)]
  library(dplyr)
  pairMods <- pairMods %>% AMPAD::removeSwappedDupKeyValueDf() %>%
    dplyr::mutate(adj=p.adjust(fisherPval,method='bonferroni')) %>%
    dplyr::filter(adj<=0.05) %>%
    dplyr::filter(from%in%enrMatTCX$ModuleNameFull & to %in% enrMatTCX$ModuleNameFull) %>%
    dplyr::mutate(weight = 1/fisherOR)

  tcx_mods <- AMPAD::buildModuleGraph(pairMods,indMods,'TCX')

  ###CBE
  enrMatCBE <- dplyr::filter(enrichmentMatrix4,brainRegion == 'CBE')
  pairMods <- synapser::synTableQuery("SELECT * FROM syn10339153 where ModuleNameFull like '%CBE' and category like '%CBE'")$asDataFrame()
  pairMods <- pairMods[,-c(1,2)]
  pairMods <- pairMods %>% AMPAD::removeSwappedDupKeyValueDf() %>%
    dplyr::mutate(adj=p.adjust(fisherPval,method='bonferroni')) %>%
    dplyr::filter(adj<=0.05) %>%
    dplyr::filter(from%in%enrMatCBE$ModuleNameFull & to %in% enrMatCBE$ModuleNameFull) %>%
    dplyr::mutate(weight = 1/fisherOR)

  cbe_mods <- AMPAD::buildModuleGraph(pairMods,indMods,'CBE')

  ####FP
  enrMatFP <- dplyr::filter(enrichmentMatrix4,brainRegion == 'FP')

  pairMods <- synapser::synTableQuery("SELECT * FROM syn10339153 where ModuleNameFull like '%FP' and category like '%FP'")$asDataFrame()
  pairMods <- pairMods[,-c(1,2)]
  #pairMods <- read.csv('fp_pairwise.csv',stringsAsFactors = F)
  #pairMods <- pairMods[,-c(1,2)]
  pairMods <- pairMods %>% AMPAD::removeSwappedDupKeyValueDf() %>%
    dplyr::mutate(adj=p.adjust(fisherPval,method='bonferroni')) %>%
    dplyr::filter(adj<=0.05) %>%
    dplyr::filter(from%in%enrMatFP$ModuleNameFull & to %in% enrMatFP$ModuleNameFull) %>%
    dplyr::mutate(weight = 1/fisherOR)

  fp_mods <- AMPAD::buildModuleGraph(pairMods,indMods,'FP')

  ####STG
  enrMatSTG <- dplyr::filter(enrichmentMatrix4,brainRegion == 'STG')

  pairMods <- synapser::synTableQuery("SELECT * FROM syn10339153 where ModuleNameFull like '%STG' and category like '%STG'")$asDataFrame()
  pairMods <- pairMods[,-c(1,2)]
  #pairMods <- read.csv('fp_pairwise.csv',stringsAsFactors = F)
  #pairMods <- pairMods[,-c(1,2)]
  pairMods <- pairMods %>% AMPAD::removeSwappedDupKeyValueDf() %>%
    dplyr::mutate(adj=p.adjust(fisherPval,method='bonferroni')) %>%
    dplyr::filter(adj<=0.05) %>%
    dplyr::filter(from%in%enrMatSTG$ModuleNameFull & to %in% enrMatSTG$ModuleNameFull) %>%
    dplyr::mutate(weight = 1/fisherOR)

  stg_mods <- AMPAD::buildModuleGraph(pairMods,indMods,'STG')

  ###IFG
  enrMatIFG <- dplyr::filter(enrichmentMatrix4,brainRegion == 'IFG')

  pairMods <- synapser::synTableQuery("SELECT * FROM syn10339153 where ModuleNameFull like '%IFG' and category like '%IFG'")$asDataFrame()
  pairMods <- pairMods[,-c(1,2)]
  #pairMods <- read.csv('fp_pairwise.csv',stringsAsFactors = F)
  #pairMods <- pairMods[,-c(1,2)]
  pairMods <- pairMods %>% AMPAD::removeSwappedDupKeyValueDf() %>%
    dplyr::mutate(adj=p.adjust(fisherPval,method='bonferroni')) %>%
    dplyr::filter(adj<=0.05) %>%
    dplyr::filter(from%in%enrMatIFG$ModuleNameFull & to %in% enrMatIFG$ModuleNameFull) %>%
    dplyr::mutate(weight = 1/fisherOR)

  ifg_mods <- AMPAD::buildModuleGraph(pairMods,indMods,'IFG')

  ####PHG
  enrMatPHG <- dplyr::filter(enrichmentMatrix4,brainRegion == 'PHG')

  pairMods <- synapser::synTableQuery("SELECT * FROM syn10339153 where ModuleNameFull like '%PHG' and category like '%PHG'")$asDataFrame()
  pairMods <- pairMods[,-c(1,2)]
  #pairMods <- read.csv('fp_pairwise.csv',stringsAsFactors = F)
  #pairMods <- pairMods[,-c(1,2)]
  pairMods <- pairMods %>% AMPAD::removeSwappedDupKeyValueDf() %>%
    dplyr::mutate(adj=p.adjust(fisherPval,method='bonferroni')) %>%
    dplyr::filter(adj<=0.05) %>%
    dplyr::filter(from%in%enrMatPHG$ModuleNameFull & to %in% enrMatPHG$ModuleNameFull) %>%
    dplyr::mutate(weight = 1/fisherOR)

  phg_mods <- AMPAD::buildModuleGraph(pairMods,indMods,'PHG')
  AggregateModuleManifest <- rbind(dlpfc_mods$df,
                                   cbe_mods$df,
                                   tcx_mods$df,
                                   ifg_mods$df,
                                   stg_mods$df,
                                   phg_mods$df,
                                   fp_mods$df)
  return(AggregateModuleManifest)
  #save(dlpfc_mods,cbe_mods,tcx_mods,ifg_mods,stg_mods,phg_mods,fp_mods,file='aggregateNonAdModules.rda')
  #rSynapseUtilities::makeTable(AggregateModuleManifest,'AMP-AD aggregate modules',projectId='syn2370594')
}
