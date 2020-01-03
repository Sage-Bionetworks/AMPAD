improvedAdRelevancePlot <- function(outputFile=FALSE){

  adGeneticsSummaryAggOld <- AMPAD::getAdGenetics(synId='syn10915669')
  adGeneticsSummaryAgg <- AMPAD::getAdGenetics(synId='syn11926100')
  adGeneticsSummaryInd <- AMPAD::getAdGenetics(synId='syn10309369')
  adGeneticsSummaryAgg2 <- AMPAD::getAdGenetics2(synId='syn11926100')
  adGeneticsSummaryAgg3 <- AMPAD::getAdGenetics(synId = 'syn11870970')
  adGeneticsSummaryTest <- AMPAD::getAdGenetics(synId = 'syn11932957')

  adList<-adGeneticsSummaryAgg2[[2]]
  amp.ad.de.geneSets <- AMPAD::makeDEGAD()

  foo2 <- synapser::synTableQuery("select distinct external_gene_name from syn10309369")$asDataFrame()
  foo2ensg <- synapser::synTableQuery("select distinct GeneID from syn10309369")$asDataFrame()

  adListensg <- lapply(adList,AMPAD::convertHgncToEnsembl)
  adListensg2 <- lapply(adListensg,function(x) unique(x$ensembl_gene_id))

  foobar <- AMPAD::outerSapplyParallel( AMPAD::fisherWrapperPval, amp.ad.de.geneSets, adListensg2,foo2ensg$GeneID)

  foobar <- data.frame(foobar,stringsAsFactors=F)
  foobar$pathway <- rownames(foobar)
  foobar2<-tidyr::gather(foobar,key='geneset',value='pval',-pathway)

  foobar3 <- AMPAD::outerSapplyParallel( AMPAD::fisherWrapperOR, amp.ad.de.geneSets, adListensg2,foo2ensg$GeneID)
  foobar3 <- data.frame(foobar3,stringsAsFactors=F)
  foobar3$pathway <- rownames(foobar)
  foobar4<-tidyr::gather(foobar3,key='geneset',value='OR',-pathway)


  mckenzieObj1 <- synapser::synGet('syn21482836')
  gaiteri_mods <- data.table::fread(mckenzieObj1$path,data.table = F)
  #gaiteri_mods <- read.csv('zhang_modules.csv',stringsAsFactors=F)
  modList <- lapply(unique(gaiteri_mods$Module), AMPAD::listify,gaiteri_mods$Gene_Symbol,gaiteri_mods$Module)
  names(modList) <- unique(gaiteri_mods$Module)
  foobarX <- AMPAD::outerSapplyParallel( AMPAD::fisherWrapperPval, modList, adList,foo2$external_gene_name)

  foobarX <- data.frame(foobarX,stringsAsFactors=F)
  foobarX$pathway <- rownames(foobarX)
  foobarX2<-tidyr::gather(foobarX,key='geneset',value='pval',-pathway)

  foobarX3 <- AMPAD::outerSapplyParallel( AMPAD::fisherWrapperOR, modList, adList,foo2$external_gene_name)
  foobarX3 <- data.frame(foobarX3,stringsAsFactors=F)
  foobarX3$pathway <- rownames(foobarX)
  foobarX4<-tidyr::gather(foobarX3,key='geneset',value='OR',-pathway)


  load(synapser::synGet('syn11914811')$path)
  all.gs <- all.gs[1:12]
  foobarY <- AMPAD::outerSapplyParallel( AMPAD::fisherWrapperPval, all.gs, adListensg2,foo2ensg$GeneID)

  foobarY <- data.frame(foobarY,stringsAsFactors=F)
  foobarY$pathway <- rownames(foobarY)
  foobarY2<-tidyr::gather(foobarY,key='geneset',value='pval',-pathway)

  foobarY3 <- AMPAD::outerSapplyParallel( AMPAD::fisherWrapperOR, all.gs, adListensg2,foo2ensg$GeneID)
  foobarY3 <- data.frame(foobarY3,stringsAsFactors=F)
  foobarY3$pathway <- rownames(foobarY)
  foobarY4<-tidyr::gather(foobarY3,key='geneset',value='OR',-pathway)

  metaAnalysisEnrichment <- dplyr::left_join(foobarY2,foobarY4)
  gaiteriEnrichment <- dplyr::left_join(foobarX2,foobarX4)
  degAnalysisEnrichment <- dplyr::left_join(foobar2,foobar4)

  #load('aggregateModules.rda')
  aggmodobj <- synapser::synGet('syn21483261')
  load(aggmodobj$path)
  combinedMatrix <- rbind(dlpfc_mods$moduleGraph,
                          cbe_mods$moduleGraph,
                          tcx_mods$moduleGraph,
                          fp_mods$moduleGraph,
                          stg_mods$moduleGraph,
                          ifg_mods$moduleGraph,
                          phg_mods$moduleGraph)
  uniqueMods <- unique(unique(combinedMatrix$from),
                       unique(combinedMatrix$to))

  adGeneticsSummaryTestDEG <- dplyr::filter(adGeneticsSummaryInd,
                                            ModuleNameFull %in% uniqueMods)


  #combine into one massive data frame that can be reorganized as necessary for effective plotttng or grouping.
  # adGeneticsSummaryAggOld <- AMPAD::getAdGenetics(synId='syn10915669')
  # adGeneticsSummaryAgg <- AMPAD::getAdGenetics(synId='syn11926100')
  # adGeneticsSummaryInd <- AMPAD::getAdGenetics(synId='syn10309369')
  # adGeneticsSummaryAgg2 <- AMPAD::getAdGenetics2(synId='syn11926100')
  # adGeneticsSummaryAgg3 <- AMPAD::getAdGenetics(synId = 'syn11870970')
  # adGeneticsSummaryTest <- AMPAD::getAdGenetics(synId = 'syn11932957')


  AggOldAD <- dplyr::select(adGeneticsSummaryAggOld,
                            ModuleNameFull,
                            GeneSetName,
                            GeneSetAssociationStatistic,
                            GeneSetEffect) %>%
    dplyr::mutate(category = rep('Original Modules',nrow(adGeneticsSummaryAggOld)))


  AggModsAD <- dplyr::select(adGeneticsSummaryAgg,
                                 ModuleNameFull,
                                 GeneSetName,
                                 GeneSetAssociationStatistic,
                                 GeneSetEffect) %>%
    dplyr::mutate(category = rep('Aggregate Brain Specific Modules',nrow(adGeneticsSummaryAgg)))

  IndModsAD <- dplyr::select(adGeneticsSummaryInd,
                             ModuleNameFull,
                             GeneSetName,
                             GeneSetAssociationStatistic,
                             GeneSetEffect) %>%
    dplyr::mutate(category = rep('Individual Modules',nrow(adGeneticsSummaryInd)))

  FinalMods <- dplyr::select(adGeneticsSummaryTest,
                             ModuleNameFull,
                             GeneSetName,
                             GeneSetAssociationStatistic,
                             GeneSetEffect) %>%
    dplyr::mutate(category = rep('Final AD Modules',nrow(adGeneticsSummaryTest)))


  # metaAnalysisEnrichment <- dplyr::left_join(foobarY2,foobarY4)
  # gaiteriEnrichment <- dplyr::left_join(foobarX2,foobarX4)
  # degAnalysisEnrichment <- dplyr::left_join(foobar2,foobar4)
  # adGeneticsSummaryTestDEG <- dplyr::filter(adGeneticsSummaryInd,
  #                                           ModuleNameFull %in% uniqueMods)

  DEGmeta <- dplyr::mutate(metaAnalysisEnrichment,category = rep('DEG Meta Analysis',nrow(metaAnalysisEnrichment)))

  cellPaper <- dplyr::mutate(gaiteriEnrichment,category = rep('Zhang et al. 2013',nrow(gaiteriEnrichment)))

  DEGbrain <- dplyr::mutate(degAnalysisEnrichment,category = rep('DEG Brain Region Specific',nrow(degAnalysisEnrichment)))

  ModsDEGEnriched <- dplyr::select(adGeneticsSummaryTestDEG,
                                      ModuleNameFull,
                                      GeneSetName,
                                      GeneSetAssociationStatistic,
                                      GeneSetEffect) %>%
    dplyr::mutate(category = rep('Individual DEG Modules',nrow(adGeneticsSummaryTestDEG)))
  colnames(DEGmeta) <- colnames(AggOldAD)
  colnames(cellPaper) <- colnames(AggOldAD)
  colnames(DEGbrain) <- colnames(AggOldAD)
  #colnames(ModsDEGEnriched) <- colnames(AggOldAD)

  fullMatrix <- rbind(AggOldAD,
                      AggModsAD,
                      IndModsAD,
                      FinalMods,
                      DEGmeta,
                      cellPaper,
                      DEGbrain,
                      ModsDEGEnriched)

  #### create a new significance column
  sigFun <- function(categoryName,fullMatrix){
    fob <- dplyr::filter(fullMatrix,category == categoryName)
    fob$adj.pval <- p.adjust(fob$GeneSetAssociationStatistic,
                             method='fdr')
    fob$significant <- fob$adj.pval <= 0.05
    return(fob)
  }

  fullMatrix2 <- lapply(unique(fullMatrix$category),
                        sigFun,
                        fullMatrix)
  fullMatrix2 <- do.call(rbind,fullMatrix2)

  summFullMatrix <- dplyr::group_by(fullMatrix2,category)
  summMatrix <- dplyr::summarise(summFullMatrix,
                                 percentSig = mean(significant),
                                 se = sd(significant)/sqrt(length(significant)))
  summMatrix <- dplyr::arrange(summMatrix,(percentSig))
  #summMatrix

  #summMatrix <- summMatrix[c(6,5,3,8,2,4,1),]
  #summMatrix <- summMatrix[7:1,]
  summMatrix <- summMatrix[-5,]
  summMatrix$category <- factor(summMatrix$category,levels = summMatrix$category)

  g <- ggplot2::ggplot(summMatrix,
                       ggplot2::aes(x=category, y=percentSig))

  g <- g + ggplot2::geom_bar(position=ggplot2::position_dodge(),
                             stat="identity")

  g <- g + ggplot2::geom_errorbar(ggplot2::aes(ymin=percentSig-se, ymax=percentSig+se),
                  width=.2,                    # Width of the error bars
                  position=ggplot2::position_dodge(.9))
  g <- g + cowplot::theme_cowplot(12)
  g <- g + ggplot2::labs(x = 'System Biology Derived AD Geneset',
                         y = 'Percent pairwise associations significant')
  #g <- g + ggplot2::theme_update(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
  #                               panel.background = ggplot2::element_rect(fill = "white", colour = "white"),
  #                               panel.border = ggplot2::element_rect(colour = "white"))


  g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))


  #
  #g
  #g2 <- cowplot::plot_grid(g)
  if(outputFile){
    g
    ggplot2::ggsave('figure1.tiff',device='tiff',units='mm',width=114,height=85,scale=2)
  }
  return(g)

}
