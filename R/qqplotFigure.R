qqplotFigure = function(){
  synapseClient::synapseLogin()
  adGeneticsSummaryAggOld <- AMPAD::getAdGenetics(synId='syn10915669')
  adGeneticsSummaryAgg <- AMPAD::getAdGenetics(synId='syn11926100')
  adGeneticsSummaryInd <- AMPAD::getAdGenetics(synId='syn10309369')
  adGeneticsSummaryAgg2 <- AMPAD::getAdGenetics2(synId='syn11926100')
  adGeneticsSummaryAgg3 <- AMPAD::getAdGenetics(synId = 'syn11870970')
  adGeneticsSummaryTest <- AMPAD::getAdGenetics(synId = 'syn11932957')

  adList<-adGeneticsSummaryAgg2[[2]]
  amp.ad.de.geneSets <- AMPAD::makeDEGAD()

  foo2 <- synapseClient::synTableQuery("select distinct external_gene_name from syn10309369")@values
  foo2ensg <- synapseClient::synTableQuery("select distinct GeneID from syn10309369")@values

  adListensg <- lapply(adList,utilityFunctions::convertHgncToEnsembl)
  adListensg2 <- lapply(adListensg,function(x) unique(x$ensembl_gene_id))

  foobar <- utilityFunctions::outerSapplyParallel( utilityFunctions::fisherWrapperPval, amp.ad.de.geneSets, adListensg2,foo2ensg$GeneID)

  foobar <- data.frame(foobar,stringsAsFactors=F)
  foobar$pathway <- rownames(foobar)
  foobar2<-tidyr::gather(foobar,key='geneset',value='pval',-pathway)

  foobar3 <- utilityFunctions::outerSapplyParallel( utilityFunctions::fisherWrapperOR, amp.ad.de.geneSets, adListensg2,foo2ensg$GeneID)
  foobar3 <- data.frame(foobar3,stringsAsFactors=F)
  foobar3$pathway <- rownames(foobar)
  foobar4<-tidyr::gather(foobar3,key='geneset',value='OR',-pathway)

  agg_eff <- adGeneticsSummaryAgg$GeneSetEffect
  ind_eff <- adGeneticsSummaryInd$GeneSetEffect
  de_eff <- foobar4$OR

  agg_eff[!is.finite(agg_eff)] <- NA
  ind_eff[!is.finite(ind_eff)] <- NA
  de_eff[!is.finite(de_eff)] <- NA
  agg_eff[agg_eff==0] <- NA
  ind_eff[ind_eff==0] <- NA
  de_eff[de_eff==0] <- NA

  gaiteri_mods <- read.csv('zhang_modules.csv',stringsAsFactors=F)
  modList <- lapply(unique(gaiteri_mods$Module), utilityFunctions::listify,gaiteri_mods$Gene_Symbol,gaiteri_mods$Module)
  names(modList) <- unique(gaiteri_mods$Module)
  foobarX <- utilityFunctions::outerSapplyParallel( utilityFunctions::fisherWrapperPval, modList, adList,foo2$external_gene_name)

  foobarX <- data.frame(foobarX,stringsAsFactors=F)
  foobarX$pathway <- rownames(foobarX)
  foobarX2<-tidyr::gather(foobarX,key='geneset',value='pval',-pathway)

  foobarX3 <- utilityFunctions::outerSapplyParallel( utilityFunctions::fisherWrapperOR, modList, adList,foo2$external_gene_name)
  foobarX3 <- data.frame(foobarX3,stringsAsFactors=F)
  foobarX3$pathway <- rownames(foobarX)
  foobarX4<-tidyr::gather(foobarX3,key='geneset',value='OR',-pathway)


  load(synapseClient::synGet('syn11914811')@filePath)
  all.gs <- all.gs[1:12]
  foobarY <- utilityFunctions::outerSapplyParallel( utilityFunctions::fisherWrapperPval, all.gs, adListensg2,foo2ensg$GeneID)

  foobarY <- data.frame(foobarY,stringsAsFactors=F)
  foobarY$pathway <- rownames(foobarY)
  foobarY2<-tidyr::gather(foobarY,key='geneset',value='pval',-pathway)

  foobarY3 <- utilityFunctions::outerSapplyParallel( utilityFunctions::fisherWrapperOR, all.gs, adListensg2,foo2ensg$GeneID)
  foobarY3 <- data.frame(foobarY3,stringsAsFactors=F)
  foobarY3$pathway <- rownames(foobarY)
  foobarY4<-tidyr::gather(foobarY3,key='geneset',value='OR',-pathway)


  load('cbe_res.rda')
  load('dlpfc_res.rda')
  load('phg_res.rda')
  load('fp_res.rda')
  load('tcx_res.rda')
  load('ifg_res.rda')
  load('stg_res.rda')
  uniqueMods <- c(CBEres$moduleGraph$from,CBEres$moduleGraph$to)%>% unique
  uniqueMods <- c(DLPFCres$moduleGraph$from,DLPFCres$moduleGraph$to) %>% unique %>% c(uniqueMods)
  uniqueMods <- c(PHGres$moduleGraph$from,PHGres$moduleGraph$to) %>% unique %>% c(uniqueMods)
  uniqueMods <- c(FPres$moduleGraph$from,FPres$moduleGraph$to) %>% unique %>% c(uniqueMods)
  uniqueMods <- c(TCXres$moduleGraph$from,TCXres$moduleGraph$to) %>% unique %>% c(uniqueMods)
  uniqueMods <- c(IFGres$moduleGraph$from,IFGres$moduleGraph$to) %>% unique %>% c(uniqueMods)
  uniqueMods <- c(STGres$moduleGraph$from,STGres$moduleGraph$to) %>% unique %>% c(uniqueMods)

  gap::qqunif(adGeneticsSummaryInd$GeneSetAssociationStatistic,xlim=c(0,5), ylim=c(0,42),pch=3)
  par(new=T)
  gap::qqunif(dplyr::filter(adGeneticsSummaryInd,ModuleNameFull %in% uniqueMods)$GeneSetAssociationStatistic,xlim=c(0,5),ylim=c(0,42),col='purple',pch=16)
  par(new=T)
  gap::qqunif(foobar2$pval,xlim=c(0,5),ylim=c(0,42),col='green',pch=17)
  par(new=T)
  gap::qqunif(foobarX2$pval,xlim=c(0,5),ylim=c(0,42),col='brown',pch=18)
  par(new=T)
  gap::qqunif(adGeneticsSummaryAgg$GeneSetAssociationStatistic,xlim=c(0,5),ylim=c(0,42),col='red',main='AD Gene-set Enrichments',pch=15)
  par(new=T)
  gap::qqunif(adGeneticsSummaryTest$GeneSetAssociationStatistic,xlim=c(0,5),ylim=c(0,42),col='darkblue')
  par(new=T)
  gap::qqunif(foobarY2$pval,xlim=c(0,5),ylim=c(0,42),col='red')
  #par(new=T)
  #gap::qqunif(adGeneticsSummaryAggOld$GeneSetAssociationStatistic,xlim=c(0,5),ylim=c(0,42),col='purple')
  legend('topleft',c('all modules','all DEG enriched modules','all DEGs','Zhang et al. 2013 modules','aggregate modules'),pch=c(3,16:18,15),col=c('blue','purple','green','brown','red'))


  #####make a series of box plots of enrichment odds ratios, also of the percentage of significant associations
  mean((adGeneticsSummaryInd$GeneSetAssociationStatistic)*nrow(adGeneticsSummaryInd)<0.05)

  ga2<-dplyr::filter(adGeneticsSummaryInd,ModuleNameFull %in% uniqueMods)
  mean(ga2$GeneSetAssociationStatistic*nrow(ga2)<0.05)
  mean(foobar2$pval*nrow(foobar2)<0.05)
  mean(foobarX2$pval*nrow(foobarX2)<0.05)
  mean((adGeneticsSummaryAgg$GeneSetAssociationStatistic)*nrow(adGeneticsSummaryAgg)<0.05)
  mean((adGeneticsSummaryAggOld$GeneSetAssociationStatistic)*nrow(adGeneticsSummaryAggOld)<0.05)
  mean((adGeneticsSummaryAgg3$GeneSetAssociationStatistic)*nrow(adGeneticsSummaryAgg3)<0.05)
  mean((adGeneticsSummaryTest$GeneSetAssociationStatistic)*nrow(adGeneticsSummaryTest)<0.05)
  mean(foobarY2$pval*nrow(foobarY2) < 0.05)
  #mean((adGeneticsSummaryInd))

  vec1 <- (adGeneticsSummaryTest$GeneSetAssociationStatistic)*nrow(adGeneticsSummaryTest)<0.05
  vec2 <- foobarY2$pval*nrow(foobarY2) < 0.05
  vec3 <- (adGeneticsSummaryInd$GeneSetAssociationStatistic)*nrow(adGeneticsSummaryInd)<0.05
  vec4 <- foobar2$pval*nrow(foobar2)<0.05
  vec5 <- foobarX2$pval*nrow(foobarX2)<0.05


  png(file='agg_mod_figure.png',width=1600,height=1200,pointsize=40)
  gap::qqunif(adGeneticsSummaryInd$GeneSetAssociationStatistic,xlim=c(0,5), ylim=c(0,42),pch=3)
  par(new=T)
  gap::qqunif(dplyr::filter(adGeneticsSummaryInd,ModuleNameFull %in% uniqueMods)$GeneSetAssociationStatistic,xlim=c(0,5),ylim=c(0,42),col='purple',pch=16)
  par(new=T)
  gap::qqunif(foobar2$pval,xlim=c(0,5),ylim=c(0,42),col='green',pch=17)
  par(new=T)
  gap::qqunif(foobarX2$pval,xlim=c(0,5),ylim=c(0,42),col='brown',pch=18)
  par(new=T)
  gap::qqunif(adGeneticsSummaryAgg$GeneSetAssociationStatistic,xlim=c(0,5),ylim=c(0,42),col='red',main='AD Gene-set Enrichments',pch=15)
  #par(new=T)
  #gap::qqunif(adGeneticsSummaryAggOld$GeneSetAssociationStatistic,xlim=c(0,5),ylim=c(0,42),col='purple')
  legend('topleft',c('all modules','all DEG enriched modules','all DEGs','Zhang et al. 2013 modules','aggregate modules'),pch=c(3,16:18,15),col=c('blue','purple','green','brown','red'))
  dev.off()
}
