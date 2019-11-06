pairwiseMatrixOfEnrichments = function(synId){
  synapser::synLogin()

  ##pull aggregate modules
  #aggregateModules <- rSynapseUtilities::loadFullTable(synId)
  aggregateModules <- synapser::synTableQuery(paste0("select * from ",synId))$asDataFrame()
  aggregateModules <- aggregateModules[,-c(1,2)]

  #View(aggregateModules)

  mats <- utilityFunctions::pairwiseMatrixOfEnrichments(key=aggregateModules$Module,
                                                        value=aggregateModules$GeneID)


  customDf <- data.frame(moduleName=c('TCXblue',
                                      'IFGyellow',
                                      'PHGyellow',
                                      'DLPFCblue',
                                      'CBEturquoise',
                                      'STGblue',
                                      'PHGturquoise',
                                      'IFGturquoise',
                                      'TCXturquoise',
                                      'FPturquoise',
                                      'IFGbrown',
                                      'STGbrown',
                                      'DLPFCyellow',
                                      'TCXgreen',
                                      'FPyellow',
                                      'CBEyellow',
                                      'PHGbrown',
                                      'DLPFCbrown',
                                      'STGyellow',
                                      'PHGgreen',
                                      'CBEbrown',
                                      'TCXyellow',
                                      'IFGblue',
                                      'FPblue',
                                      'FPbrown',
                                      'CBEblue',
                                      'DLPFCturquoise',
                                      'TCXbrown',
                                      'STGturquoise',
                                      'PHGblue'),
                         Cluster= c(rep('Consensus Cluster A',3),
                                       rep('Consensus Cluster B',7),
                                       rep('Consensus Cluster C',7),
                                       rep('Consensus Cluster D',7),
                                       rep('Consensus Cluster E',6)),
                         stringsAsFactors=F)
  rownames(customDf) <- customDf$moduleName
  customDf <- dplyr::select(customDf,Cluster)
  #get rid of underflow issues
  #red,yellow,green,blue,purple
  mats$pval <- mats$pval + 1e-300
  tiff(filename='figure2A.tiff',height=85,width=85,units='mm',res=300)
  ann_colors<-list(Cluster=c(`Consensus Cluster A`='#F8766D',
                 `Consensus Cluster B`='#A3A500',
                 `Consensus Cluster C`='#00BF7D',
                 `Consensus Cluster D`='#00B0F6',
                 `Consensus Cluster E`='#E76BF3'))
  pheatmap::pheatmap(-log10(mats$pval),
                     show_colnames = F,
                     border_color = NA,
                     fontsize=4,
                     treeheight_row=10,
                     treeheight_col=10,
                     annotation_col = customDf,
                     annotation_colors = ann_colors[1])
  dev.off()
}
