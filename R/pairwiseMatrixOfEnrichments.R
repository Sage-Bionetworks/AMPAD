pairwiseMatrixOfEnrichments = function(synId){
  synapseClient::synapseLogin()

  ##pull aggregate modules
  aggregateModules <- rSynapseUtilities::loadFullTable(synId)

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
                         Cluster= c(rep('Astrocytic',3),
                                       rep('Microglial',7),
                                       rep('Neuronal',7),
                                       rep('Oligodendroglial',7),
                                       rep('Proteostasis',6)),
                         stringsAsFactors=F)
  rownames(customDf) <- customDf$moduleName
  customDf <- dplyr::select(customDf,Cluster)
  #get rid of underflow issues
  #red,yellow,green,blue,purple
  mats$pval <- mats$pval + 1e-300
  tiff(filename='figure2A.tiff',height=85,width=85,units='mm',res=300)
  ann_colors<-list(Cluster=c(Astrocytic='#F8766D',
                 Microglial='#A3A500',
                 Neuronal='#00BF7D',
                 Oligodendroglial='#00B0F6',
                 Proteostasis='#E76BF3'))
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
