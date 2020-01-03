make_dot_plot <- function(dummyDf,xlab,ylab,cowsize=11,outputFile=FALSE,fileName='test.tiff'){
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
  dummyDf <- dplyr::left_join(dummyDf,customDf,by=c('Module'='moduleName'))
  dummyDf$Cluster <- factor(dummyDf$Cluster,levels = (c('Consensus Cluster A',
                                                        'Consensus Cluster B',
                                                        'Consensus Cluster C',
                                                        'Consensus Cluster D',
                                                        'Consensus Cluster E')))
  dummyDf$Module <- factor(dummyDf$Module,levels = (c('TCXblue',
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
                                                      'PHGblue')))



  g <- ggplot2::ggplot(dummyDf,
                       ggplot2::aes(x = Module,
                                    y = category,
                                    size = fisherOR,
                                    color = Cluster))
  g <- g + ggplot2::geom_count()

  g <- g + ggplot2::labs(y = ylab,
                         x = xlab)

  g <- g + AMPAD::cowplot_rotated(cowsize)

  if(outputFile){
    ggplot2::ggsave(fileName,device='tiff',units='mm',width=85,height=85,scale=1.8)
  }
  return(g)
}
