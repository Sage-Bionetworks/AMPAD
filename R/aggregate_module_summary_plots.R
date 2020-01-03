aggregate_module_summary_plots = function(outputFile=FALSE){
  foo <- synapser::synTableQuery("select * from syn11932957")$asDataFrame()
  foo2 <- dplyr::select(foo,GeneID,Module)
  foo2$Presence <- 1
  foo3 <- tidyr::pivot_wider(foo2,
                             id_cols = "GeneID",
                             names_from = "Module",
                             values_from = "Presence")
  foo3[is.na(foo3)] <- 0
  foo3 <- data.frame(foo3,stringsAsFactors=F)
  foo4 <- dplyr::select(foo3,GeneID,TCXblue,IFGyellow,PHGyellow)
  resu <- list()

  if(outputFile){
    tiff(filename = 'consensusClusterA.tiff', height = 4, width = 6,units='in',pointsize=14,res=300)
    UpSetR::upset(foo4,nintersects = NA,show.numbers=F)
    dev.off()
  }else{
    resu$A<-UpSetR::upset(foo4,nintersects = NA,show.numbers=F)
  }

  nUniqueGenesA <- data.frame(Module=c('TCXblue','PHGyellow','IFGyellow'),nGenes=c(979,366,127),stringsAsFactors=F)

  foo4 <- dplyr::select(foo3,GeneID,DLPFCblue,CBEturquoise,STGblue,PHGturquoise,IFGturquoise,TCXturquoise,FPturquoise)


  if(outputFile){
    tiff(filename = 'consensusClusterB.tiff', height = 4, width = 6,units='in',pointsize=14,res=300)
    UpSetR::upset(foo4,nsets=7,nintersects = NA,point.size=1,show.numbers = F)
    dev.off()
  } else{
    resu$B<-UpSetR::upset(foo4,nsets=7,nintersects = NA,point.size=1,show.numbers = F)
  }

  nUniqueGenesB <- data.frame(Module=c('CBEturquoise','DLPFCblue','IFGturquoise','PHGturquoise','STGblue','TCXturquoise','FPturquoise'),nGenes=c(593,349,275,209,163,69,40),stringsAsFactors=F)

  foo4 <- dplyr::select(foo3,GeneID,IFGbrown,STGbrown,DLPFCyellow,TCXgreen,FPyellow,CBEyellow,PHGbrown)

  if(outputFile){
    tiff(filename = 'consensusClusterC.tiff', height = 4, width = 6,units='in',pointsize=14,res=300)
    UpSetR::upset(foo4,nsets=7,nintersects = NA,point.size=1,show.numbers = F)
    dev.off()
  }else{
    resu$C<-UpSetR::upset(foo4,nsets=7,nintersects = NA,point.size=1,show.numbers = F)
  }

  nUniqueGenesC <- data.frame(Module=c('IFGbrown','FPyellow','STGbrown','DLPFCyellow','TCXgreen','PHGbrown','CBEyellow'),nGenes=c(966,641,233,178,141,139,28),stringsAsFactors=F)

  foo4 <- dplyr::select(foo3,GeneID,DLPFCbrown,STGyellow,PHGgreen,CBEbrown,TCXyellow,IFGblue,FPblue)



  if(outputFile){
    tiff(filename = 'consensusClusterD.tiff', height = 4, width = 6,units='in',pointsize=14,res=300)
    UpSetR::upset(foo4,nsets=7,nintersects = NA,point.size=1,show.numbers = F)
    dev.off()
  }else{
    resu$D<-UpSetR::upset(foo4,nsets=7,nintersects = NA,point.size=1,show.numbers = F)
  }

  nUniqueGenesD <- data.frame(Module=c('IFGblue','TCXyellow','FPblue','STGyellow','PHGgreen','DLPFCbrown','CBEbrown'),nGenes=c(1148,673,627,344,122,103,56),stringsAsFactors=F)


  foo4 <- dplyr::select(foo3,GeneID,FPbrown,CBEblue,DLPFCturquoise,TCXbrown,STGturquoise,PHGblue)


  if(outputFile){
  tiff(filename = 'consensusClusterE.tiff', height = 4, width = 6,units='in',pointsize=14,res=300)
  UpSetR::upset(foo4,nsets=6,nintersects = NA,point.size=1,show.numbers = F)
  dev.off()
  } else{
    resu$E<-UpSetR::upset(foo4,nsets=6,nintersects = NA,point.size=1,show.numbers = F)
  }

  nUniqueGenesE <- data.frame(Module=c('CBEblue','PHGblue','DLPFCturquoise','STGturquoise','TCXbrown','FPbrown'),nGenes=c(1862,951,447,423,358,201),stringsAsFactors=F)


  nUniqueGenes <- rbind(nUniqueGenesA,
                        nUniqueGenesB,
                        nUniqueGenesC,
                        nUniqueGenesD,
                        nUniqueGenesE)

  library(dplyr)
  modSize <- dplyr::group_by(foo2,Module) %>%
    dplyr::summarise(mSize=sum(Presence))

  sumMat1 <- dplyr::left_join(modSize,nUniqueGenes)
  sumMat1$percentUnique <- sumMat1$nGenes/sumMat1$mSize

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


  sumMat1 <- dplyr::left_join(sumMat1,customDf,by=c('Module'='moduleName'))

  # cat('% overlap for Consensus Clusters A-C')
  # print(summary(lm(percentUnique ~ 1,dplyr::filter(sumMat1,Cluster=="Consensus Cluster A" | Cluster=="Consensus Cluster B" | Cluster=="Consensus Cluster C"))))
  #
  # cat('% overlap for Consensus Clusters D & E\n')
  # print(summary(lm(percentUnique ~ 1,dplyr::filter(sumMat1,Cluster=="Consensus Cluster D" | Cluster=="Consensus Cluster E"))))

  return(resu)

}
