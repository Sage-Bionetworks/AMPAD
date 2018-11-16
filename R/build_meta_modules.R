build_meta_modules <- function(){
  #get aggregate modules
  synapser::synLogin()
  foo <- synapser::synTableQuery("select * from syn11932957")
  foo <- foo$asDataFrame()
  #define mapping into meta modules
  bar <- list(astrocyteModule = c('TCXblue',
                                  'IFGyellow',
                                  'PHGyellow'),
              microgliaModule = c('DLPFCblue',
                                  'CBEturquoise',
                                  'STGblue',
                                  'PHGturquoise',
                                  'IFGturquoise',
                                  'TCXturquoise',
                                  'FPturquoise'),
              neuronalModule = c('IFGbrown',
                                 'STGbrown',
                                 'DLPFCyellow',
                                 'TCXgreen',
                                 'FPyellow',
                                 'CBEyellow',
                                 'PHGbrown'),
              oligodendrocyteGlialModule = c('DLPFCbrown',
                                             'STGyellow',
                                             'PHGgreen',
                                             'CBEbrown',
                                             'TCXyellow',
                                             'IFGblue',
                                             'FPblue'),
              proteostasisModule = c('FPbrown',
                                     'CBEblue',
                                     'DLPFCturquoise',
                                     'TCXbrown',
                                     'STGturquoise',
                                     'PHGblue'))

  #create meta modules
  fxn1 <- function(modules,metamoduleName,aggmoddf){
    foo1 <- dplyr::filter(aggmoddf,Module %in% modules)
    foo2 <- foo1[!duplicated(foo1$GeneID),]
    foo2 <- foo2[,-c(1,2,4,5,6,7,8)]
    foo2$Module <- metamoduleName
    foo2$method <- 'meta'
    foo2 <- foo2$external_gene_name
    return(foo2)
  }
  res <- mapply(fxn1,bar,names(bar),MoreArgs = list(aggmoddf=foo),SIMPLIFY=F)

  #push to synapse
}

fooB <- synapser::synGet('syn12540368')
barB <- data.table::fread(fooB$path,data.table=F)

targetList <- list(AMPADTargets = barB$hgnc_symbol,
                   DrugTargets = c('ACHE',
                                   'APP',
                                   'PSEN1',
                                   'PSEN2',
                                   'MAPT',
                                   'RXRA',
                                   'GRIN1',
                                   'PLA2G7',
                                   'NGF',
                                   'BACE1',
                                   'BACE2',
                                   'GABRA1'))

#run enrichment
resA <- list()
resA$fisher <- utilityFunctions::outerSapplyParallel(utilityFunctions::fisherWrapper,
                                                    res,
                                                    targetList,
                                                    unique(unlist(res)))
resA$pval <- resA$fisher[which(1:nrow(resA$fisher)%%2==1),]
resA$OR <- resA$fisher[which(1:nrow(resA$fisher)%%2==0),]
