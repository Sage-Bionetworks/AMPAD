buildModuleGraph <- function(pairwise,allMods,tissueType){
  res <- list()
  res$moduleGraph <- pairwise
  res$moduleGraph$weight <- res$moduleGraph$weight+1e-15
  graph1 <- igraph::graph_from_data_frame(res$moduleGraph,directed=FALSE)
  set.seed(2)
  test1 <- igraph::edge.betweenness.community(graph1)
  metaGraph <- data.frame(ModuleNameFull = test1$names,
                          metaModule = test1$membership,
                          stringsAsFactors=F)
  synapser::synLogin()
  moduleSet <- synapser::synTableQuery("SELECT DISTINCT ModuleNameFull, Module, method, brainRegion from syn10338156")$asDataFrame()
  #moduleSet <- moduleSet[,-c(1,2)]
  metaGraph2<-dplyr::left_join(metaGraph,moduleSet)




  res$mods<-lapply(unique(metaGraph2$metaModule),
                   AMPAD::bootstrapFrequencies,
                   metaGraph2,
                   allMods,
                   tissueType)

  res$moduleGraphCommunities <- metaGraph2

  modLen<-sapply(res$mods,function(x) length(x$ensg))
  res$mods <- res$mods[modLen>0]

  # getMajority <- function(df){
  #   masterTableHGNC <- table(df$external_gene_name)
  #   masterTableENSG <- table(df$GeneID)
  #   #masterTable <- masterTable/len1
  #   gen <- list()
  #   masterTableHGNC <- sort(masterTableHGNC,decreasing=T)
  #   masterTableENSG <- sort(masterTableENSG,decreasing=T)
  #   if(length(masterTableHGNC)< 1000 | length(masterTableENSG) < 1000){
  #     gen$ensg <- names(masterTableENSG)
  #     gen$hgnc <- names(masterTableHGNC)
  #   }else{
  #     gen$ensg <- names(masterTableENSG)[1:1000]
  #     gen$hgnc <- names(masterTableHGNC)[1:1000]
  #   }
  #   #masterTableHGNC <- which(masterTableHGNC > 1)
  #   #masterTableENSG <- which(masterTableENSG > 1)
  #
  #   return(gen)
  # }
  #38 41 17 53  9 16 27  4
  mods <- unique(res$moduleGraphCommunities$metaModule)[modLen>0]
  # dfList <- lapply(mods,function(x,df,allMods){
  #   library(dplyr)
  #   dplyr::filter(allMods,ModuleNameFull %in% df$ModuleNameFull[df$metaModule==x]) %>%
  #     return},
  #   res$moduleGraphCommunities,
  #   allMods)
  #
  # res$mods <- lapply(dfList,getMajority)
  names(res$mods) <- paste0(tissueType,WGCNA::labels2colors(mods))
  list2Df <- function(x,module,tissueType){
    mod <- data.frame(GeneID = x$ensg,
                      Module = rep(module,length(x$ensg)),
                      method = rep('aggregate',length(x$ensg)),
                      ModuleName = paste0('aggregate',rep(module,length(x$ensg))),
                      brainRegion = rep(tissueType,length(x$ensg)),
                      ModuleNameFull = paste0('aggregate',rep(module,length(x$ensg)),tissueType),
                      stringsAsFactors=F)
    return(mod)}
  res$df<-mapply(list2Df,
                 x = res$mods,
                 module = names(res$mods),
                 MoreArgs = list(tissueType=tissueType),
                 SIMPLIFY=F)
  res$df <- do.call(rbind,res$df)
  exg <- utilityFunctions::convertEnsemblToHgnc(res$df$GeneID)
  if(sum(duplicated(exg))>0){
    exg <- exg[which(!duplicated(exg)),]
  }
  res$df <- dplyr::left_join(res$df,exg,by = c('GeneID' = 'ensembl_gene_id'))
  return(res)
}
