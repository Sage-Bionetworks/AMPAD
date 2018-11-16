makeEnrichmentPCAplot <- function(modules = c('aggregateFPbrownFP'),
                                  cached='cachedGeneSets.rda'){
  #pull on pathways of enrichment for pathways
  synapseClient::synapseLogin()
  foobar <- synapseClient::synTableQuery('select * from syn11954640')@values
  foobar <- dplyr::mutate(foobar,adj.pval = p.adjust(fisherPval,method='fdr'))
  foobar2 <- dplyr::filter(foobar,ModuleNameFull %in% modules & adj.pval <= 0.05)
  #get pathway definitions
  if(is.null(cached)){
    pathwayGeneSets <- AMPAD::collatePathways()
  }else{
    load(cached)
  }


  #convert from list 2 df
  pathwayGeneSets2 <- lapply(pathwayGeneSets,utilityFunctions::list2df)

  #add pathway source annotation
  #x: pathwayGeneSets2
  #y: name of pathway source
  pathwayGeneSets3 <- mapply(function(x,y){
    new_x <- cbind(x,rep(y,nrow(x)))
    colnames(new_x)[3] <- 'pathwaySource'
    new_x$key <- make.names(new_x$key)
    return(new_x)
  }, pathwayGeneSets2,names(pathwayGeneSets2),SIMPLIFY = F)
  pathwayGeneSets3 <- do.call(rbind,pathwayGeneSets3)
  pathwayGeneSets4 <- dplyr::mutate(pathwayGeneSets3,fullPathway = paste0(key,pathwaySource))
  #filter gene lists down
  foobar2 <- dplyr::mutate(foobar2,fullPathway = paste0(category,geneSet))
  #foobar2 <- dplyr::left_join(foobar2,pathwayGeneSets4,by=c('fullPathway'))
  rm(pathwayGeneSets3)
  rm(pathwayGeneSets2)
  rm(pathwayGeneSets)
  gc()

  #summarize

  pathwayGeneSets4 <- dplyr::filter(pathwayGeneSets4,fullPathway %in% foobar2$fullPathway)

  pathwayGeneSets5 <- dplyr::mutate(pathwayGeneSets4,isPresent = TRUE)
  pathwayGeneSets5 <- dplyr::select(pathwayGeneSets5,-key,-pathwaySource)
  pathwayGeneSets5 <- tidyr::spread(pathwayGeneSets5,fullPathway,isPresent,fill=FALSE)

  rownames(pathwayGeneSets5) <- pathwayGeneSets5$value
  adjMat <- dplyr::select(pathwayGeneSets5,-value)
  adjMat <- data.matrix(adjMat)
  #di
  adjMat <- t(adjMat)
  adjMat <- scale(adjMat)
  adjMat <- t(adjMat)
  #adjMat <- t(adjMat)
  #adjMat <- scale(adjMat)
  #fisherOR <- foobar2$fisherOR
  #names(fisherOR) <- foobar2$fullPathway
  #fisherOR[fisherOR==Inf] <- 55

  #adjMat <- t(adjMat)
  #adjMat <- adjMat * fisherOR[rownames(adjMat)]

  mat1 <- cor((adjMat),method = 'spearman')
  res_tsne <- tsne::tsne(mat1)
  colnames(res_tsne) <- c('tsne1','tsne2')
  res_tsne <- data.frame(res_tsne,stringsAsFactors=F)
  res_tsne$fullPathway <- colnames(adjMat)
  foobar2$fisherOR[foobar2$fisherOR==Inf] <- 55
  foobar3 <- dplyr::group_by(foobar2,fullPathway)
  foobar4 <- dplyr::summarise(foobar3,percentSig=length(ModuleNameFull)/length(modules),meanOR = mean(fisherOR),meanpval = mean(-log10(adj.pval)))

  res_tsne <- dplyr::left_join(res_tsne,
                               foobar4)
  set.seed(1)
  kmeans2 <- kmeans(res_tsne[,1:2],40,nstart = 1000)
  res_tsne$cluster <- kmeans2$cluster
  #res_tsne2 <- dplyr::filter(res_tsne,
  #                           percentSig>=0.5)
  #res_tsne2 <- res_tsne[grep('go_bp',res_tsne$fullPathway),]
  #res_tsne$fisherOR[res_tsne$fisherOR==Inf] <- 55
  p <- plotly::plot_ly(data=res_tsne,
                       x = ~tsne1,
                       y = ~tsne2,
                       size = ~percentSig,
                       color = ~cluster,
                       text = ~fullPathway)


  #add pathway, source, fisher, -log10 pvalue

  #unnormalized svd
  return(p)
  #reorganize based on gene name -> genes as columns, rows are pathway names
  #make a distance matrix based on similarity of the pathways


}
