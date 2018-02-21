getDataForModulePainting <- function(){
  source('dataPulling/pullExpressionAndPhenoWinsorized.R')
  source('pullADgenes.R')
  res <- list()
  names(geneExpressionForAnalysis) <- c('TCX',
                                        'CBE',
                                        'DLPFC',
                                        'FP',
                                        'STG',
                                        'PHG',
                                        'IFG')
  #####get module definition
  #synapseClient::synapseLogin()
  res$geneExpressionForAnalysis <- geneExpressionForAnalysis
  moduleSet <- synapseClient::synTableQuery("SELECT DISTINCT ModuleNameFull, Module, method, brainRegion from syn10915669")@values
  colnames(moduleSet)[c(3:4)] <- c('ModuleMethod','ModuleBrainRegion')
  res$moduleSet <- moduleSet
  adgenes <- pullADgenes()
  deggenes <- pullDegGenes()
  celltype <- pullCellTypeGenes()
  res$combined <- c(adgenes,
                    deggenes,
                    celltype)
  return(res)
}

pull_all_results <- function(moduleName,
                             annos,
                             geneExpressionForAnalysis,
                             moduleSet,
                             combined){
  #####pull expression data
  print(moduleName)
  print(annos)
  foo <- synapseClient::synTableQuery(paste0("select * from syn10915669 where ModuleNameFull =\'",moduleName,"\'"))@values
  modbr <- moduleSet$ModuleBrainRegion[moduleSet$ModuleName==moduleName]
  modIn <- function(x,y){
    return(y%in%x)
  }
  if(modbr=='CBE'){
    annos <- gsub('CBE','CER',annos)
  }
  paintMod <- sapply(combined,modIn,foo$external_gene_name)
  rownames(paintMod) <- foo$GeneID
  paintMod <- data.frame(paintMod,
                         stringsAsFactors = F)
  #paintMod <- paintMod[,annos]
  paintMod <- dplyr::select(paintMod,annos)

  paintMod$GeneID <- rownames(paintMod)
  paintMod <- dplyr::left_join(paintMod,foo)
  rownames(paintMod) <- paintMod$GeneID
  res <- list()
  exprMat <- geneExpressionForAnalysis[[modbr]][,paintMod$GeneID]
  rownames(exprMat) <- geneExpressionForAnalysis[[modbr]]$aSampleId
  res$rowAnnoMat <- geneExpressionForAnalysis[[modbr]][,c('aSampleId','logitDiagnosis')]
  annoMat <- data.matrix(paintMod[,1:length(annos)]) %>% data.frame()
  #annoMat <- apply(annoMat,2,as.factor) %>% data.frame()


  res$expr <- exprMat
  #res$anno <- annoMat
  res$anno <- paintMod

  #res$colAnno <- res$colAnno[,c(annos,'hubs')]

  bicNets <- synapseClient::synTableQuery(paste0("SELECT * FROM syn8681664 where ( (method = \'bic\') and (tissueTypeAbrv = \'",modbr,"\' )  and ( assay = \'RNAseq\'))"))@values

  load(synapseClient::synGet(bicNets$id[1])@filePath)
  library(Matrix)
  res$adjacencyMatrix <- as.matrix(bicNetworks$network[paintMod$GeneID,paintMod$GeneID])
  hubs <- data.frame(hubs = rowSums(res$adjacencyMatrix+t(res$adjacencyMatrix)),GeneID=rownames(res$adjacencyMatrix),stringsAsFactors=F)
  res$anno <- dplyr::left_join(res$anno,hubs)
  res$colAnno <- res$anno
  rownames(res$colAnno) <- res$anno$GeneID
  res$colAnno <- res$colAnno[,c(annos,'hubs')]
  return(res)
}
