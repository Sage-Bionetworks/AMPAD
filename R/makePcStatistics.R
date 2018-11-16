makePcStatistics <- function(date,synId='syn11932957'){

  aggMods <- rSynapseUtilities::loadFullTable(synId)

  geneExpressionForAnalysis <- AMPAD::pullExpressionAndPhenoWinsorized()
  names(geneExpressionForAnalysis) <- c('TCX','CBE','DLPFC','FP','STG','PHG','IFG')

  computeEigengene <- function(br,geneExp,moduleDefinitions){

    geneExp <- geneExp[[br]]

    #get modules
    mods <- dplyr::filter(moduleDefinitions,brainRegion==br)
    #convert modules into list of genes
    modsDefs <- lapply(unique(mods$ModuleNameFull),
                       utilityFunctions::listify,
                       mods$GeneID,
                       mods$ModuleNameFull)

    names(modsDefs) <- unique(mods$ModuleNameFull)

    internal <- function(mod,modsDefs,geneExp){
      geneExpMod <- dplyr::select(geneExp,modsDefs[[mod]])
      geneExpMod <- scale(geneExpMod)
      foo <- svd(geneExpMod)
      eigenGenes <- foo$u[,1:5]
      colnames(eigenGenes) <- paste0('pc',1:5)
      res <- cor(eigenGenes,geneExpMod)
      return(res)
    }

    full_res<-lapply(names(modsDefs),internal,modsDefs,geneExp)
    names(full_res) <- names(modsDefs)
    return(full_res)
  }

  fullList<-lapply(names(geneExpressionForAnalysis),computeEigengene,geneExpressionForAnalysis,aggMods)
  names(fullList) <- names(geneExpressionForAnalysis)

  #####make into a tidy data format with the following schema

  #module name
  #gene name
  #pc number
  #correlation

  #transpose
  #add column for gene
  #gather with tidyr
  #add column for modulename full
  #do.call rbind

  reformatCorrelations <- function(br,fullList){
    corMats <- fullList[[br]]
    internal <- function(modName,corMats){
      corMat <- corMats[[modName]]
      corMat <- t(corMat)
      corMat <- data.frame(corMat,stringsAsFactors = F)
      corMat$GeneID <- rownames(corMat)
      df <- tidyr::gather(corMat,
                          key = 'pc',
                          value ='correlation',
                          dplyr::starts_with('pc'))
      df$ModuleNameFull <- rep(modName,nrow(df))
      return(df)
    }
    foores<-lapply(names(corMats),internal,corMats)
    foores <- do.call(rbind,foores)
    return(foores)
  }

  foores2 <- lapply(names(fullList),reformatCorrelations,fullList)
  foores2 <- do.call(rbind,foores2)

  synapseClient::synapseLogin()
  rSynapseUtilities::makeTable(foores2,paste0('Gene Level Aggregate Module correlation with pcs ',date),'syn2370594')


  #share with Christoph
  #compute percentage variation explained

  computeEigengenePV <- function(br,geneExp,moduleDefinitions){

    geneExp <- geneExp[[br]]

    #get modules
    mods <- dplyr::filter(moduleDefinitions,brainRegion==br)
    #convert modules into list of genes
    modsDefs <- lapply(unique(mods$ModuleNameFull),
                       utilityFunctions::listify,
                       mods$GeneID,
                       mods$ModuleNameFull)

    names(modsDefs) <- unique(mods$ModuleNameFull)

    internal <- function(mod,modsDefs,geneExp){
      geneExpMod <- dplyr::select(geneExp,modsDefs[[mod]])
      geneExpMod <- scale(geneExpMod)
      foo <- svd(geneExpMod)
      pv <- foo$d^2/(sum(foo$d^2))
      names(pv) <- paste0('pc',1:length(pv))
      #eigenGenes <- foo$u[,1:5]
      #colnames(eigenGenes) <- paste0('pc',1:5)
      #res <- cor(eigenGenes,geneExpMod)
      return(pv)
    }

    full_res<-lapply(names(modsDefs),internal,modsDefs,geneExp)
    names(full_res) <- names(modsDefs)
    return(full_res)
  }

  fullList<-lapply(names(geneExpressionForAnalysis),computeEigengenePV,geneExpressionForAnalysis,aggMods)
  names(fullList) <- names(geneExpressionForAnalysis)

  cleanUpFxn <- function(br,fullList2){
    pvpcs <- fullList2[[br]]
    internal <- function(mod,pvpcs1){
      vec <- pvpcs1[[mod]]
      df <- data.frame(pc = names(vec)[1:5],
                       PV = vec[1:5],
                       ModuleNameFull = rep(mod,5),
                       stringsAsFactors=F)
      return(df)
    }
    resres <- lapply(names(pvpcs),internal,pvpcs)
    resres <- do.call(rbind,resres)
    return(resres)
  }
  resres2 <- lapply(names(fullList),cleanUpFxn,fullList)
  resres2 <- do.call(rbind,resres2)


  synapseClient::synapseLogin()
  rSynapseUtilities::makeTable(resres2,paste0('Percentage variation explained aggregate module PC analysis ',date),'syn2370594')
}
