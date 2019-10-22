bootstrapFrequencies <- function(i,df,allMods,tissueType){
  set.seed(47)
  print(i)
  bar1 <- dplyr::filter(df,metaModule==i)

  bar2 <- lapply(unique(df$method),function(x,y){
    library(dplyr)
    dplyr::filter(y,
                  method==x) %>%
      return()
  },df)
  names(bar2) <- unique(df$method)
  bar3 <- table(bar1$method)
  bar3 <- bar3[names(bar2)]
  if(sum(is.na(bar3))>0){
    bar2 <- bar2[!is.na(bar3)]
    bar3 <- bar3[!is.na(bar3)]
  }
  fxn1 <- function(modCount,
                   df){
    foo1 <- sapply(1:100,
                   function(x,dta,modcount){
                     ind <- sample(1:nrow(dta),modcount)
                     return(dta$ModuleNameFull[ind])
                   },df,modCount)
    foo1 <- as.matrix(foo1)
    if(nrow(foo1)==100){
      foo1 <- t(foo1)
    }
    return(foo1)
    #sapply(1:1000,sample(1:nrow(df),modCount))
  }
  bar4 <- mapply(fxn1,
                 bar3,
                 bar2,
                 SIMPLIFY=F)
  getFrequencies <- function(modNameMat,allMods,refGene){
    library(dplyr)
    apply(modNameMat,2,function(x,allMods,refGene){
      genes<-dplyr::filter(allMods,ModuleNameFull %in% x)
      genes <- table(genes$GeneID)
      refGene[names(genes)]<-genes
      return(refGene)
    },allMods,refGene) %>% return
  }

  refGene <- unique(dplyr::filter(allMods,brainRegion==tissueType)$GeneID)
  refGene2 <- rep(0,length(refGene))
  names(refGene2) <- refGene

  modGenes <- dplyr::filter(allMods,ModuleNameFull %in% bar1$ModuleNameFull)$GeneID
  modGenes2 <- table(modGenes)
  refGene3 <- refGene2
  refGene3[names(modGenes2)] <- modGenes2


  frequencyMatrices <- lapply(bar4,getFrequencies,allMods,refGene2)
  names(frequencyMatrices) <- names(bar4)
  masterFrequencyMatrix <- Reduce("+",frequencyMatrices)
  fivePercent <- apply(masterFrequencyMatrix,1,quantile,.95)
  gen <- list()
  gen$ensg <- names(which(refGene3>fivePercent))
  return(gen)
}
