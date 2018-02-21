getDeg <- function(synId='syn10338156'){
  moduleSet <- synapseClient::synTableQuery(paste0("SELECT DISTINCT ModuleNameFull, Module, method, brainRegion from ",synId))@values
  colnames(moduleSet)[c(3:4)] <- c('ModuleMethod','ModuleBrainRegion')
  degResObj <- synapseClient::synGet("syn10496554")
  load(degResObj@filePath)
  #source('enrichmentAnalysis/run_amp_ad_enrichment.R')
  keep <- grep('^Diagnosis',names(amp.ad.de.geneSets))
  keep <- intersect(keep,grep('AD-CONTROL',names(amp.ad.de.geneSets)))
  amp.ad.de.geneSets <- amp.ad.de.geneSets[keep]
  degResults <- AMPAD::run_amp_ad_enrichment(amp.ad.de.geneSets,
                                      "degs",
                                      hgnc = FALSE,
                                      manifestId = synId)



  parseDegName <- function(x){
    library(dplyr)
    foo1 <- strsplit(x,'\\.')[[1]]

    if(foo1[2]=='SEX' | foo1[2]=='Sex'){
      br <- foo1[3]
      cate <- paste0(foo1[c(1:2,4:(length(foo1) - 1))],collapse='_')
    }else{
      br <- foo1[2]
      cate <- paste0(foo1[c(1,3:(length(foo1) - 1))],collapse = '_')
    }

    #br <- foo1[2]
    dir <- foo1[length(foo1)]
    #cate <- foo1[2]

    #if(length(grep(paste0('.',br,'_'),cate)) > 0) {
    #  cate <- gsub(paste0('.',br,'_'),'.',cate)
    #}

    c('brainRegion' = br,
      'Direction' = dir,
      'reducedCategory' = cate,
      'Category' = x) %>% return
  }

  categoryKey <- sapply(unique(degResults$category),
                        parseDegName)
  categoryKey <- t(categoryKey)
  categoryKey <- data.frame(categoryKey,stringsAsFactors = F)

  degResults2 <- dplyr::left_join(degResults,categoryKey,by = c('category' = 'Category'))
  degResultsModified <- dplyr::select(degResults2,
                                      ModuleNameFull,
                                      reducedCategory,
                                      geneSet,
                                      fisherPval,
                                      fisherOR,
                                      brainRegion,
                                      Direction)

  colnames(degResultsModified) <- c('ModuleNameFull',
                                    'GeneSetName',
                                    'GeneSetCategoryName',
                                    'GeneSetAssociationStatistic',
                                    'GeneSetEffect',
                                    'GeneSetBrainRegion',
                                    'GeneSetDirectionAD')

  moduleSummaryDeg <- dplyr::left_join(moduleSet,degResultsModified)

  #moduleSummaryDeg$GeneSetBrainRegion[moduleSummaryDeg$GeneSetBrainRegion=='CER']<-'CBE'

  ###match brain regions for clarity sake
  moduleSummaryDeg <- dplyr::filter(moduleSummaryDeg,ModuleBrainRegion==GeneSetBrainRegion)

  ##define which ones are ad related
  #ad_related <- grep('AD',moduleSummaryDeg$GeneSetName)
  moduleSummaryDeg$GeneSetADLinked <- rep(TRUE,nrow(moduleSummaryDeg))
  #moduleSummaryDeg$GeneSetADLinked[ad_related] <- TRUE

  moduleSummaryDeg <- AMPAD::splitByBrainRegionAdjustPvalue(moduleSummaryDeg)
  moduleSummaryDeg$EvidenceClass <- rep('deg',nrow(moduleSummaryDeg))
  return(moduleSummaryDeg)
}
