run_full_deg_eigenanalysis = function(){
  #pull module definitions
  synapseClient::synapseLogin()
  agg_mods <- synapseClient::synTableQuery('SELECT * FROM syn11932957')@values

  synId <- 'syn11932957'
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
      eigenGenes <- foo$v[,1:5]
      colnames(eigenGenes) <- paste0('pc',1:5)

      #res <- cor(eigenGenes,geneExpMod)
      res <- eigenGenes
      #rownames(res) <- rownames(geneExpMod)
      return(res)
    }

    full_res<-lapply(names(modsDefs),internal,modsDefs,geneExp)
    names(full_res) <- names(modsDefs)
    return(full_res)
  }

  fullList<-lapply(names(geneExpressionForAnalysis),computeEigengene,geneExpressionForAnalysis,aggMods)
  names(fullList) <- names(geneExpressionForAnalysis)

  CBEbrown <- fullList$CBE$aggregateCBEbrownCBE
  CBEbrown <- data.frame(CBEbrown,stringsAsFactors = F)
  CBEbrown$sampleID <- geneExpressionForAnalysis$CBE$aSampleId

  #get covariates TCX
  mayoCovObj <- synapseClient::synGet('syn8466814')
  mayoCov <- data.table::fread(mayoCovObj@filePath,data.table=F)

  CBEbrown <- dplyr::left_join(CBEbrown,mayoCov,by=c('sampleID'='SampleID'))
  CBEbrown <- dplyr::filter(CBEbrown,Tissue.Diagnosis!='CBE.OTHER')


  TCXyellow$Tissue.Diagnosis <- factor(TCXyellow$Tissue.Diagnosis,levels = rev(c('TCX.AD','TCX.CONTROL')))
  TCXyellow$Sex <- factor(TCXyellow$Sex,levels = rev(c('FEMALE','MALE')))

  #pull DEG results
  #ROSMAP: syn8456721
  fob <- synapseClient::synGet('syn8456721')
  bar <- data.table::fread(fob@filePath,data.table=F)

  diagnosis <- dplyr::filter(bar,Model=='Diagnosis')
  diagnosisMale <- dplyr::filter(bar,Model=='Diagnosis.Sex' & Comparison =='AD-CONTROL.IN.MALE')
  diagnosisFemale <- dplyr::filter(bar,Model=='Diagnosis.Sex' & Comparison =='AD-CONTROL.IN.FEMALE')

  dx <- dplyr::filter(diagnosis,ensembl_gene_id%in%mods[['TCXyellow']])
  dxM <- dplyr::filter(diagnosisMale,ensembl_gene_id%in%mods[['TCXyellow']])
  dxF <- dplyr::filter(diagnosisFemale,ensembl_gene_id%in%mods[['TCXyellow']])

  #DLPFCbrown
  #STGyellow
  #PHGgreen
  #CBEbrown
  #TCXyellow
  #IFGblue
  #FPblue



  #run model similarly as Thanneer does for all genes
  #models: sex stratified
  #models: sex adjusted

}
