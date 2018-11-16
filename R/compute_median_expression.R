compute_median_expression <- function(){
  library(synapser)
  synapser::synLogin()
  rosmapExprObj <- synapser::synGet("syn8456638")
  #rosmapClinObj <- synapser::synGet("syn8456631")


  rosmapExpr <- data.table::fread(rosmapExprObj$path,
                                  data.table=F)
  rosmapExpr2 <- tidyr::gather(rosmapExpr,key='sample',value='logCPM',-1)
  library(dplyr)
  rosmapSummMatrix <- dplyr::group_by(rosmapExpr2,ensembl_gene_id) %>%
    dplyr::summarise(medianLogCPM = median(logCPM))
  rosmapSummMatrix <- data.frame(rosmapSummMatrix,stringsAsFactors = F)
  rosmapSummMatrix$tissue <- 'DLPFC'

  mayoExprObj <- synapser::synGet("syn8466816")
  mayoExpr <- data.table::fread(mayoExprObj$path,
                                data.table=F)
  tcxIds <- grep('TCX',colnames(mayoExpr))
  cbeIds <- grep('CER',colnames(mayoExpr))

  mayoTCX <- mayoExpr[,c(1,tcxIds)]

  mayoTCX2 <- tidyr::gather(mayoTCX,
                            key='sample',
                            value='logCPM',
                            -1)
  mayoTCXSummMatrix <- dplyr::group_by(mayoTCX2,ensembl_gene_id) %>%
    dplyr::summarise(medianLogCPM = median(logCPM))
  mayoTCXSummMatrix <- data.frame(mayoTCXSummMatrix,stringsAsFactors = F)
  mayoTCXSummMatrix$tissue <- 'TCX'


  mayoCBE <- mayoExpr[,c(1,cbeIds)]
  mayoCBE2 <- tidyr::gather(mayoCBE,
                            key='sample',
                            value='logCPM',
                            -1)
  mayoCBESummMatrix <- dplyr::group_by(mayoCBE2,ensembl_gene_id) %>%
    dplyr::summarise(medianLogCPM = median(logCPM))
  mayoCBESummMatrix <- data.frame(mayoCBESummMatrix,stringsAsFactors = F)
  mayoCBESummMatrix$tissue <- 'CBE'


  mssmExprObj <- synapser::synGet('syn8485017')
  mssmExpr <- data.table::fread(mssmExprObj$path)
  mssmClinObj <- synapser::synGet('syn8484996')
  mssmClin <- data.table::fread(mssmClinObj$path)

  mssmExpr2 <- tidyr::gather(mssmExpr,
                      key='sample',
                      value='logCPM',
                      -1)
  mssmClin$tissue <- sapply(mssmClin$Tissue.Diagnosis,function(x) strsplit(x,'\\.')[[1]][1])

  mssmExpr2 <- dplyr::left_join(mssmExpr2,dplyr::select(mssmClin,SampleID,tissue),by=c('sample'='SampleID'))

  mssmSTG <- dplyr::filter(mssmExpr2,tissue=='STG') %>%
    dplyr::select(-tissue)
  mssmSTGSummMatrix <- dplyr::group_by(mssmSTG,ensembl_gene_id) %>%
    dplyr::summarise(medianLogCPM = median(logCPM))
  mssmSTGSummMatrix <- data.frame(mssmSTGSummMatrix,
                                  stringsAsFactors=F)
  mssmSTGSummMatrix$tissue <- 'STG'


  mssmFP <- dplyr::filter(mssmExpr2,tissue=='FP') %>%
    dplyr::select(-tissue)
  mssmFPSummMatrix <- dplyr::group_by(mssmFP,ensembl_gene_id) %>%
    dplyr::summarise(medianLogCPM = median(logCPM))
  mssmFPSummMatrix <- data.frame(mssmFPSummMatrix,
                                  stringsAsFactors=F)
  mssmFPSummMatrix$tissue <- 'FP'


  mssmPHG <- dplyr::filter(mssmExpr2,tissue=='PHG') %>%
    dplyr::select(-tissue)
  mssmPHGSummMatrix <- dplyr::group_by(mssmPHG,ensembl_gene_id) %>%
    dplyr::summarise(medianLogCPM = median(logCPM))
  mssmPHGSummMatrix <- data.frame(mssmPHGSummMatrix,
                                  stringsAsFactors=F)
  mssmPHGSummMatrix$tissue <- 'PHG'

  mssmIFG <- dplyr::filter(mssmExpr2,tissue=='IFG') %>%
    dplyr::select(-tissue)
  mssmIFGSummMatrix <- dplyr::group_by(mssmIFG,ensembl_gene_id) %>%
    dplyr::summarise(medianLogCPM = median(logCPM))
  mssmIFGSummMatrix <- data.frame(mssmIFGSummMatrix,
                                  stringsAsFactors=F)
  mssmIFGSummMatrix$tissue <- 'IFG'


  fullDf <- rbind(rosmapSummMatrix,
                  mayoTCXSummMatrix,
                  mayoCBESummMatrix,
                  mssmSTGSummMatrix,
                  mssmFPSummMatrix,
                  mssmPHGSummMatrix,
                  mssmIFGSummMatrix)

  return(fullDf)


}
