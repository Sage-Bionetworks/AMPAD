run_magma_analysis = function(outputFile=FALSE){


  #download stage 2 data
  str <- "curl -O https://www.niagads.org/system/tdf/public_docs/Kunkle_etal_Stage2_results.txt?file=1&type=field_collection_item&id=121&force="

  system(str)


  #rename file
  str <- "mv Kunkle_etal_Stage2_results.txt\\?file\\=1 Kunkle_etal_Stage2_results.txt"
  system(str)

  #load into R, and munge into appropriate format
  pvals <- data.table::fread('Kunkle_etal_Stage2_results.txt')

  #SNP,CHR,BP,P,NOBS
  pvalsReformatted <- dplyr::select(pvals,
                                    MarkerName,
                                    Chromosome,
                                    Position,
                                    Pvalue)

  colnames(pvalsReformatted) <- c('SNP',
                                  'CHR',
                                  'BP',
                                  'P')

  pvalsReformatted$NOBS <- 82771
  pvalsReformatted$P <- as.numeric(pvalsReformatted$P)

  write.table(pvalsReformatted,sep="\t",row.names=F,file="kunkle_stage2.txt",quote=F)

  write.table(pvalsReformatted[,1:3],sep="\t",row.names=F,col.names=F,file="kunkle2snp.txt",quote=F)

  str <- ".//magma --annotate --snp-loc kunkle2snp.txt --gene-loc NCBI37.3.gene.loc --out kunkleannots2"


  system(str)

  str <- ".//magma --bfile g1000_eur --pval kunkle_stage2.txt N=82771 --gene-annot kunkleannots2.genes.annot --out kunkle2 --gene-model multi=snp-wise"
  system(str)


  #create geneset files
  #pull manifest from synapse

  synapser::synLogin()
  ampMods <- synapser::synTableQuery("select * from syn11932957")$asDataFrame()


  #convert ensembl ids to ncbi ids

  convertEnsemblToEntrez <- function(ensemblIds){

    library(biomaRt)
    ensembl=biomaRt::useMart('ENSEMBL_MART_ENSEMBL',
                             dataset = 'hsapiens_gene_ensembl',
                             host='uswest.ensembl.org')

    genes<-biomaRt::getBM(attributes = c('ensembl_gene_id','entrezgene_id'),
                          filters='ensembl_gene_id',
                          values=ensemblIds,
                          mart=ensembl)
    return(genes)
  }

  mapDf <- convertEnsemblToEntrez(unique(ampMods$GeneID))
  ampModsRed <- dplyr::select(ampMods,Module,GeneID)
  ampModsRed <- dplyr::left_join(ampModsRed,mapDf,by=c('GeneID'='ensembl_gene_id'))

  #write to file
  write.table(pvalsReformatted[,1:3],sep="\t",row.names=F,col.names=F,file="kunklesnp.txt",quote=F)
  write.table(ampModsRed[,c(1,3)],sep="\t",row.names=F,col.names=F,file='geneset.txt',quote=F)

  str <- ".//magma --gene-results kunkle2.genes.raw --set-annot geneset.txt col=2,1 --out kunkle22"

  system(str)


  res <- data.table::fread('kunkle22.gsa.out',data.table=F,skip=4)
  res$p.adjust <- p.adjust(res$P,method='fdr')

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

  res <- dplyr::left_join(res,customDf,by=c('VARIABLE'='moduleName'))
  res$VARIABLE <- factor(res$VARIABLE,levels=c("TCXblue",
                                               "IFGyellow",
                                               "PHGyellow",
                                               "DLPFCblue",
                                               "CBEturquoise",
                                               "STGblue",
                                               "PHGturquoise",
                                               "IFGturquoise",
                                               "TCXturquoise",
                                               "FPturquoise",
                                               "IFGbrown",
                                               "STGbrown",
                                               "DLPFCyellow",
                                               "TCXgreen",
                                               "FPyellow",
                                               "CBEyellow",
                                               "PHGbrown",
                                               "DLPFCbrown",
                                               "STGyellow",
                                               "PHGgreen",
                                               "CBEbrown",
                                               "TCXyellow",
                                               "IFGblue",
                                               "FPblue",
                                               "FPbrown",
                                               "CBEblue",
                                               "DLPFCturquoise",
                                               "TCXbrown",
                                               "STGturquoise",
                                               "PHGblue"))

  g <- ggplot2::ggplot(res,ggplot2::aes(x=VARIABLE,y=p.adjust,fill=Cluster))
  g <- g+ggplot2::geom_col()
  g <- g +ggplot2::geom_hline(yintercept=0.10, linetype="dashed", color = "red")
  g <- g + ggplot2::xlab("Module")
  g <- g+ ggplot2::ylab("FDR")
  g <- g + ggplot2::scale_y_continuous(trans = "log10")
  g <- g + ggplot2::scale_fill_manual(values = c('#fefd11',
                                                  '#18bebf',
                                                  '#a82828',
                                                  '#34cc37',
                                                  '#470606'))
  #g <- g + ggplot2::theme_bw()
  #g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
  g <- g + AMPAD::cowplot_rotated(11)
  if(outputFile){
    g
    ggplot2::ggsave('magma_figure.tiff',device='tiff',units='mm',width=85,height=85,scale=1.8)
  }else{
    return(g)
  }

}
