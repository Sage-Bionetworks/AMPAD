library(EWCE)
library(ggplot2)
library(cowplot)
library(limma)
library(readxl)

foo <- data.table::fread('GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt',data.table=F)

lake_cortex <- list()
lake_cortex$expr <- as.matrix(foo[,-1])
lake_cortex$expr <- apply(lake_cortex$expr,2,as.numeric)
rownames(lake_cortex$expr) <- foo$V1
colnames(lake_cortex$expr) <- colnames(foo)[-1]

fxn1 <- function(x){
  return(strsplit(x,'\\_')[[1]])
}

res <- sapply(colnames(lake_cortex$expr),fxn1)
res <- t(res)

res <- data.frame(res,stringsAsFactors = F)
res$cell_id <- colnames(lake_cortex$expr)
colnames(res)[1:3] <- c('level1class','groupNo','UMI')
lake_cortex$annot <- res

synapser::synLogin()
aggMods <- synapser::synTableQuery("select * from syn11932957")$asDataFrame()
aggMods <- aggMods[,-c(1,2)]

testMod <- dplyr::filter(aggMods,Module=='PHGblue')$external_gene_name
bckgd <- unique(aggMods$external_gene_name)

generate.celltype.data(lake_cortex$expr, list(l1=lake_cortex$annot$level1class), "lake_cortex")


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


load('CellTypeData_lake_cortex.rda')

df2 <- c()
for (i in unique(aggMods$Module)){
  testMod <- dplyr::filter(aggMods,Module==i)$external_gene_name
  full_results = bootstrap.enrichment.test(sct_data=ctd,
                                           hits=testMod,
                                           bg=bckgd,
                                           genelistSpecies = "human",
                                           sctSpecies = "human",
                                           reps=1000,
                                           annotLevel=1)
  df1 <- full_results$results
  df1$Module <- i
  df2 <- rbind(df2,df1)
}

df2 <- dplyr::left_join(df2,customDf,by=c('Module'='moduleName'))

dummyDf <- df2
dummyDf$adj.p <- p.adjust(dummyDf$p,method='bonferroni')
dummyDf$fold_change[dummyDf$adj.p > 0.05] <- NA
dummyDf$Cluster <- factor(dummyDf$Cluster,levels = (c('Consensus Cluster A',
                                                      'Consensus Cluster B',
                                                      'Consensus Cluster C',
                                                      'Consensus Cluster D',
                                                      'Consensus Cluster E')))
dummyDf$Module <- factor(dummyDf$Module,levels = (c('TCXblue',
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
                                                    'PHGblue')))



g <- ggplot2::ggplot(dummyDf,
                     ggplot2::aes(x = Module,
                                  y = CellType,
                                  size = log2(fold_change),
                                  color = Cluster))
g <- g + ggplot2::geom_count()
#g <- g + ggplot2::scale_y_log10()

#g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
#g <- g + ggplot2::theme(axis.text.x=ggplot2::element_blank(),
#        axis.ticks.x=ggplot2::element_blank())
#g <- g + ggplot2::scale_color_gradientn(colours = c(viridis::viridis(2)[2], viridis::viridis(2)[1]),values = c(0,1), breaks = c(1.5, 3,8,20,50,110),trans='log')
#g <- g + ggplot2::coord_flip()
#g <- g + ggplot2::ggtitle('Enrichment for Cell Type Specific Signatures')
g <- g + ggplot2::labs(y = 'Lake et al. Cell Type Signature',
                       x = 'AD Coexpression Module')
g <- g + AMPAD::cowplot_rotated(11)

g

ggplot2::ggsave('lake2.tiff',device='tiff',units='mm',width=85,height=85,scale=1.8)

