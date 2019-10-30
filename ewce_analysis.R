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

testMod <- dplyr::filter(aggMods,Module=='DLPFCturquoise')$external_gene_name
bckgd <- unique(aggMods$external_gene_name)

lakedata<-generate.celltype.data(lake_cortex$expr, list(l1=lake_cortex$annot$level1class), "lake_cortex")

full_results = bootstrap.enrichment.test(sct_data=lake_cortex,
                                         hits=testMod,
                                         bg=bckgd,
                                         genelistSpecies = "human",
                                         sctSpecies = "human",
                                         reps=100,
                                         annotLevel=1)

