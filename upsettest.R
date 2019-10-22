#login to synapse
synapser::synLogin()

#pull mods
foo <- synapser::synTableQuery("select * from syn11932957")$asDataFrame()
foo2 <- dplyr::select(foo,GeneID,Module)
foo2$Presence <- 1
foo3 <- tidyr::pivot_wider(foo2,
                           id_cols = "GeneID",
                           names_from = "Module",
                           values_from = "Presence")
foo3[is.na(foo3)] <- 0
foo3 <- data.frame(foo3,stringsAsFactors=F)
foo4 <- dplyr::select(foo3,GeneID,TCXblue,IFGyellow,PHGyellow)

tiff(filename = 'consensusClusterA.tiff', height = 4, width = 6,units='in',pointsize=14,res=300)
UpSetR::upset(foo4,nintersects = NA,show.numbers=F)
dev.off()



foo4 <- dplyr::select(foo3,GeneID,DLPFCblue,CBEturquoise,STGblue,PHGturquoise,IFGturquoise,TCXturquoise,FPturquoise)
tiff(filename = 'consensusClusterB.tiff', height = 4, width = 6,units='in',pointsize=14,res=300)
UpSetR::upset(foo4,nsets=7,nintersects = NA,point.size=1,show.numbers = F)
dev.off()


foo4 <- dplyr::select(foo3,GeneID,IFGbrown,STGbrown,DLPFCyellow,TCXgreen,FPyellow,CBEyellow,PHGbrown)
tiff(filename = 'consensusClusterC.tiff', height = 4, width = 6,units='in',pointsize=14,res=300)
UpSetR::upset(foo4,nsets=7,nintersects = NA,point.size=1,show.numbers = F)
dev.off()

foo4 <- dplyr::select(foo3,GeneID,DLPFCbrown,STGyellow,PHGgreen,CBEbrown,TCXyellow,IFGblue,FPblue)
tiff(filename = 'consensusClusterD.tiff', height = 4, width = 6,units='in',pointsize=14,res=300)
UpSetR::upset(foo4,nsets=7,nintersects = NA,point.size=1,show.numbers = F)
dev.off()

foo4 <- dplyr::select(foo3,GeneID,FPbrown,CBEblue,DLPFCturquoise,TCXbrown,STGturquoise,PHGblue)
tiff(filename = 'consensusClusterE.tiff', height = 4, width = 6,units='in',pointsize=14,res=300)
UpSetR::upset(foo4,nsets=6,nintersects = NA,point.size=1,show.numbers = F)
dev.off()

