compareADtoNonADmods <- function(ADaggMods,NonADaggMods,outputFile=FALSE){

  foobar <- ADaggMods

  AggregateModuleManifest <- NonADaggMods

  mod1=c(paste0('all',unique(AggregateModuleManifest$Module)),paste0('ad',unique(foobar$Module)))
  annrow <- data.frame(moduleType=c(rep('All',length(unique(AggregateModuleManifest$Module))),
                                    rep('AD',length(unique(foobar$Module)))),
                       stringsAsFactors=F)
  rownames(annrow) <- mod1


  mats <- utilityFunctions::pairwiseMatrixOfEnrichments(key=c(paste0('all',AggregateModuleManifest$Module),paste0('ad',foobar$Module)),
                                                        value=c(AggregateModuleManifest$GeneID,foobar$GeneID))
  pheatmap::pheatmap(-log10(mats$pval+1e-300),annotation_row=annrow)

  if(outputFile){
    tiff(filename='adnonadoverlap.tiff',height=85,width=85,units='mm',res=300)
    pheatmap::pheatmap(-log10(mats$pval+1e-300),
                       show_colnames = F,
                       border_color = NA,
                       fontsize=3,
                       treeheight_row=10,
                       treeheight_col=10,
                       annotation_col=annrow)
    dev.off()
  }else{
    pheatmap::pheatmap(-log10(mats$pval+1e-300),
                       show_colnames = F,
                       border_color = NA,
                       fontsize=3,
                       treeheight_row=10,
                       treeheight_col=10,
                       annotation_col=annrow)
  }

}
