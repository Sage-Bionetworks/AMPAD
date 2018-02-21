splitByBrainRegionAdjustPvalue <- function(x){
  #split by brain region and category to adjust the p-values appropriately given the multiple hypothesis testing burden
  brs <- unique(x$ModuleBrainRegion)
  #print(brs)
  cats <- unique(x$GeneSetCategoryName)
  #print(cats)
  combined <- expand.grid(brs,cats)
  fxn1 <- function(y,z,t){
    foo1 <- dplyr::filter(t,ModuleBrainRegion==y & GeneSetCategoryName==z)
    foo1 <- dplyr::mutate(foo1,GeneSetAdjustedAssociationStatistic = p.adjust(GeneSetAssociationStatistic,method='fdr'))
    return(foo1)
  }
  foo2 <- mapply(fxn1,
                 combined[,1],
                 combined[,2],
                 MoreArgs = list(t=x),
                 SIMPLIFY = FALSE)
  foo2 <- do.call(rbind,foo2)
  return(foo2)
}
