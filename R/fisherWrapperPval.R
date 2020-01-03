fisherWrapperPval <- function(moduleGenes,annotationGenes,allGenes){
  a00 <- sum(!(allGenes%in%moduleGenes) & !(allGenes%in%annotationGenes))
  a10 <- sum((allGenes%in%moduleGenes) & !(allGenes%in%annotationGenes))
  a01 <- sum(!(allGenes%in%moduleGenes) & (allGenes%in%annotationGenes))
  a11 <- sum((allGenes%in%moduleGenes) & (allGenes%in%annotationGenes))
  bar <- matrix(c(a00,a10,a01,a11),2,2)
  #print(bar)
  foo <- fisher.test(bar,alternative='greater')
  return(foo$p.value)
}
