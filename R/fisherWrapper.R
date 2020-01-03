fisherWrapper <- function(moduleGenes,annotationGenes,allGenes){
  vec1 <- allGenes%in%moduleGenes
  vec2 <- allGenes%in%annotationGenes
  a00 <- sum(!(vec1) & !(vec2))
  a10 <- sum((vec1) & !(vec2))
  a01 <- sum(!(vec1) & (vec2))
  a11 <- sum((vec1) & (vec2))
  bar <- matrix(c(a00,a10,a01,a11),2,2)
  #print(bar)
  foo <- fisher.test(bar,alternative='greater')
  return(c('pval'=foo$p.value,'OR'=as.numeric(foo$estimate)[1]))
}
