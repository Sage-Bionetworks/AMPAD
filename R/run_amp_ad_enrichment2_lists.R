run_amp_ad_enrichment2_lists <- function(geneSetList,
                                  geneSetName,
                                  refGeneSetList,
                                  hgnc = TRUE,
                                  testhgnc=FALSE){
  #INPUT:
  #geneSetList - a list of genes in hgnc or ensembl format
  #geneSetName - name of geneset (should be a single character string)
  #hgnc - boolean indicating whether the gene identifiers are hgnc (TRUE) or ensembl (FALSE)
  #manifestId - location of a synapse table with the modules definitions in terms of hgnc names


  #OUTPUT:
  #data frame with module name, gene set name, p-value, and odds ratio of enrichment, intersections, and gene set sizes
  straightHgncConversion <- function(x){
    y<-utilityFunctions::convertHgncToEnsembl(x)
    return(unique(y$ensembl_gene_id))
  }

  straightEnsemblConversion <- function(x){
    y<- utilityFunctions::convertEnsemblToHgnc(x)
    return(unique(y$external_gene_name))
  }

  if((hgnc!=testhgnc) & (!testhgnc)){
    geneSetList <- lapply(geneSetList,straightEnsemblConversion)
  }else if((hgnc!=testhgnc) & (testhgnc)){
    geneSetList <- lapply(geneSetList,straightHgncConversion)
  }


  library(dplyr)
  cat('removing genes that are not relevant from reference set...\n')
  #get unique gene keys,drop categories in both cases that are 0 in size
  uniqueModuleList <- refGeneSetList %>%
    unlist %>%
    unique

  uniqueGeneSet <- geneSetList %>%
    unlist %>%
    unique

  refGeneSet <- uniqueModuleList

  cat('running enrichments....\n')

  res <- list()
  res$fisher <- utilityFunctions::outerSapplyParallel(utilityFunctions::fisherWrapper,
                                                      refGeneSetList,
                                                      geneSetList,
                                                      refGeneSet)
  #pvalues are odd rows, odds ratios are even rows
  res$pval <- res$fisher[which(1:nrow(res$fisher)%%2==1),]
  rownames(res$pval) <- names(geneSetList)
  res$OR <- res$fisher[which(1:nrow(res$fisher)%%2==0),]
  rownames(res$OR) <- names(geneSetList)

  sizeOfInter <- function(x,y){
    return(length(intersect(x,y)))
  }

  res$inter <- utilityFunctions::outerSapplyParallel(sizeOfInter,
                                                     refGeneSetList,
                                                     geneSetList)


  cat('producing tidy data frame....\n')
  #transpose pvalues
  pval1 <- t(res$pval)
  #make into a data frame
  pval1 <- data.frame(pval1,stringsAsFactors=F)
  #create a unique key for each row for gather step
  pval1 <- dplyr::mutate(pval1,ModuleNameFull = rownames(pval1))
  #go from matrix form - module by enrichment categories - to long table form
  gatherTest1 <- tidyr::gather(pval1,category,fisherPval,-ModuleNameFull)

  #transpose odds ratios
  or1 <- t(res$OR)
  #make into a data frame
  or1 <- data.frame(or1,stringsAsFactors=F)
  #create a unique key for each row for gather step
  or1 <- dplyr::mutate(or1,ModuleNameFull = rownames(or1))
  #go from matrix form - module by enrichment categories - to long table form
  gatherTest2 <- tidyr::gather(or1,category,fisherOR,-ModuleNameFull)

  #transpose intersection sizes
  inter1 <- t(res$inter)
  #make into a data frame
  inter1 <- data.frame(inter1,stringsAsFactors=F)
  #create a unique key for each row for gather step
  inter1 <- dplyr::mutate(inter1,ModuleNameFull = rownames(inter1))
  #go from matrix form - module by enrichment categories - to long table form
  gatherTest3 <- tidyr::gather(inter1,category,nInter,-ModuleNameFull)


  ###get sizes
  aggModSize <- data.frame(ModuleNameFull = names(refGeneSetList),
                           mod_size = sapply(refGeneSetList,length),
                           stringsAsFactors = F)
  categorySize <- data.frame(category = make.names(names(geneSetList)),
                             category_size = sapply(geneSetList,length),
                             stringsAsFactors=F)



  #do a left join to combine the pvalues and odds ratios
  gatherTest <- dplyr::left_join(gatherTest1,
                                 gatherTest2)
  gatherTest <- dplyr::left_join(gatherTest,
                                 gatherTest3)
  gatherTest <- dplyr::left_join(gatherTest,
                                 aggModSize)
  gatherTest <- dplyr::left_join(gatherTest,
                                 categorySize)
  #add in the geneset names
  gatherTest <- dplyr::mutate(gatherTest,
                              geneSet = geneSetName)


  #return the data frame
  return(gatherTest)

}
