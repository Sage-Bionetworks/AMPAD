run_full_enrichment_suite <- function(){
  targetedGeneSets <- AMPAD::collateEnrichmentSets()

  targetedEnrichment <- list()
  targetedEnrichment$ad <- AMPAD::run_amp_ad_enrichment2(targetedGeneSets$ad,
                                                         'AD',
                                                         hgnc = TRUE,
                                                         manifestId = 'syn11932957')

  targetedEnrichment$cell <- AMPAD::run_amp_ad_enrichment2(targetedGeneSets$cell,
                                                           'Cell',
                                                           hgnc= TRUE,
                                                           manifestId = 'syn11932957')

  targetedEnrichment$cell2 <- AMPAD::run_amp_ad_enrichment2(targetedGeneSets$cell2,
                                                            'Mathys',
                                                            hgnc= TRUE,
                                                            manifestId = 'syn11932957')

  targetedEnrichment$scz <- AMPAD::run_amp_ad_enrichment2(targetedGeneSets$scz,
                                                          'scz',
                                                          hgnc= FALSE,
                                                          manifestId = 'syn11932957')

  targetedEnrichment$control <- AMPAD::run_amp_ad_enrichment2(targetedGeneSets$control,
                                                              'control',
                                                              hgnc= FALSE,
                                                              manifestId = 'syn11932957')

  targetedEnrichment$mssm2 <- AMPAD::run_amp_ad_enrichment2(targetedGeneSets$MSSM2,
                                                            'zhang',
                                                            hgnc= TRUE,
                                                            manifestId = 'syn11932957')

  targetedEnrichment$mito <- AMPAD::run_amp_ad_enrichment2(targetedGeneSets$mito,
                                                           'Mito',
                                                           hgnc= TRUE,
                                                           manifestId = 'syn11932957')

  targetedEnrichment$deg <- AMPAD::run_amp_ad_enrichment2(targetedGeneSets$deg,
                                                          'DEG',
                                                          hgnc=FALSE,
                                                          manifestId = 'syn11932957')

  targetedEnrichment$degMeta <- AMPAD::run_amp_ad_enrichment2(targetedGeneSets$degMeta,
                                                              'DEGmeta',
                                                              hgnc=FALSE,
                                                              manifestId = 'syn11932957')

  targetedEnrichment$targetPathways <- AMPAD::run_amp_ad_enrichment2(targetedGeneSets$targetedPathways,
                                                                     'TargetedPathways',
                                                                     hgnc=FALSE,
                                                                     manifestId = 'syn11932957')

  targetedEnrichment$ad_deg <- AMPAD::run_amp_ad_enrichment2_lists(targetedGeneSets$deg,
                                                                   'ADDEG',
                                                                   targetedGeneSets$ad,
                                                                   hgnc = T,
                                                                   testhgnc = F)

  targetedEnrichment$ad_degmeta <- AMPAD::run_amp_ad_enrichment2_lists(targetedGeneSets$degMeta,
                                                                       'ADDEGMeta',
                                                                       targetedGeneSets$ad,
                                                                       hgnc = T,
                                                                       testhgnc = F)

  targetedEnrichment$degFull <- rbind(targetedEnrichment$deg,targetedEnrichment$ad_deg)
  targetedEnrichment$degMetaFull <- rbind(targetedEnrichment$degMeta, targetedEnrichment$ad_degmeta)
  targetedEnrichment$ad$adj.pval <- p.adjust(targetedEnrichment$ad$fisherPval,method='fdr')
  targetedEnrichment$deg$adj.pval <- p.adjust(targetedEnrichment$deg$fisherPval,method='fdr')
  targetedEnrichment$cell$adj.pval <- p.adjust(targetedEnrichment$cell$fisherPval,method='fdr')
  targetedEnrichment$cell2$adj.pval <- p.adjust(targetedEnrichment$cell2$fisherPval,method='fdr')
  targetedEnrichment$scz$adj.pval <- p.adjust(targetedEnrichment$scz$fisherPval,method='fdr')
  targetedEnrichment$control$adj.pval <- p.adjust(targetedEnrichment$control$fisherPval,method='fdr')
  targetedEnrichment$mito$adj.pval <- p.adjust(targetedEnrichment$mito$fisherPval,method='fdr')
  targetedEnrichment$mssm2$adj.pval <- p.adjust(targetedEnrichment$mssm2$fisherPval,method='fdr')
  targetedEnrichment$degMeta$adj.pval <- p.adjust(targetedEnrichment$degMeta$fisherPval,method='fdr')
  targetedEnrichment$targetPathways$adj.pval <- p.adjust(targetedEnrichment$targetPathways$fisherPval,method='fdr')
  targetedEnrichment$degFull$adj.pval <- p.adjust(targetedEnrichment$degFull$fisherPval,method='fdr')
  targetedEnrichment$degMetaFull$adj.pval <- p.adjust(targetedEnrichment$degMetaFull$fisherPval,method='fdr')
  return(targetedEnrichment)
}
