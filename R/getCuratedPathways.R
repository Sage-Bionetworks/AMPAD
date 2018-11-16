getCuratedPathways <- function(wpIds = c('WP2714',
                                          'WP366',
                                          'WP2677',
                                          'WP3547',
                                          'WP3869',
                                          'WP1877',
                                          'WP2593',
                                          'WP2380',
                                          'WP382',
                                          'WP1798',
                                          'WP4096')){

  # hippo signalling - WP2714 - astrogliosis (primary reference)
  # tgfb signalling  - WP366 - neurofinlammation - (primary reference)
  # igf1 signalling - WP2677 - female specific metabolic changes (brinton paper)
  # amyloid fibril formation - WP3547 - amyloid hypothesis
  # endocanniboid signaling - WP3869 - neurofinlammation regulation
  # downregulation of g-lymphatic flow - WP1877 - reduced glymphatic flow
  # jak/stat signaling - WP2593 - astrogliosis
  # bdnf pathway - WP2380 - neuroprotection/aging
  # mapk signaling - WP382 - tau hyperphosphorylation
  # complement cascade - WP1798 - synaptic pruning/neuroinflammation


  curatedList <- lapply(wpIds,rWikiPathways::getXrefList,'En')
  names(curatedList) <- wpIds
  return(curatedList)
}
