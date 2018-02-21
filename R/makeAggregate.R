makeAggregate <- function(tableName = 'collapsed ad meta modules February 21 2018'){
  cat('building DLPFC...\n')
  DLPFCres <- AMPAD::buildTargetedModules('DLPFC')
  save(DLPFCres,file='dlpfc_res.rda')

  #rm(DLPFCres)
  #gc()
  cat('building CBE...\n')
  CBEres <- AMPAD::buildTargetedModules('CBE')
  save(CBEres,file='cbe_res.rda')
  #rm(CBEres)
  #gc()
  cat('building TCX...\n')
  TCXres <- AMPAD::buildTargetedModules('TCX')
  save(TCXres,file='tcx_res.rda')
  #rm(TCXres)
  #gc()
  cat('building IFG...\n')
  IFGres <- AMPAD::buildTargetedModules('IFG')
  save(IFGres,file='ifg_res.rda')
  #rm(IFGres)
  #gc()
  cat('building STG...\n')
  STGres <- AMPAD::buildTargetedModules('STG')
  save(STGres,file='stg_res.rda')
  #rm(STGres)
  #gc()
  cat('building PHG...\n')
  PHGres <- AMPAD::buildTargetedModules('PHG')
  save(PHGres,file='phg_res.rda')
  #rm(PHGres)
  #gc()
  cat('building FP...\n')
  FPres <- AMPAD::buildTargetedModules('FP')
  save(FPres,file='fp_res.rda')
  #rm(FPres)
  #gc()
  cat('pushing to synapsely...\n')
  AggregateModuleManifest <- rbind(DLPFCres$df,
                                   CBEres$df,
                                   TCXres$df,
                                   IFGres$df,
                                   STGres$df,
                                   PHGres$df,
                                   FPres$df)


  rSynapseUtilities::makeTable(AggregateModuleManifest,tableName,projectId='syn2370594')
}
