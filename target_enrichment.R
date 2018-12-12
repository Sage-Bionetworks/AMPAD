#enrichment of nominated targets analysis
#pull targets from synapse

synapser::synLogin()

foo <- synapser::synGet('syn12540368')
bar <- data.table::fread(foo$path,data.table=F)

targetList <- list(AMPADTargets = bar$hgnc_symbol,
                   DrugTargets = c('ACHE',
                                   'APP',
                                   'PSEN1',
                                   'PSEN2',
                                   'MAPT',
                                   'RXRA',
                                   'GRIN1',
                                   'PLA2G7',
                                   'NGF',
                                   'BACE1',
                                   'BACE2',
                                   'GABRA1'))

library(synapser);detach('package:synapser',force=T,unload=T)
library(synapseClient);detach('package:synapseClient',force=T,unload=T)

#run enrichment analysis
baz <- AMPAD::run_amp_ad_enrichment(targetList,
                                    'Targets',
                                    manifestId = 'syn11932957')

ampTargets <- dplyr::filter(baz,category == 'AMPADTargets')
ampTargets <- dplyr::mutate(ampTargets,adj.p = p.adjust(fisherPval,method='bonferroni'))

drugTargets <- dplyr::filter(baz,category == 'DrugTargets')
drugTargets <- dplyr::mutate(drugTargets,adj.p = p.adjust(fisherPval,method='bonferroni'))
