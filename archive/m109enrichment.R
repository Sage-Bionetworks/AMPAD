#synapseClient::synapseLogin()

#m109 <- synapseClient::synTableQuery('select * from syn5321231 where speakeasyModule = 109')@values$geneName

#m9 <- synapseClient::synTableQuery('select * from syn5321231 where speakeasyModule = 9')@values$geneName

#res <- AMPAD::run_amp_ad_enrichment(list(m109=m109,m9=m9),'speakeasy',hgnc=F,manifestId='syn11932957')
