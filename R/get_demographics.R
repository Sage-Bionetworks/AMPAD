get_demographics = function(){
  synapseClient::synapseLogin()
  covariateManifest <- synapseClient::synTableQuery("SELECT * FROM syn9704300 WHERE ( ( normalizationType = 'CQN' ) AND ( dataSubType = 'covariates' ) )")

  covariateList <- rSynapseUtilities::loadDelimIntoList(covariateManifest@values$id)
  covariateList$syn8484996 <- synapseClient::synGet('syn8484996',version = 8) %>%
    synapseClient::getFileLocation() %>%
    data.table::fread(data.table=F)

  covariateList$syn8456631 <- synapseClient::synGet('syn8456631',version = 15) %>%
    synapseClient::getFileLocation() %>%
    data.table::fread(data.table=F)

  covariateList$syn8466814 <- synapseClient::synGet('syn8466814',version = 12) %>%
    synapseClient::getFileLocation() %>%
    data.table::fread(data.table=F)


  # demographics
  #Study N Sex tissue type (group by) mean(AOD) sd (AOD) N_males N_females mean(PMI) sd(PMI) mean(RIN) sd(RIN) N Diagnosis

  rosmap <- covariateList$syn8456631
  uncensoredAgesObj <- synapseClient::synGet('syn7116000')
  uncensoredAges <- data.table::fread(uncensoredAgesObj@filePath,data.table=F)
  idkeyobj <- synapseClient::synGet('syn3382527')
  idkey <- data.table::fread(idkeyobj@filePath,data.table=F)
  idkey <- dplyr::select(idkey,projid,rnaseq_id)
  key <- dplyr::left_join(uncensoredAges,idkey,by='projid')
  rosmap <- dplyr::left_join(rosmap,key,by=c('SampleID'='rnaseq_id'))
  dups <- which(duplicated(rosmap$SampleID))
  rosmap <- rosmap[-dups,]

  rosmap2 <- dplyr::summarise(rosmap,
                              N=length(unique(SampleID)),
                              meanAOD = mean(age_death.y),
                              sdAOD = sd(age_death.y),
                              Nfemales = sum(1-msex),
                              Nmales = sum(msex),
                              NAD = sum(Diagnosis == 'AD'),
                              NControl = sum(Diagnosis=='CONTROL'),
                              meanPMI = mean(pmi),
                              sdPMI = sd(pmi),
                              meanRIN = mean(RINcontinuous),
                              sdRIN = sd(RINcontinuous),
                              NADfemales = sum((msex==0) & (Diagnosis == 'AD')),
                              NControlfemales = sum((msex==0) & (Diagnosis == 'CONTROL')),
                              NADmales = sum((msex==1) & (Diagnosis == 'AD')),
                              NControlmales = sum((msex == 1) & (Diagnosis == 'CONTROL')))


  #mayo
  mayo <- covariateList$syn8466814
  mayo$brainRegion <- sapply(mayo$BrainRegion.Diagnosis,function(x) strsplit(x,'\\.')[[1]][1])
  mayo$Diagnosis <- sapply(mayo$BrainRegion.Diagnosis,function(x) strsplit(x,'\\.')[[1]][2])
  mayo2 <- dplyr::group_by(mayo,brainRegion)
  mayo3 <- dplyr::summarise(mayo2,
                            NSample=length(unique(SampleID)),
                            NDonor=length(unique(Donor_ID)),
                            meanAOD = mean(AgeAtDeath),
                            sdAOD = sd(AgeAtDeath),
                            Nfemales = sum(Gender=='F'),
                            Nmales = sum(Gender=='M'),
                            NAD = sum(Diagnosis == 'AD'),
                            NControl = sum(Diagnosis == 'CONTROL'),
                            meanPMI = mean(PMI),
                            sdPMI = sd(PMI),
                            meanRIN = mean(RIN),
                            sdRIN = sd(RIN),
                            NADfemales = sum((Gender == 'F')&(Diagnosis == 'AD')),
                            NControlfemales = sum((Gender == 'F')&(Diagnosis == 'CONTROL')),
                            NADmales = sum((Gender == 'M')&(Diagnosis == 'AD')),
                            NADcontrol = sum((Gender == 'M')&(Diagnosis == 'CONTROL')))

  msbb <- covariateList$syn8484996

  msbbRealAgesObj  <- synapseClient::synGet('syn10156693')
  msbbRealAge <- data.table::fread(msbbRealAgesObj@filePath,data.table=F)

  msbbReduxObj <- synapseClient::synGet('syn6100548')
  msbbRedux <- data.table::fread(msbbReduxObj@filePath)

  msbbRedux <- dplyr::left_join(msbbRedux,msbbRealAge,by='individualIdentifier')
  msbb <- dplyr::left_join(msbb,msbbRedux,by=c('SampleID' = 'sampleIdentifier'))
  dups <- duplicated(msbb$SampleID)
  msbb <- msbb[!dups,]

  msbb$brainRegion <- sapply(msbb$BrainRegion.Diagnosis,function(x) strsplit(x,'\\.')[[1]][1])
  msbb$Diagnosis <- sapply(msbb$BrainRegion.Diagnosis,function(x) strsplit(x,'\\.')[[1]][2])
  msbb2 <- dplyr::group_by(msbb,brainRegion)
  msbb3 <- dplyr::summarise(msbb2,
                            NSample=length(unique(SampleID)),
                            NDonor=length(unique(individualIdentifier.y)),
                            meanAOD = mean(AOD.x),
                            sdAOD = sd(AOD.x),
                            Nfemales = sum(SEX=='F'),
                            Nmales = sum(SEX=='M'),
                            NAD = sum(Diagnosis == 'AD'),
                            NControl = sum(Diagnosis == 'CONTROL'),
                            meanPMI = mean(PMI),
                            sdPMI = sd(PMI),
                            meanRIN = mean(RIN.x),
                            sdRIN = sd(RIN.x),
                            NADfemales = sum((SEX == 'F')&(Diagnosis == 'AD')),
                            NControlfemales = sum((SEX == 'F')&(Diagnosis == 'CONTROL')),
                            NADmales = sum((SEX == 'M')&(Diagnosis == 'AD')),
                            NADcontrol = sum((SEX == 'M')&(Diagnosis == 'CONTROL')))

}
