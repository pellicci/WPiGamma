import FWCore.ParameterSet.Config as cms

WPiGammaAnalysis = cms.EDAnalyzer('WPiGammaAnalysis',
                                  packedPFCandidates = cms.InputTag("packedPFCandidates"),
                                  slimmedMuons       = cms.InputTag("slimmedMuons"),
                                  prunedGenParticles = cms.InputTag("prunedGenParticles"),
                                  slimmedPhotons     = cms.InputTag("slimmedPhotons"),
                                  slimmedElectrons   = cms.InputTag("slimmedElectrons"),
                                  slimmedJets        = cms.InputTag("slimmedJets"),
                                  runningOnData      = cms.bool(False),
                                  pvCollection       = cms.InputTag("offlineSlimmedPrimaryVertices"), #New Stuff 
                                  bsCollection       = cms.InputTag("offlineBeamSpot"),
                                  PileupSrc          = cms.InputTag("slimmedAddPileupInfo"),

                                  # ELE ID decisions (common to all formats)
                                  eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90"),
                                  eleTightIdMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80"),
                                  # ValueMaps with MVA results
                                  mvaValuesMap_el     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
                                  mvaCategoriesMap_el = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),

                                  # GAMMA ID decisions (common to all formats)
                                  phoMediumIdBoolMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring16-nonTrig-V1-wp90"),
                                  phoMediumIdFullInfoMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring16-nonTrig-V1-wp90"),
                                  # This is a fairly verbose mode if switched on, with full cut flow 
                                  # diagnostics for each candidate. Use it in a low event count test job.
                                  phoIdVerbose = cms.bool(False),
                                  # ValueMaps with MVA results
                                  mvaValuesMap_ph     = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring16NonTrigV1Values"),
                                  mvaCategoriesMap_ph = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring16NonTrigV1Categories")
)
