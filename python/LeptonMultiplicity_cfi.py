import FWCore.ParameterSet.Config as cms

LeptonMultiplicity = cms.EDAnalyzer('LeptonMultiplicity',
                                  packedPFCandidates = cms.InputTag("packedPFCandidates"),
                                  slimmedMuons       = cms.InputTag("slimmedMuons"),
                                  prunedGenParticles = cms.InputTag("prunedGenParticles"),
                                  slimmedPhotons     = cms.InputTag("slimmedPhotons"),
                                  slimmedElectrons   = cms.InputTag("slimmedElectrons"),
                                  slimmedJets        = cms.InputTag("slimmedJets"),
                                  runningOnData      = cms.bool(False),
                                  runningEra         = cms.int32(0), # One of the possible python types for the C++ type "int" 
                                  pvCollection       = cms.InputTag("offlineSlimmedPrimaryVertices"), #New Stuff 
                                  bsCollection       = cms.InputTag("offlineBeamSpot"),
                                  PileupSrc          = cms.InputTag("slimmedAddPileupInfo"),
                                  triggerbits        = cms.InputTag("TriggerResults","","HLT"),
                                  rho                = cms.InputTag("fixedGridRhoFastjetAll"),

                                  # ELE ID decisions (common to all formats)
                                  # eleLooseIdMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-HZZ-V1-wpLoose"),
                                  # eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90"),
                                  # eleTightIdMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80"),
                                  # ValueMaps with MVA results
                                  # mvaValuesMap_el_loose     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Values"),
                                  # mvaCategoriesMap_el_loose = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Categories"),
                                  # mvaValuesMap_el           = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
                                  # mvaCategoriesMap_el       = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),

                                  # GAMMA ID decisions (common to all formats)
                                  # phoMediumIdBoolMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring16-nonTrig-V1-wp90"),
                                  # phoMediumIdFullInfoMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring16-nonTrig-V1-wp90"),
                                  # This is a fairly verbose mode if switched on, with full cut flow 
                                  # diagnostics for each candidate. Use it in a low event count test job.
                                  phoIdVerbose = cms.bool(False),
                                  # ValueMaps with MVA results
                                  # mvaValuesMap_ph     = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring16NonTrigV1Values"),
                                  # mvaCategoriesMap_ph = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring16NonTrigV1Categories"),
                                  # Effective areas for computing PU correction for isolations
                                  effAreasConfigFile_el = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt"),
                                  effAreasConfigFile_ph = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased.txt")
                                 )

