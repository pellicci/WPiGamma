import FWCore.ParameterSet.Config as cms

WPiGammaAnalysis = cms.EDAnalyzer('WPiGammaAnalysis',
                                  packedPFCandidates = cms.InputTag("packedPFCandidates"),
                                  slimmedMuons       = cms.InputTag("slimmedMuons"),
                                  prunedGenParticles = cms.InputTag("prunedGenParticles"),
                                  slimmedPhotons     = cms.InputTag("slimmedPhotons"),
                                  slimmedElectrons   = cms.InputTag("slimmedElectrons"),
                                  slimmedJets        = cms.InputTag("slimmedJets"),
                                  slimmedMETs        = cms.InputTag("slimmedMETs"),
                                  slimmedMETsPuppi   = cms.InputTag("slimmedMETsPuppi"),
                                  runningOnData      = cms.bool(False),
                                  runningEra         = cms.int32(0), # One of the possible python types for the C++ type "int" 
                                  pvCollection       = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                  bsCollection       = cms.InputTag("offlineBeamSpot"),
                                  PileupSrc          = cms.InputTag("slimmedAddPileupInfo"),
                                  triggerBits        = cms.InputTag("TriggerResults","","HLT"),

                                  rho                = cms.InputTag("fixedGridRhoFastjetAll"),
                                  #https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
                                  Bjets_WP_loose_2016  = cms.double(0.2217), # Working point for the 2016 Bjet discriminator (loose)
                                  Bjets_WP_loose_2017  = cms.double(0.1522), # Working point for the 2017 Bjet discriminator (loose)
                                  Bjets_WP_loose_2018  = cms.double(0.1241), # Working point for the 2018 Bjet discriminator (loose)
                                  Bjets_WP_medium_2016 = cms.double(0.6321), # Working point for the 2016 Bjet discriminator (medium)
                                  Bjets_WP_medium_2017 = cms.double(0.4941), # Working point for the 2017 Bjet discriminator (medium)
                                  Bjets_WP_medium_2018 = cms.double(0.4184), # Working point for the 2018 Bjet discriminator (medium)

                                  # This is a fairly verbose mode if switched on, with full cut flow 
                                  # diagnostics for each candidate. Use it in a low event count test job.
                                  phoIdVerbose = cms.bool(False),

                                  # Effective areas for computing PU correction for isolations
                                  effAreasConfigFile_el = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt"),#The EA constants go with the ID used (Fall17), not with the release (so it is ok to use 92X if the release is 9_4_9)
                                  effAreasConfigFile_ph = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased.txt")
                                  # effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt")
                                 )
