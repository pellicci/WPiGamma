import FWCore.ParameterSet.Config as cms

WPiGammaAnalysis = cms.EDAnalyzer('WPiGammaAnalysis',
                               packedPFCandidates = cms.InputTag("packedPFCandidates"),
                               slimmedMuons = cms.InputTag("slimmedMuons"),
                               prunedGenParticles = cms.InputTag("prunedGenParticles"),
                               slimmedPhotons = cms.InputTag("slimmedPhotons"),
                               runningOnData = cms.bool(False),

                               pvCollection = cms.InputTag("offlineSlimmedPrimaryVertices"), #New Stuff 
                               bsCollection = cms.InputTag("offlineBeamSpot"),
                               PileupSrc = cms.InputTag("slimmedAddPileupInfo") # ,
)
