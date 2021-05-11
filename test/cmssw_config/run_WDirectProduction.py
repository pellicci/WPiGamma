import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
            'root://cms-xrd-global.cern.ch//store/user/rselvati/WMinusPiGamma_DirectProduction_102X_2018_v1/WMinusPiGamma_DirectProduction_GENSIM_102X_2018_v1/200603_205320/0000/WMinusPiGamma_pythia8_DirectProduction_GENSIM_2018_102.root'
                )
                            )

#Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string("WDirectProduction_output.root")
)

process.WDirectProduction = cms.EDAnalyzer('WDirectProduction'
                              )

process.p = cms.Path(process.WDirectProduction)
