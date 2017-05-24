import FWCore.ParameterSet.Config as cms
process = cms.Process("USER")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Services_cff') 
from Configuration.AlCa.GlobalTag import GlobalTag

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register('runningOnData',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "PU config flag")
options.parseArguments()

#Input source
if options.runningOnData: 
   process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v14') #which conditions to use
   print "Data Sample will be taken as input for check up of the code working "
   inputFiles="root://cms-xrd-global.cern.ch//store/data/Run2016C/BTagCSV/MINIAOD/23Sep2016-v1/70000/00B6A0EA-A783-E611-8973-02163E0165C4.root"
else:
   process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')
   print "MC Sample will be taken as input for check up of the code working "
   inputFiles='root://cms-xrd-global.cern.ch//store/user/pellicci/WPiGamma_GENSIM_80XV1/WPiGamma_MINIAODSIM_80XV1/161214_125251/0000/WPiGamma_pythia8_MINIAOD_1.root'

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring (inputFiles)
)

# Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string("WPiGammaAnalysis_output.root")
)

process.load("StandardModel.WPiGamma.WPiGammaAnalysis_cfi")
process.WPiGammaAnalysis.runningOnData = options.runningOnData

process.seq = cms.Path(process.WPiGammaAnalysis)

process.schedule = cms.Schedule(process.seq)
