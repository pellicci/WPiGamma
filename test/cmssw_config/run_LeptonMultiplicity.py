import FWCore.ParameterSet.Config as cms
process = cms.Process("USER")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Services_cff') 
from Configuration.AlCa.GlobalTag import GlobalTag

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register('runningOnData',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "PU config flag")
options.register('runningOnMuons',
                 True, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "muon trigger config flag")
options.parseArguments()

################################################################################################################
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq

#if not options.runningOnData: #Add PostReco corrections for MC
setupEgammaPostRecoSeq(process,
                       runVID=True,
                       runEnergyCorrections=False, #corrections by default are fine so no need to re-run
                       era='2016-Legacy')

################################################################################################################

#Input source
if options.runningOnData: 
   process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v10') #which conditions to use
   print "Data Sample will be taken as input for check up of the code working "
   inputFiles = "root://cms-xrd-global.cern.ch//store/data/Run2016C/BTagCSV/MINIAOD/23Sep2016-v1/70000/00B6A0EA-A783-E611-8973-02163E0165C4.root"
else:
   process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mcRun2_asymptotic_v3')
   print "MC Sample will be taken as input for check up of the code working "
   inputFiles = {"root://cms-xrd-global.cern.ch//store/user/rselvati/WMinusPiGamma_GENSIM_80XV1/WMinusPiGamma_MINIAODSIM_94XV3/190129_151107/0000/WPiGamma_pythia8_MINIAOD_9.root"}

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring (inputFiles)
)

# Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string("LeptonMultiplicity_output.root")
)

#Use the latest JECs
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
if not options.runningOnData:      #This loop is for MC
   print "Jet Energy Corrections on Monte Carlo will be applied "
   jetCorrectionsList = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')
else:
   print "Jet Energy Corrections on Data will be applied "
   jetCorrectionsList = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')# Added 'L2L3Residual'for data!
   print "Data reprocessing is already using the latest JEC"

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = jetCorrectionsList
)

process.load("StandardModel.WPiGamma.LeptonMultiplicity_cfi")
process.LeptonMultiplicity.runningOnData = options.runningOnData
if not options.runningOnData:    #Only for MC, data reprocessing already has corrections!
   process.LeptonMultiplicity.slimmedJets = cms.InputTag("updatedPatJetsUpdatedJEC")

#Add the trigger request
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt
process.trigger_filter = hlt.triggerResultsFilter.clone()
process.trigger_filter.triggerConditions = cms.vstring('HLT_IsoMu24_v*', 'HLT_IsoTkMu24_v*', 'HLT_Mu50_v*', 'HLT_Ele25_eta2p1_WPTight_Gsf_v*', 'HLT_Ele27_WPTight_Gsf_v7')
process.trigger_filter.hltResults = cms.InputTag("TriggerResults", "", "HLT")
process.trigger_filter.l1tResults = cms.InputTag("")
process.trigger_filter.throw = cms.bool( False )

if options.runningOnData and options.runningOnMuons:
   process.seq = cms.Path(process.trigger_filter_mu * (~process.trigger_filter_ele) * process.egammaPostRecoSeq * process.LeptonMultiplicity) #Excluding events when both muon and electron triggers were lit
if options.runningOnData and not options.runningOnMuons:
   process.seq = cms.Path(process.trigger_filter_ele * process.egammaPostRecoSeq * process.LeptonMultiplicity)
if not options.runningOnData:
   process.seq = cms.Path(process.trigger_filter * process.egammaPostRecoSeq * process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC * process.LeptonMultiplicity)
