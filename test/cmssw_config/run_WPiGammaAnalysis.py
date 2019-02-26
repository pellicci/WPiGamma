import FWCore.ParameterSet.Config as cms
process = cms.Process("USER")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff') 
process.load('Geometry.CaloEventSetup.CaloTowerConstituents_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff') #Could be not the coprrect one, but should contain the one without "condDBv2"
from Configuration.AlCa.GlobalTag import GlobalTag

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register('runningOnData',
                 True, #default value
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
   inputFiles = "root://cms-xrd-global.cern.ch//store/data/Run2016E/SingleMuon/MINIAOD/17Jul2018-v1/20000/3A57508D-348B-E811-8F39-008CFA197E0C.root"#"root://cms-xrd-global.cern.ch//store/data/Run2016B/SingleMuon/MINIAOD/17Jul2018_ver2-v1/00000/1EB7C205-0E8D-E811-ABD0-0CC47AC52D7A.root"
#"root://cms-xrd-global.cern.ch//store/data/Run2016E/SingleMuon/MINIAOD/17Jul2018-v1/00000/022C9A20-A78A-E811-A986-008CFA111210.root"#"root://cms-xrd-global.cern.ch//store/data/Run2016B/SingleMuon/MINIAOD/17Jul2018_ver2-v1/90000/FEADEB19-1D92-E811-BAFA-0025905C54D8.root"
else:
   process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mcRun2_asymptotic_v3')
   print "MC Sample will be taken as input for check up of the code working "
   inputFiles = {#"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/92962341-B2BE-E611-9D57-0025905B8574.root"}#,
#"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/80000/FE7324D3-87DA-E611-8CB1-0CC47A4C8E5E.root"
"root://cms-xrd-global.cern.ch//store/user/rselvati/WMinusPiGamma_GENSIM_80XV1/WMinusPiGamma_MINIAODSIM_94XV3/190129_151107/0000/WPiGamma_pythia8_MINIAOD_9.root"}


process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring (inputFiles)
)

# Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string("WPiGammaAnalysis_output.root")
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

process.load("StandardModel.WPiGamma.WPiGammaAnalysis_cfi")
process.WPiGammaAnalysis.runningOnData  = options.runningOnData # Load the boolean into the python _cfg file, so that the plugins can read it
#process.WPiGammaAnalysis.runningOnMuons = options.runningOnMuons

if not options.runningOnData:    #Only for MC, data reprocessing already has corrections!
   process.WPiGammaAnalysis.slimmedJets = cms.InputTag("updatedPatJetsUpdatedJEC")

#Add the trigger request
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt

# Trigger filters for muons
process.trigger_filter_mu = hlt.triggerResultsFilter.clone()
process.trigger_filter_mu.triggerConditions = cms.vstring('HLT_IsoMu24_v* OR HLT_IsoTkMu24_v* OR HLT_Mu50_v*') #paths for 2016 data samples SingleMu
process.trigger_filter_mu.hltResults = cms.InputTag("TriggerResults", "", "HLT")
process.trigger_filter_mu.l1tResults = cms.InputTag("")
process.trigger_filter_mu.throw = cms.bool( False )

# Trigger filters for electrons
process.trigger_filter_ele = hlt.triggerResultsFilter.clone()
process.trigger_filter_ele.triggerConditions = cms.vstring('HLT_Ele25_eta2p1_WPTight_Gsf_v* OR HLT_Ele27_WPTight_Gsf_v*') #paths for 2016 data samples SingleEle
process.trigger_filter_ele.hltResults = cms.InputTag("TriggerResults", "", "HLT")
process.trigger_filter_ele.l1tResults = cms.InputTag("")
process.trigger_filter_ele.throw = cms.bool( False )

# Trigger filters for MC samples
process.trigger_filter = hlt.triggerResultsFilter.clone()
process.trigger_filter.triggerConditions = cms.vstring('HLT_IsoMu24_v* OR HLT_IsoTkMu24_v* OR HLT_Mu50_v* OR HLT_Ele25_eta2p1_WPTight_Gsf_v* OR HLT_Ele27_WPTight_Gsf_v*') #paths for 2016 MC
process.trigger_filter.hltResults = cms.InputTag("TriggerResults", "", "HLT")
process.trigger_filter.l1tResults = cms.InputTag("")
process.trigger_filter.throw = cms.bool( False )


if options.runningOnData and options.runningOnMuons:
   process.seq = cms.Path(process.trigger_filter_mu * (~process.trigger_filter_ele) * process.egammaPostRecoSeq * process.WPiGammaAnalysis) #Excluding events when both muon and electron triggers were lit
if options.runningOnData and not options.runningOnMuons:
   process.seq = cms.Path(process.trigger_filter_ele * process.egammaPostRecoSeq * process.WPiGammaAnalysis)
if not options.runningOnData:
   process.seq = cms.Path(process.trigger_filter * process.egammaPostRecoSeq * process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC * process.WPiGammaAnalysis)

process.schedule = cms.Schedule(process.seq)
