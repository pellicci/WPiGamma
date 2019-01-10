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
    input = cms.untracked.int32(100)
)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register('runningOnData',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "PU config flag")
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
   inputFiles = "root://cms-xrd-global.cern.ch//store/data/Run2016B/SingleMuon/MINIAOD/17Jul2018_ver2-v1/90000/FEADEB19-1D92-E811-BAFA-0025905C54D8.root"
else:
   process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mcRun2_asymptotic_v3')
   print "MC Sample will be taken as input for check up of the code working "
   inputFiles = {"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/92962341-B2BE-E611-9D57-0025905B8574.root"}#,
#"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/120000/FE655795-AABF-E811-9522-44A84225D36F.root"}
                 # "root://cms-xrd-global.cern.ch//store/user/pellicci/WPiGamma_GENSIM_80XV1/WPiGamma_MINIAODSIM_80XV1/161214_125251/0000/WPiGamma_pythia8_MINIAOD_1.root",
                 # "root://cms-xrd-global.cern.ch//store/user/pellicci/WPiGamma_GENSIM_80XV1/WPiGamma_MINIAODSIM_80XV1/161214_125251/0000/WPiGamma_pythia8_MINIAOD_2.root",
                 # "root://cms-xrd-global.cern.ch//store/user/pellicci/WPiGamma_GENSIM_80XV1/WPiGamma_MINIAODSIM_80XV1/161214_125251/0000/WPiGamma_pythia8_MINIAOD_3.root",
                 # "root://cms-xrd-global.cern.ch//store/user/pellicci/WPiGamma_GENSIM_80XV1/WPiGamma_MINIAODSIM_80XV1/161214_125251/0000/WPiGamma_pythia8_MINIAOD_4.root",
                 # "root://cms-xrd-global.cern.ch//store/user/pellicci/WPiGamma_GENSIM_80XV1/WPiGamma_MINIAODSIM_80XV1/161214_125251/0000/WPiGamma_pythia8_MINIAOD_5.root",
                 # "root://cms-xrd-global.cern.ch//store/user/pellicci/WPiGamma_GENSIM_80XV1/WPiGamma_MINIAODSIM_80XV1/161214_125251/0000/WPiGamma_pythia8_MINIAOD_150.root"}

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring (inputFiles)
)

# Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string("WPiGammaAnalysis_output.root")
)

#######------------#######

# #EGAMMA ID
# #Upload the ele/gamma ID information
# from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# dataFormat = DataFormat.MiniAOD
# switchOnVIDElectronIdProducer(process, dataFormat)
# switchOnVIDPhotonIdProducer(process, dataFormat)

# # define which IDs we want to produce
# my_id_modules_el_loose = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff']
# my_id_modules_el       = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff']
# my_id_modules_ph       = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff']

# #add them to the VID producer
# for idmod in my_id_modules_el_loose:
#    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
# for idmod in my_id_modules_el:
#    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
# for idmod in my_id_modules_ph:
#    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

#######------------#######

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
process.WPiGammaAnalysis.runningOnData = options.runningOnData
if not options.runningOnData:    #Only for MC, data reprocessing already has corrections!
   process.WPiGammaAnalysis.slimmedJets = cms.InputTag("updatedPatJetsUpdatedJEC")

#Add the trigger request
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt
process.trigger_filter = hlt.triggerResultsFilter.clone()
process.trigger_filter.triggerConditions = cms.vstring('HLT_IsoMu24_v*', 'HLT_IsoTkMu24_v*', 'HLT_Mu50_v*' , 'HLT_TkMu50_v*', 'HLT_Ele25_eta2p1_WPTight_Gsf_v*', 'HLT_Ele27_WPTight_Gsf_v*')   #paths for 2016 samples
process.trigger_filter.hltResults = cms.InputTag("TriggerResults", "", "HLT")
process.trigger_filter.l1tResults = cms.InputTag("")
process.trigger_filter.throw = cms.bool( False )


# if options.runningOnData:
#    process.seq = cms.Path(process.trigger_filter * process.egmGsfElectronIDSequence * process.egmPhotonIDSequence * process.WPiGammaAnalysis)
# else:
#    process.seq = cms.Path(process.trigger_filter * process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC * process.egmGsfElectronIDSequence * process.egmPhotonIDSequence * process.WPiGammaAnalysis)

if options.runningOnData:
   process.seq = cms.Path(process.trigger_filter * process.egammaPostRecoSeq * process.WPiGammaAnalysis)
else:
   process.seq = cms.Path(process.trigger_filter * process.egammaPostRecoSeq * process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC * process.WPiGammaAnalysis)

process.schedule = cms.Schedule(process.seq)
