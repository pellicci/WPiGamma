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
    input = cms.untracked.int32(3000)
)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register('runningOnData',
                 True, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "PU config flag")
options.register('runningOnMuons',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "muon trigger config flag")
options.register('runningOn2017',
                 True, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "2016-2017 config flag")
options.parseArguments()

################################################################################################################
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq

if options.runningOn2017:
   setupEgammaPostRecoSeq(process,
                          runVID=True, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                          era='2017-Nov17ReReco')  #Re-run energy corrections for electrons, to fix scale & smearing bug. Default re-run value is True

if not options.runningOn2017: #Add PostReco corrections for MC
   setupEgammaPostRecoSeq(process,
                          runVID=True,
                          runEnergyCorrections=False, #corrections by default are fine so no need to re-run
                          era='2016-Legacy')

################################################################################################################

#Input source
if options.runningOnData: 
   if options.runningOn2017: #test 2017 data
      process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v11') #which conditions to use (Mario di usare, a parita' di nome della tag, quella con la versione piu' alta. La versione piu' alta e' quella consigliata per le JEC
      print "Data Sample (2017) will be taken as input for check up of the code working "
      inputFiles = {"root://cms-xrd-global.cern.ch//store/data/Run2017F/SingleElectron/MINIAOD/31Mar2018-v1/90002/EC099452-C938-E811-9922-0CC47A7C354C.root","root://cms-xrd-global.cern.ch//store/data/Run2017F/SingleElectron/MINIAOD/31Mar2018-v1/90002/B09F43EA-CD38-E811-9D77-0CC47A4C8EE8.root","root://cms-xrd-global.cern.ch//store/data/Run2017F/SingleElectron/MINIAOD/31Mar2018-v1/90002/A421D61D-B137-E811-9DE7-0CC47A4C8F06.root","root://cms-xrd-global.cern.ch//store/data/Run2017F/SingleElectron/MINIAOD/31Mar2018-v1/90002/A05A43B6-CA38-E811-B43A-0CC47A4D7644.root"}
   else: #test 2016 data
      process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v10') #which conditions to use
      print "Data Sample (2016) will be taken as input for check up of the code working "
      inputFiles = {"root://cms-xrd-global.cern.ch//store/data/Run2016E/SingleMuon/MINIAOD/17Jul2018-v1/20000/3A57508D-348B-E811-8F39-008CFA197E0C.root"}

else:
   if options.runningOn2017: #test 2017 MC
      process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v17') 
      print "MC Sample (2017) will be taken as input for check up of the code working "
      inputFiles = {"root://cms-xrd-global.cern.ch//store/user/rselvati/WMinusPiGamma_GENSIM_80XV1/WMinusPiGamma_MINIAODSIM_94XV3/190129_151107/0000/WPiGamma_pythia8_MINIAOD_9.root"}
   else: #test 2016 MC
      process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mcRun2_asymptotic_v3')
      print "MC Sample (2016) will be taken as input for check up of the code working "
      inputFiles = {"root://cms-xrd-global.cern.ch//store/user/rselvati/WMinusPiGamma_GENSIM_80XV1/WMinusPiGamma_MINIAODSIM_94XV3/190129_151107/0000/WPiGamma_pythia8_MINIAOD_9.root"}


process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring (inputFiles)
)

# Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string("WPiGammaAnalysis_output.root")
)


###############################################################################################################################
#                                                                                                                             #
#-------------------------------------------------------- Jet corrections ----------------------------------------------------#
#                                                                                                                             #
#                                      https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC                                     #
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections?redirectedfrom=CMS.WorkBookJetEnergyCorrections #
#                                                                                                                             #
###############################################################################################################################

#Use the latest JECs
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
if not options.runningOnData:      #This loop is for MC
   print "Jet Energy Corrections on Monte Carlo will be applied"
   jetCorrectionsList = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')# The right corrections are taken through the chosen global tag
else:
   print "Jet Energy Corrections on Data will be applied "
   jetCorrectionsList = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')# The right corrections are taken through the chosen global tag

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = jetCorrectionsList
)

process.load("StandardModel.WPiGamma.WPiGammaAnalysis_cfi")
process.WPiGammaAnalysis.runningOnData  = options.runningOnData
process.WPiGammaAnalysis.runningOn2017  = options.runningOn2017

# Apply JEC to both MC and data
process.WPiGammaAnalysis.slimmedJets = cms.InputTag("updatedPatJetsUpdatedJEC") # This tells the WPiGammaAnalysis python module (python/WPiGammaAnalysis.py) that the python object called slimmedJets should point to a list called updatedPatJetsUpdatedJEC, and not anymore to slimmedJets



###############################################
#                                             #
#------------------ Trigger ------------------#
#                                             #
###############################################


#Add the trigger request
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt

# Trigger filters for muons
process.trigger_filter_data_mu = hlt.triggerResultsFilter.clone()
if options.runningOn2017:
   process.trigger_filter_data_mu.triggerConditions = cms.vstring('HLT_IsoMu27_v* OR HLT_Mu50_v*') #paths for 2017 data samples SingleMu
else:
   process.trigger_filter_data_mu.triggerConditions = cms.vstring('HLT_IsoMu24_v* OR HLT_IsoTkMu24_v* OR HLT_Mu50_v*') #paths for 2016 data samples SingleMu
process.trigger_filter_data_mu.hltResults = cms.InputTag("TriggerResults", "", "HLT")
process.trigger_filter_data_mu.l1tResults = cms.InputTag("")
process.trigger_filter_data_mu.throw = cms.bool( False )

# Trigger filters for electrons
process.trigger_filter_data_ele = hlt.triggerResultsFilter.clone()

if options.runningOn2017:
   process.trigger_filter_data_ele.triggerConditions = cms.vstring('HLT_Ele32_WPTight_Gsf_L1DoubleEG_v* OR HLT_Ele32_WPTight_Gsf_v*') #paths for 2017 data samples SingleEle
else:
   process.trigger_filter_data_ele.triggerConditions = cms.vstring('HLT_Ele25_eta2p1_WPTight_Gsf_v* OR HLT_Ele27_WPTight_Gsf_v*') #paths for 2016 data samples SingleEle

process.trigger_filter_data_ele.hltResults = cms.InputTag("TriggerResults", "", "HLT")
process.trigger_filter_data_ele.l1tResults = cms.InputTag("")
process.trigger_filter_data_ele.throw = cms.bool( False )

# Trigger filters for MC samples
process.trigger_filter_MC = hlt.triggerResultsFilter.clone()
if options.runningOn2017:
   process.trigger_filter_MC.triggerConditions = cms.vstring('HLT_IsoMu27_v* OR HLT_Mu50_v* OR HLT_Ele32_WPTight_Gsf_L1DoubleEG_v*') #paths for 2017 MC
else:
   process.trigger_filter_MC.triggerConditions = cms.vstring('HLT_IsoMu24_v* OR HLT_IsoTkMu24_v* OR HLT_Mu50_v* OR HLT_Ele25_eta2p1_WPTight_Gsf_v* OR HLT_Ele27_WPTight_Gsf_v*') #paths for 2016 MC
process.trigger_filter_MC.hltResults = cms.InputTag("TriggerResults", "", "HLT")
process.trigger_filter_MC.l1tResults = cms.InputTag("")
process.trigger_filter_MC.throw = cms.bool( False )


###############################################
#                                             #
#----------------- Sequence ------------------#
#                                             #
###############################################

if options.runningOnData and options.runningOnMuons and options.runningOn2017: # Data, Muons 2017
   process.seq = cms.Path(process.trigger_filter_data_mu  * process.egammaPostRecoSeq * process.WPiGammaAnalysis) #Excluding events when both muon and electron triggers were lit

if options.runningOnData and not options.runningOnMuons and options.runningOn2017: # Data, Electrons 2017
   process.seq = cms.Path(process.trigger_filter_data_ele * (~process.trigger_filter_data_mu) * process.egammaPostRecoSeq * process.WPiGammaAnalysis)

if options.runningOnData and options.runningOnMuons and not options.runningOn2017: # Data, Muons 2016
   process.seq = cms.Path(process.trigger_filter_data_mu  * process.egammaPostRecoSeq * process.WPiGammaAnalysis) #Excluding events when both muon and electron triggers were lit

if options.runningOnData and not options.runningOnMuons and not options.runningOn2017: # Data, Electrons 2016
   process.seq = cms.Path(process.trigger_filter_data_ele * (~process.trigger_filter_data_mu) * process.egammaPostRecoSeq * process.WPiGammaAnalysis)

if not options.runningOnData and options.runningOn2017: # MC 2017
   process.seq = cms.Path(process.trigger_filter_MC * process.egammaPostRecoSeq * process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC * process.WPiGammaAnalysis)

if not options.runningOnData and not options.runningOn2017: # MC 2016
   process.seq = cms.Path(process.trigger_filter_MC * process.egammaPostRecoSeq * process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC * process.WPiGammaAnalysis)

process.schedule = cms.Schedule(process.seq)
