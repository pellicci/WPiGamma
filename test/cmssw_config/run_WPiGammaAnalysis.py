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
    input = cms.untracked.int32(20000)
)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register('runningOnData',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "PU config flag")
options.register('runningOnMuons',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "muon trigger config flag")
options.register('runningEra',
                 1, #default value. 0 is 2016, 1 is 2017, 2 is 2018
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "2016-2017-2018 config flag")
options.parseArguments()

################################################################################################################
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq

if options.runningEra == 0: #Add PostReco corrections for MC
   setupEgammaPostRecoSeq(process,
                          runVID=True,
                          runEnergyCorrections=True, #corrections by default are fine so no need to re-run
                          era='2016-Legacy')

   if not options.runningOnData:
      #ECAL prefiring corrections https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe
      process.prefiringweight = cms.EDProducer("L1ECALPrefiringWeightProducer",
                                               ThePhotons = cms.InputTag("slimmedPhotons"),
                                               TheJets = cms.InputTag("slimmedJets"),
                                               L1Maps = cms.string("L1PrefiringMaps_new.root"), # update this line with the location of this file
                                               DataEra = cms.string("2016BtoH"), #2016
                                               UseJetEMPt = cms.bool(False), #can be set to true to use jet prefiring maps parametrized vs pt(em) instead of pt
                                               PrefiringRateSystematicUncty = cms.double(0.2) #Minimum relative prefiring uncty per object
                                            )

if options.runningEra == 1: #running on 2017
   setupEgammaPostRecoSeq(process,
                          runVID=True, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                          era='2017-Nov17ReReco')  #Re-run energy corrections for electrons, to fix scale & smearing bug. Default re-run value is True

   if not options.runningOnData:
      #ECAL prefiring corrections https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe
      process.prefiringweight = cms.EDProducer("L1ECALPrefiringWeightProducer",
                                               ThePhotons = cms.InputTag("slimmedPhotons"),
                                               TheJets = cms.InputTag("slimmedJets"),
                                               L1Maps = cms.string("L1PrefiringMaps_new.root"), # update this line with the location of this file
                                               DataEra = cms.string("2017BtoF"), #2017
                                               UseJetEMPt = cms.bool(False), #can be set to true to use jet prefiring maps parametrized vs pt(em) instead of pt
                                               PrefiringRateSystematicUncty = cms.double(0.2) #Minimum relative prefiring uncty per object
                                            )

################################################################################################################

#Input source
if options.runningOnData:

   if options.runningEra == 0: #test 2016 data
      process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v10') #which conditions to use
      print "Data Sample (2016) will be taken as input for check up of the code working "
      inputFiles = {#"root://cms-xrd-global.cern.ch//store/data/Run2016E/SingleMuon/MINIAOD/17Jul2018-v1/20000/3A57508D-348B-E811-8F39-008CFA197E0C.root"}
         "root://cms-xrd-global.cern.ch//store/data/Run2016B/SingleElectron/MINIAOD/17Jul2018_ver2-v1/80000/FEAB7D8A-048C-E811-A513-AC1F6B1AEFFC.root"}

   if options.runningEra == 1: #test 2017 data
      process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v11') #which conditions to use (Mario dice di usare, a parita' di nome della tag, quella con la versione piu' alta. La versione piu' alta e' quella consigliata per le JEC)
      print "Data Sample (2017) will be taken as input for check up of the code working "
      inputFiles = {"root://cms-xrd-global.cern.ch//store/data/Run2017F/SingleElectron/MINIAOD/31Mar2018-v1/90002/EC099452-C938-E811-9922-0CC47A7C354C.root","root://cms-xrd-global.cern.ch//store/data/Run2017F/SingleElectron/MINIAOD/31Mar2018-v1/90002/B09F43EA-CD38-E811-9D77-0CC47A4C8EE8.root","root://cms-xrd-global.cern.ch//store/data/Run2017F/SingleElectron/MINIAOD/31Mar2018-v1/90002/A421D61D-B137-E811-9DE7-0CC47A4C8F06.root","root://cms-xrd-global.cern.ch//store/data/Run2017F/SingleElectron/MINIAOD/31Mar2018-v1/90002/A05A43B6-CA38-E811-B43A-0CC47A4D7644.root"}
         #"/store/data/Run2017B/SingleMuon/MINIAOD/31Mar2018-v1/100000/481AD1C9-D538-E811-9B11-0025905B8600.root"}

else:

   if options.runningEra == 0: #test 2016 MC
      process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mcRun2_asymptotic_v3')
      print "MC Sample (2016) will be taken as input for check up of the code working "
      inputFiles = {"root://cms-xrd-global.cern.ch//store/user/rselvati/WMinusPiGamma_GENSIM_80XV1/WMinusPiGamma_MINIAODSIM_94XV3/190129_151107/0000/WPiGamma_pythia8_MINIAOD_9.root","root://cms-xrd-global.cern.ch//store/user/rselvati/WPlusPiGamma_GENSIM_80XV1/WPlusPiGamma_MINIAODSIM_94XV3/190129_150618/0000/WPiGamma_pythia8_MINIAOD_80.root"}
         #"/store/mc/RunIISummer16MiniAODv2/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/CE408A70-56C1-E611-8187-FA163EDC366E.root"}

   if options.runningEra == 1: #test 2017 MC
      process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v17') 
      print "MC Sample (2017) will be taken as input for check up of the code working "
      inputFiles = {#"root://cms-xrd-global.cern.ch//store/user/rselvati/WMinusPiGamma_GENSIM_80XV1/WMinusPiGamma_MINIAODSIM_94XV3/190129_151107/0000/WPiGamma_pythia8_MINIAOD_9.root"}
         #"/store/mc/RunIIFall17MiniAODv2/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/60000/FCC2AFA9-4BBB-E811-B35F-0CC47AFB7D48.root"}
         "/store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/20000/10928512-C04D-E811-AD45-F01FAFE5E4F8.root"}


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
process.WPiGammaAnalysis.runningOnData = options.runningOnData
process.WPiGammaAnalysis.runningEra    = options.runningEra

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

if options.runningEra == 0:
   process.trigger_filter_data_mu.triggerConditions = cms.vstring('HLT_IsoMu24_v* OR HLT_IsoTkMu24_v* OR HLT_Mu50_v*') #paths for 2016 data samples SingleMu

if options.runningEra == 1:
   process.trigger_filter_data_mu.triggerConditions = cms.vstring('HLT_IsoMu27_v* OR HLT_Mu50_v*') #paths for 2017 data samples SingleMu. IsoMu27 is the lowest pT unprescaled for 2017

process.trigger_filter_data_mu.hltResults = cms.InputTag("TriggerResults", "", "HLT")
process.trigger_filter_data_mu.l1tResults = cms.InputTag("")
process.trigger_filter_data_mu.throw = cms.bool( False )

# Trigger filters for electrons
process.trigger_filter_data_ele = hlt.triggerResultsFilter.clone()

if options.runningEra == 0:
   process.trigger_filter_data_ele.triggerConditions = cms.vstring('HLT_Ele25_eta2p1_WPTight_Gsf_v* OR HLT_Ele27_WPTight_Gsf_v*') #paths for 2016 data samples SingleEle

if options.runningEra == 1:
   process.trigger_filter_data_ele.triggerConditions = cms.vstring('HLT_Ele32_WPTight_Gsf_L1DoubleEG_v* OR HLT_Ele32_WPTight_Gsf_v*') #paths for 2017 data samples SingleEle. The trigger WITHOUT L1DoubleEG is active only for around 27/fb. When both the triggers are on, L1DoubleEG is on ~2% times more than simple WPTight

process.trigger_filter_data_ele.hltResults = cms.InputTag("TriggerResults", "", "HLT")
process.trigger_filter_data_ele.l1tResults = cms.InputTag("")
process.trigger_filter_data_ele.throw = cms.bool( False )


# Trigger filters for MC samples
process.trigger_filter_MC = hlt.triggerResultsFilter.clone()

if options.runningEra == 0:
   process.trigger_filter_MC.triggerConditions = cms.vstring('HLT_IsoMu24_v* OR HLT_IsoTkMu24_v* OR HLT_Mu50_v* OR HLT_Ele25_eta2p1_WPTight_Gsf_v* OR HLT_Ele27_WPTight_Gsf_v*') #paths for 2016 MC

if options.runningEra == 1:
   process.trigger_filter_MC.triggerConditions = cms.vstring('HLT_IsoMu27_v* OR HLT_Mu50_v* OR HLT_Ele32_WPTight_Gsf_L1DoubleEG_v* OR HLT_Ele32_WPTight_Gsf_v*') #paths for 2017 MC

process.trigger_filter_MC.hltResults = cms.InputTag("TriggerResults", "", "HLT")
process.trigger_filter_MC.l1tResults = cms.InputTag("")
process.trigger_filter_MC.throw = cms.bool( False )


###############################################
#                                             #
#----------------- Sequence ------------------#
#                                             #
###############################################

if options.runningOnData and options.runningOnMuons: # Data, Muons
   process.seq = cms.Path(process.trigger_filter_data_mu  * process.egammaPostRecoSeq * process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC * process.WPiGammaAnalysis) 

if options.runningOnData and not options.runningOnMuons: # Data, Electrons
   process.seq = cms.Path(process.trigger_filter_data_ele * (~process.trigger_filter_data_mu) * process.egammaPostRecoSeq * process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC * process.WPiGammaAnalysis) #Excluding events when both muon and electron triggers were lit

if not options.runningOnData: # MC
   process.seq = cms.Path(process.trigger_filter_MC * process.egammaPostRecoSeq * process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC * process.prefiringweight * process.WPiGammaAnalysis)

process.schedule = cms.Schedule(process.seq)
