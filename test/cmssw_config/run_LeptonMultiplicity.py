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
    input = cms.untracked.int32(2000)
)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register('runningOnData',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "PU config flag")
options.register('runningEra',
                 0, #default value. 0 is 2016, 1 is 2017, 2 is 2018
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

if options.runningEra == 1: #running on 2017
   setupEgammaPostRecoSeq(process,
                          runVID=True, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                          era='2017-Nov17ReReco')  #Re-run energy corrections for electrons, to fix scale & smearing bug. Default re-run value is True

if options.runningEra == 2: #running on 2018
    setupEgammaPostRecoSeq(process,
                           era='2018-Prompt')

################################################################################################################

#Input source
if options.runningEra == 0: #test 2016 MC
   process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mcRun2_asymptotic_v3')
   print "MC Sample (2016) will be taken as input for check up of the code working "
   inputFiles = {"root://cms-xrd-global.cern.ch//store/user/rselvati/WPlusPiGamma_GENSIM_80X_v2/WPlusPiGamma_MINIAODSIM_94X_2016_v2/190824_230046/0002/WPiGamma_pythia8_MINIAOD_2797.root"}

if options.runningEra == 1: #test 2017 MC
   process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v17') 
   print "MC Sample (2017) will be taken as input for check up of the code working "
   inputFiles = {"root://cms-xrd-global.cern.ch//store/user/rselvati/WMinusPiGamma_GENSIM_80XV1/WMinusPiGamma_MINIAODSIM_94XV3/190129_151107/0000/WPiGamma_pythia8_MINIAOD_9.root"}
        
if options.runningEra == 2: #test 2018 MC
   process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15') 
   print "MC Sample (2018) will be taken as input for check up of the code working "
   inputFiles = {"root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18MiniAOD/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/00000/38CB48DD-A494-4044-B9CD-64241707E25F.root"}


process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring (inputFiles)
)

# Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string("LeptonMultiplicity_output.root")
)


###############################################################################################################################
#                                                                                                                             #
#-------------------------------------------------------- Jet corrections ----------------------------------------------------#
#                                                                                                                             #
#                                      https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC                                     #
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections?redirectedfrom=CMS.WorkBookJetEnergyCorrections #
#                                                                                                                             #
###############################################################################################################################

process.load("StandardModel.WPiGamma.LeptonMultiplicity_cfi")
process.LeptonMultiplicity.runningOnData = options.runningOnData
process.LeptonMultiplicity.runningEra    = options.runningEra


###############################################
#                                             #
#------------------ Trigger ------------------#
#                                             #
###############################################

#Add the trigger request
import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt

# Trigger filters for MC samples
process.trigger_filter_MC = hlt.triggerResultsFilter.clone()

if options.runningEra == 0:
   process.trigger_filter_MC.triggerConditions = cms.vstring('HLT_IsoMu24_v* OR HLT_IsoTkMu24_v* OR HLT_Mu50_v* OR HLT_Ele25_eta2p1_WPTight_Gsf_v* OR HLT_Ele27_WPTight_Gsf_v*') #paths for 2016 MC

if options.runningEra == 1:
   process.trigger_filter_MC.triggerConditions = cms.vstring('HLT_IsoMu27_v* OR HLT_Mu50_v* OR HLT_Ele32_WPTight_Gsf_L1DoubleEG_v* OR HLT_Ele32_WPTight_Gsf_v*') #paths for 2017 MC

if options.runningEra == 2:
   process.trigger_filter_MC.triggerConditions = cms.vstring('HLT_IsoMu24_v* OR HLT_Mu50_v* OR HLT_Ele32_WPTight_Gsf_v*') #paths for 2018 MC

process.trigger_filter_MC.hltResults = cms.InputTag("TriggerResults", "", "HLT")
process.trigger_filter_MC.l1tResults = cms.InputTag("")
process.trigger_filter_MC.throw = cms.bool( False )


###############################################
#                                             #
#----------------- Sequence ------------------#
#                                             #
###############################################

if not options.runningOnData:
   process.seq = cms.Path(process.trigger_filter_MC * process.egammaPostRecoSeq * process.LeptonMultiplicity)

process.schedule = cms.Schedule(process.seq)
