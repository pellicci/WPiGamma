from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = Configuration()

config.section_('General')
config.General.transferOutputs = True

runningEra = 1 # 0 = 2016, 1 = 2017, 2 = 2018

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/run_WPiGammaAnalysis.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['WPiGammaAnalysis_output.root']
config.JobType.pyCfgParams = ['runningOnData=False']


if runningEra == 0:
    config.General.workArea = 'crab_projects/samples_MC_2016/'
    config.JobType.inputFiles = ['PU/MCpileUp_2016_25ns_Moriond17MC_PoissonOOTPU.root','PU/MyDataPileupHistogram_2016.root'] #MC and data files for PileUp reweighting (2016)

if runningEra == 1:
    config.General.workArea = 'crab_projects/samples_MC_2017/'
    config.JobType.inputFiles = ['PU/MCpileUp_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU.root','PU/MyDataPileupHistogram_2017.root'] #MC and data files for PileUp reweighting (2017)


config.section_('Data')
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException
    from multiprocessing import Process

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    #################################################
    #                                               #
    #--------------- Running 2016 MC ---------------#
    #                                               #
    #################################################


    if runningEra == 0:

        config.JobType.pyCfgParams = ['runningEra=0'] # Configure 2016 MC jobs 

        config.General.requestName = 'WPiGammaAnalysis_ttbar'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_ttbarlnu'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_ttbarWQQ'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_ttbarWlnu_1'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v3/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_ttbarWlnu_2' ################### NEW
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_ttbarZQQ'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_ttbarZlnu_1'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_ttbarZlnu_2' ################## NEW
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_ttbarZlnu_3' ################## NEW
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext3-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_SingleTop_tW'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_SingleAntiTop_tW'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_WJetsToLNu'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_DY_10_50_1' ################################ NEW
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_DY_10_50_2'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_DY_50_1'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_DY_50_2' ################################ NEW
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_HT100to200'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_HT200to300_1'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join() 
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_HT200to300_2'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join() 
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_HT300to500_1'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_HT300to500_2'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_HT500to700_1'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_HT500to700_2'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_HT700to1000_1'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_HT700to1000_2'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_HT1000to1500_1'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_HT1000to1500_2'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join() 
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_HT1500to2000_1'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_HT1500to2000_2'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_HT2000toInf_1'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_HT2000toInf_2'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        # config.General.requestName = 'WPiGammaAnalysis_ZZ'
        # config.Data.unitsPerJob = 5
        # config.Data.inputDataset = '/ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        # p = Process(target=submit, args=(config,))
        # p.start()
        # p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_WW'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/WWTo4Q_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_WZ_1'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_WZ_2' ########################### NEW
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_GammaJets_20_40'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISummer16MiniAODv2-PUMoriond17_backup_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_GammaJets_40_Inf'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_GammaJets_20_Inf'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_WGToLNuG_1'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_WGToLNuG_2'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_WGToLNuG_3'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext3-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_DoubleEMEnriched_30to40'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_DoubleEMEnriched_30toInf'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_DoubleEMEnriched_40toInf'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_TTGJets_1'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_TTGJets_2'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_ZGTo2LG_1' ####################### NEW
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_ZGTo2LG_2' ####################### NEW
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()


    #################################################
    #                                               #
    #--------------- Running 2017 MC ---------------#
    #                                               #
    #################################################


    if runningEra == 1:

        config.JobType.pyCfgParams = ['runningEra=1'] # Configure 2017 MC jobs 

        config.General.requestName = 'WPiGammaAnalysis_ttbar_ToHadronic'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = 'WPiGammaAnalysis_ttbar_ToSemiLeptonic'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = 'WPiGammaAnalysis_ttbarlnu'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = 'WPiGammaAnalysis_ttbarWQQ'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = 'WPiGammaAnalysis_ttbarWlnu'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = 'WPiGammaAnalysis_ttbarZQQ'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_ttbarZlnu'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = 'WPiGammaAnalysis_SingleTop_tW'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_SingleAntiTop_tW'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = 'WPiGammaAnalysis_WJetsToLNu_1'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = 'WPiGammaAnalysis_WJetsToLNu_2'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_DY_10_50'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_DY_50_1'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = 'WPiGammaAnalysis_DY_50_2'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = 'WPiGammaAnalysis_QCD_HT100to200'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT100to200_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_HT200to300'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join() 
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_HT300to500'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_HT500to700'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = 'WPiGammaAnalysis_QCD_HT700to1000'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = 'WPiGammaAnalysis_QCD_HT1000to1500'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = 'WPiGammaAnalysis_QCD_HT1000to1500'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = 'WPiGammaAnalysis_QCD_HT1500to2000'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = 'WPiGammaAnalysis_QCD_HT2000toInf'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_WW'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/WWTo4Q_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_WZ' # Non sembra esserci un campione con le giuste condizioni di PU
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/WZ_TuneCP5_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' 
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_GammaJets_20_40' # Non sembra esserci un campione con le giuste condizioni di PU
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_GammaJets_40_Inf' # Non sembra esserci un campione con le giuste condizioni di PU (nonostante ci sia "v2")
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_GammaJets_20_Inf' # Non sembra esserci un campione con le giuste condizioni di PU
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_WGToLNuG' # Nonostante il "v1", ci sono le condizioni di PU corrette
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_DoubleEMEnriched_30to40' # Non sembra esserci un campione con le giuste condizioni di PU
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_DoubleEMEnriched_30toInf' # Non sembra esserci un campione con le giuste condizioni di PU
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_QCD_DoubleEMEnriched_40toInf' # Non sembra esserci un campione con le giuste condizioni di PU
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_TTGJets_1'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_TTGJets_2'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = 'WPiGammaAnalysis_TTGJets_3'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = 'WPiGammaAnalysis_ZGTo2LG'
        config.Data.unitsPerJob = 5
        config.Data.inputDataset = '/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/MINIAODSIM'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()



        


    
