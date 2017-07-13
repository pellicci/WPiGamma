from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects/samples_Medium/'

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/run_WPiGammaAnalysis.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['MCpileUp_25ns_Recent2016.root','pileUpHistogramFromjson_Nominal.root' ] #data files for PileUp reweighting
config.JobType.outputFiles = ['WPiGammaAnalysis_output.root']
config.JobType.pyCfgParams = ['runningOnData=False']

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

    config.General.requestName = 'WPiGammaAnalysis_ttbarWlnu'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v3/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'WPiGammaAnalysis_ttbarZQQ'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'WPiGammaAnalysis_ttbarZlnu'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
    config.General.requestName = 'WPiGammaAnalysis_SingleTop_tW'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    p.join()

    config.General.requestName = 'WPiGammaAnalysis_SingleAntiTop_tW'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'WPiGammaAnalysis_WJetsToLNu'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'WPiGammaAnalysis_DY_10_50'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'WPiGammaAnalysis_DY_50'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM'
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

    config.General.requestName = 'WPiGammaAnalysis_ZZ'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'WPiGammaAnalysis_WW'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/WWTo4Q_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'WPiGammaAnalysis_WZ'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
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



    
