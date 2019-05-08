from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = Configuration()

config.section_('General')
config.General.transferOutputs = True

runningEra = 0 # 0 = 2016, 1 = 2017, 2 = 2018

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/run_WPiGammaAnalysis.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['WPiGammaAnalysis_output.root']
config.section_('Data')


if runningEra == 0:
    config.General.workArea = 'crab_projects/samples_data_2016/'
    #config.JobType.inputFiles = ['PU/MCpileUp_2016_25ns_Moriond17MC_PoissonOOTPU.root','PU/MyDataPileupHistogram_2016.root'] #MC and data files for PileUp reweighting (2016)
    config.Data.lumiMask = 'json/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'

if runningEra == 1:
    config.General.workArea = 'crab_projects/samples_data_2017/'
    #config.JobType.inputFiles = ['PU/MCpileUp_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU.root','PU/MyDataPileupHistogram_2017.root'] #MC and data files for PileUp reweighting (2017)
    config.Data.lumiMask = 'json/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'


config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
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

    #First the muon datasets

    if runningEra == 0:

        config.JobType.pyCfgParams = ['runningOnData=True','runningOnMuons=True','runningEra=0'] # Configure 2016 data jobs - muons
    
        config.General.requestName = '2016_WPiGammaAnalysis_SingleMu_B'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2016B-17Jul2018_ver2-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleMu_C'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2016C-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleMu_D'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2016D-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleMu_E'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2016E-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleMu_F'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2016F-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleMu_G'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2016G-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleMu_H'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2016H-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

    if runningEra == 1:

        config.JobType.pyCfgParams = ['runningOnData=True','runningOnMuons=True','runningEra=1'] # Configure 2017 data jobs - muons

        config.General.requestName = '2017_WPiGammaAnalysis_SingleMu_B'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2017B-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2017_WPiGammaAnalysis_SingleMu_C'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2017C-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2017_WPiGammaAnalysis_SingleMu_D'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2017D-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2017_WPiGammaAnalysis_SingleMu_E'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2017E-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2017_WPiGammaAnalysis_SingleMu_F'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2017F-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()


    #Now the electron datasets    

    if runningEra == 0:

        config.JobType.pyCfgParams = ['runningOnData=True','runningOnMuons=False','runningEra=0'] # Configure 2016 data jobs - electrons
    
        config.General.requestName = '2016_WPiGammaAnalysis_SingleEle_B'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2016B-17Jul2018_ver2-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleEle_C'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2016C-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleEle_D'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2016D-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleEle_E'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2016E-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleEle_F'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2016F-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleEle_G'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2016G-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleEle_H'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2016H-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

    if runningEra == 1:

        config.JobType.pyCfgParams = ['runningOnData=True','runningOnMuons=False','runningEra=1'] # Configure 2017 data jobs - electrons
        
        config.General.requestName = '2017_WPiGammaAnalysis_SingleEle_B'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2017B-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2017_WPiGammaAnalysis_SingleEle_C'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2017C-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2017_WPiGammaAnalysis_SingleEle_D'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2017D-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2017_WPiGammaAnalysis_SingleEle_E'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2017E-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2017_WPiGammaAnalysis_SingleEle_F'
        config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2017F-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
