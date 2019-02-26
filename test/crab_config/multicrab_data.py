from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = Configuration()

config.section_('General')
config.General.transferOutputs = True

config.General.workArea = 'crab_projects/dataprocess/'

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/run_WPiGammaAnalysis.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['MCpileUp_2016_25ns_Moriond17MC_PoissonOOTPU.root','MyDataPileupHistogram.root'] #MC and data files for PileUp reweighting
config.JobType.outputFiles = ['WPiGammaAnalysis_output.root']

config.JobType.pyCfgParams = ['runningOnData=True']

config.section_('Data')
config.Data.lumiMask = 'json/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
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

    config.JobType.pyCfgParams = ['runningOnMuons=True']
    
    config.General.requestName = 'WPiGammaAnalysis_SingleMu_B'
    config.Data.unitsPerJob = 50
    config.Data.inputDataset = '/SingleMuon/Run2016B-17Jul2018_ver2-v1/MINIAOD'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
    config.General.requestName = 'WPiGammaAnalysis_SingleMu_C'
    config.Data.unitsPerJob = 50
    config.Data.inputDataset = '/SingleMuon/Run2016C-17Jul2018-v1/MINIAOD'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'WPiGammaAnalysis_SingleMu_D'
    config.Data.unitsPerJob = 50
    config.Data.inputDataset = '/SingleMuon/Run2016D-17Jul2018-v1/MINIAOD'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
    config.General.requestName = 'WPiGammaAnalysis_SingleMu_E'
    config.Data.unitsPerJob = 50
    config.Data.inputDataset = '/SingleMuon/Run2016E-17Jul2018-v1/MINIAOD'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'WPiGammaAnalysis_SingleMu_F'
    config.Data.unitsPerJob = 50
    config.Data.inputDataset = '/SingleMuon/Run2016F-17Jul2018-v1/MINIAOD'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
    config.General.requestName = 'WPiGammaAnalysis_SingleMu_G'
    config.Data.unitsPerJob = 50
    config.Data.inputDataset = '/SingleMuon/Run2016G-17Jul2018-v1/MINIAOD'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
    config.General.requestName = 'WPiGammaAnalysis_SingleMu_H'
    config.Data.unitsPerJob = 50
    config.Data.inputDataset = '/SingleMuon/Run2016H-17Jul2018-v1/MINIAOD'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()


    #Now the electron datasets    
    config.JobType.pyCfgParams = ['runningOnMuons=False']
    
    config.General.requestName = 'WPiGammaAnalysis_SingleEle_B'
    config.Data.unitsPerJob = 50
    config.Data.inputDataset = '/SingleElectron/Run2016B-17Jul2018_ver2-v1/MINIAOD'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'WPiGammaAnalysis_SingleEle_C'
    config.Data.unitsPerJob = 50
    config.Data.inputDataset = '/SingleElectron/Run2016C-17Jul2018-v1/MINIAOD'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
    config.General.requestName = 'WPiGammaAnalysis_SingleEle_D'
    config.Data.unitsPerJob = 50
    config.Data.inputDataset = '/SingleElectron/Run2016D-17Jul2018-v1/MINIAOD'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
    config.General.requestName = 'WPiGammaAnalysis_SingleEle_E'
    config.Data.unitsPerJob = 50
    config.Data.inputDataset = '/SingleElectron/Run2016E-17Jul2018-v1/MINIAOD'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
    config.General.requestName = 'WPiGammaAnalysis_SingleEle_F'
    config.Data.unitsPerJob = 50
    config.Data.inputDataset = '/SingleElectron/Run2016F-17Jul2018-v1/MINIAOD'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
    config.General.requestName = 'WPiGammaAnalysis_SingleEle_G'
    config.Data.unitsPerJob = 50
    config.Data.inputDataset = '/SingleElectron/Run2016G-17Jul2018-v1/MINIAOD'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    
    config.General.requestName = 'WPiGammaAnalysis_SingleEle_H'
    config.Data.unitsPerJob = 50
    config.Data.inputDataset = '/SingleElectron/Run2016H-17Jul2018-v1/MINIAOD'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

