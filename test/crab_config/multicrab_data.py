from CRABAPI.RawCommand import crabCommand
from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config

config = Configuration()

config.section_('General')
config.General.transferOutputs = True

runningEra = 2 # 0 = 2016, 1 = 2017, 2 = 2018

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/run_WPiGammaAnalysis.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['WPiGammaAnalysis_output.root']
config.JobType.allowUndistributedCMSSW = True #Otherwise get an error for incompatibility of architecture(slc7)/release(9_4_10). It is safe according to https://hypernews.cern.ch/HyperNews/CMS/get/computing-tools/4935/2.html

config.section_('Data')


if runningEra == 0:
    config.General.workArea = 'crab_projects/samples_data_2016/'
    config.Data.lumiMask = 'json/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
    config.JobType.inputFiles = ['DeepCSV_2016LegacySF_WP_V1.csv','bTagEff_2016.root']

if runningEra == 1:
    config.General.workArea = 'crab_projects/samples_data_2017/'
    config.Data.lumiMask = 'json/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
    config.JobType.inputFiles = ['DeepCSV_94XSF_WP_V4_B_F.csv','bTagEff_2017.root']

if runningEra == 2:
    config.General.workArea = 'crab_projects/samples_data_2018/'
    config.Data.lumiMask = 'json/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
    config.JobType.inputFiles = ['DeepCSV_102XSF_WP_V1.csv','bTagEff_2018.root']


config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
config.Data.outLFNDirBase = '/store/user/rselvati/'
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


    ################################################
    #                                              #
    #----------------- Muons 2016 -----------------#
    #                                              #
    ################################################

    if runningEra == 0:

        config.JobType.pyCfgParams = ['runningOnData=True','runningOnMuons=True','runningEra=0','runningOn2018D=False'] # Configure 2016 data jobs - muons
    
        config.General.requestName = '2016_WPiGammaAnalysis_SingleMu_B'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2016B-17Jul2018_ver2-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleMu_C'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2016C-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleMu_D'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2016D-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleMu_E'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2016E-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleMu_F'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2016F-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleMu_G'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2016G-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleMu_H'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2016H-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()


    ################################################
    #                                              #
    #----------------- Muons 2017 -----------------#
    #                                              #
    ################################################

    if runningEra == 1:

        config.JobType.pyCfgParams = ['runningOnData=True','runningOnMuons=True','runningEra=1','runningOn2018D=False'] # Configure 2017 data jobs - muons

        config.General.requestName = '2017_WPiGammaAnalysis_SingleMu_B'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2017B-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2017_WPiGammaAnalysis_SingleMu_C'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2017C-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2017_WPiGammaAnalysis_SingleMu_D'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2017D-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2017_WPiGammaAnalysis_SingleMu_E'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2017E-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2017_WPiGammaAnalysis_SingleMu_F'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2017F-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()


    ################################################
    #                                              #
    #----------------- Muons 2018 -----------------#
    #                                              #
    ################################################

    if runningEra == 2: 

        config.JobType.pyCfgParams = ['runningOnData=True','runningOnMuons=True','runningEra=2','runningOn2018D=False'] # Configure 2018 data jobs - muons - eras A,B,C

        config.General.requestName = '2018_WPiGammaAnalysis_SingleMu_A'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2018A-17Sep2018-v2/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = '2018_WPiGammaAnalysis_SingleMu_B'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2018B-17Sep2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = '2018_WPiGammaAnalysis_SingleMu_C'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2018C-17Sep2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        #Muons 2018 - Era D
        config.JobType.pyCfgParams = ['runningOnData=True','runningOnMuons=True','runningEra=2','runningOn2018D=True'] # Configure 2018 data jobs - muons - era D

        config.General.requestName = '2018_WPiGammaAnalysis_SingleMu_D'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleMuon/Run2018D-22Jan2019-v2/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()


    ################################################
    #                                              #
    #--------------- Electrons 2016 ---------------#
    #                                              #
    ################################################

    if runningEra == 0:

        config.JobType.pyCfgParams = ['runningOnData=True','runningOnMuons=False','runningEra=0','runningOn2018D=False'] # Configure 2016 data jobs - electrons
    
        config.General.requestName = '2016_WPiGammaAnalysis_SingleEle_B'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2016B-17Jul2018_ver2-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleEle_C'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2016C-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleEle_D'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2016D-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleEle_E'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2016E-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleEle_F'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2016F-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleEle_G'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2016G-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_SingleEle_H'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2016H-17Jul2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()


    ################################################
    #                                              #
    #--------------- Electrons 2017 ---------------#
    #                                              #
    ################################################

    if runningEra == 1:

        config.JobType.pyCfgParams = ['runningOnData=True','runningOnMuons=False','runningEra=1','runningOn2018D=False'] # Configure 2017 data jobs - electrons
        
        config.General.requestName = '2017_WPiGammaAnalysis_SingleEle_B'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2017B-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2017_WPiGammaAnalysis_SingleEle_C'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2017C-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2017_WPiGammaAnalysis_SingleEle_D'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2017D-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2017_WPiGammaAnalysis_SingleEle_E'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2017E-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2017_WPiGammaAnalysis_SingleEle_F'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/SingleElectron/Run2017F-31Mar2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()


    ################################################
    #                                              #
    #--------------- Electrons 2018 ---------------#
    #                                              #
    ################################################

    if runningEra == 2: 

        config.JobType.pyCfgParams = ['runningOnData=True','runningOnMuons=False','runningEra=2','runningOn2018D=False'] # Configure 2018 data jobs - electrons - eras A,B,C

        config.General.requestName = '2018_WPiGammaAnalysis_SingleEle_A'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/EGamma/Run2018A-17Sep2018-v2/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = '2018_WPiGammaAnalysis_SingleEle_B'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/EGamma/Run2018B-17Sep2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        config.General.requestName = '2018_WPiGammaAnalysis_SingleEle_C'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/EGamma/Run2018C-17Sep2018-v1/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

        #Electrons 2018 - Era D
        config.JobType.pyCfgParams = ['runningOnData=True','runningOnMuons=False','runningEra=2','runningOn2018D=True'] # Configure 2018 data jobs - electrons - era D

        config.General.requestName = '2018_WPiGammaAnalysis_SingleEle_D'
        #config.Data.unitsPerJob = 50
        config.Data.inputDataset = '/EGamma/Run2018D-22Jan2019-v2/MINIAOD'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
