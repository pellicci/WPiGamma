from CRABAPI.RawCommand import crabCommand
from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config
config = Configuration()

config.section_('General')
config.General.transferOutputs = True

runningEra = 2 # 0 = 2016, 1 = 2017, 2 = 2018

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/run_WPiGammaAnalysis.py'
#config.JobType.psetName = 'cmssw_config/run_LeptonMultiplicity.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['WPiGammaAnalysis_output.root']
#config.JobType.outputFiles = ['LeptonMultiplicity_output.root']
config.JobType.allowUndistributedCMSSW = True #Otherwise get an error for incompatibility of architecture(slc7)/release(9_4_10). It is safe according to https://hypernews.cern.ch/HyperNews/CMS/get/computing-tools/4935/2.html

if runningEra == 0:
    config.General.workArea = 'crab_projects/samples_MC_2016/'
    #config.General.workArea = 'crab_projects/samples_MC_2016_bTagEfficiency_medium/'
    #config.General.workArea = 'crab_projects/samples_LeptonStudy_2016/'
    config.JobType.inputFiles = ['MCpileUp_2016_25ns_Moriond17MC_PoissonOOTPU.root','MyDataPileupHistogram_2016.root','DeepCSV_2016LegacySF_WP_V1.csv','bTagEff_2016.root'] #MC and data files for PileUp reweighting (2016)

if runningEra == 1:
    config.General.workArea = 'crab_projects/samples_MC_2017/'
    #config.General.workArea = 'crab_projects/samples_MC_2017_bTagEfficiency_medium/'
    #config.General.workArea = 'crab_projects/samples_LeptonStudy_2017/'
    config.JobType.inputFiles = ['MCpileUp_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU.root','MyDataPileupHistogram_2017.root','DeepCSV_94XSF_WP_V4_B_F.csv','bTagEff_2017.root'] #MC and data files for PileUp reweighting (2017)

if runningEra == 2:
    config.General.workArea = 'crab_projects/samples_MC_2018_3referee/'
    #config.General.workArea = 'crab_projects/samples_MC_2018_bTagEfficiency_medium/'
    #config.General.workArea = 'crab_projects/samples_LeptonStudy_2018/'
    config.JobType.inputFiles = ['MCpileUp_2018_25ns_JuneProjectionFull18_PoissonOOTPU.root','MyDataPileupHistogram_2018.root','DeepCSV_102XSF_WP_V1.csv','bTagEff_2018.root'] #MC and data files for PileUp reweighting (2017)


config.section_('Data')
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'Automatic'
#config.Data.splitting = 'LumiBased'
#config.Data.unitsPerJob = 5
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

    if runningEra == 0: #2016

        config.JobType.pyCfgParams = ['runningOnData=False','runningEra=0'] # Configure 2016 MC signal jobs 

        config.General.requestName = '2016_WPiGammaAnalysis_Signal_WPlus'
        #config.Data.inputDataset = '/WPlusPiGamma_GENSIM_80X_v2/rselvati-WPlusPiGamma_MINIAODSIM_94X_2016_v2-793fab47b42ceb2f28d6e316553fc80e/USER'
        config.Data.inputDataset = '/TTToSemilepWPlusPiGamma_GENSIM_80X_2016_v1/pellicci-WPlusPiGamma_MINIAODSIM_102X_2016_v1-793fab47b42ceb2f28d6e316553fc80e/USER'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_Signal_WMinus'
        #config.Data.inputDataset = '/WMinusPiGamma_GENSIM_80X_v2/rselvati-WMinusPiGamma_MINIAODSIM_94X_2016_v2-793fab47b42ceb2f28d6e316553fc80e/USER'
        config.Data.inputDataset = '/TTToSemilepWMinusPiGamma_GENSIM_80X_2016_v1/pellicci-WMinusPiGamma_MINIAODSIM_102X_2016_v1-793fab47b42ceb2f28d6e316553fc80e/USER'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

    if runningEra == 1: #2017

        config.JobType.pyCfgParams = ['runningOnData=False','runningEra=1'] # Configure 2017 MC signal jobs 

        config.General.requestName = '2017_WPiGammaAnalysis_Signal_WPlus'
        #config.Data.inputDataset = '/WPlusPiGamma_GENSIM_94X_2017_v7/rselvati-WPlusPiGamma_MINIAODSIM_94X_2017_v7-4fa7a8e0675466927bc37686a776c034/USER'
        config.Data.inputDataset = '/TTToSemilepWPlusPiGamma_GENSIM_94X_2017_v1/pellicci-WPlusPiGamma_MINIAODSIM_102X_2017_v1-93b25218fa2b50cbb76b8d1c9edf4a08/USER'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2017_WPiGammaAnalysis_Signal_WMinus'
        #config.Data.inputDataset = '/WMinusPiGamma_GENSIM_94X_2017_v7/rselvati-WMinusPiGamma_MINIAODSIM_94X_2017_v7-4fa7a8e0675466927bc37686a776c034/USER'
        config.Data.inputDataset = '/TTToSemilepWMinusPiGamma_GENSIM_94X_2017_v1/pellicci-WMinusPiGamma_MINIAODSIM_102X_2017_v1-93b25218fa2b50cbb76b8d1c9edf4a08/USER'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

    if runningEra == 2: #2018

        config.JobType.pyCfgParams = ['runningOnData=False','runningEra=2'] # Configure 2018 MC signal jobs 

        config.General.requestName = '2018_WPiGammaAnalysis_Signal_WPlus'
        #config.Data.inputDataset = '/WPlusPiGamma_102X_2018/rselvati-WPlusPiGamma_MINIAODSIM_102X_2018_v1-9caaa56992b50a7cf6d420a9959af8e5/USER'
        config.Data.inputDataset = '/TTToSemilepWPlusPiGamma_102X_2018_v1/pellicci-WPlusPiGamma_MINIAODSIM_102X_2018_v1-9caaa56992b50a7cf6d420a9959af8e5/USER'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2018_WPiGammaAnalysis_Signal_WMinus'
        #config.Data.inputDataset = '/WMinusPiGamma_102X_2018/rselvati-WMinusPiGamma_MINIAODSIM_102X_2018_v1-9caaa56992b50a7cf6d420a9959af8e5/USER'
        config.Data.inputDataset = '/TTToSemilepWMinusPiGamma_102X_2018_v1/pellicci-WMinusPiGamma_MINIAODSIM_102X_2018_v1-9caaa56992b50a7cf6d420a9959af8e5/USER'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
