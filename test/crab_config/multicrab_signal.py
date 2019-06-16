from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = Configuration()

config.section_('General')
config.General.transferOutputs = True

runningEra = 1 # 0 = 2016, 1 = 2017, 2 = 2018

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/run_WPiGammaAnalysis.py'
#config.JobType.psetName = 'cmssw_config/run_LeptonMultiplicity.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['WPiGammaAnalysis_output.root']
#config.JobType.outputFiles = ['LeptonMultiplicity_output.root']

if runningEra == 0:
    config.General.workArea = 'crab_projects/samples_MC_2016/'
    #config.General.workArea = 'crab_projects/samples_LeptonStudy_2016/'
    config.JobType.inputFiles = ['MCpileUp_2016_25ns_Moriond17MC_PoissonOOTPU.root','MyDataPileupHistogram_2016.root'] #MC and data files for PileUp reweighting (2016)

if runningEra == 1:
    config.General.workArea = 'crab_projects/samples_MC_2017/'
    #config.General.workArea = 'crab_projects/samples_LeptonStudy_2017/'
    config.JobType.inputFiles = ['MCpileUp_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU.root','MyDataPileupHistogram_2017.root'] #MC and data files for PileUp reweighting (2017)


config.section_('Data')
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'Automatic'
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

    if runningEra == 0: #2016

        config.JobType.pyCfgParams = ['runningOnData=False','runningEra=0'] # Configure 2016 MC signal jobs 

        config.General.requestName = '2016_WPiGammaAnalysis_Signal_WPlus'
        config.Data.inputDataset = '/WPlusPiGamma_GENSIM_80XV1/rselvati-WPlusPiGamma_MINIAODSIM_94XV3-9a5eaf7e0651dc0135fee9d652526a74/USER'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2016_WPiGammaAnalysis_Signal_WMinus'
        config.Data.inputDataset = '/WMinusPiGamma_GENSIM_80XV1/rselvati-WMinusPiGamma_MINIAODSIM_94XV3-9a5eaf7e0651dc0135fee9d652526a74/USER'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

    if runningEra == 1: #2017

        config.JobType.pyCfgParams = ['runningOnData=False','runningEra=1'] # Configure 2017 MC signal jobs 

        config.General.requestName = '2017_WPiGammaAnalysis_Signal_WPlus'
        config.Data.inputDataset = '/WPlusPiGamma_GENSIM_94X_2017_v5/rselvati-WPlusPiGamma_MINIAODSIM_94X_2017_v5-10a67796329f4238191918f07d0f7633/USER'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
        
        config.General.requestName = '2017_WPiGammaAnalysis_Signal_WMinus'
        config.Data.inputDataset = '/WMinusPiGamma_GENSIM_94X_2017_v3/rselvati-WMinusPiGamma_MINIAODSIM_94X_2017_v3-10a67796329f4238191918f07d0f7633/USER'
        p = Process(target=submit, args=(config,))
        p.start()
        p.join()
