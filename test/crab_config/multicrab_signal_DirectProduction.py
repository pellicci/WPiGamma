from CRABAPI.RawCommand import crabCommand
from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config
config = Configuration()

config.section_('General')
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/run_WPiGammaAnalysis.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['WPiGammaAnalysis_output.root']
config.JobType.allowUndistributedCMSSW = True #Otherwise get an error for incompatibility of architecture(slc7)/release(9_4_10). It is safe according to https://hypernews.cern.ch/HyperNews/CMS/get/computing-tools/4935/2.html

config.General.workArea = 'crab_projects/samples_MC_2018/'
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

    config.JobType.pyCfgParams = ['runningOnData=False','runningEra=2'] # Configure 2018 MC signal jobs 

    config.General.requestName = '2018_WPiGammaAnalysis_WPlusPiGamma_DirectProduction'
    config.Data.inputDataset = '/WPlusPiGamma_DirectProduction_102X_2018_v1/rselvati-WPlusPiGamma_DirectProduction_MINIAOD_102X_2018_v1-9caaa56992b50a7cf6d420a9959af8e5/USER'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
        
    config.General.requestName = '2018_WPiGammaAnalysis_WMinusPiGamma_DirectProduction'
    config.Data.inputDataset = '/WMinusPiGamma_DirectProduction_102X_2018_v1/rselvati-WMinusPiGamma_DirectProduction_MINIAOD_102X_2018_v1-9caaa56992b50a7cf6d420a9959af8e5/USER'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
