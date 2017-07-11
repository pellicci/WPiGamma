from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
#config.General.workArea = 'crab_projects/samples/'
config.General.workArea = 'crab_projects/samples_LeptonStudy/'

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/run_WPiGammaAnalysis.py'
#config.JobType.psetName = 'cmssw_config/run_LeptonMultiplicity.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['MCpileUp_25ns_Recent2016.root','pileUpHistogramFromjson_Nominal.root' ] #data files for PileUp reweighting
config.JobType.outputFiles = ['WPiGammaAnalysis_output.root']
#config.JobType.outputFiles = ['LeptonMultiplicity_output.root']
config.JobType.pyCfgParams = ['runningOnData=False']

config.section_('Data')
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
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

    config.General.requestName = 'WPiGammaAnalysis_Signal_WPlus'
    config.Data.inputDataset = '/WPlusPiGamma_GENSIM_80XV1/pellicci-WPlusPiGamma_MINIAODSIM_80XV1-1ccbf35f8ddfb81066cfc2789c02310c/USER'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'WPiGammaAnalysis_Signal_WMinus'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/WMinusPiGamma_GENSIM_80XV1/pellicci-WMinusPiGamma_MINIAODSIM_80XV1-1ccbf35f8ddfb81066cfc2789c02310c/USER'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
