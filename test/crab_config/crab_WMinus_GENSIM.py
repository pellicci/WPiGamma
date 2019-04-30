from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

runningOn2017 = True # Decide whether to run on 2016 or 2017 data

config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects'
config.section_('JobType')
config.JobType.pluginName = 'PrivateMC'
config.section_('Data')

if runningOn2017:
    config.General.requestName = 'WMinusPiGamma_Pythia8_GENSIM_94X'
    config.JobType.psetName = 'cmssw_config/WMinusPiGamma_13TeV_pythia8_GENSIM_2017_cfg.py'
    config.JobType.outputFiles = ['WPiGamma_pythia8_GENSIM_2017.root']
    config.Data.outputPrimaryDataset = 'WMinusPiGamma_GENSIM_94X'
else:
    config.General.requestName = 'WMinusPiGamma_Pythia8_GENSIM_80XV1'
    config.JobType.psetName = 'cmssw_config/WMinusPiGamma_13TeV_pythia8_GENSIM_cfg.py'
    config.JobType.outputFiles = ['WPiGamma_pythia8_GENSIM.root']
    config.Data.outputPrimaryDataset = 'WMinusPiGamma_GENSIM_80XV1'

config.Data.splitting = 'EventBased' # Can only be set to EventBased if pluginName = PrivateMC 
config.Data.unitsPerJob = 10
NJOBS = 6000
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True

if runningOn2017:
    config.Data.outputDatasetTag = 'WMinusPiGamma_GENSIM_94X'
else:
    config.Data.outputDatasetTag = 'WMinusPiGamma_GENSIM_80XV1'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'
