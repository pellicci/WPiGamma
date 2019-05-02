from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

runningEra = 1 # 0 = 2016, 1 = 2017, 2 = 2018

config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects'
config.section_('JobType')
config.JobType.pluginName = 'PrivateMC'
config.section_('Data')

if runningEra == 0:
    config.General.requestName = 'WPlusPiGamma_Pythia8_GENSIM_80XV1'
    config.JobType.psetName = 'cmssw_config/WPlusPiGamma_13TeV_pythia8_GENSIM_cfg.py'
    config.JobType.outputFiles = ['WPiGamma_pythia8_GENSIM.root']
    config.Data.outputPrimaryDataset = 'WPlusPiGamma_GENSIM_80XV1'

if runningEra == 1:
    config.General.requestName = 'WPlusPiGamma_Pythia8_GENSIM_94X'
    config.JobType.psetName = 'cmssw_config/WPlusPiGamma_13TeV_pythia8_GENSIM_2017_cfg.py'
    config.JobType.outputFiles = ['WPiGamma_pythia8_GENSIM_2017.root']
    config.Data.outputPrimaryDataset = 'WPlusPiGamma_GENSIM_94X'

config.Data.splitting = 'EventBased' # Can only be set to EventBased if pluginName = PrivateMC
config.Data.unitsPerJob = 10
NJOBS = 6000
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True

if runningEra == 0:
    config.Data.outputDatasetTag = 'WPlusPiGamma_GENSIM_80XV1'
if runningEra == 1:
    config.Data.outputDatasetTag = 'WPlusPiGamma_GENSIM_94X'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'
