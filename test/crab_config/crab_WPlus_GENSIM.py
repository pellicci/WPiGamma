from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

runningEra = 2 # 0 = 2016, 1 = 2017, 2 = 2018

config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects'
config.section_('JobType')
config.JobType.pluginName = 'PrivateMC'

config.JobType.allowUndistributedCMSSW = True #Otherwise get an error for incompatibility of architecture(slc7)/release(10_2_10). It is safe according to https://hypernews.cern.ch/HyperNews/CMS/get/computing-tools/4935/2.html

config.section_('Data')

if runningEra == 0:
    config.General.requestName = 'WPlusPiGamma_Pythia8_GENSIM_80X'
    config.JobType.psetName = 'cmssw_config/WPlusPiGamma_13TeV_pythia8_GENSIM_2016_cfg.py'
    config.JobType.outputFiles = ['WPiGamma_pythia8_GENSIM.root']
    config.Data.outputPrimaryDataset = 'WPlusPiGamma_80X_2016'

if runningEra == 1:
    config.General.requestName = 'WPlusPiGamma_Pythia8_GENSIM_94X_2017_v7'
    config.JobType.psetName = 'cmssw_config/WPlusPiGamma_13TeV_pythia8_GENSIM_2017_cfg.py'
    config.JobType.outputFiles = ['WPiGamma_pythia8_GENSIM_2017.root']
    config.Data.outputPrimaryDataset = 'WPlusPiGamma_94X_2017_v7'

if runningEra == 2:
    config.General.requestName = 'WPlusPiGamma_Pythia8_GENSIM_102X_2018_v1'
    config.JobType.psetName = 'cmssw_config/WPlusPiGamma_13TeV_pythia8_GENSIM_2018_cfg.py'
    config.JobType.outputFiles = ['WPiGamma_pythia8_GENSIM_2018.root']
    config.Data.outputPrimaryDataset = 'WPlusPiGamma_102X_2018' #That is the name of the first "block": WPlusPiGamma_102X_2018/*/*

Config.Data.splitting = 'EventBased' # Can only be set to EventBased if pluginName = PrivateMC
config.Data.unitsPerJob = 5
NJOBS = 8000
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True

if runningEra == 0:
    config.Data.outputDatasetTag = 'WPlusPiGamma_GENSIM_80X_2016'
if runningEra == 1:
    config.Data.outputDatasetTag = 'WPlusPiGamma_GENSIM_94X_2017_v7'
if runningEra == 2:
    config.Data.outputDatasetTag = 'WPlusPiGamma_GENSIM_102X_2018_v1'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'
