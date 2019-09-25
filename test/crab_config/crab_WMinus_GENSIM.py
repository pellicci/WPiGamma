from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'WMinusPiGamma_Pythia8_GENSIM_80X_v2'
config.General.workArea = 'crab_projects'
config.JobType.allowUndistributedCMSSW = True #Otherwise get an error for incompatibility of architecture(slc7)/release(8_0_28). It is safe according to https://hypernews.cern.ch/HyperNews/CMS/get/computing-tools/4935/2.html

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/WMinusPiGamma_13TeV_pythia8_GENSIM_2016_cfg.py'
config.JobType.pluginName = 'PrivateMC'
config.JobType.outputFiles = ['WPiGamma_pythia8_GENSIM.root']

config.section_('Data')
config.Data.outputPrimaryDataset = 'WMinusPiGamma_GENSIM_80X_v2'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 5
NJOBS = 10000 #Do not increase: maximum number of jobs per task is 10k
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'WMinusPiGamma_GENSIM_80X_v2'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'
