from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'WPiGamma_Pythia8_GENSIM_80XV1'
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = 'WPiGamma_13TeV_pythia8_GENSIM_cfg.py'
config.JobType.pluginName = 'PrivateMC'
config.JobType.outputFiles = ['WPiGamma_pythia8_GENSIM.root']

config.section_('Data')
config.Data.outputPrimaryDataset = 'WPiGamma_GENSIM_80XV1'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 10
NJOBS = 5000
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'WPiGamma_GENSIM_80XV1'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'
