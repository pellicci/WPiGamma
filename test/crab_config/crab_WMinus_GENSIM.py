from CRABClient.UserUtilities import config
config = config()

runningEra = 0 # 0 = 2016, 1 = 2017, 2 = 2018

config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects'
config.section_('JobType')
config.JobType.pluginName = 'PrivateMC'

config.JobType.allowUndistributedCMSSW = True

config.section_('Data')

if runningEra == 0:
    config.General.requestName = 'WMinusPiGamma_Pythia8_GENSIM_80X_2016_v1'
    config.JobType.psetName = 'cmssw_config/WMinusPiGamma_13TeV_pythia8_GENSIM_2016_cfg.py'
    config.JobType.outputFiles = ['WPiGamma_pythia8_GENSIM.root']
    config.Data.outputPrimaryDataset = 'TTToSemilepWMinusPiGamma_GENSIM_80X_2016_v1'

if runningEra == 1:
    config.General.requestName = 'WMinusPiGamma_Pythia8_GENSIM_94X_2017_v1'
    config.JobType.psetName = 'cmssw_config/WMinusPiGamma_13TeV_pythia8_GENSIM_2017_cfg.py'
    config.JobType.outputFiles = ['WPiGamma_pythia8_GENSIM_2017.root']
    config.Data.outputPrimaryDataset = 'TTToSemilepWMinusPiGamma_GENSIM_94X_2017_v1'

if runningEra == 2:
    config.General.requestName = 'WMinusPiGamma_Pythia8_GENSIM_102X_2018_v1'
    config.JobType.psetName = 'cmssw_config/WMinusPiGamma_13TeV_pythia8_GENSIM_2018_cfg.py'
    config.JobType.outputFiles = ['WPiGamma_pythia8_GENSIM_2018.root']
    config.Data.outputPrimaryDataset = 'TTToSemilepWMinusPiGamma_102X_2018_v1'

config.Data.splitting = 'EventBased' # Can only be set to EventBased if pluginName = PrivateMC 
config.Data.unitsPerJob = 5
NJOBS = 8000
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = True

if runningEra == 0:
    config.Data.outputDatasetTag = 'TTToSemilepWMinusPiGamma_GENSIM_80X_2016_v1'
if runningEra == 1:
    config.Data.outputDatasetTag = 'TTToSemilepWMinusPiGamma_GENSIM_94X_2017_v1'
if runningEra == 2:
    config.Data.outputDatasetTag = 'TTToSemilepWMinusPiGamma_GENSIM_102X_2018_v1'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'
