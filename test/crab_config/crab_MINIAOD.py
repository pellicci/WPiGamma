from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

doPlus = True
runningEra = 1 # 0 = 2016, 1 = 2017, 2 = 2018 

config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects'
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.section_('Data')

if runningEra == 0:

    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_MINIAOD_cfg.py'

    if doPlus:
        config.General.requestName = 'WPlusPiGamma_Pythia8_MINIAODSIM_94XV3'
        config.Data.inputDataset = '/WPlusPiGamma_GENSIM_80XV1/rselvati-WPlusPiGamma_RECOSIM_80XV1-c536d85e5d9fce8caa236321c5af92c3/USER'
    else:
        config.General.requestName = 'WMinusPiGamma_Pythia8_MINIAODSIM_94XV3'
        config.Data.inputDataset = '/WMinusPiGamma_GENSIM_80XV1/rselvati-WMinusPiGamma_RECOSIM_80XV1-c536d85e5d9fce8caa236321c5af92c3/USER'


if runningEra == 1:

    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_MINIAOD_2017_cfg.py'

    if doPlus:
        config.General.requestName = 'WPlusPiGamma_Pythia8_MINIAODSIM_94X_2017_v5'
        config.Data.inputDataset = '/WPlusPiGamma_GENSIM_94X_2017_v5/rselvati-WPlusPiGamma_RECOSIM_94X_2017_v5-c3d6de13a4792afb4dd0c4ab58e49a3d/USER'
    else:
        config.General.requestName = 'WMinusPiGamma_Pythia8_MINIAODSIM_94X_2017_v3'
        config.Data.inputDataset = '/WMinusPiGamma_GENSIM_94X_2017_v3/rselvati-WMinusPiGamma_RECOSIM_94X_2017_v3-c3d6de13a4792afb4dd0c4ab58e49a3d/USER'


config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True

if runningEra == 0:

    if doPlus:
        config.Data.outputDatasetTag = 'WPlusPiGamma_MINIAODSIM_94X'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_MINIAODSIM_94X'

if runningEra == 1:

    if doPlus:
        config.Data.outputDatasetTag = 'WPlusPiGamma_MINIAODSIM_94X_2017_v5'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_MINIAODSIM_94X_2017_v3'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'
