from CRABClient.UserUtilities import config
config = config()

doPlus = True
runningEra = 2 # 0 = 2016, 1 = 2017, 2 = 2018 

config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects'
config.section_('JobType')
config.JobType.pluginName = 'Analysis'

config.JobType.allowUndistributedCMSSW = True

config.section_('Data')

if runningEra == 0:

    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_MINIAOD_2016_cfg.py'

    if doPlus:
        config.General.requestName = 'WPlusPiGamma_Pythia8_MINIAODSIM_102X_2016_v1'
        config.Data.inputDataset = '/TTToSemilepWPlusPiGamma_GENSIM_80X_2016_v1/pellicci-WPlusPiGamma_RECOSIM_80X_2016_v1-8388a14248b3bf0d0a82d70d71bce005/USER'
    else:
        config.General.requestName = 'WMinusPiGamma_Pythia8_MINIAODSIM_102X_2016_v1'
        config.Data.inputDataset = '/TTToSemilepWMinusPiGamma_GENSIM_80X_2016_v1/pellicci-WMinusPiGamma_RECOSIM_80X_2016_v1-8388a14248b3bf0d0a82d70d71bce005/USER'

if runningEra == 1:

    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_MINIAOD_2017_cfg.py'

    if doPlus:
        config.General.requestName = 'WPlusPiGamma_Pythia8_MINIAODSIM_102X_2017_v1'
        config.Data.inputDataset = '/TTToSemilepWPlusPiGamma_GENSIM_94X_2017_v1/pellicci-WPlusPiGamma_RECOSIM_94X_2017_v1-fcfc615a65be9fb627e3afc83a7469ff/USER'
    else:
        config.General.requestName = 'WMinusPiGamma_Pythia8_MINIAODSIM_102X_2017_v1'
        config.Data.inputDataset = '/TTToSemilepWMinusPiGamma_GENSIM_94X_2017_v1/pellicci-WMinusPiGamma_RECOSIM_94X_2017_v1-fcfc615a65be9fb627e3afc83a7469ff/USER'

if runningEra == 2:

    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_MINIAOD_2018_cfg.py'

    if doPlus:
        config.General.requestName = 'WPlusPiGamma_Pythia8_MINIAODSIM_102X_2018_v1'
        config.Data.inputDataset = '/TTToSemilepWPlusPiGamma_102X_2018_v1/pellicci-WPlusPiGamma_RECOSIM_102X_2018_v1-f640515c3dc7ddab735258ddeb1359c7/USER'
    else:
        config.General.requestName = 'WMinusPiGamma_Pythia8_MINIAODSIM_102X_2018_v1'
        config.Data.inputDataset = '/TTToSemilepWMinusPiGamma_102X_2018_v1/pellicci-WMinusPiGamma_RECOSIM_102X_2018_v1-f640515c3dc7ddab735258ddeb1359c7/USER'


config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
config.Data.publication = True

if runningEra == 0:

    if doPlus:
        config.Data.outputDatasetTag = 'WPlusPiGamma_MINIAODSIM_102X_2016_v1'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_MINIAODSIM_102X_2016_v1'

if runningEra == 1:

    if doPlus:
        config.Data.outputDatasetTag = 'WPlusPiGamma_MINIAODSIM_102X_2017_v1'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_MINIAODSIM_102X_2017_v1'

if runningEra == 2:

    if doPlus:
        config.Data.outputDatasetTag = 'WPlusPiGamma_MINIAODSIM_102X_2018_v1'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_MINIAODSIM_102X_2018_v1'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'
