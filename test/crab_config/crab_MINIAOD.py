from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

doPlus = True
runningEra = 1 # 0 = 2016, 1 = 2017, 2 = 2018 

config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects'
config.section_('JobType')
config.JobType.pluginName = 'Analysis'

config.JobType.allowUndistributedCMSSW = True #Otherwise get an error for incompatibility of architecture(slc7)/release(9_4_10). It is safe according to https://hypernews.cern.ch/HyperNews/CMS/get/computing-tools/4935/2.html

config.section_('Data')

if runningEra == 0:

    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_MINIAOD_2016_cfg.py'

    if doPlus:
        config.General.requestName = 'WPlusPiGamma_Pythia8_MINIAODSIM_94X_2016_v2'
        config.Data.inputDataset = '/WPlusPiGamma_GENSIM_80X_v2/rselvati-WPlusPiGamma_RECOSIM_80X_v2-8388a14248b3bf0d0a82d70d71bce005/USER'
    else:
        config.General.requestName = 'WMinusPiGamma_Pythia8_MINIAODSIM_94X_2016_v2'
        config.Data.inputDataset = '/WMinusPiGamma_GENSIM_80X_v2/rselvati-WMinusPiGamma_RECOSIM_80X_v2-8388a14248b3bf0d0a82d70d71bce005/USER'


if runningEra == 1:

    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_MINIAOD_2017_cfg.py'

    if doPlus:
        config.General.requestName = 'WPlusPiGamma_Pythia8_MINIAODSIM_94X_2017_v7'
        config.Data.inputDataset = '/WPlusPiGamma_GENSIM_94X_2017_v7/rselvati-WPlusPiGamma_RECOSIM_94X_2017_v7-fcfc615a65be9fb627e3afc83a7469ff/USER'
    else:
        config.General.requestName = 'WMinusPiGamma_Pythia8_MINIAODSIM_94X_2017_v7'
        config.Data.inputDataset = '/WMinusPiGamma_GENSIM_94X_2017_v7/rselvati-WMinusPiGamma_RECOSIM_94X_2017_v7-fcfc615a65be9fb627e3afc83a7469ff/USER'


config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True

if runningEra == 0:

    if doPlus:
        config.Data.outputDatasetTag = 'WPlusPiGamma_MINIAODSIM_94X_2016_v2'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_MINIAODSIM_94X_2016_v2'

if runningEra == 1:

    if doPlus:
        config.Data.outputDatasetTag = 'WPlusPiGamma_MINIAODSIM_94X_2017_v7'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_MINIAODSIM_94X_2017_v7'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'
