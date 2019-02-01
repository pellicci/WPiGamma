from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

doPlus = False

config.section_('General')
config.General.transferOutputs = True
if doPlus:
    config.General.requestName = 'WPlusPiGamma_Pythia8_MINIAODSIM_94XV3'
else:
    config.General.requestName = 'WMinusPiGamma_Pythia8_MINIAODSIM_94XV3'

config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_MINIAOD_cfg.py'
config.JobType.pluginName = 'Analysis'

config.section_('Data')

if doPlus:
    config.Data.inputDataset = '/WPlusPiGamma_GENSIM_80XV1/rselvati-WPlusPiGamma_RECOSIM_80XV1-c536d85e5d9fce8caa236321c5af92c3/USER'
else:
    config.Data.inputDataset = '/WMinusPiGamma_GENSIM_80XV1/rselvati-WMinusPiGamma_RECOSIM_80XV1-c536d85e5d9fce8caa236321c5af92c3/USER'

config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True

if doPlus:
    config.Data.outputDatasetTag = 'WPlusPiGamma_MINIAODSIM_94XV3'
else:
    config.Data.outputDatasetTag = 'WMinusPiGamma_MINIAODSIM_94XV3'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'
