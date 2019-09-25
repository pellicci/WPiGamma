from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()
 
doPlus = False
runningEra = 0 # 0 = 2016, 1 = 2017, 2 = 2018 

config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects'
config.section_('JobType')
config.JobType.pluginName = 'Analysis'

config.JobType.allowUndistributedCMSSW = True #Otherwise get an error for incompatibility of architecture(slc7)/release(9_4_10) or (8_0_28). It is safe according to https://hypernews.cern.ch/HyperNews/CMS/get/computing-tools/4935/2.html

config.section_('Data')

if runningEra == 0:

    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_RAW2DIGIRECO_2016_cfg.py'

    if doPlus: 
        config.General.requestName = 'WPlusPiGamma_Pythia8_RECOSIM_80X_v2'
        config.Data.inputDataset = '/WPlusPiGamma_GENSIM_80X_v2/rselvati-WPlusPiGamma_DIGIL1HLT_80X_v2-6868a286d3522be326f21d858c30a51e/USER'
    else:
        config.General.requestName = 'WMinusPiGamma_Pythia8_RECOSIM_80X_v2'
        config.Data.inputDataset = '/WMinusPiGamma_GENSIM_80X_v2/rselvati-WMinusPiGamma_DIGIL1HLT_80X_v2-6868a286d3522be326f21d858c30a51e/USER'

if runningEra == 1:

    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_RAW2DIGIRECO_2017_cfg.py'

    if doPlus: 
        config.General.requestName = 'WPlusPiGamma_Pythia8_RECOSIM_94X_2017_v5'
        config.Data.inputDataset = '/WPlusPiGamma_GENSIM_94X_2017_v5/rselvati-WPlusPiGamma_DIGIHLT_94X_2017_v5-5b9cd2c7eef36524de7af1c8e43b0ebc/USER'
    else:
        config.General.requestName = 'WMinusPiGamma_Pythia8_RECOSIM_94X_2017_v3'
        config.Data.inputDataset = '/WMinusPiGamma_GENSIM_94X_2017_v3/rselvati-WMinusPiGamma_DIGIHLT_94X_2017_v3-5b9cd2c7eef36524de7af1c8e43b0ebc/USER'


config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True

if runningEra == 0:

    if doPlus:
        config.Data.outputDatasetTag = 'WPlusPiGamma_RECOSIM_80X_v2'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_RECOSIM_80X_v2'

if runningEra == 1:

    if doPlus:
        config.Data.outputDatasetTag = 'WPlusPiGamma_RECOSIM_94X_2017_v5'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_RECOSIM_94X_2017_v3'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'
