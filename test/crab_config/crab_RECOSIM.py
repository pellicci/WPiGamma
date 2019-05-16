from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()
 
doPlus = False
runningEra = 1 # 0 = 2016, 1 = 2017, 2 = 2018 

config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects'
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.section_('Data')

if runningEra == 0:

    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_RECO_cfg.py'

    if doPlus: 
        config.General.requestName = 'WPlusPiGamma_Pythia8_RECOSIM_80XV1'
        config.Data.inputDataset = '/WPlusPiGamma_GENSIM_80XV1/pellicci-WPlusPiGamma_GENSIM_80XV1-8d45dc557bfa1dd3d3626668174f8fd0/USER'
    else:
        config.General.requestName = 'WMinusPiGamma_Pythia8_RECOSIM_80XV1'
        config.Data.inputDataset = '/WMinusPiGamma_GENSIM_80XV1/pellicci-WMinusPiGamma_GENSIM_80XV1-f51ba7e770fc51cae8263bad5b777b64/USER'

if runningEra == 1:

    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_RAW2DIGIRECO_2017_cfg.py'

    if doPlus: 
        config.General.requestName = 'WPlusPiGamma_Pythia8_RECOSIM_94X'
        config.Data.inputDataset = '/WPlusPiGamma_GENSIM_80XV1/pellicci-WPlusPiGamma_GENSIM_80XV1-8d45dc557bfa1dd3d3626668174f8fd0/USER'
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
        config.Data.outputDatasetTag = 'WPlusPiGamma_RECOSIM_80XV1'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_RECOSIM_80XV1'

if runningEra == 1:

    if doPlus:
        config.Data.outputDatasetTag = 'WPlusPiGamma_RECOSIM_94X_2017'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_RECOSIM_94X_2017_v3'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'
