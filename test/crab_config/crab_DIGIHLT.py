from CRABClient.UserUtilities import config
config = config()
 
doPlus = False
runningEra = 2 # 0 = 2016, 1 = 2017, 2 = 2018
 
config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects'
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.section_('Data')

config.JobType.allowUndistributedCMSSW = True

if runningEra == 0:

    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_DIGIL1HLT_2016_cfg.py'

    if doPlus: 
        config.General.requestName = 'WPlusPiGamma_Pythia8_DIGIL1HLT_80X_2016_v1'
        config.Data.inputDataset = '/TTToSemilepWPlusPiGamma_GENSIM_80X_2016_v1/pellicci-TTToSemilepWPlusPiGamma_GENSIM_80X_2016_v1-1945d29c337cf9cb6c6f6a9c6060d5c4/USER'
    else:
        config.General.requestName = 'WMinusPiGamma_Pythia8_DIGIL1HLT_80X_2016_v1'
        config.Data.inputDataset = '/TTToSemilepWMinusPiGamma_GENSIM_80X_2016_v1/pellicci-TTToSemilepWMinusPiGamma_GENSIM_80X_2016_v1-6d49b0bcd77835c8ca05583901fc82de/USER'

if runningEra == 1:

    config.JobType.inputFiles = ['WPiGamma_MixingModule_2017.py'] #Import correct mixing module file
    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_DIGIL1HLT_2017_cfg.py'

    if doPlus: 
        config.General.requestName = 'WPlusPiGamma_Pythia8_DIGIL1HLT_94X_2017_v1'
        config.Data.inputDataset = '/TTToSemilepWPlusPiGamma_GENSIM_94X_2017_v1/pellicci-TTToSemilepWPlusPiGamma_GENSIM_94X_2017_v1-62827fe3595cf28cf32f246d8a66ca3d/USER'
    else:
        config.General.requestName = 'WMinusPiGamma_Pythia8_DIGIL1HLT_94X_2017_v1'
        config.Data.inputDataset = '/TTToSemilepWMinusPiGamma_GENSIM_94X_2017_v1/pellicci-TTToSemilepWMinusPiGamma_GENSIM_94X_2017_v1-7a3edcfd8f83a535f7f2253939f3a07e/USER'

if runningEra == 2:

    config.JobType.inputFiles = ['WPiGamma_MixingModule_2018.py'] #Import correct mixing module file
    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_DIGIL1HLT_2018_cfg.py'

    if doPlus: 
        config.General.requestName = 'WPlusPiGamma_Pythia8_DIGIL1HLT_102X_2018_v1'
        config.Data.inputDataset = '/TTToSemilepWPlusPiGamma_102X_2018_v1/pellicci-TTToSemilepWPlusPiGamma_GENSIM_102X_2018_v1-3e4e7186e38635d81b8617fb25a86806/USER'
    else:
        config.General.requestName = 'WMinusPiGamma_Pythia8_DIGIL1HLT_102X_2018_v1'
        config.Data.inputDataset = '/TTToSemilepWMinusPiGamma_102X_2018_v1/pellicci-TTToSemilepWMinusPiGamma_GENSIM_102X_2018_v1-9675e0ce03dbef0bfdc505437f2a0f8a/USER'


config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.publication = True

if runningEra == 0:

    if doPlus:
        config.Data.outputDatasetTag = 'WPlusPiGamma_DIGIL1HLT_80X_2016_v1'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_DIGIL1HLT_80X_2016_v1'

if runningEra == 1:

    if doPlus:
        config.Data.outputDatasetTag = 'WPlusPiGamma_DIGIHLT_94X_2017_v1'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_DIGIHLT_94X_2017_v1'

if runningEra == 2:

    if doPlus:
        config.Data.outputDatasetTag = 'WPlusPiGamma_DIGIHLT_102X_2018_v1'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_DIGIHLT_102X_2018_v1'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'
