from CRABClient.UserUtilities import config
config = config()
 
doPlus = False
runningEra = 2 # 0 = 2016, 1 = 2017, 2 = 2018 

config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects'
config.section_('JobType')
config.JobType.pluginName = 'Analysis'

config.JobType.allowUndistributedCMSSW = True

config.section_('Data')

if runningEra == 0:

    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_RECO_cfg.py'

    if doPlus: 
        config.General.requestName = 'WPlusPiGamma_Pythia8_RECOSIM_80X_2016_v1'
        config.Data.inputDataset = '/TTToSemilepWPlusPiGamma_GENSIM_80X_2016_v1/pellicci-WPlusPiGamma_DIGIL1HLT_80X_2016_v1-6868a286d3522be326f21d858c30a51e/USER'
    else:
        config.General.requestName = 'WMinusPiGamma_Pythia8_RECOSIM_80X_2016_v1'
        config.Data.inputDataset = '/TTToSemilepWMinusPiGamma_GENSIM_80X_2016_v1/pellicci-WMinusPiGamma_DIGIL1HLT_80X_2016_v1-6868a286d3522be326f21d858c30a51e/USER'

if runningEra == 1:

    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_RAW2DIGIRECO_2017_cfg.py'

    if doPlus: 
        config.General.requestName = 'WPlusPiGamma_Pythia8_RECOSIM_94X_2017_v1'
        config.Data.inputDataset = '/TTToSemilepWPlusPiGamma_GENSIM_94X_2017_v1/pellicci-WPlusPiGamma_DIGIHLT_94X_2017_v1-80ee49574e9e2ea92e8268e286fe8a92/USER'
    else:
        config.General.requestName = 'WMinusPiGamma_Pythia8_RECOSIM_94X_2017_v1'
        config.Data.inputDataset = '/TTToSemilepWMinusPiGamma_GENSIM_94X_2017_v1/pellicci-WMinusPiGamma_DIGIHLT_94X_2017_v1-80ee49574e9e2ea92e8268e286fe8a92/USER'

if runningEra == 2:

    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_RAW2DIGIRECO_2018_cfg.py'

    if doPlus: 
        config.General.requestName = 'WPlusPiGamma_Pythia8_RECOSIM_102X_2018_v1'
        config.Data.inputDataset = '/TTToSemilepWPlusPiGamma_102X_2018_v1/pellicci-WPlusPiGamma_DIGIHLT_102X_2018_v1-561697a7e1ccc674784d5e0d3e6ef789/USER'
    else:
        config.General.requestName = 'WMinusPiGamma_Pythia8_RECOSIM_102X_2018_v1'
        config.Data.inputDataset = '/TTToSemilepWMinusPiGamma_102X_2018_v1/pellicci-WMinusPiGamma_DIGIHLT_102X_2018_v1-561697a7e1ccc674784d5e0d3e6ef789/USER'


config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.publication = True

if runningEra == 0:

    if doPlus:
        config.Data.outputDatasetTag = 'WPlusPiGamma_RECOSIM_80X_2016_v1'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_RECOSIM_80X_2016_v1'

if runningEra == 1:

    if doPlus:
        config.Data.outputDatasetTag = 'WPlusPiGamma_RECOSIM_94X_2017_v1'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_RECOSIM_94X_2017_v1'

if runningEra == 2:

    if doPlus:
        config.Data.outputDatasetTag = 'WPlusPiGamma_RECOSIM_102X_2018_v1'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_RECOSIM_102X_2018_v1'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'
