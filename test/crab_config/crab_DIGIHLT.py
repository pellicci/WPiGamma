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

config.JobType.allowUndistributedCMSSW = True #Otherwise get an error for incompatibility of architecture(slc7)/release(9_4_10). It is safe according to https://hypernews.cern.ch/HyperNews/CMS/get/computing-tools/4935/2.html

if runningEra == 0:

    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_DIGIL1HLT_2016_cfg.py'

    if doPlus: 
        config.General.requestName = 'WPlusPiGamma_Pythia8_RECOSIM_80XV1'
        config.Data.inputDataset = '/WPlusPiGamma_GENSIM_80XV1/pellicci-WPlusPiGamma_GENSIM_80XV1-8d45dc557bfa1dd3d3626668174f8fd0/USER'
    else:
        config.General.requestName = 'WMinusPiGamma_Pythia8_RECOSIM_80XV1'
        config.Data.inputDataset = '/WMinusPiGamma_GENSIM_80XV1/pellicci-WMinusPiGamma_GENSIM_80XV1-f51ba7e770fc51cae8263bad5b777b64/USER'

if runningEra == 1:

    config.JobType.inputFiles = ['WPiGamma_MixingModule_2017.py'] #Import correct mixing module file
    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_DIGIL1HLT_2017_cfg.py'

    if doPlus: 
        config.General.requestName = 'WPlusPiGamma_Pythia8_DIGIL1HLT_94X_2017_v7'
        config.Data.inputDataset = '/WPlusPiGamma_GENSIM_94X_2017_v7/rselvati-WPlusPiGamma_GENSIM_94X_2017_v7-5da65549ab31696241a4a516c7cdf1f6/USER'
    else:
        config.General.requestName = 'WMinusPiGamma_Pythia8_DIGIL1HLT_94X_2017_v7'
        config.Data.inputDataset = '/WMinusPiGamma_GENSIM_94X_2017_v7/rselvati-WMinusPiGamma_GENSIM_94X_2017_v7-7ca3ebea458deee460a6fa996e37a2ff/USER'


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
        config.Data.outputDatasetTag = 'WPlusPiGamma_DIGIHLT_94X_2017_v7'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_DIGIHLT_94X_2017_v7'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'
