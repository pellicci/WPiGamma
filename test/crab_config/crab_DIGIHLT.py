from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()
 
doPlus = False
runningEra = 2 # 0 = 2016, 1 = 2017, 2 = 2018
 
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

if runningEra == 2:

    config.JobType.inputFiles = ['WPiGamma_MixingModule_2018.py'] #Import correct mixing module file
    config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_DIGIL1HLT_2018_cfg.py'

    if doPlus: 
        config.General.requestName = 'WPlusPiGamma_Pythia8_DIGIL1HLT_102X_2018_v1'
        config.Data.inputDataset = '/WPlusPiGamma_102X_2018/rselvati-WPlusPiGamma_102X_2018-716cd3531e67f611a8522846236281ff/USER'
    else:
        config.General.requestName = 'WMinusPiGamma_Pythia8_DIGIL1HLT_102X_2018_v2'
        config.Data.inputDataset = '/WMinusPiGamma_102X_2018/rselvati-WMinusPiGamma_102X_2018-1729d1fddfbbf7a02185b2ad26a98ae1/USER'


config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True

if runningEra == 0:

    if doPlus:
        config.Data.outputDatasetTag = 'WPlusPiGamma_DIGIHLT_80X_2016'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_DIGIHLT_80X_2016'

if runningEra == 1:

    if doPlus:
        config.Data.outputDatasetTag = 'WPlusPiGamma_DIGIHLT_94X_2017_v7'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_DIGIHLT_94X_2017_v7'

if runningEra == 2:

    if doPlus:
        config.Data.outputDatasetTag = 'WPlusPiGamma_DIGIHLT_102X_2018_v1'
    else:
        config.Data.outputDatasetTag = 'WMinusPiGamma_DIGIHLT_102X_2018_v2'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'
