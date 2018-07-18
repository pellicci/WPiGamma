from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()
 
doPlus = False
 
config.section_('General')
config.General.transferOutputs = True
if doPlus: 
    config.General.requestName = 'WPlusPiGamma_Pythia8_RECOSIM_80XV1'
else:
    config.General.requestName = 'WMinusPiGamma_Pythia8_RECOSIM_80XV1'
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/WPiGamma_13TeV_pythia8_RECO_cfg.py'
config.JobType.pluginName = 'Analysis'

config.section_('Data')
if doPlus:
    config.Data.inputDataset = '/WPlusPiGamma_GENSIM_80XV1/pellicci-WPlusPiGamma_GENSIM_80XV1-8d45dc557bfa1dd3d3626668174f8fd0/USER'
else:
    config.Data.inputDataset = '/WMinusPiGamma_GENSIM_80XV1/pellicci-WMinusPiGamma_GENSIM_80XV1-f51ba7e770fc51cae8263bad5b777b64/USER'

config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True

if doPlus:
    config.Data.outputDatasetTag = 'WPlusPiGamma_RECOSIM_80XV1'
else:
    config.Data.outputDatasetTag = 'WMinusPiGamma_RECOSIM_80XV1'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'
