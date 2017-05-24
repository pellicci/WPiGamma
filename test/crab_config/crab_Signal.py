from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'WPiGamma_Pythia8_MINIAODSIM_signal'
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/run_WPiGammaAnalysis.py'
config.JobType.pluginName = 'Analysis'

config.section_('Data')
config.Data.inputDataset = '/WPiGamma_GENSIM_80XV1/pellicci-WPiGamma_MINIAODSIM_80XV1-0e6df83a66ae3f10341eebd4053e4881/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.JobType.outputFiles = ['WPiGammaAnalysis_output.root']

config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'
