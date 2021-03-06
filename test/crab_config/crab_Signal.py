from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'WPiGammaAnalysis_Signal'
config.General.workArea = 'crab_projects/samples'

config.section_('JobType')
config.JobType.psetName = 'cmssw_config/run_WPiGammaAnalysis.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['MCpileUp_25ns_Recent2016.root','pileUpHistogramFromjson_Nominal.root' ] #data files for PileUp reweighting
config.JobType.outputFiles = ['WPiGammaAnalysis_output.root']

config.JobType.allowUndistributedCMSSW = True #Otherwise get an error for incompatibility of architecture(slc7)/release(9_4_10). It is safe according to https://hypernews.cern.ch/HyperNews/CMS/get/computing-tools/4935/2.html

config.section_('Data')
config.Data.inputDataset = '/WPiGamma_GENSIM_80XV1/pellicci-WPiGamma_MINIAODSIM_80XV1-0e6df83a66ae3f10341eebd4053e4881/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False

config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'
