# WPiGamma

- Create a new CMSSW release
   
   cmsrel CMSSW_10_2_10
   
   cd CMSSW_10_2_10/src
   
   cmsenv

- Get the EGamma recommended tags (https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2)

   git cms-init
   
   git cms-merge-topic cms-egamma:PhotonIDValueMapSpeedup1029

   git cms-merge-topic cms-egamma:slava77-btvDictFix_10210

   git cms-addpkg EgammaAnalysis/ElectronTools

   rm EgammaAnalysis/ElectronTools/data -rf

   git clone https://github.com/cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data

- In src, get the package

   git clone https://github.com/pellicci/WPiGamma.git StandardModel/WPiGamma


- Compile
   
   scram b -j 8

