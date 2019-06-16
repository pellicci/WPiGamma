# WPiGamma

- Create a new CMSSW release
   
   cmsrel CMSSW_9_4_10
   
   cd CMSSW_9_4_10/src
   
   cmsenv

- Get the EGamma recommended tags

   git cms-merge-topic cms-egamma:EgammaPostRecoTools_940

- In src, get the package

   git clone https://github.com/pellicci/WPiGamma.git StandardModel/WPiGamma


- Compile
   
   scram b -j 4

