# WPiGamma

- Create a new CMSSW release
   
   cmsrel CMSSW_8_0_28
   
   cd CMSSW_8_0_28/src
   
   cmsenv

- Get the EGamma recommended tags

   git cms-merge-topic ikrav:egm_id_80X_v3_photons

- In src, get the package

   git clone -b WPiGamma_80X https://github.com/pellicci/WPiGamma.git StandardModel/WPiGamma


- Compile
   
   scram b -j 4

