//ROOT includes
#include <TH1F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "Math/VectorUtil.h"
#include <stdlib.h>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h" 
 
#include "WDirectProduction.h"
 
// constructors and destructor
WDirectProduction::WDirectProduction(const edm::ParameterSet& iConfig)
{
  genParticlesToken_            = consumes<std::vector<reco::GenParticle> >(edm::InputTag("genParticles"));

  create_trees();
}

WDirectProduction::~WDirectProduction()
{
}

// ------------ method called for each event  ------------
void WDirectProduction::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<reco::GenParticle>  > genParticles;
  iEvent.getByToken(genParticlesToken_, genParticles);

  int genW_ID    = -999;
  float genW_pT  = -999.;
  float genW_eta = -999.;
  float genW_phi = -999.;
  float genW_E   = -999.;

  int genPh_ID    = -999;
  float genPh_pT  = -999.;
  float genPh_eta = -999.;
  float genPh_phi = -999.;
  float genPh_E   = -999.;

  int genPi_ID    = -999;
  float genPi_pT  = -999.;
  float genPi_eta = -999.;
  float genPi_phi = -999.;
  float genPi_E   = -999.;

  genW_ID_tree  = 0;
  genW_pT_tree  = 0.;
  genW_eta_tree = 0.;
  genW_phi_tree = 0.;
  genW_E_tree   = 0.;

  genPh_ID_tree  = 0;
  genPh_pT_tree  = 0.;
  genPh_eta_tree = 0.;
  genPh_phi_tree = 0.;
  genPh_E_tree   = 0.;

  genPi_ID_tree  = 0;
  genPi_pT_tree  = 0.;
  genPi_eta_tree = 0.;
  genPi_phi_tree = 0.;
  genPi_E_tree   = 0.;


  for(auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){


    if(fabs(gen->pdgId()) == 24 && gen->numberOfDaughters() == 2){


      for(int i = 0; i < 2; i++){
	
	if(!(fabs(gen->daughter(i)->pdgId()) == 22 || fabs(gen->daughter(i)->pdgId()) == 211)) continue;
  	genW_ID  = gen->pdgId();
  	genW_pT  = gen->pt();
  	genW_eta = gen->eta();
  	genW_phi = gen->phi();
  	genW_E   = gen->energy();
	
  	if(fabs(gen->daughter(i)->pdgId()) == 22){
  	  genPh_ID  = gen->daughter(i)->pdgId();
  	  genPh_pT  = gen->daughter(i)->pt();
  	  genPh_eta = gen->daughter(i)->eta();
  	  genPh_phi = gen->daughter(i)->phi();
  	  genPh_E   = gen->daughter(i)->energy();
  	} else {
  	  genPi_ID  = gen->daughter(i)->pdgId();
  	  genPi_pT  = gen->daughter(i)->pt();
  	  genPi_eta = gen->daughter(i)->eta();
  	  genPi_phi = gen->daughter(i)->phi();
  	  genPi_E   = gen->daughter(i)->energy();
  	}	  
      }
    }
    
    // if(fabs(gen->pdgId()) == 24 && gen->numberOfDaughters() == 1 && fabs(gen->daughter(0)->pdgId()) == 24){
      
    //   for(int i = 0; i < 2; i++){
	
    // 	if(!(fabs(gen->daughter(0)->daughter(i)->pdgId()) == 22 || fabs(gen->daughter(0)->daughter(i)->pdgId()) == 211)) continue;
    // 	genW_ID  = gen->pdgId();
    // 	genW_pT  = gen->pt();
    // 	genW_eta = gen->eta();
    // 	genW_phi = gen->phi();
    // 	genW_E   = gen->energy();
    // 	std::cout << "ciao" << std::endl;
    // 	if(fabs(gen->daughter(0)->daughter(i)->pdgId()) == 22){
    // 	  genPh_ID  = gen->daughter(0)->daughter(i)->pdgId();
    // 	  genPh_pT  = gen->daughter(0)->daughter(i)->pt();
    // 	  genPh_eta = gen->daughter(0)->daughter(i)->eta();
    // 	  genPh_phi = gen->daughter(0)->daughter(i)->phi();
    // 	  genPh_E   = gen->daughter(0)->daughter(i)->energy();
    // 	} else {
    // 	  genPi_ID  = gen->daughter(0)->daughter(i)->pdgId();
    // 	  genPi_pT  = gen->daughter(0)->daughter(i)->pt();
    // 	  genPi_eta = gen->daughter(0)->daughter(i)->eta();
    // 	  genPi_phi = gen->daughter(0)->daughter(i)->phi();
    // 	  genPi_E   = gen->daughter(0)->daughter(i)->energy();
    // 	}	  
    //   }
    // }
  }
  
  genW_ID_tree  = genW_ID;
  genW_pT_tree  = genW_pT;
  genW_eta_tree = genW_eta;
  genW_phi_tree = genW_phi;
  genW_E_tree   = genW_E;
  
  genPh_ID_tree  = genPh_ID;
  genPh_pT_tree  = genPh_pT;
  genPh_eta_tree = genPh_eta;
  genPh_phi_tree = genPh_phi;
  genPh_E_tree   = genPh_E;
  
  genPi_ID_tree  = genPi_ID;
  genPi_pT_tree  = genPi_pT;
  genPi_eta_tree = genPi_eta;
  genPi_phi_tree = genPi_phi;
  genPi_E_tree   = genPi_E;
  
  mytree->Fill();
}

//*************************************************************//
//                                                             //
//---------------------- Create the tree ----------------------//
//                                                             //
//*************************************************************//

void WDirectProduction::create_trees()
{
  mytree = fs->make<TTree>("mytree", "Tree containing gen info");

  mytree->Branch("genW_ID",&genW_ID_tree);
  mytree->Branch("genW_pT",&genW_pT_tree);
  mytree->Branch("genW_eta",&genW_eta_tree);
  mytree->Branch("genW_phi",&genW_phi_tree);
  mytree->Branch("genW_E",&genW_E_tree);
  
  mytree->Branch("genPh_ID",&genPh_ID_tree);
  mytree->Branch("genPh_pT",&genPh_pT_tree);
  mytree->Branch("genPh_eta",&genPh_eta_tree);
  mytree->Branch("genPh_phi",&genPh_phi_tree);
  mytree->Branch("genPh_E",&genPh_E_tree);
  
  mytree->Branch("genPi_ID",&genPi_ID_tree);
  mytree->Branch("genPi_pT",&genPi_pT_tree);
  mytree->Branch("genPi_eta",&genPi_eta_tree);
  mytree->Branch("genPi_phi",&genPi_phi_tree);
  mytree->Branch("genPi_E",&genPi_E_tree);
    
}

void WDirectProduction::beginJob()
{
}

void WDirectProduction::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(WDirectProduction);
