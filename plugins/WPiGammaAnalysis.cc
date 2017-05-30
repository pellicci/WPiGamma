//ROOT includes
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
 
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
  
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
 
// vertex inclusions
#include "DataFormats/VertexReco/interface/Vertex.h" 
#include "DataFormats/BeamSpot/interface/BeamSpot.h" 
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
 
typedef math::XYZTLorentzVector LorentzVector;
 
#include "WPiGammaAnalysis.h"
 
// constructors and destructor
WPiGammaAnalysis::WPiGammaAnalysis(const edm::ParameterSet& iConfig) :
  packedPFCandidates_(iConfig.getParameter<edm::InputTag>("packedPFCandidates")),
  slimmedMuons_(iConfig.getParameter<edm::InputTag>("slimmedMuons")), 
  prunedGenParticles_(iConfig.getParameter<edm::InputTag>("prunedGenParticles")),
  slimmedPhotons_(iConfig.getParameter<edm::InputTag>("slimmedPhotons")),
  slimmedElectrons_(iConfig.getParameter<edm::InputTag>("slimmedElectrons")),
  slimmedJets_(iConfig.getParameter<edm::InputTag>("slimmedJets")),
  runningOnData_(iConfig.getParameter<bool>("runningOnData")),
  pvCollection_(iConfig.getParameter<edm::InputTag>("pvCollection")),   
  bsCollection_(iConfig.getParameter<edm::InputTag>("bsCollection")),  
  PileupSrc_(iConfig.getParameter<edm::InputTag>("PileupSrc"))
{
  packedPFCandidatestoken_  = consumes<std::vector<pat::PackedCandidate> >(packedPFCandidates_); 
  slimmedMuonstoken_        = consumes<std::vector<pat::Muon> >(slimmedMuons_);
  prunedGenParticlestoken_  = consumes<std::vector<reco::GenParticle> >(prunedGenParticles_);
  slimmedPhotonstoken_      = consumes<std::vector<pat::Photon> >(slimmedPhotons_);
  slimmedElectronstoken_    = consumes<std::vector<pat::Electron> >(slimmedElectrons_);
  slimmedJetstoken_         = consumes<std::vector<pat::Jet> >(slimmedJets_);
  tok_Vertex_               = consumes<std::vector<reco::Vertex> > (pvCollection_);  
  tok_beamspot_             = consumes<reco::BeamSpot> (edm::InputTag(bsCollection_));
  pileupSummaryToken_       = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag(PileupSrc_));

  nevent = 0;

  cand_total = 0;
  cand_passing_selection = 0;
  events_least_one_pi = 0; //number of events with at least one reconstructed pi candidate
  single_pi_counter = 0;
  pi_from_w_correction = 0;
  mu_selection = 0;
  mu_selection_event = 0;
  ph_total = 0;
  ph_passing_selection = 0;
  events_least_one_ph = 0;
  pi_from_w = 0;
  photon_from_w = 0;
  pi_and_photon_from_w = 0;
  mu_per_event = 0;
  el_per_event = 0;

  //WPiGammaAnalysis_output = fs->make<TFile>("WPiGammaAnalysis_output.root","recreate");
  inv_mass_1 = fs->make<TH1F>("Mw - no match with MC Truth", "Mw no match", 200,0,120);
  inv_mass_2 = fs->make<TH1F>("Mw - match with MC Truth", "Mw match", 200,0,120);
  track_iso_hist = fs->make<TH1F>("Track iso", "Track isolation", 200,0,60);
  ecal_iso_hist = fs->make<TH1F>("Ecal iso", "Ecal isolation", 200,0,40);
  hcal_iso_hist = fs->make<TH1F>("Hcal iso", "Hcal isolation", 200,0,25);
  calo_iso_hist = fs->make<TH1F>("Calo iso", "Calo isolation", 200,0,140);
  iso_sum_hist = fs->make<TH1F>("Sum iso", "Sum isolation", 150,0,3);

  create_trees();
}

WPiGammaAnalysis::~WPiGammaAnalysis()
{
}

// ------------ method called for each event  ------------
void WPiGammaAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<pat::PackedCandidate>  > PFCandidates;
  iEvent.getByLabel(packedPFCandidates_, PFCandidates);

  edm::Handle<std::vector<pat::Muon>  > slimmedMuons;
  iEvent.getByLabel(slimmedMuons_, slimmedMuons);

  edm::Handle<std::vector<reco::GenParticle>  > genParticles;
  if(!runningOnData_)iEvent.getByLabel(prunedGenParticles_, genParticles);

  edm::Handle<std::vector<pat::Photon>  > slimmedPhotons;
  iEvent.getByLabel(slimmedPhotons_, slimmedPhotons);

  edm::Handle<std::vector<pat::Electron > > slimmedElectrons;
  iEvent.getByLabel(slimmedElectrons_, slimmedElectrons);

  edm::Handle<std::vector<pat::Jet > > slimmedJets;
  iEvent.getByLabel(slimmedJets_, slimmedJets);

  edm::Handle<std::vector<reco::Vertex > > slimmedPV;
  iEvent.getByLabel(pvCollection_, slimmedPV);


  nevent += 1;

  float pTmuMin = -1000.;
  float pTpiMax = -1000.;
  float eTphMax = -1000.;

  mu_per_event_tree = 0;
  el_per_event_tree = 0;
  mu_per_event = 0;
  el_per_event = 0;

  float mu_pT = 0.;
  float mu_eta = 0.;
  float mu_phi = 0.;
  int   mu_ID = 0;

  float el_pT = 0.;
  float el_eta = 0.;
  float el_phi = 0.;
  int   el_ID = 0;

  float pi_pT = 0.;
  float pi_eta = 0.;
  float pi_phi = 0.;
  LorentzVector pi_p4;

  float ph_pT = 0.;
  float ph_eta = 0.;
  float ph_phi = 0.;
  LorentzVector ph_p4;

  bool tag_lepton_found = false;
  bool is_muon = false;
  bool is_signal_mu = false;
  bool is_signal = false;

  //These variables will go in the tree
  lepton_pT_tree = 0.;
  lepton_eta_tree = 0.;
  lepton_phi_tree = 0.;

  pi_pT_tree = 0.;
  pi_eta_tree = 0.;
  pi_phi_tree = 0.;

  photon_eT_tree = 0.;
  photon_eta_tree = 0.;
  photon_phi_tree = 0.;

  is_signal_mu_tree = false;
  is_signal_tree = false;
  is_muon_tree = false;
  electron_over_right_mu_tree = false;

  //--------- Asking MC truth if tag muon or tag electron. Then finding mu and electron multiplicity --------
  /* if(!runningOnData_){
     for(auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){
     if((gen->pdgId()==13 || gen->pdgId()==-13) && (gen->mother()->pdgId()==24 || gen->mother()->pdgId()==-24)){
     is_signal_mu = true;}
     if(gen->numberOfDaughters()==2){
     if((gen->pdgId()==24 || gen->pdgId()==-24) && ((gen->daughter(0)->pdgId()==211 && gen->daughter(1)->pdgId()==22) || (gen->daughter(1)->pdgId()==211 && gen->daughter(0)->pdgId()==22) || (gen->daughter(0)->pdgId()==-211 && gen->daughter(1)->pdgId()==22) || (gen->daughter(1)->pdgId()==-211 && gen->daughter(0)->pdgId()==22))){
     is_signal = true;}//actually useless, since there's a WPiGamma decay for each event
     }
     }
     }

     for (auto mu = slimmedMuons->begin(); mu != slimmedMuons->end(); ++mu){
     if (mu->pt()>24. && mu->pt() > pTmuMin && mu->isMediumMuon() ){
     pTmuMin = mu->pt();
     //std::cout << "mu pT :" << mu->pt() << "Eta: " << mu->eta() << "phi:" << mu->phi() << std::endl;

     mu_ID = mu->pdgId();
     is_muon = true;
     mu_selection += 1; //checking how many mu over total pass the selection
     mu_per_event += 1;
     //mu_eta = mu->eta();
     //mu_phi = mu->phi();
     //mu_pT = mu->pt();
     if (tag_lepton_found == false){//checking how many events contain at least one mu passing the selection
     mu_selection_event += 1;}
     tag_lepton_found = true;

     lepton_pT_tree = mu->pt();
     lepton_eta_tree = mu->eta();
     lepton_phi_tree = mu->phi();

     }
     }


     for (auto el = slimmedElectrons->begin(); el != slimmedElectrons->end(); ++el){
     if (el->pt()>24){
     el_ID = el->pdgId();
     is_muon = false;
     //mu_selection += 1; //checking how many mu over total pass the selection
     //if (checker == false){ //checking how many events contain mu passing the if selection
     //mu_selection_event += 1;}
     //tag_lepton_found = true;
     el_per_event +=1;
     //tag_el_hist->Fill(el);
     lepton_pT_tree = el->pt();
     lepton_eta_tree = el->eta();
     lepton_phi_tree = el->phi();
     //std::cout << "pT del leptone" << lepton_pT_tree << std::endl;
     }
     }

     mu_per_event_tree = mu_per_event;
     el_per_event_tree = el_per_event;
     if(is_signal_mu==true){is_signal_mu_tree = true;}
     if(is_signal==true){is_signal_tree = true;}*/


  //----------- Event reconstruction -----------

  for(auto mu = slimmedMuons->begin(); mu != slimmedMuons->end(); ++mu){
    if(mu->pt() < 25. || mu->pt() < pTmuMin || !mu->isTightMuon(slimmedPV->at(0))) continue;
    if( (mu->chargedHadronIso() + std::max(0., mu->neutralHadronIso() + mu->photonIso() - 0.5*mu->puChargedHadronIso())/mu->pt()) > 0.2) continue;
    pTmuMin = mu->pt();
    //std::cout << "mu pT :" << mu->pt() << "Eta: " << mu->eta() << "phi:" << mu->phi() << std::endl;

    mu_ID         = mu->pdgId();
    is_muon       = true;
    mu_selection += 1; //checking how many mu over total pass the selection
    mu_per_event += 1;
    mu_eta        = mu->eta();
    mu_phi        = mu->phi();
    mu_pT         = mu->pt();
    if(tag_lepton_found == false) mu_selection_event += 1; //checking how many events contain at least one mu passing the selection
    tag_lepton_found = true;
  }

  //---------- Comparison between genparticles and reco mu -----------
  /* 
     if (tag_lepton_found == true){
     if(!runningOnData_){
     for(auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){
     deltaR = TMath::Sqrt((mu_eta-gen->eta())*(mu_eta-gen->eta())+(mu_phi-gen->phi())*(mu_phi-gen->phi()));
     deltapT = TMath::Abs(mu_pT-gen->pt());

     if (deltaR <= deltaRMax && deltapT < deltapTMax){
     deltapTMax = deltapT;
     gen_ID = gen->pdgId();
     gen_mother = gen->mother()->pdgId();
     }
     }
     if((gen_ID == 13 || gen_ID == -13) && (gen_mother == 24 || gen_mother == -24)){
     is_signal_mu = true;
	  
     }
     }
     }*/


  for(auto el = slimmedElectrons->begin(); el != slimmedElectrons->end(); ++el){
    if(el->pt() < 26. || el->pt() < pTmuMin || el->trackIso() > 3. || el->caloIso() > 5.) continue;
    el_ID   = el->pdgId();
    is_muon = false;
    el_ID   = el->pdgId();
    el_eta  = el->eta();
    el_phi  = el->phi();
    el_pT   = el->pt();
    el_per_event += 1;
    tag_lepton_found      = true;
  }

  if(tag_lepton_found && is_muon){
    lepton_pT_tree  = mu_pT;
    lepton_eta_tree = mu_eta;
    lepton_phi_tree = mu_phi;
  }

  if(tag_lepton_found && !is_muon){
    lepton_pT_tree  = el_pT;
    lepton_eta_tree = el_eta;
    lepton_phi_tree = el_phi;
  }

  if(is_muon) is_muon_tree = true;
  else is_muon_tree = false;

  //Do NOT continue if you didn't find either a muon or an electron
  if(!tag_lepton_found) return;

  // mu_per_event_tree = mu_per_event;
  //el_per_event_tree = el_per_event;
  //if(is_signal_mu==true){is_signal_mu_tree = true;}
  //if(is_signal==true){is_signal_tree = true;}
  //if(is_signal_mu==true && in_electron_selection==true){electron_over_right_mu_tree=true;}

  //----------- Starting to search for pi and gamma -------------

  cand_passing_selection = 0;
  bool cand_pion_found  = false;
  bool is_pi_a_pi       = false; //together with is_photon_a_photon, used to calculate how many times both reconstructed pi and gamma come from generated particles whose mother is a W
  bool is_last_pi_a_pi  = false; //used if only the last pi to pass the selection is actually a pi in generation too
  bool is_bad_single_pi = false; //used to avoid to count the case in which only a single reco-pion, not matching with a gen-pion, passes the selection
  single_pi_counter = 0;

  for(auto cand = PFCandidates->begin(); cand != PFCandidates->end(); ++cand){
    if(cand->pdgId()*mu_ID < 0 && cand->pdgId()*el_ID < 0) continue;    
    if(!tag_lepton_found || cand->pt() < 20. || !cand->trackHighPurity() || cand->fromPV() != 3) continue;
    if(cand->pt() < pTpiMax) continue;

    pTpiMax = cand->pt();
    cand_pion_found = true;

    pi_pT  = cand->pt();
    pi_eta = cand->eta();
    pi_phi = cand->phi();
    pi_p4  = cand->p4();

    float deltapTMax = 10000.;
    const float deltaRMax  = 0.3;
    int   gen_mother = 0;
    int   gen_ID = 0;

    //gen_ID = cand.pdgId()#so if it doesn't enter the following loop-and-if, it doesn't display "reconstructed particle is different..." either 
    if(!runningOnData_){
      for (auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){//matching candidate for W reconstruction with MC truth
	float deltaR = sqrt((pi_eta-gen->eta())*(pi_eta-gen->eta())+(pi_phi-gen->phi())*(pi_phi-gen->phi()));
	float deltapT = fabs(pi_pT-gen->pt());

	if(deltaR > deltaRMax || deltapT > deltapTMax) continue;
	deltapTMax = deltapT;
	gen_ID = gen->pdgId();
	gen_mother = gen->mother()->pdgId();
      }
    }

    if(fabs(gen_ID) != 211 && single_pi_counter == 0){
      //std::cout << "|||corresponding generated particle is not a pion, but: " << gen_ID << std::endl;
      is_bad_single_pi = true; //avoiding to count the case in which only a single reco-pion, not matching with a gen-pion, passes the selection
    } 

    if(fabs(gen_ID) == 211){
      //std::cout << "\\identity of generated PION's mother: " << gen_mother << std::endl;
      //std::cout << "reco px: " << pxpi << "reco py: " << pypi << "reco pz: " << pzpi << "reco pT: " << pTpi << std::endl;
      //std::cout << "gen px: " << gen_px << "gen py: " << gen_py << "gen pz: " << gen_pz << "gen pT: " << gen_pt << std::endl;
      if(fabs(gen_mother) == 24){
	pi_from_w += 1;
	is_pi_a_pi = true;
	gen_mother = 0;
	gen_ID = 0;
	single_pi_counter += 1;
	if(cand_passing_selection >= 1){
	  pi_from_w += 1;
	  is_last_pi_a_pi = true;} //activates in case the first reco particle to fulfill the conditions is not a pion from a W, while the second is
      }
    }
      
    if (single_pi_counter >= 1 && cand_passing_selection >= 1){
      pi_from_w = pi_from_w-1; //without this, a case in which the first particle to pass the selection is a pion but the second and final is not, is counted as a good reconstruction 
      is_pi_a_pi = false;}
    if(is_last_pi_a_pi) is_pi_a_pi = true;
                
    cand_passing_selection += 1;
  }

  if(cand_passing_selection != 0) events_least_one_pi += 1;

  //Do NOT continue if you didn't find a pion
  if(!cand_pion_found) return;

  cand_total += cand_passing_selection;

  pi_pT_tree  = pi_pT;
  pi_eta_tree = pi_eta;
  pi_phi_tree = pi_phi;

  ph_passing_selection = 0;
  bool cand_photon_found = false;
  bool is_photon_a_photon = false; //together with is_pi_a_pi, used to calculate how many times both reconstructed pi and gamma come from generated particles whose mother is a W

  for (auto photon = slimmedPhotons->begin(); photon != slimmedPhotons->end(); ++photon){

    if(photon->et() < 20. || photon->et() < eTphMax) continue;

    eTphMax = photon->et();
    //std::cout << "px: " << photon->px() << "py: " << photon->py() << "pz: " << photon->pz() << "energy: " << photon->energy() << std::endl;

    ph_pT  = photon->pt();
    ph_eta = photon->eta();
    ph_phi = photon->phi();
    ph_p4  = photon->p4();

    cand_photon_found = true;
    ph_passing_selection +=1;

    float track_iso = photon->trackIso();
    float ecal_iso = photon->ecalIso();
    float hcal_iso = photon->hcalIso();
    float calo_iso = photon->caloIso();
    float iso_sum = (track_iso+calo_iso)/photon->et();

    track_iso_hist->Fill(track_iso);
    ecal_iso_hist->Fill(ecal_iso);
    hcal_iso_hist->Fill(hcal_iso);
    calo_iso_hist->Fill(calo_iso);
    iso_sum_hist->Fill(iso_sum);

    float deltapTMax = 10000.;
    const float deltaRMax = 0.3;
    int   gen_mother = 0;
    int   gen_ID = 0;

    //gen_ID = cand.pdgId()#so if it doesn't enter the following loop-and-if, it doesn't display "reconstructed particle is different..." either 
    if(!runningOnData_){
      for (auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){
	float deltaR = sqrt((ph_eta-gen->eta())*(ph_eta-gen->eta())+(ph_phi-gen->phi())*(ph_phi-gen->phi()));
	float deltapT = fabs(ph_pT-gen->pt());

	if(deltaR > deltaRMax || deltapT > deltapTMax) continue;
	deltapTMax = deltapT;
	gen_ID = gen->pdgId();
	gen_mother = gen->mother()->pdgId();
      }
    }
                
    if(gen_ID != 22){
      //std::cout << "|||corresponding generated particle is not a photon, but: " << gen_ID << std::endl;
      //std::cout <<  "______reco px: " << pxph << "reco py: " << pyph << "reco pz: " << pzph << "reco pT: "<< pTph << std::endl;
      //std::cout << "gen px: " << gen_px << "gen py: " << gen_py << "gen pz: " << gen_pz << "gen pT: " << gen_pt << std::endl;
    }
            
    if(gen_ID == 22){
      //std::cout << "\\identity of generated photon's mother: " << gen_mother << std::endl;
      //std::cout << "______reco px: " << pxph << "reco py: " << pyph << "reco pz: " << pzph << "reco pT: " << pTph << std::endl;
      //std::cout << "gen px: " << gen_px << "gen py: " << gen_py << "gen pz: " << gen_pz << "gen pT: " << gen_pt << std::endl;
      if( fabs(gen_mother) == 24){
	photon_from_w += 1;
	is_photon_a_photon = true;
      }
    }
  }

  if(cand_pion_found && !cand_photon_found && !is_bad_single_pi) pi_from_w_correction += 1;

  if(ph_passing_selection != 0) events_least_one_ph += 1; //incrementing the number of events with at least one photon

  ph_total += ph_passing_selection;

  photon_eT_tree =  ph_pT;
  photon_eta_tree = ph_eta;
  photon_phi_tree = ph_phi;

  //Do not continue if there's no pions or kaons    
  if(!cand_pion_found || !cand_photon_found) return;
  //std::cout << "INV MASS: " << (candidate_ph + candidate_pi).M() << "costheta: " << (pxpi*pxph+pypi*pyph+pzpi*pzph)/(TMath::Sqrt(pxpi*pxpi+pypi*pypi+pzpi*pzpi)*TMath::Sqrt(pxph*pxph+pyph*pyph+pzph*pzph));
  if (!is_pi_a_pi || !is_photon_a_photon){
    inv_mass_1->SetLineColor(3);
    inv_mass_1->Fill((pi_p4 + ph_p4).M());
  }
  if(is_pi_a_pi && is_photon_a_photon){
    pi_and_photon_from_w += 1;
    inv_mass_2->SetLineColor(2);
    inv_mass_2->Fill((pi_p4 + ph_p4).M());
  }

  n_bjets = 0;
  for (auto jet = slimmedJets->begin(); jet != slimmedJets->end(); ++jet){
    if(jet->pt() < 30.) continue;
    if(jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") < 0.46) continue;   //0.46 = loose
    n_bjets++;
  }

  mytree->Fill();

  //std::cout << "n fotoni: " << events_least_one_ph << std::endl;
}

void WPiGammaAnalysis::create_trees(){
  mytree = fs->make<TTree>("mytree", "Tree containing gen&reco");
  mytree->Branch("tag_lepton_pT",&lepton_pT_tree);
  mytree->Branch("tag_lepton_eta",&lepton_eta_tree);
  mytree->Branch("tag_lepton_phi",&lepton_phi_tree);
  mytree->Branch("is_muon",&is_muon_tree);
  mytree->Branch("cand_pi_pT",&pi_pT_tree);
  mytree->Branch("cand_pi_eta",&pi_eta_tree);
  mytree->Branch("cand_pi_phi",&pi_phi_tree);
  mytree->Branch("photon_eT",&photon_eT_tree);
  mytree->Branch("photon_eta",&photon_eta_tree);
  mytree->Branch("photon_phi",&photon_phi_tree);
  mytree->Branch("n_bjets",&n_bjets);
}

//define this as a plug-in
DEFINE_FWK_MODULE(WPiGammaAnalysis);
