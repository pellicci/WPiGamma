// system include files
#include <memory>
#include <algorithm>

//ROOT includes
#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TMath.h>

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
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
  runningOnData_(iConfig.getParameter<bool>("runningOnData")),
  pvCollection_(iConfig.getParameter<edm::InputTag>("pvCollection")),   
  bsCollection_(iConfig.getParameter<edm::InputTag>("bsCollection")),  
  PileupSrc_(iConfig.getParameter<edm::InputTag>("PileupSrc"))
{
  packedPFCandidatestoken_  = consumes<std::vector<pat::PackedCandidate> >(packedPFCandidates_); 
  slimmedMuonstoken_  = consumes<std::vector<pat::Muon> >(slimmedMuons_);
  prunedGenParticlestoken_  = consumes<std::vector<reco::GenParticle> >(prunedGenParticles_);
  slimmedPhotonstoken_  = consumes<std::vector<pat::Photon> >(slimmedPhotons_);
  slimmedElectronstoken_ = consumes<std::vector<pat::Electron> >(slimmedElectrons_);
  tok_Vertex_         = consumes<reco::VertexCollection>(pvCollection_);  
  tok_beamspot_       = consumes<reco::BeamSpot>(edm::InputTag(bsCollection_));
  pileupSummaryToken_ = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag(PileupSrc_));

  cand_total = 0;
  cand_passing_selection = 0;
  events_least_one_pi = 0; //number of events with at least one reconstructed pi candidate
  single_pi_counter = 0;
  pi_from_w_correction = 0;
  mu_selection = 0;
  mu_selection_event = 0;
  mu_ID = 0;
  el_ID = 0;
  gen_ID = 0;
  ph_total = 0;
  ph_passing_selection = 0;
  events_least_one_ph = 0;
  gen_mother = 0;
  pi_from_w = 0;
  photon_from_w = 0;
  pi_and_photon_from_w = 0;


  pTpi = 0.;
  pTpiMax = -1000.;
  pTgamma = 0.;
  pTe = 0.; 
  pTmu = 0.; 
  deta = 0.;
  dphi = 0.;
  etapi = 0.;
  phipi = 0.;
  deltaRMax = 10000.;
  deltapTMax = 10000.;
  pTmuMin = -1000.;
  eTphMin = -1000.;
  pTpiMin = -1000.;
  candidate_eta = 0.;
  candidate_phi = 0.;
  photon_eta = 0.;
  photon_phi = 0.;
  pxpi = 0.;
  pypi = 0.;
  pzpi = 0.;
  pxph = 0.;
  pyph = 0.;
  pzph = 0.;
  pTph = 0.;
  track_iso = 0.; 
  ecal_iso = 0.; 
  hcal_iso = 0.;
  calo_iso = 0.; 
  iso_sum = 0.;
  gen_px = 0.;
  gen_py = 0.;
  gen_pz = 0.;
  gen_pt = 0.;
  deltaR = 0.;
  deltapT = 0.;

  //useful bools
  tag_lepton_found = false;
  cand_pion_found = false;
  cand_photon_found = false;
  is_pi_a_pi = false;
  is_photon_a_photon = false;

  //bools used for a correct implementation of counters 
  is_last_pi_a_pi = false; //used if only the last pi to pass the selection is actually a pi in generation too
  is_bad_single_pi = false; //used to avoid to count the case in which only a single reco-pion, not matching with a gen-pion, passes the selection

  //tree variables
  lepton_pT_tree = 0.;
  lepton_eta_tree = 0.;
  lepton_phi_tree = 0.;
  is_muon = false;

  //WPiGammaAnalysis_output = fs->make<TFile>("WPiGammaAnalysis_output.root","recreate");
  inv_mass_1 = fs->make<TH1F>("Mw - no match with MC Truth", "Mw no match", 200,0,120);
  inv_mass_2 = fs->make<TH1F>("Mw - match with MC Truth", "Mw match", 200,0,120);
  track_iso_hist = fs->make<TH1F>("Track iso", "Track isolation", 200,0,60);
  ecal_iso_hist = fs->make<TH1F>("Ecal iso", "Ecal isolation", 200,0,40);
  hcal_iso_hist = fs->make<TH1F>("Hcal iso", "Hcal isolation", 200,0,25);
  calo_iso_hist = fs->make<TH1F>("Calo iso", "Calo isolation", 200,0,140);
  iso_sum_hist = fs->make<TH1F>("Sum iso", "Sum isolation", 150,0,3);

  create_trees();

  /*Create the histograms and let TFileService handle them
    create_Histos_and_Trees();*/
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

<<<<<<< HEAD
  
  pTmuMin = -1000;
  tag_lepton_found = false;
=======
  //std::cout<< "Event n: " << nevents << "/" << maxEvents << std::endl;
  
  pTmuMin = -1000.;
  checker = false;
>>>>>>> 8fd4d7d7eec1bebc49964549b341f61e8cac8642
  lepton_pT_tree = 0.;
  lepton_eta_tree = 0.;
  lepton_phi_tree = 0.;
  mu_ID = 0;
  el_ID = 0;

  for (auto mu = slimmedMuons->begin(); mu != slimmedMuons->end(); ++mu){
        if (mu->pt()>24. && mu->pt() > pTmuMin && mu->isMediumMuon() ){
          pTmuMin = mu->pt();
	  //std::cout << "mu pT :" << mu->pt() << "Eta: " << mu->eta() << "phi:" << mu->phi() << std::endl;

	  mu_ID = mu->pdgId();
	  is_muon = true;
	  mu_selection += 1; //checking how many mu over total pass the selection
	  if (tag_lepton_found == false){ //checking how many events contain mu passing the if selection
	    mu_selection_event += 1;}
	  tag_lepton_found = true;

	lepton_pT_tree = mu->pt();
	lepton_eta_tree = mu->eta();
	lepton_phi_tree = mu->phi();
	//std::cout << "pT del leptone" << lepton_pT_tree << std::endl;
	}
  }

  for (auto el = slimmedElectrons->begin(); el != slimmedElectrons->end(); ++el){
        if (el->pt()>24 && el->pt() > pTmuMin){
          pTmuMin = el->pt();
	  //std::cout << "mu pT :" << mu->pt() << "Eta: " << mu->eta() << "phi:" << mu->phi() << std::endl;

	  el_ID = el->pdgId();
	  is_muon = false;
	  //mu_selection += 1; //checking how many mu over total pass the selection
	  //if (checker == false){ //checking how many events contain mu passing the if selection
	  //mu_selection_event += 1;}
	  tag_lepton_found = true;

	lepton_pT_tree = el->pt();
	lepton_eta_tree = el->eta();
	lepton_phi_tree = el->phi();
	//std::cout << "pT del leptone" << lepton_pT_tree << std::endl;
	}
  }
  
<<<<<<< HEAD
  cand_passing_selection = 0;
  pTpiMax = -1000;
  cand_pion_found = false;
  is_pi_a_pi = false; //together with is_photon_a_photon, used to calculate how many times both reconstructed pi and gamma come from generated particles whose mother is a W
  is_last_pi_a_pi = false;
  is_bad_single_pi = false;
  single_pi_counter = 0;
=======
  cont1 = 0;
  pTpiMax = -1000.;
  checker1 = false;
  checker3 = false; //together with checker4, used to calculate how many times both reconstructed pi and gamma come from generated particles whose mother is a W
  checker5 = false;
  checker6 = false;
  counter = 0;
>>>>>>> 8fd4d7d7eec1bebc49964549b341f61e8cac8642
  for (auto cand = PFCandidates->begin(); cand != PFCandidates->end(); ++cand){
    if ((cand->pdgId()*mu_ID > 0 || cand->pdgId()*el_ID > 0) && tag_lepton_found == true && cand->pt()>=20 && cand->trackHighPurity()==true && cand->pt()>pTpiMax && cand->fromPV()==3){

      pTpiMax = cand->pt();
      candidate_pi = cand->p4();
      cand_pion_found = true;

      pxpi = cand->px();
      pypi = cand->py();
      pzpi = cand->pz();
      pTpi = cand->pt();
      candidate_eta = cand->eta();
      candidate_phi = cand->phi();
            
      deltapTMax = 10000;
      deltaRMax = 0.3;
      gen_mother = 0;
      gen_px = 0;
      gen_py = 0;
      gen_pz = 0;
      gen_pt = 0;

    //gen_ID = cand.pdgId()#so if it doesn't enter the following loop-and-if, it doesn't display "reconstructed particle is different..." either 
    if(!runningOnData_){
    for (auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){//matching candidate for W reconstruction with MC truth
      deltaR = TMath::Sqrt((candidate_eta-gen->eta())*(candidate_eta-gen->eta())+(candidate_phi-gen->phi())*(candidate_phi-gen->phi()));
      deltapT = TMath::Abs(cand->pt()-gen->pt());

      if (deltaR <= deltaRMax && deltapT < deltapTMax){
	  deltapTMax = deltapT;
	  gen_ID = gen->pdgId();
	  gen_mother = gen->mother()->pdgId();
	  gen_px = gen->px();
	  gen_py = gen->py();
	  gen_pz = gen->pz();
	  gen_pt = gen->pt();
      }
    }}

      if (gen_ID != 211 && gen_ID != -211){
	//std::cout << "|||corresponding generated particle is not a pion, but: " << gen_ID << std::endl;
	if (single_pi_counter == 0){
	  is_bad_single_pi = true;}} //avoiding to count the case in which only a single reco-pion, not matching with a gen-pion, passes the selection

      if (gen_ID == 211 || gen_ID == -211){
	//std::cout << "\\identity of generated PION's mother: " << gen_mother << std::endl;
	//std::cout << "reco px: " << pxpi << "reco py: " << pypi << "reco pz: " << pzpi << "reco pT: " << pTpi << std::endl;
	//std::cout << "gen px: " << gen_px << "gen py: " << gen_py << "gen pz: " << gen_pz << "gen pT: " << gen_pt << std::endl;
	if (gen_mother == 24 || gen_mother == -24){
	  pi_from_w += 1;
	  is_pi_a_pi = true;
	  gen_mother = 0;
	  gen_ID = 0;
	  single_pi_counter += 1;
	  if (cand_passing_selection >= 1){
	    pi_from_w += 1;
	    is_last_pi_a_pi = true;} //activates in case the first reco particle to fulfill the conditions is not a pion from a W, while the second is
	}}
      
      if (counter >= 1 && cand_passing_selection >= 1){
	pi_from_w = pi_from_w-1; //without this, a case in which the first particle to pass the selection is a pion but the second and final is not, is counted as a good reconstruction 
	is_pi_a_pi = false;}
      if (is_last_pi_a_pi == true){
	is_pi_a_pi = true;}
                
      cand_passing_selection += 1;
 
    }}

  if (cand_passing_selection != 0){
		events_least_one_pi += 1;}

  cand_total += cand_passing_selection;


  ph_passing_selection = 0;
  eTphMin = -1000;
  cand_photon_found = false;
  is_photon_a_photon = false; //together with is_pi_a_pi, used to calculate how many times both reconstructed pi and gamma come from generated particles whose mother is a W
    for (auto photon = slimmedPhotons->begin(); photon != slimmedPhotons->end(); ++photon){
      if (photon->et()>20 && photon->et()>eTphMin && cand_pion_found == true){
      eTphMin = photon->et();
      std::cout << "px: " << photon->px() << "py: " << photon->py() << "pz: " << photon->pz() << "energy: " << photon->energy() << std::endl;
      candidate_ph = photon->p4();
      cand_photon_found = true;
      ph_passing_selection +=1;
            
      pxph = photon->px();
      pyph = photon->py();
      pzph = photon->pz();
      pTph = photon->pt();
      photon_eta = photon->eta();
      photon_phi = photon->phi();
      track_iso = photon->trackIso();
      ecal_iso = photon->ecalIso();
      hcal_iso = photon->hcalIso();
      calo_iso = photon->caloIso();
      iso_sum = (track_iso+calo_iso)/photon->et();

      track_iso_hist->Fill(track_iso);
      ecal_iso_hist->Fill(ecal_iso);
      hcal_iso_hist->Fill(hcal_iso);
      calo_iso_hist->Fill(calo_iso);
      iso_sum_hist->Fill(iso_sum);
            
            
      gen_px = 0;
      gen_py = 0;
      gen_pz = 0;
      gen_pt = 0;
      deltapTMax = 10000;
      deltaRMax = 0.3;

      //gen_ID = cand.pdgId()#so if it doesn't enter the following loop-and-if, it doesn't display "reconstructed particle is different..." either 
      if(!runningOnData_){
      for (auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){
	deltaR = TMath::Sqrt((photon_eta-gen->eta())*(photon_eta-gen->eta())+(photon_phi-gen->phi())*(photon_phi-gen->phi()));
	deltapT = TMath::Abs(photon->pt()-gen->pt());

	if (deltaR <= deltaRMax && deltapT < deltapTMax){
	  deltapTMax = deltapT;
	  gen_ID = gen->pdgId();
	  gen_mother = gen->mother()->pdgId();
	  gen_px = gen->px();
	  gen_py = gen->py();
	  gen_pz = gen->pz();
	  gen_pt = gen->pt();
	}
      }}
                
      if (gen_ID != 22){
	std::cout << "|||corresponding generated particle is not a photon, but: " << gen_ID << std::endl;
	std::cout <<  "______reco px: " << pxph << "reco py: " << pyph << "reco pz: " << pzph << "reco pT: "<< pTph << std::endl;
	std::cout << "gen px: " << gen_px << "gen py: " << gen_py << "gen pz: " << gen_pz << "gen pT: " << gen_pt << std::endl;
      }
            
      if (gen_ID == 22){
	std::cout << "\\identity of generated photon's mother: " << gen_mother << std::endl;
	std::cout << "______reco px: " << pxph << "reco py: " << pyph << "reco pz: " << pzph << "reco pT: " << pTph << std::endl;
	std::cout << "gen px: " << gen_px << "gen py: " << gen_py << "gen pz: " << gen_pz << "gen pT: " << gen_pt << std::endl;
	if (gen_mother == 24 || gen_mother == -24){
	  photon_from_w += 1;
	is_photon_a_photon = true;}
      }}}

      if (cand_pion_found == true && cand_photon_found == false && is_bad_single_pi == false){
        pi_from_w_correction += 1;}

      if (ph_passing_selection != 0){
        events_least_one_ph += 1;} //incrementing the number of events with at least one photon

      ph_total += ph_passing_selection;

    
      if (cand_pion_found == true && cand_photon_found == true){
	std::cout << "INV MASS: " << (candidate_ph + candidate_pi).M() << "costheta: " << (pxpi*pxph+pypi*pyph+pzpi*pzph)/(TMath::Sqrt(pxpi*pxpi+pypi*pypi+pzpi*pzpi)*TMath::Sqrt(pxph*pxph+pyph*pyph+pzph*pzph));
        if (is_pi_a_pi == false || is_photon_a_photon == false){
	  inv_mass_1->SetLineColor(3);
	  inv_mass_1->Fill((candidate_ph + candidate_pi).M());}
        if (is_pi_a_pi == true && is_photon_a_photon == true){
	  pi_and_photon_from_w += 1;
	  inv_mass_2->SetLineColor(2);
	  inv_mass_2->Fill((candidate_ph + candidate_pi).M());}
      }
      std::cout << "n fotoni: " << events_least_one_ph << std::endl;
      mytree->Fill();
}

void WPiGammaAnalysis::create_trees(){
  mytree = fs->make<TTree>("mytree", "Tree containing gen&reco");
  mytree->Branch("tag_lepton_pT",&lepton_pT_tree);
  mytree->Branch("tag_lepton_eta",&lepton_eta_tree);
  mytree->Branch("tag_lepton_phi",&lepton_phi_tree);
  mytree->Branch("is_muon",&is_muon);
}

//define this as a plug-in
DEFINE_FWK_MODULE(WPiGammaAnalysis);
