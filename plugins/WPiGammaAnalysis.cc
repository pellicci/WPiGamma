//ROOT includes
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "Math/VectorUtil.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
  
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

// electron ID stuff
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

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
  PileupSrc_(iConfig.getParameter<edm::InputTag>("PileupSrc")),
  eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  mvaValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap"))),
  mvaCategoriesMapToken_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap")))
{
  packedPFCandidatestoken_  = consumes<std::vector<pat::PackedCandidate> >(packedPFCandidates_); 
  slimmedMuonstoken_        = consumes<std::vector<pat::Muon> >(slimmedMuons_);
  prunedGenParticlestoken_  = consumes<std::vector<reco::GenParticle> >(prunedGenParticles_);
  slimmedPhotonstoken_      = consumes<std::vector<pat::Photon> >(slimmedPhotons_);
  electronsMiniAODToken_    = mayConsume<edm::View<reco::GsfElectron> > (slimmedElectrons_);
  slimmedJetstoken_         = consumes<std::vector<pat::Jet> >(slimmedJets_);
  tok_Vertex_               = consumes<std::vector<reco::Vertex> > (pvCollection_);  
  tok_beamspot_             = consumes<reco::BeamSpot> (edm::InputTag(bsCollection_));
  pileupSummaryToken_       = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag(PileupSrc_));

  h_Events = fs->make<TH1F>("h_Events", "Event counting in different steps", 7, 0., 7.);
  _Nevents_processed  = 0;
  _Nevents_isMuon     = 0;
  _Nevents_isElectron = 0;
  _Nevents_isLepton   = 0;
  _Nevents_isPion     = 0;
  _Nevents_isPhotons  = 0;

  inv_mass_1 = fs->make<TH1F>("Mw - no match with MC Truth", "Mw no match", 200,0,120);
  inv_mass_2 = fs->make<TH1F>("Mw - match with MC Truth", "Mw match", 200,0,120);

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

  edm::Handle<edm::View<reco::GsfElectron> > slimmedElectrons;
  iEvent.getByToken(electronsMiniAODToken_,slimmedElectrons);

  edm::Handle<std::vector<pat::Jet > > slimmedJets;
  iEvent.getByLabel(slimmedJets_, slimmedJets);

  edm::Handle<std::vector<reco::Vertex > > slimmedPV;
  iEvent.getByLabel(pvCollection_, slimmedPV);

  //Get the electron ID data from the event stream.
  // Note: this implies that the VID ID modules have been run upstream.
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions; 
  iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(eleTightIdMapToken_,tight_id_decisions);

  // Get MVA values and categories (optional)
  edm::Handle<edm::ValueMap<float> > mvaValues;
  edm::Handle<edm::ValueMap<int> > mvaCategories;
  iEvent.getByToken(mvaValuesMapToken_,mvaValues);
  iEvent.getByToken(mvaCategoriesMapToken_,mvaCategories);

  _Nevents_processed++;

  nMuons     = 0;
  nElectrons = 0;
  nPions     = 0;
  nPhotons   = 0;
  nBjets     = 0;

  is_pi_a_pi         = false;
  is_pi_matched      = false;
  is_photon_a_photon = false;
  is_photon_matched  = false;

  float pTmuMax  = -1000.;
  float pTeleMax = -1000.;
  float pTpiMax  = -1000.;
  float eTphMax  = -1000.;

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

  _Wmass = 0.;

  bool is_muon = false;
  bool is_ele  = false;

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

  is_muon_tree = false;

  //Loop over muons
  for(auto mu = slimmedMuons->begin(); mu != slimmedMuons->end(); ++mu){
    if(mu->pt() < 25. || mu->pt() < pTmuMax || !mu->isTightMuon(slimmedPV->at(0))) continue;
    if( (mu->chargedHadronIso() + std::max(0., mu->neutralHadronIso() + mu->photonIso() - 0.5*mu->puChargedHadronIso())/mu->pt()) > 0.2) continue;
    pTmuMax = mu->pt();
    //std::cout << "mu pT :" << mu->pt() << "Eta: " << mu->eta() << "phi:" << mu->phi() << std::endl;

    mu_ID    = mu->pdgId();
    is_muon  = true;
    mu_eta   = mu->eta();
    mu_phi   = mu->phi();
    mu_pT    = mu->pt();
    nMuons++;
  }

  //Loop over electrons
  for (size_t i = 0; i < slimmedElectrons->size(); ++i){
    const auto el = slimmedElectrons->ptrAt(i);
    if(el->pt() < 26. || el->pt() < pTeleMax) continue;

    bool isPassMedium = (*medium_id_decisions)[el];
    bool isPassTight = (*tight_id_decisions)[el];

    std::cout << isPassMedium << " " << isPassTight << std::endl;

    el_ID       = el->pdgId();
    is_ele      = true;
    el_ID       = el->pdgId();
    el_eta      = el->eta();
    el_phi      = el->phi();
    el_pT       = el->pt();
    nElectrons++;
  }

  if(is_muon){
    lepton_pT_tree  = mu_pT;
    lepton_eta_tree = mu_eta;
    lepton_phi_tree = mu_phi;
    _Nevents_isMuon++;
  }

  if(!is_muon && is_ele){
    lepton_pT_tree  = el_pT;
    lepton_eta_tree = el_eta;
    lepton_phi_tree = el_phi;
    _Nevents_isElectron++;
  }

  //Do NOT continue if you didn't find either a muon or an electron
  if(!is_muon && !is_ele) return;

  if(is_muon) is_muon_tree = true;
  else is_muon_tree = false;

  _Nevents_isLepton++;

  //In signal, identify if there's a real mu or ele from W
  is_Ele_signal = false;
  is_Mu_signal = false;
  if(!runningOnData_){
    for (auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){
      if(fabs(gen->pdgId() ) == 11 && fabs(gen->mother()->pdgId()) == 24) is_Ele_signal = true;
      if(fabs(gen->pdgId() ) == 13 && fabs(gen->mother()->pdgId()) == 24) is_Mu_signal = true;
    }
  }

  //----------- Starting to search for pi and gamma -------------
  bool cand_pion_found = false;

  for(auto cand = PFCandidates->begin(); cand != PFCandidates->end(); ++cand){
    if(cand->pdgId()*mu_ID < 0 && cand->pdgId()*el_ID < 0) continue;    
    if(cand->pt() < 20. || !cand->trackHighPurity() || cand->fromPV() != 3) continue;
    if(cand->pt() < pTpiMax) continue;

    pTpiMax = cand->pt();
    cand_pion_found = true;
    nPions++;

    pi_pT  = cand->pt();
    pi_eta = cand->eta();
    pi_phi = cand->phi();
    pi_p4  = cand->p4();

    float deltapTMax = 10000.;
    const float deltaRMax  = 0.3;
    int   gen_mother = 0;
    int   gen_ID = 0;

    is_pi_a_pi = false;
    is_pi_matched = false;

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

      if(fabs(gen_ID) == 211) is_pi_a_pi = true;
      if(fabs(gen_mother) == 24) is_pi_matched = true;

      //std::cout << "\\identity of generated PION's mother: " << gen_mother << std::endl;
      //std::cout << "reco px: " << pxpi << "reco py: " << pypi << "reco pz: " << pzpi << "reco pT: " << pTpi << std::endl;
      //std::cout << "gen px: " << gen_px << "gen py: " << gen_py << "gen pz: " << gen_pz << "gen pT: " << gen_pt << std::endl;
    }
  }

  //Do NOT continue if you didn't find a pion
  if(!cand_pion_found) return;

  _Nevents_isPion++;

  pi_pT_tree  = pi_pT;
  pi_eta_tree = pi_eta;
  pi_phi_tree = pi_phi;

  bool cand_photon_found = false;

  for (auto photon = slimmedPhotons->begin(); photon != slimmedPhotons->end(); ++photon){

    ph_pT  = photon->pt();
    ph_eta = photon->eta();
    ph_phi = photon->phi();
    ph_p4  = photon->p4();

    if(is_ele && fabs(ph_pT - lepton_pT_tree) < 1.) continue;
    if(photon->et() < 20. || photon->et() < eTphMax) continue;
    eTphMax = photon->et();

    cand_photon_found = true;
    nPhotons++;

    float deltapTMax = 10000.;
    const float deltaRMax = 0.3;
    int   gen_mother = 0;
    int   gen_ID = 0;

    is_photon_a_photon = false;
    is_photon_matched = false;

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
            
      if(gen_ID == 22) is_photon_a_photon = true;
      if(gen_ID != 22) std::cout << "ph gen ID = " << gen_ID << std::endl;
      if(gen_ID != 22 && fabs(gen_mother) == 24) std::cout << "ph gen ID when matched = " << gen_ID << std::endl;
      if(fabs(gen_mother) == 24) is_photon_matched = true;
    }

  }

  photon_eT_tree =  ph_pT;
  photon_eta_tree = ph_eta;
  photon_phi_tree = ph_phi;

  //Do not continue if there's no photons
  if(!cand_photon_found) return;

  _Nevents_isPhotons++;

  _Wmass = (pi_p4 + ph_p4).M();

  if (!is_pi_a_pi || !is_photon_a_photon){
    inv_mass_1->SetLineColor(3);
    inv_mass_1->Fill(_Wmass);
  }
  if(is_pi_a_pi && is_photon_a_photon && is_pi_matched && is_photon_matched){
    inv_mass_2->SetLineColor(2);
    inv_mass_2->Fill(_Wmass);
  }

  nBjets = 0;
  for (auto jet = slimmedJets->begin(); jet != slimmedJets->end(); ++jet){
    if(jet->pt() < 30.) continue;
    if(jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") < 0.46) continue;   //0.46 = loose
    nBjets++;
  }

  mytree->Fill();

  //std::cout << "n fotoni: " << events_least_one_ph << std::endl;
}

void WPiGammaAnalysis::create_trees(){
  mytree = fs->make<TTree>("mytree", "Tree containing gen&reco");
  mytree->Branch("lepton_pT",&lepton_pT_tree);
  mytree->Branch("lepton_eta",&lepton_eta_tree);
  mytree->Branch("lepton_phi",&lepton_phi_tree);
  mytree->Branch("is_muon",&is_muon_tree);
  mytree->Branch("pi_pT",&pi_pT_tree);
  mytree->Branch("pi_eta",&pi_eta_tree);
  mytree->Branch("pi_phi",&pi_phi_tree);
  mytree->Branch("photon_eT",&photon_eT_tree);
  mytree->Branch("photon_eta",&photon_eta_tree);
  mytree->Branch("photon_phi",&photon_phi_tree);

  mytree->Branch("Wmass",&_Wmass);

  mytree->Branch("nMuons",&nMuons);
  mytree->Branch("nElectrons",&nElectrons);
  mytree->Branch("nPions",&nPions);
  mytree->Branch("nPhotons",&nPhotons);
  mytree->Branch("nBjets",&nBjets);

  //Save MC truth
  if(!runningOnData_){
    mytree->Branch("isMuonSignal",&is_Mu_signal);
    mytree->Branch("isEleSignal",&is_Ele_signal);
    mytree->Branch("isPionTrue",&is_pi_a_pi);
    mytree->Branch("isPionMatched",&is_pi_matched);
    mytree->Branch("isPhotonTrue",&is_photon_a_photon);
    mytree->Branch("isPhotonMatched",&is_photon_matched);
  }

}

void WPiGammaAnalysis::endJob() 
{
  h_Events->Fill(0.5,_Nevents_processed);
  h_Events->Fill(1.5,_Nevents_isMuon);
  h_Events->Fill(2.5,_Nevents_isElectron);
  h_Events->Fill(3.5,_Nevents_isLepton);
  h_Events->Fill(4.5,_Nevents_isPion);
  h_Events->Fill(5.5,_Nevents_isPhotons);
}


//define this as a plug-in
DEFINE_FWK_MODULE(WPiGammaAnalysis);
