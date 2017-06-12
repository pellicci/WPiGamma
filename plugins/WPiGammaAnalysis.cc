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

//Photon ID stuff
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

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
  mvaValuesMapToken_el_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap_el"))),
  mvaCategoriesMapToken_el_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap_el"))),
  phoMediumIdBoolMapToken_(consumes<edm::ValueMap<bool> > (iConfig.getParameter<edm::InputTag>("phoMediumIdBoolMap"))),
  phoMediumIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> > (iConfig.getParameter<edm::InputTag>("phoMediumIdFullInfoMap"))),
  mvaValuesMapToken_ph_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap_ph"))),
  mvaCategoriesMapToken_ph_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap_ph"))),
  verboseIdFlag_(iConfig.getParameter<bool>("phoIdVerbose"))

{
  packedPFCandidatestoken_  = consumes<std::vector<pat::PackedCandidate> >(packedPFCandidates_); 
  slimmedMuonstoken_        = consumes<std::vector<pat::Muon> >(slimmedMuons_);
  prunedGenParticlestoken_  = consumes<std::vector<reco::GenParticle> >(prunedGenParticles_);
  photonsMiniAODToken_      = mayConsume<edm::View<reco::Photon> > (slimmedPhotons_);
  electronsMiniAODToken_    = mayConsume<edm::View<reco::GsfElectron> > (slimmedElectrons_);
  slimmedJetstoken_         = consumes<std::vector<pat::Jet> >(slimmedJets_);
  tok_Vertex_               = consumes<std::vector<reco::Vertex> > (pvCollection_);  
  tok_beamspot_             = consumes<reco::BeamSpot> (edm::InputTag(bsCollection_));
  pileupSummaryToken_       = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag(PileupSrc_));

  h_Events = fs->make<TH1F>("h_Events", "Event counting in different steps", 8, 0., 8.);
  _Nevents_processed  = 0;
  _Nevents_muVeto     = 0;
  _Nevents_eleVeto    = 0;
  _Nevents_isLepton   = 0;
  _Nevents_TwoLepton  = 0;
  _Nevents_isPion     = 0;
  _Nevents_isPhotons  = 0;
  _Nevents_isWmass    = 0;

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

  edm::Handle<edm::View<reco::Photon> > slimmedPhotons;
  iEvent.getByToken(photonsMiniAODToken_,slimmedPhotons);

  edm::Handle<edm::View<reco::GsfElectron> > slimmedElectrons;
  iEvent.getByToken(electronsMiniAODToken_,slimmedElectrons);

  edm::Handle<std::vector<pat::Jet > > slimmedJets;
  iEvent.getByLabel(slimmedJets_, slimmedJets);

  edm::Handle<std::vector<reco::Vertex > > slimmedPV;
  iEvent.getByLabel(pvCollection_, slimmedPV);

  //Get the electron ID data from the event stream.
  // Note: this implies that the VID ID modules have been run upstream.
  edm::Handle<edm::ValueMap<bool> > el_medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > el_tight_id_decisions; 
  iEvent.getByToken(eleMediumIdMapToken_,el_medium_id_decisions);
  iEvent.getByToken(eleTightIdMapToken_,el_tight_id_decisions);

  // Get MVA values and categories (optional)
  edm::Handle<edm::ValueMap<float> > el_mvaValues;
  edm::Handle<edm::ValueMap<int> > el_mvaCategories;
  iEvent.getByToken(mvaValuesMapToken_el_,el_mvaValues);
  iEvent.getByToken(mvaCategoriesMapToken_el_,el_mvaCategories);

  // Get the photon ID data from the event stream.
  // Note: this implies that the VID ID modules have been run upstream.
  // The first map simply has pass/fail for each particle
  edm::Handle<edm::ValueMap<bool> > ph_medium_id_decisions;
  iEvent.getByToken(phoMediumIdBoolMapToken_,ph_medium_id_decisions);
  //
  // The second map has the full info about the cut flow
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > ph_medium_id_cutflow_data;
  iEvent.getByToken(phoMediumIdFullInfoMapToken_,ph_medium_id_cutflow_data);

  // Get MVA values and categories (optional)
  edm::Handle<edm::ValueMap<float> > ph_mvaValues;
  edm::Handle<edm::ValueMap<int> > ph_mvaCategories;
  iEvent.getByToken(mvaValuesMapToken_ph_,ph_mvaValues);
  iEvent.getByToken(mvaCategoriesMapToken_ph_,ph_mvaCategories);

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

  float deltaphi_lep_pi = 0.;

  //These variables will go in the tree
  pi_pT = 0.;
  pi_eta = 0.;
  pi_phi = 0.;
  LorentzVector pi_p4;

  ph_pT = 0.;
  ph_eta = 0.;
  ph_phi = 0.;
  LorentzVector ph_p4;

  _Wmass = 0.;

  is_muon = false;
  bool is_ele  = false;

  lepton_pT_tree = 0.;
  lepton_eta_tree = 0.;
  lepton_phi_tree = 0.;

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

  if(nMuons > 1) return;
  _Nevents_muVeto++;

  //Loop over electrons
  for (size_t i = 0; i < slimmedElectrons->size(); ++i){
    const auto el = slimmedElectrons->ptrAt(i);
    if(el->pt() < 26. || el->pt() < pTeleMax) continue;

    bool isPassTight = (*el_tight_id_decisions)[el];
    if(!isPassTight) continue;

    el_ID       = el->pdgId();
    is_ele      = true;
    el_ID       = el->pdgId();
    el_eta      = el->eta();
    el_phi      = el->phi();
    el_pT       = el->pt();

    pTeleMax = el_pT;
    nElectrons++;
  }

  if(nElectrons > 1) return;
  _Nevents_eleVeto++;

  if(is_muon){
    lepton_pT_tree  = mu_pT;
    lepton_eta_tree = mu_eta;
    lepton_phi_tree = mu_phi;
  }

  if(!is_muon && is_ele){
    lepton_pT_tree  = el_pT;
    lepton_eta_tree = el_eta;
    lepton_phi_tree = el_phi;
  }

  //Do NOT continue if you didn't find either a muon or an electron
  if(!is_muon && !is_ele) return;
  _Nevents_isLepton++;

  if(is_muon && is_ele) return;
  _Nevents_TwoLepton++;

  //In signal, identify if there's a real mu or ele from W
  bool is_Wplus_from_t = false;
  bool is_Wminus_from_tbar = false;
  bool is_Wplus_in_lep = false;
  bool is_Wminus_in_lep = false;
  is_ttbar_lnu = false;
  is_Ele_signal = false;
  is_Mu_signal = false;
  if(!runningOnData_){
    for (auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){
      if(fabs(gen->pdgId() ) == 11 && fabs(gen->mother()->pdgId()) == 24) is_Ele_signal = true;
      if(fabs(gen->pdgId() ) == 13 && fabs(gen->mother()->pdgId()) == 24) is_Mu_signal = true;
      if(gen->pdgId() == 24 && gen->mother()->pdgId() == 6) is_Wplus_from_t = true;
      if(gen->pdgId() == -24 && gen->mother()->pdgId() == -6) is_Wminus_from_tbar = true;
      if(fabs(gen->pdgId() != 24) || gen->numberOfDaughters() != 2) continue;
      for (int i=0; i<2; i++){
	if(gen->daughter(i)->pdgId() == 11 || gen->daughter(i)->pdgId() == 13 || gen->daughter(i)->pdgId() == 15) is_Wminus_in_lep = true;
	if(gen->daughter(i)->pdgId() == -11 || gen->daughter(i)->pdgId() == -13 || gen->daughter(i)->pdgId() == -15) is_Wplus_in_lep = true;
      }
    }
    if(is_Wplus_from_t == true && is_Wminus_from_tbar == true && is_Wminus_in_lep == true && is_Wplus_in_lep == true) is_ttbar_lnu = true;
  }
  
  //----------- Starting to search for pi and gamma -------------
  bool cand_pion_found = false;

  for(auto cand = PFCandidates->begin(); cand != PFCandidates->end(); ++cand){
    if(cand->pdgId()*mu_ID < 0 && cand->pdgId()*el_ID < 0) continue;    
    if(cand->pt() < 20. || !cand->trackHighPurity() || cand->fromPV() != 3) continue;
    if(cand->pt() < pTpiMax) continue;

    if(is_muon){
      deltaphi_lep_pi = fabs(mu_phi-cand->phi());
      if(deltaphi_lep_pi > 3.14) deltaphi_lep_pi = 6.28-deltaphi_lep_pi;
    } 
    if(!is_muon && is_ele){
      deltaphi_lep_pi = fabs(el_phi-cand->phi());
      if(deltaphi_lep_pi > 3.14) deltaphi_lep_pi = 6.28-deltaphi_lep_pi;
    } 

    if(deltaphi_lep_pi < 0.00005) continue;

    pTpiMax = cand->pt();
    cand_pion_found = true;
    nPions++;

    pi_pT  = cand->pt();
    pi_eta = cand->eta();
    pi_phi = cand->phi();
    pi_p4  = cand->p4();

    is_pi_a_pi = false;
    is_pi_matched = false;

    float deltapTMax = 10000.;
    const float deltaRMax  = 0.3;
    int   gen_mother = 0;
    int   gen_ID = 0;

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

    }
  }

  //Do NOT continue if you didn't find a pion
  if(!cand_pion_found) return;
  _Nevents_isPion++;

  bool cand_photon_found = false;

  for (size_t i = 0; i < slimmedPhotons->size(); ++i){
    const auto photon = slimmedPhotons->ptrAt(i);

    if(photon->et() < 20. || photon->et() < eTphMax) continue;
    if(photon->hasPixelSeed()) continue;   //electron veto

    // The minimal info
    bool isPassMedium = (*ph_medium_id_decisions)[photon];
    if(!isPassMedium) continue;

    eTphMax = photon->et();

    ph_pT  = photon->pt();
    ph_eta = photon->eta();
    ph_phi = photon->phi();
    ph_p4  = photon->p4();

    cand_photon_found = true;
    nPhotons++;

    float deltapTMax = 10000.;
    const float deltaRMax = 0.3;
    int   gen_mother = 0;
    int   gen_ID = 0;

    is_photon_a_photon = false;
    is_photon_matched = false;

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
      //if(gen_ID != 22) std::cout << "ph gen ID = " << gen_ID << std::endl;
      //if(gen_ID != 22 && fabs(gen_mother) == 24) std::cout << "ph gen ID when matched = " << gen_ID << std::endl;
      if(fabs(gen_mother) == 24) is_photon_matched = true;
    }

  }

  //Do not continue if there's no photons
  if(!cand_photon_found) return;
  _Nevents_isPhotons++;

  _Wmass = (pi_p4 + ph_p4).M();

  //Only save events in a certain range
  if(_Wmass < 20. || _Wmass > 120.) return;
  _Nevents_isWmass++;


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
  mytree->Branch("is_muon",&is_muon);
  mytree->Branch("pi_pT",&pi_pT);
  mytree->Branch("pi_eta",&pi_eta);
  mytree->Branch("pi_phi",&pi_phi);
  mytree->Branch("photon_eT",&ph_pT);
  mytree->Branch("photon_eta",&ph_eta);
  mytree->Branch("photon_phi",&ph_phi);

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
    mytree->Branch("isttbarlnu",&is_ttbar_lnu);
  }

}

void WPiGammaAnalysis::endJob() 
{
  h_Events->Fill(0.5,_Nevents_processed);
  h_Events->Fill(1.5,_Nevents_muVeto);
  h_Events->Fill(2.5,_Nevents_eleVeto);
  h_Events->Fill(3.5,_Nevents_isLepton);
  h_Events->Fill(4.5,_Nevents_TwoLepton);
  h_Events->Fill(5.5,_Nevents_isPion);
  h_Events->Fill(6.5,_Nevents_isPhotons);
  h_Events->Fill(7.5,_Nevents_isWmass);

  h_Events->GetXaxis()->SetBinLabel(1,"Events triggered");
  h_Events->GetXaxis()->SetBinLabel(2,"After mu veto");
  h_Events->GetXaxis()->SetBinLabel(3,"After ele veto");
  h_Events->GetXaxis()->SetBinLabel(4,"Events with lept");
  h_Events->GetXaxis()->SetBinLabel(5,"Two lep veto");
  h_Events->GetXaxis()->SetBinLabel(6,"Events with pi");
  h_Events->GetXaxis()->SetBinLabel(7,"Events with gam");
  h_Events->GetXaxis()->SetBinLabel(8,"Events in W mass");

}


//define this as a plug-in
DEFINE_FWK_MODULE(WPiGammaAnalysis);
