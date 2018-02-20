//ROOT includes
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "Math/VectorUtil.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
  
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
 
//Vertex inclusions
#include "DataFormats/VertexReco/interface/Vertex.h" 
#include "DataFormats/BeamSpot/interface/BeamSpot.h" 
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//Electron ID stuff
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

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
  slimmedMETs_(iConfig.getParameter<edm::InputTag>("slimmedMETs")),
  runningOnData_(iConfig.getParameter<bool>("runningOnData")),
  pvCollection_(iConfig.getParameter<edm::InputTag>("pvCollection")),   
  bsCollection_(iConfig.getParameter<edm::InputTag>("bsCollection")),  
  PileupSrc_(iConfig.getParameter<edm::InputTag>("PileupSrc")),
  triggerBits_(consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("triggerbits"))),
  eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  mvaValuesMapToken_el_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap_el"))),
  mvaCategoriesMapToken_el_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap_el"))),
  phoMediumIdBoolMapToken_(consumes<edm::ValueMap<bool> > (iConfig.getParameter<edm::InputTag>("phoMediumIdBoolMap"))),
  phoMediumIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> > (iConfig.getParameter<edm::InputTag>("phoMediumIdFullInfoMap"))),
  mvaValuesMapToken_ph_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap_ph"))),
  mvaCategoriesMapToken_ph_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap_ph"))),
  verboseIdFlag_(iConfig.getParameter<bool>("phoIdVerbose")),
  effectiveAreas_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath() )

{
  packedPFCandidatestoken_  = consumes<std::vector<pat::PackedCandidate> >(packedPFCandidates_); 
  slimmedMuonstoken_        = consumes<std::vector<pat::Muon> >(slimmedMuons_);
  prunedGenParticlestoken_  = consumes<std::vector<reco::GenParticle> >(prunedGenParticles_);
  photonsMiniAODToken_      = mayConsume<edm::View<reco::Photon> > (slimmedPhotons_);
  electronsMiniAODToken_    = mayConsume<edm::View<reco::GsfElectron> > (slimmedElectrons_);
  slimmedJetstoken_         = consumes<std::vector<pat::Jet> >(slimmedJets_);
  slimmedMETstoken_         = consumes<std::vector<pat::MET> >(slimmedMETs_);
  tok_Vertex_               = consumes<std::vector<reco::Vertex> > (pvCollection_);  
  tok_beamspot_             = consumes<reco::BeamSpot> (edm::InputTag(bsCollection_));
  pileupSummaryToken_       = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag(PileupSrc_));
  rhoToken_                 = consumes<double> (iConfig.getParameter <edm::InputTag>("rho"));

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

  edm::Handle<std::vector<pat::MET > > slimmedMETs;
  iEvent.getByLabel(slimmedMETs_, slimmedMETs);

  edm::Handle<std::vector<reco::Vertex > > slimmedPV;
  iEvent.getByLabel(pvCollection_, slimmedPV);

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);

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


  //Retrieve the run number
  if(runningOnData_){
    run_number = iEvent.id().run();
  }

  //Count the number of vertices
  nPV = -1;
  if(slimmedPV->size()<=0) return;
  for(reco::VertexCollection::const_iterator vtx=slimmedPV->begin();vtx!=slimmedPV->end();++vtx) {
    // check that the primary vertex is not a fake one, that is the beamspot (it happens when no primary vertex is reconstructed)
    if(!vtx->isFake()) {
      nPV++;
    }
  } 

  //PileUp code for examining the Pileup information
  PU_Weight = 1.;
  float npT = -1.;

  if(!runningOnData_){
    edm::Handle<std::vector< PileupSummaryInfo>>  PupInfo;
    iEvent.getByLabel(PileupSrc_, PupInfo);
  
    std::vector<PileupSummaryInfo>::const_iterator PVI; 
 
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      const int BX = PVI->getBunchCrossing();
      if(BX == 0) {
	npT  = PVI->getTrueNumInteractions();
      }
    }

    // calculate weight using above code
    PU_Weight = Lumiweights_.weight(npT);
  }

  //Examine the trigger information
  isSingleMuTrigger_24 = false;
  isSingleMuTrigger_50 = false;
  isSingleEleTrigger = false;

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i){
    if(!triggerBits->accept(i)) continue;
    std::string tmp_triggername = names.triggerName(i);
    if( tmp_triggername.find("HLT_IsoMu24_v") != std::string::npos ||
	tmp_triggername.find("HLT_IsoTkMu24_v") != std::string::npos){
      isSingleMuTrigger_24 = true;
    }
    if( tmp_triggername.find("HLT_Mu50_v") != std::string::npos ||
	tmp_triggername.find("HLT_TkMu50_v") != std::string::npos){
      isSingleMuTrigger_50 = true;
    }
    if( tmp_triggername.find("HLT_Ele25_eta2p1_WPTight_Gsf_v") != std::string::npos ||
	tmp_triggername.find("HLT_Ele27_WPTight_Gsf_v") != std::string::npos){
      isSingleEleTrigger = true;
    }
  }

  nMuons     = 0;
  nElectrons = 0;
  nPions     = 0;
  nPhotons   = 0;
  nBjets     = 0;
  nBjets_25  = 0;

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
  lepton_iso = 0.;
  ph_iso_ChargedHadron = 0.;
  ph_iso_NeutralHadron = 0.;
  ph_iso_Photon = 0.;
  //ph_iso_Track = 0.;
  ph_iso_eArho = 0.;

  pi_pT = 0.;
  pi_eta = 0.;
  pi_phi = 0.;
  pi_energy = 0.;
  LorentzVector pi_p4;

  ph_eT = 0.;
  ph_eta = 0.;
  ph_phi = 0.;
  ph_energy = 0.;
  LorentzVector ph_p4;

  _Wmass = 0.;

  met_pT = 0.;

  is_muon = false;
  bool is_ele  = false;

  lepton_pT_tree = 0.;
  lepton_eta_tree = 0.;
  lepton_phi_tree = 0.;


  for(auto met = slimmedMETs->begin(); met != slimmedMETs->end(); ++met){
    met_pT = met->pt();
  }

  //Loop over muons
  for(auto mu = slimmedMuons->begin(); mu != slimmedMuons->end(); ++mu){
    if(mu->pt() < 24. || mu->pt() < pTmuMax || !mu->isMediumMuon() || abs(mu->eta()) > 2.4) continue;
    lepton_iso = (mu->chargedHadronIso() + std::max(0., mu->neutralHadronIso() + mu->photonIso() - 0.5*mu->puChargedHadronIso()))/mu->pt();
    if(lepton_iso > 0.25) continue;
    pTmuMax = mu->pt();
    //std::cout << "mu pT :" << mu->pt() << "Eta: " << mu->eta() << "phi:" << mu->phi() << std::endl;

    mu_ID   = mu->pdgId();
    is_muon = true;
    mu_eta  = mu->eta();
    mu_phi  = mu->phi();
    mu_pT   = mu->pt();
    nMuons++;
  }

  if(nMuons > 1) return;
  _Nevents_muVeto++;

  // Get rho value
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;

  //Loop over electrons
  for (size_t i = 0; i < slimmedElectrons->size(); ++i){
    const auto el = slimmedElectrons->ptrAt(i);
    if(el->pt() < 26. || el->pt() < pTeleMax) continue;

    //bool isPassTight = (*el_tight_id_decisions)[el];
    //if(!isPassTight) continue;

    bool isPassMedium = (*el_medium_id_decisions)[el];
    if(!isPassMedium) continue;

    //PflowIsolationVariables pfIso = el->pfIsolationVariables();
    float abseta =  abs(el->superCluster()->eta());
    float eA = effectiveAreas_.getEffectiveArea(abseta);
    lepton_iso = (el->pfIsolationVariables().sumChargedHadronPt + std::max( 0.0f, el->pfIsolationVariables().sumNeutralHadronEt + el->pfIsolationVariables().sumPhotonEt - eA*rho_))/el->pt();
    if(lepton_iso > 0.4) continue;

    el_ID    = el->pdgId();
    is_ele   = true;
    el_ID    = el->pdgId();
    el_eta   = el->eta();
    el_phi   = el->phi();
    el_pT    = el->pt();

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
    if(is_Wplus_from_t && is_Wminus_from_tbar && is_Wminus_in_lep && is_Wplus_in_lep) is_ttbar_lnu = true;
  }
  
  //----------- Starting to search for pi and gamma -------------
  bool cand_pion_found = false;
  sum_pT_03 = 0.;
  sum_pT_05 = 0.;

  for(auto cand = PFCandidates->begin(); cand != PFCandidates->end(); ++cand){
    if(cand->pdgId()*mu_ID < 0 && cand->pdgId()*el_ID < 0) continue;    
    if(cand->pt() < 20. || !cand->trackHighPurity() || cand->fromPV() != 3) continue;
    if(cand->pt() < pTpiMax) continue;
    //if(cand->trackIso() > 5) continue;
    //std::cout << cand->pfIsolationVariables().sumChargedHadronPt << std::endl;
    
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
    pi_energy = cand->energy();
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

  //Computing pion isolation
  for(auto cand_iso = PFCandidates->begin(); cand_iso != PFCandidates->end(); ++cand_iso){
    float deltaR = sqrt((pi_eta-cand_iso->eta())*(pi_eta-cand_iso->eta())+(pi_phi-cand_iso->phi())*(pi_phi-cand_iso->phi()));
    if(deltaR <= 0.3 && deltaR >= 0.02) sum_pT_03 += cand_iso->pt();
    if(deltaR <= 0.5 && deltaR >= 0.02) sum_pT_05 += cand_iso->pt();
  }


  bool cand_photon_found = false;

  for (size_t i = 0; i < slimmedPhotons->size(); ++i){
    const auto photon = slimmedPhotons->ptrAt(i);

    if(photon->et() < 20. || abs(photon->eta()) > 2.5 || photon->et() < eTphMax) continue;
    if(photon->hasPixelSeed()) continue;   //electron veto

    // The minimal info
    bool isPassMedium = (*ph_medium_id_decisions)[photon];
    if(!isPassMedium) continue;

    float abseta = abs(photon->superCluster()->eta());
    float eA = effectiveAreas_.getEffectiveArea(abseta);
    //photon_iso = (pfIso.sumChargedHadronPt + std::max( 0.0f, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - eA*rho_))/photon->et();
    if(photon->chargedHadronIso()/photon->et() > 0.3 || photon->photonIso() > 4.) continue; //|| photon->trackIso() > 6
    ph_iso_ChargedHadron = photon->chargedHadronIso();
    ph_iso_NeutralHadron = photon->neutralHadronIso();
    ph_iso_Photon        = photon->photonIso();
    //ph_iso_Track         = photon->trackIso();
    ph_iso_eArho         = eA*rho_;

    eTphMax = photon->et();

    ph_eT     = photon->et();
    ph_eta    = photon->eta();
    ph_phi    = photon->phi();
    ph_energy = photon->energy();
    ph_p4     = photon->p4();

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
	float deltapT = fabs(ph_eT-gen->pt());

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
    if(jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") < 0.46) continue;   //0.46 = loose
    if(jet->pt() < 25.) continue;
    nBjets_25++;
    if(jet->pt() < 30.) continue;
    nBjets++;
  }

  mytree->Fill();

}

void WPiGammaAnalysis::create_trees()
{
  mytree = fs->make<TTree>("mytree", "Tree containing gen&reco");

  mytree->Branch("nPV",&nPV);
  mytree->Branch("isSingleMuTrigger_24",&isSingleMuTrigger_24);
  mytree->Branch("isSingleMuTrigger_50",&isSingleMuTrigger_50);
  mytree->Branch("isSingleEleTrigger",&isSingleEleTrigger);

  //Save run number info when running on data
  if(runningOnData_){
    mytree->Branch("run_number",&run_number);
  }

  mytree->Branch("lepton_pT",&lepton_pT_tree);
  mytree->Branch("lepton_eta",&lepton_eta_tree);
  mytree->Branch("lepton_phi",&lepton_phi_tree);
  mytree->Branch("lepton_iso",&lepton_iso);
  mytree->Branch("is_muon",&is_muon);
  mytree->Branch("pi_pT",&pi_pT);
  mytree->Branch("pi_eta",&pi_eta);
  mytree->Branch("pi_phi",&pi_phi);
  mytree->Branch("pi_energy",&pi_energy);
  mytree->Branch("sum_pT_03",&sum_pT_03);
  mytree->Branch("sum_pT_05",&sum_pT_05);
  mytree->Branch("photon_eT",&ph_eT);
  mytree->Branch("photon_eta",&ph_eta);
  mytree->Branch("photon_phi",&ph_phi);
  mytree->Branch("photon_energy",&ph_energy);
  mytree->Branch("photon_iso_ChargedHadron",&ph_iso_ChargedHadron);
  mytree->Branch("photon_iso_NeutralHadron",&ph_iso_NeutralHadron);
  mytree->Branch("photon_iso_Photon",&ph_iso_Photon);
  //mytree->Branch("photon_iso_Track",&ph_iso_Track);
  mytree->Branch("photon_iso_eArho",&ph_iso_eArho);

  mytree->Branch("Wmass",&_Wmass);

  mytree->Branch("nMuons",&nMuons);
  mytree->Branch("nElectrons",&nElectrons);
  mytree->Branch("nPions",&nPions);
  mytree->Branch("nPhotons",&nPhotons);
  mytree->Branch("nBjets",&nBjets);
  mytree->Branch("nBjets_25",&nBjets_25);
  mytree->Branch("met_pT",&met_pT);

  //Save MC info
  if(!runningOnData_){
    mytree->Branch("PU_Weight",&PU_Weight);

    mytree->Branch("isMuonSignal",&is_Mu_signal);
    mytree->Branch("isEleSignal",&is_Ele_signal);
    mytree->Branch("isPionTrue",&is_pi_a_pi);
    mytree->Branch("isPionMatched",&is_pi_matched);
    mytree->Branch("isPhotonTrue",&is_photon_a_photon);
    mytree->Branch("isPhotonMatched",&is_photon_matched);
    mytree->Branch("isttbarlnu",&is_ttbar_lnu);
  }

}

void WPiGammaAnalysis::beginJob()
{
  //Flag for PileUp reweighting
  if (!runningOnData_){
   Lumiweights_ = edm::LumiReWeighting("pileUpHistogramFromjson_Nominal.root","MCpileUp_25ns_Recent2016.root", "pileup", "pileup");
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
