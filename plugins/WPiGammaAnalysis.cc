//ROOT includes
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "Math/VectorUtil.h"
#include <stdlib.h>

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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

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
//#include "RecoEgamma/PhotonIdentification/interface/PhotonMVAEstimator.h"

typedef math::XYZTLorentzVector LorentzVector;
 
#include "WPiGammaAnalysis.h"
 
// constructors and destructor
WPiGammaAnalysis::WPiGammaAnalysis(const edm::ParameterSet& iConfig) :
  runningOnData_(iConfig.getParameter<bool>("runningOnData")),
  runningOn2017_(iConfig.getParameter<bool>("runningOn2017")),
  verboseIdFlag_(iConfig.getParameter<bool>("phoIdVerbose")),
  effectiveAreas_el_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile_el")).fullPath() ),
  effectiveAreas_ph_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile_ph")).fullPath() ),
  Bjets_WP_(iConfig.getParameter<double>("Bjets_WP"))

{
  packedPFCandidatesToken_            = consumes<std::vector<pat::PackedCandidate> >(edm::InputTag("packedPFCandidates")); 
  slimmedMuonsToken_                  = consumes<std::vector<pat::Muon> >(edm::InputTag("slimmedMuons"));
  prunedGenParticlesToken_            = consumes<std::vector<reco::GenParticle> >(edm::InputTag("prunedGenParticles"));
  photonsMiniAODToken_                = consumes<std::vector<pat::Photon> > (edm::InputTag("slimmedPhotons"));
  electronsMiniAODToken_              = consumes<std::vector<pat::Electron> > (edm::InputTag("slimmedElectrons"));
  slimmedJetsToken_                   = consumes<std::vector<pat::Jet> >(edm::InputTag("slimmedJets"));
  slimmedMETsToken_                   = consumes<std::vector<pat::MET> >(edm::InputTag("slimmedMETs"));
  slimmedMETsPuppiToken_              = consumes<std::vector<pat::MET> >(edm::InputTag("slimmedMETsPuppi"));
  offlineSlimmedPrimaryVerticesToken_ = consumes<std::vector<reco::Vertex> > (edm::InputTag("offlineSlimmedPrimaryVertices"));  
  offlineBeamSpotToken_               = consumes<reco::BeamSpot> (edm::InputTag("offlineBeamSpot"));
  pileupSummaryToken_                 = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo"));
  GenInfoToken_                       = consumes<GenEventInfoProduct> (edm::InputTag("generator"));
  triggerBitsToken_                   = consumes<edm::TriggerResults> (edm::InputTag("TriggerResults","","HLT"));
  rhoToken_                           = consumes<double> (iConfig.getParameter <edm::InputTag>("rho"));

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
  h_pileup   = fs->make<TH1F>("pileup", "pileup", 75,0,75);

  create_trees();
}

WPiGammaAnalysis::~WPiGammaAnalysis()
{
}

// ------------ method called for each event  ------------
void WPiGammaAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<pat::PackedCandidate>  > PFCandidates;
  iEvent.getByToken(packedPFCandidatesToken_, PFCandidates);

  edm::Handle<std::vector<pat::Muon>  > slimmedMuons;
  iEvent.getByToken(slimmedMuonsToken_, slimmedMuons);

  edm::Handle<std::vector<reco::GenParticle>  > genParticles;
  if(!runningOnData_)iEvent.getByToken(prunedGenParticlesToken_, genParticles);

  edm::Handle<std::vector<pat::Photon> > slimmedPhotons;
  iEvent.getByToken(photonsMiniAODToken_,slimmedPhotons);

  edm::Handle<std::vector<pat::Electron> > slimmedElectrons;
  iEvent.getByToken(electronsMiniAODToken_,slimmedElectrons);

  edm::Handle<std::vector<pat::Jet > > slimmedJets;
  iEvent.getByToken(slimmedJetsToken_, slimmedJets);

  edm::Handle<std::vector<pat::MET > > slimmedMETs;
  iEvent.getByToken(slimmedMETsToken_, slimmedMETs);

  edm::Handle<std::vector<pat::MET > > slimmedMETsPuppi;
  iEvent.getByToken(slimmedMETsPuppiToken_, slimmedMETsPuppi);

  edm::Handle<std::vector<reco::Vertex > > slimmedPV;
  iEvent.getByToken(offlineSlimmedPrimaryVerticesToken_, slimmedPV);

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsToken_, triggerBits);


  _Nevents_processed++;


  //Retrieve the run number
  if(runningOnData_){
    run_number = iEvent.id().run();
  }

  //*************************************************************//
  //                                                             //
  //-------------------------- Vertices -------------------------//
  //                                                             //
  //*************************************************************//

  //Count the number of vertices
  nPV = -1;

  if(slimmedPV->size()<=0) return;
  for(reco::VertexCollection::const_iterator vtx=slimmedPV->begin();vtx!=slimmedPV->end();++vtx) {
    // check that the primary vertex is not a fake one, that is the beamspot (it happens when no primary vertex is reconstructed)
    if(!vtx->isFake()) {
      nPV++;
    }
  } 
  // std::cout << "slimmedPV size: " << slimmedPV->size() << "   PV: " << &(slimmedPV->at(0))  << std::endl;

  //*************************************************************//
  //                                                             //
  //--------------------------- Pile Up -------------------------//
  //                                                             //
  //*************************************************************//

  PU_Weight = -1.;
  float npT = -1.;

  if(!runningOnData_){
    edm::Handle<std::vector< PileupSummaryInfo>>  PupInfo;
    iEvent.getByToken(pileupSummaryToken_, PupInfo);
  
    std::vector<PileupSummaryInfo>::const_iterator PVI; 
 
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      const int BX = PVI->getBunchCrossing();
      if(BX == 0) {
	npT  = PVI->getTrueNumInteractions();
      }
    }

    if(npT == -1) {
      std::cout << "!!!! npT = -1 !!!!" << std::endl;
      abort();
    }
    // Calculate weight using above code
    PU_Weight = Lumiweights_.weight(npT);

    // Fill histogram with PU distribution
    h_pileup->Fill(npT);
  }

  //*************************************************************//
  //                                                             //
  //-------------------------- MC Weight ------------------------//
  //                                                             //
  //*************************************************************//

  MC_Weight = -10000000.;

  if(!runningOnData_){
    edm::Handle<GenEventInfoProduct> GenInfo;
    iEvent.getByToken(GenInfoToken_, GenInfo);
    
    float _aMCatNLOweight = GenInfo->weight();
    MC_Weight = _aMCatNLOweight;

    if(MC_Weight == -10000000.) {
      std::cout << "!!!! MC_Weight = -10000000 !!!!" << std::endl;
      abort();
    }
  }



  //*************************************************************//
  //                                                             //
  //--------------------------- Trigger -------------------------//
  //                                                             //
  //*************************************************************//

  //Examine the trigger information
  isSingleMuTrigger_24 = false;
  isSingleMuTrigger_27 = false;
  isSingleMuTrigger_50 = false;
  isSingleEleTrigger_25 = false;
  isSingleEleTrigger_27 = false;
  isSingleEleTrigger_32_DoubleEG = false;
  isSingleEleTrigger_32 = false;

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i){
    if(!triggerBits->accept(i)) continue;
    std::string tmp_triggername = names.triggerName(i);
    if( tmp_triggername.find("HLT_IsoMu24_v") != std::string::npos ||
	tmp_triggername.find("HLT_IsoTkMu24_v") != std::string::npos){
      isSingleMuTrigger_24 = true;
    }
    if( tmp_triggername.find("HLT_IsoMu27_v") != std::string::npos ){
      isSingleMuTrigger_27 = true;
    }
    if( tmp_triggername.find("HLT_Mu50_v") != std::string::npos ){
      isSingleMuTrigger_50 = true;
    }
    if( tmp_triggername.find("HLT_Ele25_eta2p1_WPTight_Gsf_v") != std::string::npos){
      isSingleEleTrigger_25 = true;
    }
    if( tmp_triggername.find("HLT_Ele27_WPTight_Gsf_v") != std::string::npos){
      isSingleEleTrigger_27 = true;
    }
    if( tmp_triggername.find("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v") != std::string::npos){
      isSingleEleTrigger_32_DoubleEG = true;
    }
    if( tmp_triggername.find("HLT_Ele32_WPTight_Gsf_v") != std::string::npos){
      isSingleEleTrigger_32 = true;
    }
  }


  //*************************************************************//
  //                                                             //
  //------------------ Variable initialization ------------------//
  //                                                             //
  //*************************************************************//

  nMuons          = 0;
  nElectrons      = 0;
  nElectronsLoose = 0;
  nPions          = 0;
  nPhotons        = 0;
  nBjets          = 0;
  nBjets_25       = 0;

  is_pi_a_pi         = false;
  is_pi_matched      = false;
  is_photon_a_photon = false;
  is_photon_matched  = false;

  pTmuMax  = -1000.;
  pTeleMax = -1000.;
  eTphMax  = -1000.;

  mu_pT  = 0.;
  mu_eta = 0.;
  mu_phi = 0.;
  mu_ID  = 0;
  mu_dxy = 0.;
  mu_dz  = 0.;
  mu_iso = 0.;
  best_mu_iso = 0.;

  el_pT    = 0.;
  el_eta   = 0.;
  el_etaSC = 0.;
  el_phi   = 0.;
  el_ID    = 0;
  el_dxy   = 0.;
  el_dz    = 0.;
  el_iso   = 0.;
  best_el_iso = 0.;
  LorentzVector el_p4;

  deltaphi_lep_pi = 0.;

  //These variables will go in the tree
  lepton_iso_tree = 0.;
  ph_iso_ChargedHadron = 0.;
  ph_iso_NeutralHadron = 0.;
  ph_iso_Photon = 0.;
  //ph_iso_Track = 0.;
  ph_iso_eArho = 0.;

  pi_pT     = 0.;
  pi_eta    = 0.;
  pi_phi    = 0.;
  pi_energy = 0.;
  pi_dxy    = 0.;
  pi_dz     = 0.;
  LorentzVector pi_p4;

  ph_eT     = 0.;
  ph_eta    = 0.;
  ph_etaSC  = 0.;
  ph_phi    = 0.;
  ph_energy = 0.;
  LorentzVector ph_p4;

  _Wmass = 0.;

  met_pT = 0.;

  is_muon = false;
  is_ele  = false;

  lepton_pT_tree  = 0.;
  lepton_eta_tree = 0.;
  lepton_etaSC_tree = 0.;
  lepton_phi_tree = 0.;
  lepton_dxy_tree = 0.;
  lepton_dz_tree  = 0.;

  gen_ph_pT_tree     = 0.;
  gen_ph_mother_tree = 0;


  //*************************************************************//
  //                                                             //
  //----------------------------- MET ---------------------------//
  //                                                             //
  //*************************************************************//

  for(auto met = slimmedMETs->begin(); met != slimmedMETs->end(); ++met){
    met_pT = met->pt();
  }

  for(auto metpuppi = slimmedMETsPuppi->begin(); metpuppi != slimmedMETsPuppi->end(); ++metpuppi){
    metpuppi_pT = metpuppi->pt();
  }

  //*************************************************************//
  //                                                             //
  //---------------------------- Muons --------------------------//
  //                                                             //
  //*************************************************************//

  for(auto mu = slimmedMuons->begin(); mu != slimmedMuons->end(); ++mu){
    if(mu->pt() < 20. || !mu->CutBasedIdMedium || fabs(mu->eta()) > 2.4 || fabs(mu->muonBestTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(mu->muonBestTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;
    //mu_iso = (mu->chargedHadronIso() + std::max(0., mu->neutralHadronIso() + mu->photonIso() - 0.5*mu->puChargedHadronIso()))/mu->pt();
    if(!mu->PFIsoLoose) continue;

    is_muon = true;
    nMuons++;
    
    if(mu->pt() > pTmuMax){

      mu_ID   = mu->pdgId();
      mu_eta  = mu->eta();
      mu_phi  = mu->phi();
      mu_pT   = mu->pt();
      mu_dxy  = mu->muonBestTrack()->dxy((&slimmedPV->at(0))->position());
      mu_dz   = mu->muonBestTrack()->dz((&slimmedPV->at(0))->position());
      best_mu_iso = mu_iso; //Save the value of mu_iso of the best candidate (highest pT) muon passing the selection
      
      pTmuMax = mu_pT;
    }
  }

  if(nMuons > 1) return;
  _Nevents_muVeto++;


  //*************************************************************//
  //                                                             //
  //-------------------------- Electrons ------------------------//
  //                                                             //
  //*************************************************************//

  // Get rho value
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;

  float corr_pt = 0.;

  for(auto el = slimmedElectrons->begin(); el != slimmedElectrons->end(); ++el){
    //Calculate electron p4, correct it with the Scale&Smearing correction and extract the pT
    el_p4 = el->p4() * el->userFloat("ecalTrkEnergyPostCorr")/el->energy();
    corr_pt = el_p4.pt();

      //if(el->pt() < 20. || fabs(el->eta()) > 2.5 || fabs(el->gsfTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(el->gsfTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;
      if(corr_pt < 20. || fabs(el->eta()) > 2.5 || fabs(el->gsfTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(el->gsfTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;

    //PflowIsolationVariables pfIso = el->pfIsolationVariables();
    float abseta = fabs(el->superCluster()->eta());
    float eA     = effectiveAreas_el_.getEffectiveArea(abseta);

    el_iso   = (el->pfIsolationVariables().sumChargedHadronPt + std::max( 0.0f, el->pfIsolationVariables().sumNeutralHadronEt + el->pfIsolationVariables().sumPhotonEt - eA*rho_))/corr_pt;

    if(el_iso > 0.35) continue;

    //-------------Conditions on loose/medium MVA electron ID-------------//
    if(!is_muon){
      if(el->electronID("mvaEleID-Fall17-iso-V2-wpLoose") == 0) continue;
      nElectronsLoose++;
    }

    if(el->electronID("mvaEleID-Fall17-iso-V2-wp90") == 0) continue;
    
    is_ele = true;
    nElectrons++;
    
    if(corr_pt > pTeleMax){
      
      el_ID    = el->pdgId();
      el_eta   = el->eta();
      el_etaSC = el->superCluster()->eta();
      el_phi   = el->phi();
      el_pT    = corr_pt;
      el_dxy   = el->gsfTrack()->dxy((&slimmedPV->at(0))->position());
      el_dz    = el->gsfTrack()->dz((&slimmedPV->at(0))->position());
      best_el_iso = el_iso; //Save the value of el_iso of the best candidate (highest pT) electron passing the selection
      pTeleMax = el_pT;
    } 
  }


  //*************************************************************//
  //                                                             //
  //------------------------- Lepton VETO -----------------------//
  //                                                             //
  //*************************************************************//
  

  if((!is_muon && nElectronsLoose > 1) || (!is_muon && nElectronsLoose==1 && nElectrons!=1)) return;
     _Nevents_eleVeto++;

  //Do NOT continue if you didn't find either a muon or an electron
  if(!is_muon && !is_ele) return;
  _Nevents_isLepton++;

  if(is_muon && is_ele) return;
  _Nevents_TwoLepton++;

  if(is_muon){
    lepton_pT_tree  = mu_pT;
    lepton_eta_tree = mu_eta;
    lepton_phi_tree = mu_phi;
    lepton_dxy_tree = mu_dxy;
    lepton_dz_tree  = mu_dz;
    lepton_iso_tree = best_mu_iso;
  }

  if(!is_muon && is_ele){
    lepton_pT_tree    = el_pT;
    lepton_eta_tree   = el_eta;
    lepton_etaSC_tree = el_etaSC;
    lepton_phi_tree   = el_phi;
    lepton_dxy_tree   = el_dxy;
    lepton_dz_tree    = el_dz;
    lepton_iso_tree   = best_el_iso;
  }


  //*************************************************************//
  //                                                             //
  //----------------------- Access MC Truth ---------------------//
  //                                                             //
  //*************************************************************//

  //In signal, identify if there's a real mu or ele from W
  is_Ele_signal       = false;
  is_Mu_signal        = false;
  is_Wplus_from_t     = false;
  is_Wminus_from_tbar = false;
  is_Wminus_in_lep    = false;
  is_Wplus_in_lep     = false;
  is_signal_Wplus     = false;
  is_signal_Wminus    = false;
  is_ttbar_lnu        = false;
  Wplus_pT = -999.;
  Wminus_pT = -999.;

  if(!runningOnData_){
    for(auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){
      if(fabs(gen->pdgId() ) == 11 && fabs(gen->mother()->pdgId()) == 24) is_Ele_signal = true;
      if(fabs(gen->pdgId() ) == 13 && fabs(gen->mother()->pdgId()) == 24) is_Mu_signal = true;
      if(gen->numberOfDaughters() != 2) continue; //Avoid Ws <<decaying>> in another W (especially in Signal)
      if(gen->pdgId() == 24){// && gen->mother()->pdgId() == 6){
	is_Wplus_from_t = true;
	Wplus_pT = gen->pt();
      }
      if(gen->pdgId() == -24){// && gen->mother()->pdgId() == -6){
	is_Wminus_from_tbar = true;
	Wminus_pT = gen->pt();
      }
      //if(fabs(gen->pdgId() != 24) || gen->numberOfDaughters() != 2) continue;
      for (int i=0; i<2; i++){
      	if(gen->daughter(i)->pdgId() == 11 || gen->daughter(i)->pdgId() == 13 || gen->daughter(i)->pdgId() == 15) is_Wminus_in_lep = true;
      	if(gen->daughter(i)->pdgId() == -11 || gen->daughter(i)->pdgId() == -13 || gen->daughter(i)->pdgId() == -15) is_Wplus_in_lep = true;
      	if(gen->daughter(i)->pdgId() == 22 && is_Wplus_from_t && !is_signal_Wminus) is_signal_Wplus = true;
      	if(gen->daughter(i)->pdgId() == 22 && is_Wminus_from_tbar && !is_signal_Wplus) is_signal_Wminus = true;
      }
    }
    if(is_Wplus_from_t && is_Wminus_from_tbar && is_Wminus_in_lep && is_Wplus_in_lep) is_ttbar_lnu = true;
  }
  

  //*************************************************************//
  //                                                             //
  //---------------------------- Pions --------------------------//
  //                                                             //
  //*************************************************************//

  bool cand_pion_found = false;
  are_lep_pi_opposite_charge = false;
  sum_pT_03    = 0.;
  sum_pT_05    = 0.;
  sum_pT_05_ch = 0.;

  pTpiMax  = -1000.;

  for(auto cand = PFCandidates->begin(); cand != PFCandidates->end(); ++cand){
    //if(cand->pdgId()*mu_ID < 0 || cand->pdgId()*el_ID < 0) continue; // WARNING: this condition works only if paired with muon/electron veto
    
    if(cand->pt() < 20. || !cand->trackHighPurity() || fabs(cand->dxy()) >= 0.2 || fabs(cand->dz()) >= 0.5 ) continue;

    if(is_muon){
      deltaphi_lep_pi = fabs(mu_phi-cand->phi());
      if(deltaphi_lep_pi > 3.14) deltaphi_lep_pi = 6.28-deltaphi_lep_pi;
    } 
    if(!is_muon && is_ele){
      deltaphi_lep_pi = fabs(el_phi-cand->phi());
      if(deltaphi_lep_pi > 3.14) deltaphi_lep_pi = 6.28-deltaphi_lep_pi;
    } 

    if(deltaphi_lep_pi < 0.00005) continue;

    if(cand->pt() < pTpiMax) continue;

    pTpiMax = cand->pt();
    cand_pion_found = true;
    nPions++;

    pi_pT     = cand->pt();
    pi_eta    = cand->eta();
    pi_phi    = cand->phi();
    pi_energy = cand->energy();
    pi_p4     = cand->p4();
    pi_dxy    = cand->dxy();
    pi_dz     = cand->dz();

    if(cand->pdgId()*mu_ID > 0 || cand->pdgId()*el_ID > 0){are_lep_pi_opposite_charge = true;}

    is_pi_a_pi = false;
    is_pi_matched = false;

    float deltapTMax = 10000.;
    const float deltaRMax  = 0.3;
    int   gen_mother = 0;
    int   gen_ID = 0;

    if(!runningOnData_){
      for (auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){ // Matching candidate for W reconstruction with MC truth
	
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


  //*************************************************************//
  //                                                             //
  //----------------------- Pion isolation ----------------------//
  //                                                             //
  //*************************************************************//

  for(auto cand_iso = PFCandidates->begin(); cand_iso != PFCandidates->end(); ++cand_iso){
    float deltaR = sqrt((pi_eta-cand_iso->eta())*(pi_eta-cand_iso->eta())+(pi_phi-cand_iso->phi())*(pi_phi-cand_iso->phi()));
    if(deltaR <= 0.3 && deltaR >= 0.02) sum_pT_03 += cand_iso->pt();
    if(deltaR <= 0.5 && deltaR >= 0.02) sum_pT_05 += cand_iso->pt();
    if(cand_iso->charge() != 0 && (fabs(cand_iso->dxy()) >= 0.2 || fabs(cand_iso->dz()) >= 0.5) ) continue; // Requesting charged particles to come from PV
    if(deltaR <= 0.5 && deltaR >= 0.02) sum_pT_05_ch += cand_iso->pt();
  }


  //*************************************************************//
  //                                                             //
  //--------------------------- Photons -------------------------//
  //                                                             //
  //*************************************************************//

  bool cand_photon_found = false;
  float corr_et = 0.;

  for(auto photon = slimmedPhotons->begin(); photon != slimmedPhotons->end(); ++photon){

    corr_et = photon->et() * photon->userFloat("ecalEnergyPostCorr")/photon->energy();

    // Apply energy scale corrections to MC (not available for 2016)
    // if(!runningOnData_){
    //   corr_et = photon->et() * photon->userFloat("ecalEnergyPostCorr")/photon->energy();
    // }

    if(corr_et < 20. || fabs(photon->eta()) > 2.5 || corr_et < eTphMax) continue;
    if(photon->hasPixelSeed()) continue;   //electron veto

    if(photon->photonID("mvaPhoID-RunIIFall17-v2-wp90") == 0) continue;

    float abseta = fabs(photon->superCluster()->eta());
    float eA = effectiveAreas_ph_.getEffectiveArea(abseta);
    //photon_iso = (pfIso.sumChargedHadronPt + std::max( 0.0f, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - eA*rho_))/photon->et();
    if(photon->chargedHadronIso()/corr_et > 0.3 || photon->photonIso() > 4.) continue; //|| photon->trackIso() > 6
    ph_iso_ChargedHadron = photon->chargedHadronIso();
    ph_iso_NeutralHadron = photon->neutralHadronIso();
    ph_iso_Photon        = photon->photonIso();
    //ph_iso_Track         = photon->trackIso();
    ph_iso_eArho         = eA*rho_;

    eTphMax = corr_et;

    ph_eT     = corr_et;
    ph_eta    = photon->eta();
    ph_etaSC  = photon->superCluster()->eta();
    ph_phi    = photon->phi();
    //ph_energy = photon->energy();
    //ph_p4     = photon->p4();

    // Apply energy scale corrections to MC
    ph_energy = photon->userFloat("ecalEnergyPostCorr");
    ph_p4     = photon->p4() * photon->userFloat("ecalEnergyPostCorr")/photon->energy();
    

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


  //*************************************************************//
  //                                                             //
  //--------------------- Photons in MC Truth -------------------//
  //                                                             //
  //*************************************************************//

  float gen_ph_pT_min = -9999.;
  float gen_ph_pT = 0.;
  int gen_ph_mother = 0;
  is_gen_ph = false;

  if(!runningOnData_){
    for (auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){
      if(gen->pdgId() != 22 || gen->pt() < 20. || fabs(gen->mother()->pdgId()) == 11) continue;
      if(gen->pt() < gen_ph_pT_min) continue;
      gen_ph_pT_min = gen->pt();
      gen_ph_pT = gen->pt();
      gen_ph_mother = gen->mother()->pdgId();
      is_gen_ph = true;
    }
  }

  if(is_gen_ph){
    gen_ph_pT_tree = gen_ph_pT;
    gen_ph_mother_tree = gen_ph_mother;
  }


  //*************************************************************//
  //                                                             //
  //--------------------------- b-jets --------------------------//
  //                                                             //
  //*************************************************************//

  nBjets = 0;
  for (auto jet = slimmedJets->begin(); jet != slimmedJets->end(); ++jet){
    if(jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") < Bjets_WP_) continue;   //0.46 = loose
    if(jet->pt() < 25.) continue;
    nBjets_25++;
    if(jet->pt() < 30.) continue;
    nBjets++;
  }

  mytree->Fill();

}


//*************************************************************//
//                                                             //
//---------------------- Create the tree ----------------------//
//                                                             //
//*************************************************************//

void WPiGammaAnalysis::create_trees()
{
  mytree = fs->make<TTree>("mytree", "Tree containing gen&reco");

  mytree->Branch("nPV",&nPV);
  mytree->Branch("isSingleMuTrigger_24",&isSingleMuTrigger_24);
  mytree->Branch("isSingleMuTrigger_27",&isSingleMuTrigger_27);
  mytree->Branch("isSingleMuTrigger_50",&isSingleMuTrigger_50);
  mytree->Branch("isSingleEleTrigger_25",&isSingleEleTrigger_25);
  mytree->Branch("isSingleEleTrigger_27",&isSingleEleTrigger_27);
  mytree->Branch("isSingleEleTrigger_32_DoubleEG",&isSingleEleTrigger_32_DoubleEG);
  mytree->Branch("isSingleEleTrigger_32",&isSingleEleTrigger_32);

  //Save run number info when running on data
  if(runningOnData_){
    mytree->Branch("run_number",&run_number);
  }
  
  mytree->Branch("LepPiOppositeCharge",&are_lep_pi_opposite_charge);
  mytree->Branch("lepton_pT",&lepton_pT_tree);
  mytree->Branch("lepton_eta",&lepton_eta_tree);
  mytree->Branch("lepton_etaSC",&lepton_etaSC_tree);
  mytree->Branch("lepton_phi",&lepton_phi_tree);
  mytree->Branch("lepton_dxy",&lepton_dxy_tree);
  mytree->Branch("lepton_dz",&lepton_dz_tree);
  mytree->Branch("lepton_iso",&lepton_iso_tree);
  mytree->Branch("is_muon",&is_muon);
  mytree->Branch("pi_pT",&pi_pT);
  mytree->Branch("pi_eta",&pi_eta);
  mytree->Branch("pi_phi",&pi_phi);
  mytree->Branch("pi_energy",&pi_energy);
  mytree->Branch("pi_dxy",&pi_dxy);
  mytree->Branch("pi_dz",&pi_dz);
  mytree->Branch("sum_pT_03",&sum_pT_03);
  mytree->Branch("sum_pT_05",&sum_pT_05);
  mytree->Branch("sum_pT_05_ch",&sum_pT_05_ch);
  mytree->Branch("photon_eT",&ph_eT);
  mytree->Branch("photon_eta",&ph_eta);
  mytree->Branch("photon_etaSC",&ph_etaSC);
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
  mytree->Branch("metpuppi_pT",&metpuppi_pT);

  //Save MC info
  if(!runningOnData_){
    mytree->Branch("PU_Weight",&PU_Weight);
    mytree->Branch("MC_Weight",&MC_Weight);

    mytree->Branch("isMuonSignal",&is_Mu_signal);
    mytree->Branch("isEleSignal",&is_Ele_signal);
    mytree->Branch("isPionTrue",&is_pi_a_pi);
    mytree->Branch("isPionMatched",&is_pi_matched);
    mytree->Branch("isPhotonTrue",&is_photon_a_photon);
    mytree->Branch("isPhotonMatched",&is_photon_matched);
    mytree->Branch("isttbarlnu",&is_ttbar_lnu);
    mytree->Branch("Wplus_pT",&Wplus_pT);
    mytree->Branch("Wminus_pT",&Wminus_pT);
    mytree->Branch("is_signal_Wplus",&is_signal_Wplus);
    mytree->Branch("is_signal_Wminus",&is_signal_Wminus);
    mytree->Branch("is_gen_ph",&is_gen_ph);
    mytree->Branch("gen_ph_pT",&gen_ph_pT_tree);
    mytree->Branch("gen_ph_mother",&gen_ph_mother_tree);
  }

}

void WPiGammaAnalysis::beginJob()
{
  //Flag for PileUp reweighting
  if (!runningOnData_ && !runningOn2017_){ // PU reweighting for 2016
   Lumiweights_ = edm::LumiReWeighting("PU/MCpileUp_2016_25ns_Moriond17MC_PoissonOOTPU.root", "PU/MyDataPileupHistogram_2016.root", "pileup", "pileup");
  }
  if (!runningOnData_ && runningOn2017_){ // PU reweighting for 2017
   Lumiweights_ = edm::LumiReWeighting("PU/MCpileUp_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU.root", "PU/MyDataPileupHistogram_2017.root", "pileup", "pileup");
  }
}


//*************************************************************//
//                                                             //
//------------------- Fill event loss histos ------------------//
//                                                             //
//*************************************************************//

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
