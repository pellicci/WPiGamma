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
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
  
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
  runningEra_(iConfig.getParameter<int>("runningEra")),
  verboseIdFlag_(iConfig.getParameter<bool>("phoIdVerbose")),
  effectiveAreas_el_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile_el")).fullPath() ),
  effectiveAreas_ph_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile_ph")).fullPath() ),
  Bjets_WP_2016_(iConfig.getParameter<double>("Bjets_WP_2016")),
  Bjets_WP_2017_(iConfig.getParameter<double>("Bjets_WP_2017")),
  Bjets_WP_2018_(iConfig.getParameter<double>("Bjets_WP_2018"))

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
  triggerObjectsTokenMC2016_          = consumes<std::vector<pat::TriggerObjectStandAlone> > (edm::InputTag("slimmedPatTrigger","","PAT"));
  triggerObjectsTokenData2016_        = consumes<std::vector<pat::TriggerObjectStandAlone> > (edm::InputTag("slimmedPatTrigger","","DQM"));
  triggerObjectsToken2017_            = consumes<std::vector<pat::TriggerObjectStandAlone> > (edm::InputTag("slimmedPatTrigger","","PAT"));
  triggerObjectsTokenMC2018_          = consumes<std::vector<pat::TriggerObjectStandAlone> > (edm::InputTag("slimmedPatTrigger","","PAT"));
  triggerObjectsTokenData2018_        = consumes<std::vector<pat::TriggerObjectStandAlone> > (edm::InputTag("slimmedPatTrigger","","RECO"));
  rhoToken_                           = consumes<double> (iConfig.getParameter <edm::InputTag>("rho"));
  PrefiringWeightToken_               = consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"));

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

  h2_BTaggingEff_Num_b   = fs->make<TH2D>("h2_BTaggingEff_Num_b", ";p_{T} [GeV];#eta", 75, 25., 1000., 50, -2.5, 2.5);
  h2_BTaggingEff_Denom_b = fs->make<TH2D>("h2_BTaggingEff_Denom_b", ";p_{T} [GeV];#eta", 75, 25., 1000., 50, -2.5, 2.5);


  create_trees();
}

WPiGammaAnalysis::~WPiGammaAnalysis()
{
}

//--------- match reco electron with trigger object (for 2017 trigger) ---------//
namespace{
  std::vector<const pat::TriggerObjectStandAlone*> getMatchedObjs(const float eta,const float phi,const std::vector<pat::TriggerObjectStandAlone>& trigObjs,const float maxDeltaR=0.1)
  {
    std::vector<const pat::TriggerObjectStandAlone*> matchedObjs;
    const float maxDR2 = maxDeltaR*maxDeltaR;
    for(auto& trigObj : trigObjs){
      const float dR2 = reco::deltaR2(eta,phi,trigObj.eta(),trigObj.phi());
      if(dR2<maxDR2) matchedObjs.push_back(&trigObj);
    }
    return matchedObjs;
  }
}
//------------------------------------------------------------------------------//

// ------------ method called for each event  ------------
void WPiGammaAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<pat::PackedCandidate>  > PFCandidates;
  iEvent.getByToken(packedPFCandidatesToken_, PFCandidates);

  edm::Handle<std::vector<pat::Muon>  > slimmedMuons;
  iEvent.getByToken(slimmedMuonsToken_, slimmedMuons);

  edm::Handle<std::vector<reco::GenParticle>  > genParticles;
  if(!runningOnData_){iEvent.getByToken(prunedGenParticlesToken_, genParticles);}

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

  edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
  if(!runningOnData_ && runningEra_ == 0){
    iEvent.getByToken(triggerObjectsTokenMC2016_, triggerObjects);
  }
  if(runningOnData_ && runningEra_ == 0){
    iEvent.getByToken(triggerObjectsTokenData2016_, triggerObjects);
  }
  std::vector<pat::TriggerObjectStandAlone> unpackedTrigObjs;
  if(runningEra_ == 1){
    iEvent.getByToken(triggerObjectsToken2017_, triggerObjects);
    //Create unpacked filter names to convert Ele32_DoubleEG trigger into Ele32_WPTight
    for(auto& trigObj : *triggerObjects){
      unpackedTrigObjs.push_back(trigObj);
      unpackedTrigObjs.back().unpackFilterLabels(iEvent,*triggerBits);
    }
  }
  if(!runningOnData_ && runningEra_ == 2){
    iEvent.getByToken(triggerObjectsTokenMC2018_, triggerObjects);
  }
  if(runningOnData_ && runningEra_ == 2){
    iEvent.getByToken(triggerObjectsTokenData2018_, triggerObjects);
  }

  Prefiring_Weight = -10000;
  edm::Handle<double> PrefiringWeight;
  if(!runningOnData_ && (runningEra_ == 0 || runningEra_ == 1)){
    iEvent.getByToken(PrefiringWeightToken_, PrefiringWeight);
    Prefiring_Weight = (*PrefiringWeight);
  }


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
    // if( tmp_triggername.find("HLT_Ele25_eta2p1_WPTight_Gsf_v") != std::string::npos){
    //   isSingleEleTrigger_25 = true;
    // }
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


  BTagCalibration calib("DeepCSV", "DeepCSV_2016LegacySF_WP_V1.csv");
  BTagCalibrationReader reader(BTagEntry::OP_LOOSE,  // operating point
			       "central",            // central sys type
			       {"up", "down"});      // other sys types
    
  reader.load(calib,              // calibration instance
	       BTagEntry::FLAV_B,  // btag flavour
	       "comb");            // measurement type

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
  nBjets_scaled   = 0;

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
  //----------- Highest pT mu, ele, gamma in MC Truth -----------//
  //                                                             //
  //*************************************************************//
  
  MCT_HpT_mu_pT_Max  = -1000.;
  MCT_HpT_mu_pT      = -1000.;
  MCT_HpT_mu_eta     = -1000.;
  MCT_HpT_mu_phi     = -1000.;
  MCT_HpT_ele_pT_Max = -1000.;
  MCT_HpT_ele_pT     = -1000.;
  MCT_HpT_ele_eta    = -1000.;
  MCT_HpT_ele_phi    = -1000.;
  MCT_HeT_ph_eT_Max  = -1000.;
  MCT_HeT_ph_eT      = -1000.;
  MCT_HeT_ph_eta     = -1000.;
  MCT_HeT_ph_phi     = -1000.;

  if(!runningOnData_){
    for(auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){
      if(fabs(gen->pdgId()) == 13 && gen->status() == 1 && gen->pt() > MCT_HpT_mu_pT_Max){
	MCT_HpT_mu_pT  = gen->pt();
	MCT_HpT_mu_eta = gen->eta();
	MCT_HpT_mu_phi = gen->phi();
	MCT_HpT_mu_pT_Max = MCT_HpT_mu_pT;
      }
      if(fabs(gen->pdgId()) == 11 && gen->status() == 1 && gen->pt() > MCT_HpT_ele_pT_Max){
	MCT_HpT_ele_pT  = gen->pt();
	MCT_HpT_ele_eta = gen->eta();
	MCT_HpT_ele_phi = gen->phi();
	MCT_HpT_ele_pT_Max = MCT_HpT_ele_pT;
      }
      if(fabs(gen->pdgId()) == 22 && gen->status() == 1 && gen->et() > MCT_HeT_ph_eT_Max){
	MCT_HeT_ph_eT  = gen->et();
	MCT_HeT_ph_eta = gen->eta();
	MCT_HeT_ph_phi = gen->phi();
	MCT_HeT_ph_eT_Max = MCT_HeT_ph_eT;
      }
    }
  }

  //*************************************************************//
  //                                                             //
  //----------------------------- MET ---------------------------//
  //                                                             //
  //*************************************************************//
  met_pT_scaled = 0.;

  for(auto met = slimmedMETs->begin(); met != slimmedMETs->end(); ++met){

    if(met->phi() > -1.57 && met->phi() < 0.87 && met->eta() > -2.5 && met->eta() < -1.3){
      met_pT_scaled = 0.8*met->pt();
    }
    else if(met->phi() > -1.57 && met->phi() < 0.87 && met->eta() > -3.0 && met->eta() < -2.5){
      met_pT_scaled = 0.65*met->pt();
    }
    else if( (met->phi() > -1.57 && met->phi() < 0.87 && (met->eta() <= -3.0 || met->eta() >= -1.3)) || met->phi() <= -1.57 || met->phi() >= 0.87){
      met_pT_scaled = met->pt();
    }
    
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

      // // //Now check muon trigger matching
      // // if(runningEra_ == 0 && (mu->triggered("HLT_IsoMu24_v*") || mu->triggered("HLT_IsoTkMu24_v*") || mu->triggered("HLT_Mu50_v*"))) isTriggerMatched = true;
      // // if(runningEra_ == 1 && (mu->triggered("HLT_IsoMu27_v*") || mu->triggered("HLT_Mu50_v*"))) isTriggerMatched = true;
      // // if(runningEra_ == 2 && (mu->triggered("HLT_IsoMu24_v*") || mu->triggered("HLT_Mu50_v*"))) isTriggerMatched = true;
      
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

      if(corr_pt < 20. || fabs(el->eta()) > 2.5 || fabs(el->gsfTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(el->gsfTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;

    //PflowIsolationVariables pfIso = el->pfIsolationVariables();
    float abseta = fabs(el->superCluster()->eta());
    float eA     = effectiveAreas_el_.getEffectiveArea(abseta);

    el_iso   = (el->pfIsolationVariables().sumChargedHadronPt + std::max( 0.0f, el->pfIsolationVariables().sumNeutralHadronEt + el->pfIsolationVariables().sumPhotonEt - eA*rho_))/corr_pt;

    //if(el_iso > 0.35) continue;

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
      el_phiSC = el->superCluster()->phi();
      el_phi   = el->phi();
      el_pT    = corr_pt;
      el_dxy   = el->gsfTrack()->dxy((&slimmedPV->at(0))->position());
      el_dz    = el->gsfTrack()->dz((&slimmedPV->at(0))->position());
      best_el_iso = el_iso; //Save the value of el_iso of the best candidate (highest pT) electron passing the selection
      pTeleMax = el_pT; //el_pT has assumed the value of corr_pt a few lines above
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
    lepton_phiSC_tree = el_phiSC;
    lepton_phi_tree   = el_phi;
    lepton_dxy_tree   = el_dxy;
    lepton_dz_tree    = el_dz;
    lepton_iso_tree   = best_el_iso;
  }

  //Cut events with lepton pT lower than 25 GeV, since the lowest offline cut will be 25 GeV anyway
  if(lepton_pT_tree < 25.) return;

  //**********************************************************************//
  //                                                                      //
  //---------------------- Electron Trigger Matching ---------------------//
  //                                                                      //
  //**********************************************************************//
  //https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#Trigger

  float deltaR_lep_trigger_min = 0.1;
  bool isEle32_WPTight_equivalent = false;
  bool isTriggerMatched = false;

  const edm::TriggerNames &TNames = iEvent.triggerNames(*triggerBits);

  for (pat::TriggerObjectStandAlone obj : *triggerObjects){ // note: not "const &" since we want to call unpackPathNames
    
    bool isAcceptedPath = false;
    obj.unpackPathNames(TNames);
    
    std::vector pathNamesAll = obj.pathNames(false);
    
    for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
      // Record also if the object is associated to a 'l3' filter (always true for the definition used
      // in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which means
      // that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
      bool isSuccessfulTrigger = obj.hasPathName( pathNamesAll[h], true, true );
      //std::cout << "   " << pathNamesAll[h];
      if(!isSuccessfulTrigger) continue;
      //std::cout << pathNamesAll[h] << "(L,3)" << std::endl; 

      if((is_muon && runningEra_ == 0) && (pathNamesAll[h].find("HLT_IsoMu24_v") != std::string::npos || pathNamesAll[h].find("HLT_IsoTkMu24_v") != std::string::npos || pathNamesAll[h].find("HLT_Mu50_v") != std::string::npos )) {
	isAcceptedPath = true;
	continue;
      }

      if((is_muon && runningEra_ == 1) && (pathNamesAll[h].find("HLT_IsoMu27_v") != std::string::npos || pathNamesAll[h].find("HLT_Mu50_v") != std::string::npos)) {
	isAcceptedPath = true;
	continue;
      }

      if((is_muon && runningEra_ == 2) && (pathNamesAll[h].find("HLT_IsoMu24_v") != std::string::npos || pathNamesAll[h].find("HLT_Mu50_v") != std::string::npos)) {
	isAcceptedPath = true;
	continue;
      } 
      
      if((!is_muon && runningEra_ == 0) && pathNamesAll[h].find("HLT_Ele27_WPTight_Gsf_v") != std::string::npos) {
	isAcceptedPath = true;
	continue;
      }
      
      if((!is_muon && runningEra_ == 1) && pathNamesAll[h].find("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v") != std::string::npos) {
	//now match ALL objects in a cone of DR<0.1
	//it is important to match all objects as there are different ways to reconstruct the same electron
	//eg, L1 seeded, unseeded, as a jet etc
	//and so you want to be sure you get all possible objects
	std::vector<const pat::TriggerObjectStandAlone*> matchedTrigObjs = getMatchedObjs(lepton_etaSC_tree,lepton_phiSC_tree,unpackedTrigObjs,0.1);
	for(const auto trigObj : matchedTrigObjs){
	  //now just check if it passes the two filters
	  if(trigObj->hasFilterLabel("hltEle32L1DoubleEGWPTightGsfTrackIsoFilter") && trigObj->hasFilterLabel("hltEGL1SingleEGOrFilter") ) {
	    isEle32_WPTight_equivalent = true;
	    break;
	  }
	}
	if(isEle32_WPTight_equivalent){
	  //std::cout << "Ele32 EMULATED" << std::endl;
	  isAcceptedPath = true;
	  continue;
	}
	else{
	  return;
	}
      }
      
      if((!is_muon && runningEra_ == 2) &&  pathNamesAll[h].find("HLT_Ele32_WPTight_Gsf_v") != std::string::npos) {
	isAcceptedPath = true;
	continue;
      }
    }
    //std::cout << std::endl;
    
    if(!isAcceptedPath) continue;
    
    float deltaEta_lep_trigger;
    float deltaPhi_lep_trigger;

    if(is_muon){
    deltaEta_lep_trigger = fabs(lepton_eta_tree - obj.eta());
    deltaPhi_lep_trigger = fabs(lepton_phi_tree - obj.phi());
    }
    else{//Use SC eta and phi for electrons
    deltaEta_lep_trigger = fabs(lepton_etaSC_tree - obj.eta());
    deltaPhi_lep_trigger = fabs(lepton_phiSC_tree - obj.phi());
    }
    
    if(deltaPhi_lep_trigger > 3.14){
      deltaPhi_lep_trigger = 6.28 - deltaPhi_lep_trigger;
    }
    
    float deltaR_lep_trigger = sqrt(deltaEta_lep_trigger*deltaEta_lep_trigger + deltaPhi_lep_trigger*deltaPhi_lep_trigger);
    //float deltapTrel_lep_trigger = fabs(lepton_pT_tree - obj.pt())/obj.pt();
    if(deltaR_lep_trigger <= deltaR_lep_trigger_min){
      isTriggerMatched = true;
      //std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
      break;
    }
  }
  //if(!isTriggerMatched) return;
  isTriggerMatched_tree = isTriggerMatched;


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
      if(gen->numberOfDaughters() != 2) continue; //Avoid Ws <<decaying>> into another W (especially in Signal)
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

    if(cand->pt() < pTpiMax) continue;

    deltaphi_lep_pi = fabs(lepton_phi_tree-cand->phi());
    if(deltaphi_lep_pi > 3.14) deltaphi_lep_pi = 6.28 - deltaphi_lep_pi;
    if(deltaphi_lep_pi < 0.00005) continue;

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
	float deltaPhi = fabs(pi_phi - gen->phi());
	if(deltaPhi > 3.14){
	  deltaPhi = 6.28 - deltaPhi;
	}
	float deltaR = sqrt((pi_eta - gen->eta())*(pi_eta - gen->eta()) + deltaPhi*deltaPhi);
	float deltapT = fabs(pi_pT - gen->pt());

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
    float deltaPhi = fabs(pi_phi - cand_iso->phi());
    if(deltaPhi > 3.14){
      deltaPhi = 6.28 - deltaPhi;
    }
    float deltaR = sqrt((pi_eta - cand_iso->eta())*(pi_eta - cand_iso->eta()) + deltaPhi*deltaPhi);
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

    corr_et = photon->et()*photon->userFloat("ecalEnergyPostCorr")/photon->energy();

    if(corr_et < 20. || fabs(photon->eta()) > 2.5 || corr_et < eTphMax) continue;
    if(photon->hasPixelSeed()) continue;   //electron veto

    if(photon->photonID("mvaPhoID-RunIIFall17-v2-wp90") == 0) continue;

    float abseta = fabs(photon->superCluster()->eta());
    float eA = effectiveAreas_ph_.getEffectiveArea(abseta);
    //float photon_iso = (photon->IsolationVariables().sumChargedHadronPt + std::max( 0.0f, photon->IsolationVariables().sumNeutralHadronEt + photon->IsolationVariables().sumPhotonEt - eA*rho_))/photon->et();
    //if(photon->chargedHadronIso()/corr_et > 0.3 || photon->photonIso() > 4.) continue; //|| photon->trackIso() > 6
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
    ph_p4     = photon->p4()*photon->userFloat("ecalEnergyPostCorr")/photon->energy();
    

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
	float deltaPhi = fabs(ph_phi - gen->phi());
	if(deltaPhi > 3.14){
	  deltaPhi = 6.28 - deltaPhi;
	}
	float deltaR = sqrt((ph_eta - gen->eta())*(ph_eta - gen->eta()) + deltaPhi*deltaPhi);
	float deltapT = fabs(ph_eT - gen->pt());

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
      if(gen->pdgId() != 22 || gen->pt() < 20. || fabs(gen->mother()->pdgId()) == 11 || !gen->isPromptFinalState()) continue;
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

  nBjets = 0; //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
  float jet_pT_temp = 0.;

  for (auto jet = slimmedJets->begin(); jet != slimmedJets->end(); ++jet){

    if(runningEra_ == 1 && jet->eta() > 2.65 && jet->eta() < 3.139 && jet->pt() < 50.) continue; //Exclude noisy region in EE in 2017


    //---------------Fill b-tagging efficiency histograms---------------//
    if(!runningOnData_){
      int hadronFlavour = jet->hadronFlavour();
      
      if(fabs(hadronFlavour) == 5){
    	h2_BTaggingEff_Denom_b->Fill(jet->pt(), jet->eta());
      
      if(runningEra_ == 0 && (jet->bDiscriminator("pfDeepCSVJetTags:probb") + jet->bDiscriminator("pfDeepCSVJetTags:probbb")) >= Bjets_WP_2016_ && jet->pt() >= 25.) h2_BTaggingEff_Num_b->Fill(jet->pt(), jet->eta());
      if(runningEra_ == 1 && (jet->bDiscriminator("pfDeepCSVJetTags:probb") + jet->bDiscriminator("pfDeepCSVJetTags:probbb")) >= Bjets_WP_2017_ && jet->pt() >= 25.) h2_BTaggingEff_Num_b->Fill(jet->pt(), jet->eta());
      if(runningEra_ == 2 && (jet->bDiscriminator("pfDeepCSVJetTags:probb") + jet->bDiscriminator("pfDeepCSVJetTags:probbb")) >= Bjets_WP_2018_ && jet->pt() >= 25.) h2_BTaggingEff_Num_b->Fill(jet->pt(), jet->eta());
      }
    }
    //-----------------------------------------------------------------//

    if(runningEra_ == 0 && (jet->bDiscriminator("pfDeepCSVJetTags:probb") + jet->bDiscriminator("pfDeepCSVJetTags:probbb")) < Bjets_WP_2016_) continue; //loose 2016
    if(runningEra_ == 1 && (jet->bDiscriminator("pfDeepCSVJetTags:probb") + jet->bDiscriminator("pfDeepCSVJetTags:probbb")) < Bjets_WP_2017_) continue; //loose 2017
    if(runningEra_ == 2 && (jet->bDiscriminator("pfDeepCSVJetTags:probb") + jet->bDiscriminator("pfDeepCSVJetTags:probbb")) < Bjets_WP_2018_) continue; //loose 2018

    //----Scale down jet pT for 2018 HEM16/17 issue----//
    if(jet->phi() > -1.57 && jet->phi() < 0.87 && jet->eta() > -2.5 && jet->eta() < -1.3){
      jet_pT_temp = 0.8*jet->pt();
      if(jet_pT_temp >= 25.){
	nBjets_scaled++;
      }
    }
    else if(jet->phi() > -1.57 && jet->phi() < 0.87 && jet->eta() > -3.0 && jet->eta() < -2.5){
      jet_pT_temp = 0.65*jet->pt();
      if(jet_pT_temp >= 25.){
	nBjets_scaled++;
      }
    }
    else if( (jet->phi() > -1.57 && jet->phi() < 0.87 && (jet->eta() <= -3.0 || jet->eta() >= -1.3)) || jet->phi() <= -1.57 || jet->phi() >= 0.87){
      nBjets_scaled++;
    }
    //------------------------------------------------//

    // double jet_SF = reader.eval_auto_bounds("central", 
    // 					    BTagEntry::FLAV_B, 
    // 					    fabs(jet->eta()), // absolute value of eta
    // 					    jet->pt()); 

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

  mytree->Branch("runningEra",&runningEra_);
  mytree->Branch("nPV",&nPV);
  mytree->Branch("isSingleMuTrigger_24",&isSingleMuTrigger_24);
  mytree->Branch("isSingleMuTrigger_27",&isSingleMuTrigger_27);
  mytree->Branch("isSingleMuTrigger_50",&isSingleMuTrigger_50);
  mytree->Branch("isSingleEleTrigger_25",&isSingleEleTrigger_25);
  mytree->Branch("isSingleEleTrigger_27",&isSingleEleTrigger_27);
  mytree->Branch("isSingleEleTrigger_32_DoubleEG",&isSingleEleTrigger_32_DoubleEG);
  mytree->Branch("isSingleEleTrigger_32",&isSingleEleTrigger_32);
  mytree->Branch("isTriggerMatched",&isTriggerMatched_tree);

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
  // mytree->Branch("photon_iso_ChargedHadron",&ph_iso_ChargedHadron);
  // mytree->Branch("photon_iso_NeutralHadron",&ph_iso_NeutralHadron);
  // mytree->Branch("photon_iso_Photon",&ph_iso_Photon);
  // mytree->Branch("photon_iso_Track",&ph_iso_Track);
  // mytree->Branch("photon_iso_eArho",&ph_iso_eArho);

  mytree->Branch("Wmass",&_Wmass);

  mytree->Branch("nMuons",&nMuons);
  mytree->Branch("nElectrons",&nElectrons);
  mytree->Branch("nPions",&nPions);
  mytree->Branch("nPhotons",&nPhotons);
  mytree->Branch("nBjets",&nBjets);
  mytree->Branch("nBjets_25",&nBjets_25);
  mytree->Branch("nBjets_scaled",&nBjets_scaled);
  mytree->Branch("met_pT",&met_pT);
  mytree->Branch("met_pT_scaled",&met_pT_scaled);
  mytree->Branch("metpuppi_pT",&metpuppi_pT);

  //Save MC info
  if(!runningOnData_){
    mytree->Branch("PU_Weight",&PU_Weight);
    mytree->Branch("MC_Weight",&MC_Weight);
    if(runningEra_ == 0 || runningEra_ == 1){
      mytree->Branch("Prefiring_Weight",&Prefiring_Weight);
    }
    mytree->Branch("MCT_HpT_mu_pT",&MCT_HpT_mu_pT);
    mytree->Branch("MCT_HpT_mu_eta",&MCT_HpT_mu_eta);
    mytree->Branch("MCT_HpT_mu_phi",&MCT_HpT_mu_phi);
    mytree->Branch("MCT_HpT_ele_pT",&MCT_HpT_ele_pT);
    mytree->Branch("MCT_HpT_ele_eta",&MCT_HpT_ele_eta);
    mytree->Branch("MCT_HpT_ele_phi",&MCT_HpT_ele_phi);
    mytree->Branch("MCT_HeT_ph_eT",&MCT_HeT_ph_eT);
    mytree->Branch("MCT_HeT_ph_eta",&MCT_HeT_ph_eta);
    mytree->Branch("MCT_HeT_ph_phi",&MCT_HeT_ph_phi);
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
  if (!runningOnData_ && runningEra_ == 0){ // PU reweighting and b-tagging SF for 2016

   Lumiweights_ = edm::LumiReWeighting("MCpileUp_2016_25ns_Moriond17MC_PoissonOOTPU.root", "MyDataPileupHistogram_2016.root", "pileup", "pileup");

   // BTagCalibration calib("DeepCSV", "DeepCSV_2016LegacySF_WP_V1.csv");
   // BTagCalibrationReader reader(BTagEntry::OP_LOOSE,  // operating point
   // 				"central",            // central sys type
   // 				{"up", "down"});      // other sys types
   
   // reader.load(calib,              // calibration instance
   // 	       BTagEntry::FLAV_B,  // btag flavour
   // 	       "comb");            // measurement type
  }
  if (!runningOnData_ && runningEra_ == 1){ // PU reweighting and b-tagging SF for 2017

   Lumiweights_ = edm::LumiReWeighting("MCpileUp_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU.root", "MyDataPileupHistogram_2017.root", "pileup", "pileup");

   // BTagCalibration calib("DeepCSV", "DeepCSV_94XSF_WP_V4_B_F.csv");
   // BTagCalibrationReader reader(BTagEntry::OP_LOOSE,  // operating point
   // 				"central",            // central sys type
   // 				{"up", "down"});      // other sys types
   
   // reader.load(calib,             // calibration instance
   // 	       BTagEntry::FLAV_B, // btag flavour
   // 	       "comb");           // measurement type
  }
  if (!runningOnData_ && runningEra_ == 2){ // PU reweighting and b-tagging SF for 2018

   Lumiweights_ = edm::LumiReWeighting("MCpileUp_2018_25ns_JuneProjectionFull18_PoissonOOTPU.root", "MyDataPileupHistogram_2018.root", "pileup", "pileup");

   // BTagCalibration calib("DeepCSV", "DeepCSV_102XSF_WP_V1.csv");
   // BTagCalibrationReader reader(BTagEntry::OP_LOOSE,  // operating point
   // 				"central",            // central sys type
   // 				{"up", "down"});      // other sys types
   
   // reader.load(calib,             // calibration instance
   // 	       BTagEntry::FLAV_B, // btag flavour
   // 	       "comb");           // measurement type
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
