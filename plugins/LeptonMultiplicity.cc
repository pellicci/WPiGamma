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
 
#include "LeptonMultiplicity.h"
 
// constructors and destructor
LeptonMultiplicity::LeptonMultiplicity(const edm::ParameterSet& iConfig) :
  runningOnData_(iConfig.getParameter<bool>("runningOnData")),
  verboseIdFlag_(iConfig.getParameter<bool>("phoIdVerbose")),
  effectiveAreas_el_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile_el")).fullPath() ),
  effectiveAreas_ph_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile_ph")).fullPath() )
{
  packedPFCandidatesToken_            = consumes<std::vector<pat::PackedCandidate> >(edm::InputTag("packedPFCandidates")); 
  slimmedMuonsToken_                  = consumes<std::vector<pat::Muon> >(edm::InputTag("slimmedMuons"));
  prunedGenParticlesToken_            = consumes<std::vector<reco::GenParticle> >(edm::InputTag("prunedGenParticles"));
  electronsMiniAODToken_              = consumes<std::vector<pat::Electron> > (edm::InputTag("slimmedElectrons"));
  offlineSlimmedPrimaryVerticesToken_ = consumes<std::vector<reco::Vertex> > (edm::InputTag("offlineSlimmedPrimaryVertices"));  
  pileupSummaryToken_                 = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo"));
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

  create_trees();
}

LeptonMultiplicity::~LeptonMultiplicity()
{
}

// ------------ method called for each event  ------------
void LeptonMultiplicity::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<pat::PackedCandidate>  > PFCandidates;
  iEvent.getByToken(packedPFCandidatesToken_, PFCandidates);

  edm::Handle<std::vector<pat::Muon>  > slimmedMuons;
  iEvent.getByToken(slimmedMuonsToken_, slimmedMuons);

  edm::Handle<std::vector<reco::GenParticle>  > genParticles;
  if(!runningOnData_)iEvent.getByToken(prunedGenParticlesToken_, genParticles);

  edm::Handle<std::vector<pat::Electron> > slimmedElectrons;
  iEvent.getByToken(electronsMiniAODToken_,slimmedElectrons);

  edm::Handle<std::vector<reco::Vertex > > slimmedPV;
  iEvent.getByToken(offlineSlimmedPrimaryVerticesToken_, slimmedPV);

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsToken_, triggerBits);

  _Nevents_processed++;

  //Count the number of vertices
  nPV = -1;
  if(slimmedPV->size()<=0) return;
  for(reco::VertexCollection::const_iterator vtx=slimmedPV->begin();vtx!=slimmedPV->end();++vtx) {
    // check that the primary vertex is not a fake one, that is the beamspot (it happens when no primary vertex is reconstructed)
    if(!vtx->isFake()) {
      nPV++;
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

  nMuons          = 0;
  nElectrons      = 0;
  nElectronsLoose = 0;

  //float pTmuMax  = -1000.;
  //float pTeleMax = -1000.;

  mu_pT = 0.;
  mu_eta = 0.;
  mu_phi = 0.;
  //int   mu_ID = 0;

  el_pT = 0.;
  el_eta = 0.;
  el_phi = 0.;
  //int   el_ID = 0;

  mu_iso = 0.;
  el_iso = 0.;

  is_muon = false;
  is_ele  = false;

  //Loop over muons
  for(auto mu = slimmedMuons->begin(); mu != slimmedMuons->end(); ++mu){
    if(mu->pt() < 20. || !mu->CutBasedIdMedium || abs(mu->eta()) > 2.4 || fabs(mu->muonBestTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(mu->muonBestTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;
    if(!mu->PFIsoLoose) continue;
 
    is_muon  = true;
    mu_eta   = mu->eta();
    mu_phi   = mu->phi();
    mu_pT    = mu->pt();
    nMuons++;

  }


  // Get rho value
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;

  //Loop over electrons
  for(auto el = slimmedElectrons->begin(); el != slimmedElectrons->end(); ++el){

    if(el->pt() < 20. || fabs(el->eta()) > 2.5 || fabs(el->gsfTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(el->gsfTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;
    if(el->electronID("mvaEleID-Fall17-iso-V2-wp90") == 0) continue;

    // float abseta =  abs(el->superCluster()->eta());
    // float eA = effectiveAreas_el_.getEffectiveArea(abseta);
    // el_iso = (el->pfIsolationVariables().sumChargedHadronPt + std::max( 0.0f, el->pfIsolationVariables().sumNeutralHadronEt + el->pfIsolationVariables().sumPhotonEt - eA*rho_))/el->pt();
    // if(el_iso > 0.35) continue;

    is_ele      = true;
    el_eta      = el->eta();
    el_phi      = el->phi();
    el_pT       = el->pt();

    nElectrons++;

  }

  for(auto el = slimmedElectrons->begin(); el != slimmedElectrons->end(); ++el){
 
    if(el->pt() < 20. || fabs(el->eta()) > 2.5 || fabs(el->gsfTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(el->gsfTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;
 
    if(el->electronID("mvaEleID-Fall17-iso-V2-wpLoose") == 0) continue;

    // float abseta =  abs(el->superCluster()->eta());
    // float eA = effectiveAreas_el_.getEffectiveArea(abseta);
    // if((el->pfIsolationVariables().sumChargedHadronPt + std::max( 0.0f, el->pfIsolationVariables().sumNeutralHadronEt + el->pfIsolationVariables().sumPhotonEt - eA*rho_))/el->pt() > 0.35) continue;

    nElectronsLoose++;

  }


  //Do NOT continue if you didn't find either a muon or an electron
  if(!is_muon && !is_ele) return;

  is_Ele_signal = false;
  is_Mu_signal = false;
  if(!runningOnData_){
    for (auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){
      if(fabs(gen->pdgId() ) == 11 && fabs(gen->mother()->pdgId()) == 24) is_Ele_signal = true;
      if(fabs(gen->pdgId() ) == 13 && fabs(gen->mother()->pdgId()) == 24) is_Mu_signal = true;
  }

  mytree->Fill();

  }
}

void LeptonMultiplicity::create_trees()
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

  mytree->Branch("mu_pT",&mu_pT);
  mytree->Branch("mu_eta",&mu_eta);
  mytree->Branch("mu_phi",&mu_phi);
  mytree->Branch("ele_pT",&el_pT);
  mytree->Branch("ele_eta",&el_eta);
  mytree->Branch("ele_phi",&el_phi);
  mytree->Branch("is_muon",&is_muon);
  mytree->Branch("mu_iso",&mu_iso);
  mytree->Branch("ele_iso",&el_iso);
  mytree->Branch("nMuons",&nMuons);
  mytree->Branch("nElectrons",&nElectrons);
  mytree->Branch("nElectronsLoose",&nElectronsLoose);

  //Save MC info
  if(!runningOnData_){
    mytree->Branch("PU_Weight",&PU_Weight);

    mytree->Branch("isMuonSignal",&is_Mu_signal);
    mytree->Branch("isEleSignal",&is_Ele_signal);
    //mytree->Branch("isttbarlnu",&is_ttbar_lnu);
  }

}

// void LeptonMultiplicity::beginJob()
// {
// }


//define this as a plug-in
DEFINE_FWK_MODULE(LeptonMultiplicity);
