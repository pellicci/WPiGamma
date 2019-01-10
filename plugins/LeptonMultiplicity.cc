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
  triggerBits_(consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("triggerbits"))),
  // eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
  // eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  // eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  // mvaValuesMapToken_el_loose_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap_el_loose"))),
  // mvaCategoriesMapToken_el_loose_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap_el_loose"))),
  // mvaValuesMapToken_el_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap_el"))),
  // mvaCategoriesMapToken_el_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap_el"))),
  // phoMediumIdBoolMapToken_(consumes<edm::ValueMap<bool> > (iConfig.getParameter<edm::InputTag>("phoMediumIdBoolMap"))),
  // phoMediumIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> > (iConfig.getParameter<edm::InputTag>("phoMediumIdFullInfoMap"))),
  // mvaValuesMapToken_ph_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap_ph"))),
  // mvaCategoriesMapToken_ph_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap_ph"))),
  verboseIdFlag_(iConfig.getParameter<bool>("phoIdVerbose")),
  effectiveAreas_el_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile_el")).fullPath() )

{
  packedPFCandidatestoken_  = consumes<std::vector<pat::PackedCandidate> >(packedPFCandidates_); 
  slimmedMuonstoken_        = consumes<std::vector<pat::Muon> >(slimmedMuons_);
  prunedGenParticlestoken_  = consumes<std::vector<reco::GenParticle> >(prunedGenParticles_);
  // photonsMiniAODToken_      = mayConsume<edm::View<reco::Photon> > (slimmedPhotons_);
  // electronsMiniAODToken_    = mayConsume<edm::View<reco::GsfElectron> > (slimmedElectrons_);
  slimmedJetstoken_         = consumes<std::vector<pat::Jet> >(slimmedJets_);
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

LeptonMultiplicity::~LeptonMultiplicity()
{
}

// ------------ method called for each event  ------------
void LeptonMultiplicity::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<pat::PackedCandidate>  > PFCandidates;
  iEvent.getByLabel(packedPFCandidates_, PFCandidates);

  edm::Handle<std::vector<pat::Muon>  > slimmedMuons;
  iEvent.getByLabel(slimmedMuons_, slimmedMuons);

  edm::Handle<std::vector<reco::GenParticle>  > genParticles;
  if(!runningOnData_)iEvent.getByLabel(prunedGenParticles_, genParticles);

  //edm::Handle<edm::View<reco::Photon> > slimmedPhotons;
  edm::Handle<std::vector<pat::Photon> > slimmedPhotons;
  iEvent.getByToken(photonsMiniAODToken_,slimmedPhotons);

  // edm::Handle<edm::View<reco::GsfElectron> > slimmedElectrons;
  edm::Handle<std::vector<pat::Electron> > slimmedElectrons;
  iEvent.getByToken(electronsMiniAODToken_,slimmedElectrons);

  edm::Handle<std::vector<pat::Jet > > slimmedJets;
  iEvent.getByLabel(slimmedJets_, slimmedJets);

  edm::Handle<std::vector<reco::Vertex > > slimmedPV;
  iEvent.getByLabel(pvCollection_, slimmedPV);

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);

  //Get the electron ID data from the event stream.
  // Note: this implies that the VID ID modules have been run upstream.
  // edm::Handle<edm::ValueMap<bool> > el_loose_id_decisions; 
  // edm::Handle<edm::ValueMap<bool> > el_medium_id_decisions;
  // edm::Handle<edm::ValueMap<bool> > el_tight_id_decisions; 
  // iEvent.getByToken(eleLooseIdMapToken_,el_loose_id_decisions);
  // iEvent.getByToken(eleMediumIdMapToken_,el_medium_id_decisions);
  // iEvent.getByToken(eleTightIdMapToken_,el_tight_id_decisions);

  // Get MVA values and categories (optional)
  // edm::Handle<edm::ValueMap<float> > el_loose_mvaValues;
  // edm::Handle<edm::ValueMap<int> > el_loose_mvaCategories;
  // edm::Handle<edm::ValueMap<float> > el_mvaValues;
  // edm::Handle<edm::ValueMap<int> > el_mvaCategories;
  // iEvent.getByToken(mvaValuesMapToken_el_loose_,el_loose_mvaValues);
  // iEvent.getByToken(mvaCategoriesMapToken_el_loose_,el_loose_mvaCategories);
  // iEvent.getByToken(mvaValuesMapToken_el_,el_mvaValues);
  // iEvent.getByToken(mvaCategoriesMapToken_el_,el_mvaCategories);

  // Get the photon ID data from the event stream.
  // Note: this implies that the VID ID modules have been run upstream.
  // The first map simply has pass/fail for each particle
  // edm::Handle<edm::ValueMap<bool> > ph_medium_id_decisions;
  // iEvent.getByToken(phoMediumIdBoolMapToken_,ph_medium_id_decisions);
  
  // The second map has the full info about the cut flow
  // edm::Handle<edm::ValueMap<vid::CutFlowResult> > ph_medium_id_cutflow_data;
  // iEvent.getByToken(phoMediumIdFullInfoMapToken_,ph_medium_id_cutflow_data);

  // Get MVA values and categories (optional)
  // edm::Handle<edm::ValueMap<float> > ph_mvaValues;
  // edm::Handle<edm::ValueMap<int> > ph_mvaCategories;
  // iEvent.getByToken(mvaValuesMapToken_ph_,ph_mvaValues);
  // iEvent.getByToken(mvaCategoriesMapToken_ph_,ph_mvaCategories);

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
  isSingleEleTrigger   = false;

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
    //if(mu->pt() < 25. || mu->pt() < pTmuMax || !mu->isTightMuon(slimmedPV->at(0))) continue;
    //if(mu->pt() < 25. || mu->pt() < pTmuMax || !mu->isMediumMuon()) continue;
    if(mu->pt() < 20. || !mu->isMediumMuon() || abs(mu->eta()) > 2.4 || fabs(mu->muonBestTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(mu->muonBestTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;
    mu_iso = (mu->chargedHadronIso() + std::max(0., mu->neutralHadronIso() + mu->photonIso() - 0.5*mu->puChargedHadronIso()))/mu->pt();
    if(mu_iso > 0.25) continue;
    //pTmuMax = mu->pt();


    //mu_ID    = mu->pdgId();
    is_muon  = true;
    mu_eta   = mu->eta();
    mu_phi   = mu->phi();
    mu_pT    = mu->pt();
    nMuons++;
    //std::cout << "mu pT :" << mu->pt() << "Eta: " << mu->eta() << "phi:" << mu->phi() << "iso: " << mu_iso << "number: " << nMuons << std::endl;
  }


  // Get rho value
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;

  //Loop over electrons
  for(auto el = slimmedElectrons->begin(); el != slimmedElectrons->end(); ++el){
    //for (size_t i = 0; i < slimmedElectrons->size(); ++i){
    //const auto el = slimmedElectrons->ptrAt(i);

    if(el->pt() < 20. || fabs(el->gsfTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(el->gsfTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;

    // bool isPassMedium = (*el_medium_id_decisions)[el];
    // if(!isPassMedium) continue;
    if(el->electronID("mvaEleID-Fall17-iso-V2-wp90") == 0) continue;

    //PflowIsolationVariables pfIso = el->pfIsolationVariables();
    float abseta =  abs(el->superCluster()->eta());
    float eA = effectiveAreas_el_.getEffectiveArea(abseta);
    el_iso = (el->pfIsolationVariables().sumChargedHadronPt + std::max( 0.0f, el->pfIsolationVariables().sumNeutralHadronEt + el->pfIsolationVariables().sumPhotonEt - eA*rho_))/el->pt();
    if(el_iso > 0.35) continue;

    //el_ID       = el->pdgId();
    is_ele      = true;
    el_eta      = el->eta();
    el_phi      = el->phi();
    el_pT       = el->pt();

    //pTeleMax = el_pT;
    nElectrons++;
    //std::cout << "el pT :" << el_pT << "Eta: " << el_eta << "phi:" << el_phi << "iso: " << el_iso << "number: " << nElectrons << std::endl;
  }

  for(auto el = slimmedElectrons->begin(); el != slimmedElectrons->end(); ++el){
    //for (size_t i = 0; i < slimmedElectrons->size(); ++i){
    //const auto el = slimmedElectrons->ptrAt(i);
    // if(el->pt() < 26. || el->pt() < pTeleMax) continue;
    if(el->pt() < 20. || fabs(el->gsfTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(el->gsfTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;
    //float deltaR = sqrt((el_eta-el->eta())*(el_eta-el->eta())+(el_phi-el->phi())*(el_phi-el->phi()));
    //float deltapT = el_pT - el->pt();
    //if (deltaR < 0.01 && deltapT < 0.01 ) continue;

    // bool isPassLoose = (*el_loose_id_decisions)[el];
    // if(!isPassLoose) continue;
    if(el->electronID("mvaEleID-Fall17-iso-V2-wpLoose") == 0) continue;

    float abseta =  abs(el->superCluster()->eta());
    float eA = effectiveAreas_el_.getEffectiveArea(abseta);
    if((el->pfIsolationVariables().sumChargedHadronPt + std::max( 0.0f, el->pfIsolationVariables().sumNeutralHadronEt + el->pfIsolationVariables().sumPhotonEt - eA*rho_))/el->pt() > 0.35) continue;

    nElectronsLoose++;
    //std::cout << "el pT :" << el_pT << "Eta: " << el_eta << "phi:" << el_phi << "iso: " << el_iso << "number: " << nElectrons << std::endl;
  }


  //Do NOT continue if you didn't find either a muon or an electron
  if(!is_muon && !is_ele) return;
  //_Nevents_isLepton++;

  //if(is_muon && is_ele) return;
  //_Nevents_TwoLepton++;

  //In signal, identify if there's a real mu or ele from W
  /*bool is_Wplus_from_t = false;
  bool is_Wminus_from_tbar = false;
  bool is_Wplus_in_lep = false;
  bool is_Wminus_in_lep = false;
  is_ttbar_lnu = false;*/
  is_Ele_signal = false;
  is_Mu_signal = false;
  if(!runningOnData_){
    for (auto gen = genParticles->begin(); gen != genParticles->end(); ++gen){
      if(fabs(gen->pdgId() ) == 11 && fabs(gen->mother()->pdgId()) == 24) is_Ele_signal = true;
      if(fabs(gen->pdgId() ) == 13 && fabs(gen->mother()->pdgId()) == 24) is_Mu_signal = true;
      //---------- Define events of ttbar->WbWb->lnulnu -----------
      /*if(gen->pdgId() == 24 && gen->mother()->pdgId() == 6) is_Wplus_from_t = true;
      if(gen->pdgId() == -24 && gen->mother()->pdgId() == -6) is_Wminus_from_tbar = true;
      if(fabs(gen->pdgId() != 24) || gen->numberOfDaughters() != 2) continue;
      for (int i=0; i<2; i++){
	if(gen->daughter(i)->pdgId() == 11 || gen->daughter(i)->pdgId() == 13 || gen->daughter(i)->pdgId() == 15) is_Wminus_in_lep = true;
	if(gen->daughter(i)->pdgId() == -11 || gen->daughter(i)->pdgId() == -13 || gen->daughter(i)->pdgId() == -15) is_Wplus_in_lep = true;
      }
    }
    if(is_Wplus_from_t && is_Wminus_from_tbar && is_Wminus_in_lep && is_Wplus_in_lep) is_ttbar_lnu = true;*/
  }

  mytree->Fill();

  //std::cout << "n fotoni: " << events_least_one_ph << std::endl;
  }
}

void LeptonMultiplicity::create_trees()
{
  mytree = fs->make<TTree>("mytree", "Tree containing gen&reco");

  mytree->Branch("nPV",&nPV);
  mytree->Branch("isSingleMuTrigger_24",&isSingleMuTrigger_24);
  mytree->Branch("isSingleMuTrigger_50",&isSingleMuTrigger_50);
  mytree->Branch("isSingleEleTrigger",&isSingleEleTrigger);

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

void LeptonMultiplicity::beginJob()
{
  //Flag for PileUp reweighting
  if (!runningOnData_){
   Lumiweights_ = edm::LumiReWeighting("MCpileUp_2016_25ns_Moriond17MC_PoissonOOTPU.root", "MyDataPileupHistogram.root", "pileup", "pileup");
  }
}


/*void LeptonMultiplicity::endJob() 
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

  }*/

//define this as a plug-in
DEFINE_FWK_MODULE(LeptonMultiplicity);
