//---------- class declaration----------

class WPiGammaAnalysis : public edm::EDAnalyzer {
public:
  explicit WPiGammaAnalysis(const edm::ParameterSet&);
  ~WPiGammaAnalysis();

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() override;

  bool runningOnData_;
  int runningEra_;
  const edm::InputTag packedPFCandidates_;
  const edm::InputTag slimmedMuons_; 
  const edm::InputTag prunedGenParticles_;
  const edm::InputTag slimmedPhotons_;
  const edm::InputTag slimmedElectrons_;
  const edm::InputTag slimmedJets_;
  const edm::InputTag slimmedMETs_;
  const edm::InputTag slimmedMETsPuppi_;
  const edm::InputTag pvCollection_;  
  const edm::InputTag bsCollection_;  
  const edm::InputTag PileupSrc_;
  const edm::InputTag GenInfo_;


  edm::LumiReWeighting Lumiweights_;

  edm::Service<TFileService> fs;

  void create_trees();

  // ----------member data ---------------------------
  TH1F* h_Events;

  TH1F* inv_mass_1;
  TH1F* inv_mass_2;

  TH1F* h_pileup;

  TH2D* h2_BTaggingEff_Num_b;
  TH2D* h2_BTaggingEff_Denom_b;

  //Counters
  int nPV;

  int nMuons;
  int nElectrons;
  int nElectronsLoose;
  int nPions;
  int nPhotons;
  int njets;
  int nBjets_30;
  int nBjets_25;
  int nBjets_scaled;
  int nBjets_medium_30;
  int nBjets_medium_25;
  int nBjets_medium_scaled;

  int _Nevents_processed;
  int _Nevents_muVeto;
  int _Nevents_eleVeto;
  int _Nevents_isLepton;
  int _Nevents_TwoLepton;
  int _Nevents_isPion;
  int _Nevents_isPhotons;
  int _Nevents_isWmass;

  bool is_one_bJet_found;
  bool bJet_outside_eta_bounds;
  std::shared_ptr<TH2> h_bTagEfficiency_;
  BTagCalibrationReader reader;
  double jetSF;
  float jet_pT_temp;
  float prod_1_minus_eff;
  float prod_1_minus_eff_SF;
  float prod_eff;
  float prod_eff_SF;
  float bTag_Weight;

  //TTree and TTree variables
  TTree *mytree;

  int run_number;

  bool are_lep_pi_opposite_charge;
  float lepton_pT_tree;
  float lepton_eta_tree;
  float lepton_etaSC_tree;
  float lepton_phiSC_tree;
  float lepton_phi_tree;
  float lepton_dxy_tree;
  float lepton_dz_tree;
  float lepton_iso_tree;
  bool isTriggerMatched_tree;
  float Wplus_pT;
  float Wminus_pT;

  float mu_pT;
  float mu_eta;
  float mu_phi;
  int   mu_ID ;
  float mu_dxy;
  float mu_dz;
  float mu_iso;
  float best_mu_iso;

  float el_pT;
  float el_eta;
  float el_etaSC;
  float el_phiSC;
  float el_phi;
  int   el_ID;
  float el_dxy;
  float el_dz;
  float el_iso;
  float best_el_iso;

  float pi_pT;
  float pi_eta;
  float pi_phi;
  float pi_energy;
  float pi_dxy;
  float pi_dz;
  float sum_pT_03;
  float sum_pT_05;
  float sum_pT_05_ch;
  int nTracks_in_piIso;

  float ph_eT;
  float ph_eta;
  float ph_etaSC;
  float ph_phi;
  float ph_energy;
  float ph_iso_ChargedHadron;
  float ph_iso_NeutralHadron;
  float ph_iso_Photon;
  //float ph_iso_Track;
  float ph_iso_eArho;

  float pTmuMax;
  float pTeleMax;
  float pTpiMax;
  float eTphMax;

  float deltaphi_lep_pi;

  float _Wmass;

  float met_pT;
  float met_phi;
  float met_pT_scaled;
  float metpuppi_pT;
  float metpuppi_phi;
  
  bool is_muon;
  bool is_ele;
  bool isSingleMuTrigger_24;
  bool isSingleMuTrigger_27;
  bool isSingleMuTrigger_50;
  bool isSingleEleTrigger_25;
  bool isSingleEleTrigger_27;
  bool isSingleEleTrigger_32_DoubleEG;
  bool isSingleEleTrigger_32;

  //MC truth
  float PU_Weight;
  float MC_Weight;
  float MCT_HpT_mu_pT_Max;
  float MCT_HpT_mu_pT;
  float MCT_HpT_mu_eta;
  float MCT_HpT_mu_phi;
  float MCT_HpT_ele_pT_Max;
  float MCT_HpT_ele_pT;
  float MCT_HpT_ele_eta;
  float MCT_HpT_ele_phi;
  float MCT_HeT_ph_eT_Max;
  float MCT_HeT_ph_eT;
  float MCT_HeT_ph_eta;
  float MCT_HeT_ph_phi;

  bool is_signal_Wplus ;
  bool is_signal_Wminus;
  bool is_Wplus_from_t;
  bool is_Wminus_from_tbar;
  bool is_Wplus_in_lep;
  bool is_Wminus_in_lep;
  bool is_Mu_signal;
  bool is_Ele_signal;
  bool is_Mu_signal_fromTau;
  bool is_Ele_signal_fromTau;

  bool is_pi_a_pi;
  bool is_pi_matched;

  bool is_photon_a_photon;
  bool is_photon_matched;
  bool is_gen_pi;
  bool is_gen_ph;
  float gen_pi_pT_tree;
  float gen_pi_eta_tree;
  float gen_pi_phi_tree;
  float gen_pi_energy_tree;
  int gen_pi_mother_tree;
  int gen_pi_grandmother_tree;
  float gen_pi_grandmother_pT_tree;
  float gen_pi_grandmother_eta_tree;
  float gen_pi_grandmother_phi_tree;
  float gen_pi_grandmother_energy_tree;

  float gen_ph_pT_tree;
  float gen_ph_eta_tree;
  float gen_ph_phi_tree;
  float gen_ph_energy_tree;
  int gen_ph_mother_tree;

  int pi_gen_ID_tree;
  int pi_gen_mother_ID_tree;
  int pi_gen_nDaughters_tree;

  bool is_ttbar_lnu;

  //rho for isolation
  float rho_;

  double Prefiring_Weight;

  //Tokens
  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > packedPFCandidatesToken_; 
  edm::EDGetTokenT<std::vector<pat::Muon> > slimmedMuonsToken_; 
  edm::EDGetTokenT<std::vector<reco::GenParticle> > prunedGenParticlesToken_; 
  edm::EDGetTokenT<std::vector<pat::Jet> > slimmedJetsToken_;
  edm::EDGetTokenT<std::vector<pat::MET> > slimmedMETsToken_;
  edm::EDGetTokenT<std::vector<pat::MET> > slimmedMETsPuppiToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > offlineSlimmedPrimaryVerticesToken_; 
  edm::EDGetTokenT<reco::BeamSpot> offlineBeamSpotToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_;
  edm::EDGetTokenT<GenEventInfoProduct> GenInfoToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > triggerObjectsTokenMC_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > triggerObjectsTokenData2016_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > triggerObjectsTokenData2017_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > triggerObjectsTokenData2018_;

  //Ele ID decisions objects
  edm::EDGetToken electronsMiniAODToken_;

  //Photon ID decisions
  edm::EDGetToken photonsMiniAODToken_;


  //rho (PU energy density)
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<double> PrefiringWeightToken_;


  bool verboseIdFlag_;

  //Effective areas for isolation
  EffectiveAreas   effectiveAreas_el_;
  EffectiveAreas   effectiveAreas_ph_;
  double Bjets_WP_loose;
  double Bjets_WP_medium;

};
