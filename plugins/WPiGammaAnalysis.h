
//---------- class declaration----------

class WPiGammaAnalysis : public edm::EDAnalyzer {
public:
  explicit WPiGammaAnalysis(const edm::ParameterSet&);
  ~WPiGammaAnalysis();

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() override;

  const edm::InputTag packedPFCandidates_;
  const edm::InputTag slimmedMuons_; 
  const edm::InputTag prunedGenParticles_;
  const edm::InputTag slimmedPhotons_;
  const edm::InputTag slimmedElectrons_;
  const edm::InputTag slimmedJets_;
  const edm::InputTag slimmedMETs_;
  const edm::InputTag slimmedMETsPuppi_;
  bool runningOnData_;
  const edm::InputTag pvCollection_;  
  const edm::InputTag bsCollection_;  
  const edm::InputTag PileupSrc_;


  edm::LumiReWeighting Lumiweights_;

  edm::Service<TFileService> fs;

  void create_trees();

  // ----------member data ---------------------------
  TH1F* h_Events;

  TH1F* inv_mass_1;
  TH1F* inv_mass_2;

  //Counters
  int nPV;

  int nMuons;
  int nElectrons;
  int nElectronsLoose;
  int nPions;
  int nPhotons;
  int nBjets;
  int nBjets_25;

  int _Nevents_processed;
  int _Nevents_muVeto;
  int _Nevents_eleVeto;
  int _Nevents_isLepton;
  int _Nevents_TwoLepton;
  int _Nevents_isPion;
  int _Nevents_isPhotons;
  int _Nevents_isWmass;

  //TTree and TTree variables
  TTree *mytree;

  int run_number;

  bool are_lep_pi_opposite_charge;
  float lepton_pT_tree;
  float lepton_eta_tree;
  float lepton_etaSC_tree;
  float lepton_phi_tree;
  float lepton_dxy_tree;
  float lepton_dz_tree;
  float lepton_iso_tree;
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
  float metpuppi_pT;
  
  bool is_muon;
  bool is_ele;
  bool isSingleMuTrigger_24;
  bool isSingleMuTrigger_50;
  bool isSingleEleTrigger;

  //MC truth
  float PU_Weight;

  bool is_signal_Wplus ;
  bool is_signal_Wminus;
  bool is_Wplus_from_t;
  bool is_Wminus_from_tbar;
  bool is_Wplus_in_lep;
  bool is_Wminus_in_lep;
  bool is_Mu_signal;
  bool is_Ele_signal;

  bool is_pi_a_pi;
  bool is_pi_matched;

  bool is_photon_a_photon;
  bool is_photon_matched;
  bool is_gen_ph;
  float gen_ph_pT_tree;
  int gen_ph_mother_tree;

  bool is_ttbar_lnu;

  //rho for isolation

  float rho_;

  //Tokens
  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > packedPFCandidatestoken_; 
  edm::EDGetTokenT<std::vector<pat::Muon> > slimmedMuonstoken_; 
  edm::EDGetTokenT<std::vector<reco::GenParticle> > prunedGenParticlestoken_; 
  edm::EDGetTokenT<std::vector<pat::Jet> > slimmedJetstoken_;
  edm::EDGetTokenT<std::vector<pat::MET> > slimmedMETstoken_;
  edm::EDGetTokenT<std::vector<pat::MET> > slimmedMETsPuppitoken_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > tok_Vertex_; 
  edm::EDGetTokenT<reco::BeamSpot> tok_beamspot_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;

  //Ele ID decisions objects
  edm::EDGetToken electronsMiniAODToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
  // MVA values and categories (optional)
  edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_el_loose_;
  edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_el_loose_;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_el_;
  edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_el_;

  //Photon ID decisions
  edm::EDGetToken photonsMiniAODToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > phoMediumIdBoolMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > phoMediumIdFullInfoMapToken_;

  //rho (PU energy density)
  edm::EDGetTokenT<double> rhoToken_;

  // MVA values and categories (optional)
  edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_ph_;
  edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_ph_;

  bool verboseIdFlag_;

  //Effective areas for isolation
  EffectiveAreas   effectiveAreas_;
  double Bjets_WP_;

};
