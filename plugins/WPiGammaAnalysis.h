
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
  int nMuons;
  int nElectrons;
  int nPions;
  int nPhotons;
  int nBjets;

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

  float lepton_pT_tree;
  float lepton_eta_tree;
  float lepton_phi_tree;

  float pi_pT;
  float pi_eta;
  float pi_phi;

  float ph_pT;
  float ph_eta;
  float ph_phi;

  float _Wmass;
  
  bool is_muon;
  bool isSingleMuTrigger;
  bool isSingleEleTrigger;

  //MC truth
  float PU_Weight;

  bool is_Mu_signal;
  bool is_Ele_signal;

  bool is_pi_a_pi;
  bool is_pi_matched;

  bool is_photon_a_photon;
  bool is_photon_matched;

  bool is_ttbar_lnu;

  //Tokens
  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > packedPFCandidatestoken_; 
  edm::EDGetTokenT<std::vector<pat::Muon> > slimmedMuonstoken_; 
  edm::EDGetTokenT<std::vector<reco::GenParticle> > prunedGenParticlestoken_; 
  edm::EDGetTokenT<std::vector<pat::Jet> > slimmedJetstoken_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > tok_Vertex_; 
  edm::EDGetTokenT<reco::BeamSpot> tok_beamspot_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;

  //Ele ID decisions objects
  edm::EDGetToken electronsMiniAODToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
  // MVA values and categories (optional)
  edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_el_;
  edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_el_;

  //Photon ID decisions
  edm::EDGetToken photonsMiniAODToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > phoMediumIdBoolMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > phoMediumIdFullInfoMapToken_;
  // MVA values and categories (optional)
  edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_ph_;
  edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_ph_;
  bool verboseIdFlag_;

};
