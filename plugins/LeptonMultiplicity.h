
//---------- class declaration----------

class LeptonMultiplicity : public edm::EDAnalyzer {
public:
  explicit LeptonMultiplicity(const edm::ParameterSet&);
  ~LeptonMultiplicity();

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  //virtual void endJob() override;

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
  int nPV;

  int nMuons;
  int nElectrons;
  int nElectronsLoose;
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

  bool is_muon;
  bool is_ele;
  bool isSingleMuTrigger_24;
  bool isSingleMuTrigger_50;
  bool isSingleEleTrigger;
  float mu_pT = 0.;
  float mu_eta = 0.;
  float mu_phi = 0.;
  float el_pT = 0.;
  float el_eta = 0.;
  float el_phi = 0.;
  float mu_iso = 0.;
  float el_iso = 0.;

  //MC truth
  float PU_Weight;

  bool is_Mu_signal;
  bool is_Ele_signal;

  //bool is_ttbar_lnu;

  //rho for isolation

  float rho_;

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
  // edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
  // edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
  // edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
  // MVA values and categories (optional)
  // edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_el_loose_;
  // edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_el_loose_;
  // edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_el_;
  // edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_el_;

  //Photon ID decisions
  edm::EDGetToken photonsMiniAODToken_;
  // edm::EDGetTokenT<edm::ValueMap<bool> > phoMediumIdBoolMapToken_;
  // edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > phoMediumIdFullInfoMapToken_;

  //rho (PU energy density)
  edm::EDGetTokenT<double> rhoToken_;

  // MVA values and categories (optional)
  // edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_ph_;
  // edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_ph_;
  bool verboseIdFlag_;

  //Effective areas for isolation
  EffectiveAreas   effectiveAreas_el_;
};
