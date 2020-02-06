
//---------- class declaration----------

class LeptonMultiplicity : public edm::EDAnalyzer {
public:
  explicit LeptonMultiplicity(const edm::ParameterSet&);
  ~LeptonMultiplicity();

private:
  //virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  //virtual void endJob() override;

  bool runningOnData_;
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
  bool isSingleMuTrigger_27;
  bool isSingleMuTrigger_50;
  bool isSingleEleTrigger_25;
  bool isSingleEleTrigger_27;
  bool isSingleEleTrigger_32_DoubleEG;
  bool isSingleEleTrigger_32;

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
  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > packedPFCandidatesToken_; 
  edm::EDGetTokenT<std::vector<pat::Muon> > slimmedMuonsToken_; 
  edm::EDGetTokenT<std::vector<reco::GenParticle> > prunedGenParticlesToken_; 
  edm::EDGetTokenT<std::vector<reco::Vertex> > offlineSlimmedPrimaryVerticesToken_; 
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;

  //Ele ID decisions objects
  edm::EDGetToken electronsMiniAODToken_;

  //Photon ID decisions
  edm::EDGetToken photonsMiniAODToken_;


  //rho (PU energy density)
  edm::EDGetTokenT<double> rhoToken_;

  bool verboseIdFlag_;

  //Effective areas for isolation
  EffectiveAreas   effectiveAreas_el_;
  EffectiveAreas   effectiveAreas_ph_;
};
