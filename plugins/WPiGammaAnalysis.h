//---------- class declaration----------

class WPiGammaAnalysis : public edm::EDAnalyzer {
public:
  explicit WPiGammaAnalysis(const edm::ParameterSet&);
  ~WPiGammaAnalysis();

private:
  //virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void multiplicity(const edm::Event&, const edm::EventSetup&);
  //virtual void endJob() override;

  const edm::InputTag packedPFCandidates_;
  const edm::InputTag slimmedMuons_; 
  const edm::InputTag prunedGenParticles_;
  const edm::InputTag slimmedPhotons_;
  const edm::InputTag slimmedElectrons_;
  bool runningOnData_;
  const edm::InputTag pvCollection_;  
  const edm::InputTag bsCollection_;  
  const edm::InputTag PileupSrc_;

  edm::Service<TFileService> fs;

  void create_trees();

  // ----------member data ---------------------------
  //TFile *WPiGammaAnalysis_output;
  TCanvas* canvas1;
  TH1F* inv_mass_1;
  TH1F* inv_mass_2;

  TH1F* track_iso_hist;
  TH1F* ecal_iso_hist;
  TH1F* hcal_iso_hist;
  TH1F* calo_iso_hist;
  TH1F* iso_sum_hist;
  TH2F* tag_mu_hist;
  TH2F* tag_el_hist;

  LorentzVector candidate_pi;
  LorentzVector candidate_ph;

  Int_t nevent;
  Int_t cand_total;
  Int_t cand_passing_selection;
  Int_t events_least_one_pi;
  Int_t single_pi_counter;
  Int_t pi_from_w_correction;
  Int_t mu_selection;
  Int_t mu_selection_event;
  Int_t mu_ID;
  Int_t el_ID;
  Int_t gen_ID;
  Int_t ph_total;
  Int_t ph_passing_selection;
  Int_t events_least_one_ph;
  Int_t gen_mother;
  Int_t pi_from_w;
  Int_t photon_from_w;
  Int_t pi_and_photon_from_w;
  Int_t mu_per_event;
  Int_t el_per_event;

  float pTpi;
  float pTpiMax;
  float pTgamma;
  float pTe; 
  float pTmu; 
  float deta;
  float dphi;
  float etapi;
  float phipi;
  float deltaRMax;
  float deltapTMax;
  float pTmuMin;
  float eTphMin;
  float pTpiMin;
  float candidate_eta;
  float candidate_phi;
  float photon_eta;
  float photon_phi;
  float pxpi;
  float pypi;
  float pzpi;
  float pxph;
  float pyph;
  float pzph;
  float pTph;
  float track_iso; 
  float ecal_iso; 
  float hcal_iso;
  float calo_iso; 
  float iso_sum;
  float gen_px;
  float gen_py;
  float gen_pz;
  float gen_pt;
  float deltaR;
  float deltapT;

  bool tag_lepton_found;
  bool cand_pion_found;
  bool cand_photon_found;
  bool is_pi_a_pi;
  bool is_photon_a_photon;
  bool is_last_pi_a_pi;
  bool is_bad_single_pi;
  bool in_electron_selection;

  //TTree and TTree variables
  TTree *mytree;
  float lepton_pT_tree;
  float lepton_eta_tree;
  float lepton_phi_tree;
  bool is_muon;

  //Tokens
  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > packedPFCandidatestoken_; 
  edm::EDGetTokenT<std::vector<pat::Muon> > slimmedMuonstoken_; 
  edm::EDGetTokenT<std::vector<reco::GenParticle> > prunedGenParticlestoken_; 
  edm::EDGetTokenT<std::vector<pat::Photon> > slimmedPhotonstoken_;
  edm::EDGetTokenT<std::vector<pat::Electron> > slimmedElectronstoken_;
  edm::EDGetTokenT<reco::VertexCollection> tok_Vertex_; 
  edm::EDGetTokenT<reco::BeamSpot>         tok_beamspot_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupSummaryToken_;
};
