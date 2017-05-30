
//---------- class declaration----------

class WPiGammaAnalysis : public edm::EDAnalyzer {
public:
  explicit WPiGammaAnalysis(const edm::ParameterSet&);
  ~WPiGammaAnalysis();

private:
  //virtual void beginJob() override;
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

  edm::Service<TFileService> fs;

  void create_trees();

  // ----------member data ---------------------------
  //TFile *WPiGammaAnalysis_output;
  TH1F* inv_mass_1;
  TH1F* inv_mass_2;

  TH1F* track_iso_hist;
  TH1F* ecal_iso_hist;
  TH1F* hcal_iso_hist;
  TH1F* calo_iso_hist;
  TH1F* iso_sum_hist;

  Int_t nevent;
  Int_t cand_total;
  Int_t cand_passing_selection;
  Int_t events_least_one_pi;
  Int_t single_pi_counter;
  Int_t pi_from_w_correction;
  Int_t mu_selection;
  Int_t mu_selection_event;
  Int_t ph_total;
  Int_t ph_passing_selection;
  Int_t events_least_one_ph;
  Int_t pi_from_w;
  Int_t photon_from_w;
  Int_t pi_and_photon_from_w;
  Int_t mu_per_event;
  Int_t el_per_event;

  //TTree and TTree variables
  TTree *mytree;
  float mu_per_event_tree;
  float el_per_event_tree;

  float lepton_pT_tree;
  float lepton_eta_tree;
  float lepton_phi_tree;

  float pi_pT_tree;
  float pi_eta_tree;
  float pi_phi_tree;

  float photon_eT_tree;
  float photon_eta_tree;
  float photon_phi_tree;
  
  int n_bjets;

  bool is_signal_mu_tree;
  bool is_signal_tree;
  bool is_muon_tree;
  bool electron_over_right_mu_tree;

  //Tokens
  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > packedPFCandidatestoken_; 
  edm::EDGetTokenT<std::vector<pat::Muon> > slimmedMuonstoken_; 
  edm::EDGetTokenT<std::vector<reco::GenParticle> > prunedGenParticlestoken_; 
  edm::EDGetTokenT<std::vector<pat::Photon> > slimmedPhotonstoken_;
  edm::EDGetTokenT<std::vector<pat::Electron> > slimmedElectronstoken_;
  edm::EDGetTokenT<std::vector<pat::Jet> > slimmedJetstoken_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > tok_Vertex_; 
  edm::EDGetTokenT<reco::BeamSpot> tok_beamspot_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_;
};
