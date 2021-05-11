//---------- class declaration----------

class WDirectProduction : public edm::EDAnalyzer {

public:
  explicit WDirectProduction(const edm::ParameterSet&);
  ~WDirectProduction();

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() override;

  const edm::InputTag genParticles_;

  edm::Service<TFileService> fs;

  void create_trees();

  // ---------- member data ----------- //

  TTree *mytree;

  int genW_ID_tree;
  float genW_pT_tree;
  float genW_eta_tree;
  float genW_phi_tree;
  float genW_E_tree;

  int genPh_ID_tree;
  float genPh_pT_tree;
  float genPh_eta_tree;
  float genPh_phi_tree;
  float genPh_E_tree;

  int genPi_ID_tree;
  float genPi_pT_tree;
  float genPi_eta_tree;
  float genPi_phi_tree;
  float genPi_E_tree;

  edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticlesToken_; 
  //edm::EDGetTokenT<GenEventInfoProduct> GenInfoToken_;
};
