import ROOT
import os

#---------------------------------#

isMuon = True  # Switch from muon to electron channel
data_sidebands = False  # Switch from data sidebands to MC for training (background)

#---------------------------------#

if isMuon and not data_sidebands:
    fIn_bkg = ROOT.TFile("Tree_MC_Background_mu.root")
    tree_bkg = fIn_bkg.Get("minitree_background_mu")
    fIn_sig = ROOT.TFile("Tree_MC_Signal_mu.root")
    tree_sig = fIn_sig.Get("minitree_signal_mu")
    # fOut = ROOT.TFile("outputs/Nominal_training_mu_Wmass.root","RECREATE")
    fOut = ROOT.TFile("outputs/Nominal_training_mu.root","RECREATE")
    # fOut = ROOT.TFile("outputs/Nominal_training_mu_puppi.root","RECREATE")
if not isMuon and not data_sidebands:
    fIn_bkg = ROOT.TFile("Tree_MC_Background_ele.root")
    tree_bkg = fIn_bkg.Get("minitree_background_ele")
    fIn_sig = ROOT.TFile("Tree_MC_Signal_ele.root")
    tree_sig = fIn_sig.Get("minitree_signal_ele")
    # fOut = ROOT.TFile("outputs/Nominal_training_ele_Wmass.root","RECREATE")
    fOut = ROOT.TFile("outputs/Nominal_training_ele.root","RECREATE")
    # fOut = ROOT.TFile("outputs/Nominal_training_ele_puppi.root","RECREATE")

if isMuon and data_sidebands:
    fIn_bkg_DATA = ROOT.TFile("Tree_MC_Background_mu_DATA.root")
    tree_bkg_DATA = fIn_bkg_DATA.Get("minitree_background_mu_DATA")
    fIn_bkg = ROOT.TFile("Tree_MC_Background_mu.root")
    tree_bkg = fIn_bkg.Get("minitree_background_mu")
    fIn_sig_training = ROOT.TFile("Tree_MC_Signal_mu_training.root")
    tree_sig_training = fIn_sig_training.Get("minitree_signal_mu_training")
    fIn_sig_test = ROOT.TFile("Tree_MC_Signal_mu_test.root")
    tree_sig_test = fIn_sig_test.Get("minitree_signal_mu_test")
    fOut = ROOT.TFile("outputs/Nominal_training_mu_sidebands.root","RECREATE")

if not isMuon and data_sidebands:
    fIn_bkg_DATA = ROOT.TFile("Tree_MC_Background_ele_DATA.root")
    tree_bkg_DATA = fIn_bkg_DATA.Get("minitree_background_ele_DATA")
    fIn_bkg = ROOT.TFile("Tree_MC_Background_ele.root")
    tree_bkg = fIn_bkg.Get("minitree_background_ele")
    fIn_sig_training = ROOT.TFile("Tree_MC_Signal_ele_training.root")
    tree_sig_training = fIn_sig_training.Get("minitree_signal_ele_training")
    fIn_sig_test = ROOT.TFile("Tree_MC_Signal_ele_test.root")
    tree_sig_test = fIn_sig_test.Get("minitree_signal_ele_test")
    fOut = ROOT.TFile("outputs/Nominal_training_ele_sidebands.root","RECREATE")

ROOT.TMVA.Tools.Instance()

factory = ROOT.TMVA.Factory("TMVAClassification", fOut,":".join(["!V","Transformations=I;D;P;G,D","AnalysisType=Classification"]))


factory.AddVariable("pi_pT","F")
factory.AddVariable("gamma_eT","F")
# factory.AddVariable("pi_pT/Wmass","F")
# factory.AddVariable("gamma_eT/Wmass","F")
factory.AddVariable("nBjets_25","I")
# factory.AddVariable("deltaphi_lep_pi","F")
factory.AddVariable("lep_pT","F")
factory.AddVariable("piRelIso_05_ch","F")
#factory.AddVariable("pi_dxy","F")
# factory.AddVariable("ele_gamma_InvMass")
factory.AddVariable("MET","F")
#factory.AddVariable("MET_puppi","F")


if data_sidebands:
    factory.AddSignalTree(tree_sig_training, ROOT.TMVA.Types.kTraining)
    factory.AddSignalTree(tree_sig_test, ROOT.TMVA.Types.kTesting)
    factory.AddBackgroundTree(tree_bkg_DATA, ROOT.TMVA.Types.kTraining)
    factory.AddBackgroundTree(tree_bkg, ROOT.TMVA.Types.kTesting)

else:
    factory.AddSignalTree(tree_sig)
    factory.AddBackgroundTree(tree_bkg)

factory.SetWeightExpression("weight")

mycuts = ROOT.TCut("weight > 0.")
mycutb = ROOT.TCut("weight > 0.")

if isMuon:
    factory.PrepareTrainingAndTestTree(mycuts, mycutb, ":".join(["!V","nTrain_Signal=7905:nTrain_Background=88204:nTest_Signal=0:nTest_Background=0"]) )
else:
    factory.PrepareTrainingAndTestTree(mycuts, mycutb, ":".join(["!V","nTrain_Signal=6350:nTrain_Background=77225:nTest_Signal=0:nTest_Background=0"]) )

if isMuon:
    method_btd  = factory.BookMethod(ROOT.TMVA.Types.kBDT, "BDT", ":".join(["H","!V","NTrees=800", "MinNodeSize=2.5%","MaxDepth=3","BoostType=AdaBoost","AdaBoostBeta=0.25","nCuts=20"]))
else:
    method_btd  = factory.BookMethod(ROOT.TMVA.Types.kBDT, "BDT", ":".join(["H","!V","NTrees=1100", "MinNodeSize=2.5%","MaxDepth=2","BoostType=AdaBoost","AdaBoostBeta=0.25","nCuts=40"]))


factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

fOut.Close()

weightfile_dir = "weights/TMVAClassification_BDT.weights.xml"

if data_sidebands:
    weightfile_mu  = "weights/TMVAClassification_BDT.weights_mu_DATA.xml"
    weightfile_ele = "weights/TMVAClassification_BDT.weights_ele_DATA.xml"
else:
    weightfile_mu  = "weights/TMVAClassification_BDT.weights_mu.xml"
    weightfile_ele = "weights/TMVAClassification_BDT.weights_ele.xml"
    # weightfile_mu  = "weights/TMVAClassification_BDT.weights_mu_puppi.xml"
    # weightfile_ele = "weights/TMVAClassification_BDT.weights_ele_puppi.xml"
    # weightfile_mu  = "weights/TMVAClassification_BDT.weights_mu_Wmass.xml"
    # weightfile_ele = "weights/TMVAClassification_BDT.weights_ele_Wmass.xml"

if isMuon:
    rename_weightfile = "mv " + weightfile_dir + " " + weightfile_mu
    os.system(rename_weightfile)
else:
    rename_weightfile = "mv " + weightfile_dir + " " + weightfile_ele
    os.system(rename_weightfile)
