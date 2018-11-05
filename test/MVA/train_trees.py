import ROOT
import os
import argparse

#---------------------------------#
p = argparse.ArgumentParser(description='Select whether the MVA will be performed on the muon or the electron sample')
p.add_argument('isMuon_option', help='Type <<muon>> or <<electron>>')
args = p.parse_args()

# Switch from muon to electron channel
if args.isMuon_option == "muon":
    isMuon = True
if args.isMuon_option == "electron":
    isMuon = False

data_sidebands = False  # Switch from data sidebands to MC for training (background)

#---------------------------------#

if isMuon and not data_sidebands:
    fIn_bkg = ROOT.TFile("Tree_MC_Background_mu.root")
    tree_bkg = fIn_bkg.Get("minitree_background_mu")
    fIn_sig = ROOT.TFile("Tree_MC_Signal_mu.root")
    tree_sig = fIn_sig.Get("minitree_signal_mu")
    fOut = ROOT.TFile("outputs/Nominal_training_mu.root","RECREATE")
    # fOut = ROOT.TFile("outputs/Nominal_training_mu_Wmass.root","RECREATE")

if not isMuon and not data_sidebands:
    fIn_bkg = ROOT.TFile("Tree_MC_Background_ele.root")
    tree_bkg = fIn_bkg.Get("minitree_background_ele")
    fIn_sig = ROOT.TFile("Tree_MC_Signal_ele.root")
    tree_sig = fIn_sig.Get("minitree_signal_ele")
    fOut = ROOT.TFile("outputs/Nominal_training_ele.root","RECREATE")
    # fOut = ROOT.TFile("outputs/Nominal_training_ele_Wmass.root","RECREATE")

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


factory.AddVariable("pi_pT","D")
factory.AddVariable("gamma_eT","D")
# factory.AddVariable("pi_pT/Wmass","D")
# factory.AddVariable("gamma_eT/Wmass","D")
factory.AddVariable("nBjets_25","I")
factory.AddVariable("lep_pT","D")
factory.AddVariable("piRelIso_05_ch","D")
factory.AddVariable("MET","D")


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
    factory.PrepareTrainingAndTestTree(mycuts, mycutb, ":".join(["!V","nTrain_Signal=9758:nTrain_Background=230843:nTest_Signal=0:nTest_Background=0"]) ) # To be set with 75/25 ratio of training and testing events
else:
    factory.PrepareTrainingAndTestTree(mycuts, mycutb, ":".join(["!V","nTrain_Signal=7145:nTrain_Background=174851:nTest_Signal=0:nTest_Background=0"]) ) # To be set with 75/25 ratio of training and testing events

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
    # weightfile_mu  = "weights/TMVAClassification_BDT.weights_mu_Wmass.xml"
    # weightfile_ele = "weights/TMVAClassification_BDT.weights_ele_Wmass.xml"

if isMuon:
    rename_weightfile = "mv " + weightfile_dir + " " + weightfile_mu
    os.system(rename_weightfile)
else:
    rename_weightfile = "mv " + weightfile_dir + " " + weightfile_ele
    os.system(rename_weightfile)
