import ROOT
import os
import argparse

#---------------------------------#
p = argparse.ArgumentParser(description='Select whether the MVA will be performed on the muon or the electron sample')
p.add_argument('isMuon_option', help='Type <<muon>> or <<electron>>')
#p.add_argument('year_option', help='Type <<2016>> or <<2017>> or <<2018>>')
args = p.parse_args()

# Switch from muon to electron channel
if args.isMuon_option == "muon":
    isMuon = True
if args.isMuon_option == "electron":
    isMuon = False

# year = args.year_option

data_sidebands = False  # Switch from data sidebands to MC for training (background)

#---------------------------------#

if isMuon and not data_sidebands:
    fIn_bkg = ROOT.TFile("Tree_MC_Background_mu.root")
    tree_bkg = fIn_bkg.Get("minitree_background_mu")
    fIn_sig = ROOT.TFile("Tree_MC_Signal_mu.root")
    tree_sig = fIn_sig.Get("minitree_signal_mu")
    fOut = ROOT.TFile("outputs/Nominal_training_mu.root","RECREATE")
    #fOut = ROOT.TFile("outputs/Nominal_training_mu_Wmass.root","RECREATE")

if not isMuon and not data_sidebands:
    fIn_bkg = ROOT.TFile("Tree_MC_Background_ele.root")
    tree_bkg = fIn_bkg.Get("minitree_background_ele")
    fIn_sig = ROOT.TFile("Tree_MC_Signal_ele.root")
    tree_sig = fIn_sig.Get("minitree_signal_ele")
    fOut = ROOT.TFile("outputs/Nominal_training_ele.root","RECREATE")
    #fOut = ROOT.TFile("outputs/Nominal_training_ele_Wmass.root","RECREATE")

if isMuon and data_sidebands:
    fIn_bkg_DATA = ROOT.TFile("Tree_MC_Background_mu_DATA.root")
    tree_bkg_DATA = fIn_bkg_DATA.Get("minitree_background_mu_DATA")
    fIn_bkg = ROOT.TFile("Tree_MC_Background_mu.root")
    tree_bkg = fIn_bkg.Get("minitree_background_mu")
    fIn_sig_training = ROOT.TFile("Tree_MC_Signal_mu_training.root")
    tree_sig_training = fIn_sig_training.Get("minitree_signal_mu_training")
    fIn_sig_test = ROOT.TFile("Tree_MC_Signal_mu_test.root")
    tree_sig_test = fIn_sig_test.Get("minitree_signal_mu_test")
    fOut = ROOT.TFile("outputs/" + year + "/Nominal_training_mu_sidebands.root","RECREATE")

if not isMuon and data_sidebands:
    fIn_bkg_DATA = ROOT.TFile("Tree_MC_Background_ele_DATA.root")
    tree_bkg_DATA = fIn_bkg_DATA.Get("minitree_background_ele_DATA")
    fIn_bkg = ROOT.TFile("Tree_MC_Background_ele.root")
    tree_bkg = fIn_bkg.Get("minitree_background_ele")
    fIn_sig_training = ROOT.TFile("Tree_MC_Signal_ele_training.root")
    tree_sig_training = fIn_sig_training.Get("minitree_signal_ele_training")
    fIn_sig_test = ROOT.TFile("Tree_MC_Signal_ele_test.root")
    tree_sig_test = fIn_sig_test.Get("minitree_signal_ele_test")
    fOut = ROOT.TFile("outputs/" + year + "/Nominal_training_ele_sidebands.root","RECREATE")

ROOT.TMVA.Tools.Instance()

factory = ROOT.TMVA.Factory("TMVAClassification", fOut,":".join(["!V","Transformations=I;D;P;G,D","AnalysisType=Classification"]))

dataloader = ROOT.TMVA.DataLoader()

dataloader.AddVariable("pi_pT","F") # Both Float and Double variable types must be indicated as F
dataloader.AddVariable("gamma_eT","F")
#dataloader.AddVariable("pi_pT/Wmass","F")
#dataloader.AddVariable("gamma_eT/Wmass","F")
dataloader.AddVariable("nBjets_25","I")
dataloader.AddVariable("lep_pT","F")
dataloader.AddVariable("piRelIso_05_ch","F")
dataloader.AddVariable("MET","F")
#dataloader.AddVariable("deltaphi_lep_pi","F")
#dataloader.AddVariable("deltaphi_lep_gamma","F")

sig_weight = 1.
bkg_weight = 1.

if data_sidebands:
    dataloader.AddSignalTree(tree_sig_training, sig_weight, "Training")
    dataloader.AddSignalTree(tree_sig_test, sig_weight, "Testing")
    dataloader.AddBackgroundTree(tree_bkg_DATA, bkg_weight, "Training")
    dataloader.AddBackgroundTree(tree_bkg, bkg_weight, "Testing")

else:
    dataloader.AddSignalTree(tree_sig, sig_weight)
    dataloader.AddBackgroundTree(tree_bkg, bkg_weight)

dataloader.SetWeightExpression("weight")

mycutSig = ROOT.TCut("")
mycutBkg = ROOT.TCut("")

if isMuon:
    dataloader.PrepareTrainingAndTestTree(mycutSig, mycutBkg, ":".join(["!V","nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0"]))
    method_btd  = factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "BDT", ":".join(["H","!V","NTrees=800", "MinNodeSize=2.5%","MaxDepth=3","BoostType=AdaBoost","AdaBoostBeta=0.25","nCuts=20","NegWeightTreatment=IgnoreNegWeightsInTraining"]))
else:
    dataloader.PrepareTrainingAndTestTree(mycutSig, mycutBkg, ":".join(["!V","nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0"])) 
    method_btd  = factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "BDT", ":".join(["H","!V","NTrees=800", "MinNodeSize=2.5%","MaxDepth=3","BoostType=AdaBoost","AdaBoostBeta=0.25","nCuts=20","NegWeightTreatment=IgnoreNegWeightsInTraining"]))

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

fOut.Close()

weightfile_dir = "default/weights/TMVAClassification_BDT.weights.xml"

if data_sidebands:
    weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_DATA.xml"
    weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_DATA.xml"
else:
    weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu.xml"
    weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele.xml"
    #weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_Wmass.xml"
    #weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_Wmass.xml"

if isMuon:
    rename_weightfile = "mv " + weightfile_dir + " " + weightfile_mu
    os.system(rename_weightfile)
else:
    rename_weightfile = "mv " + weightfile_dir + " " + weightfile_ele
    os.system(rename_weightfile)
