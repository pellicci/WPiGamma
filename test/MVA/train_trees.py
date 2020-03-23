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

evaluate_BDT_systematic = False  # Train on shifted variables, test on nominal case

#---------------------------------#

if isMuon and not evaluate_BDT_systematic:
    fIn_bkg = ROOT.TFile("Tree_MC_Background_mu.root")
    tree_bkg = fIn_bkg.Get("minitree_background_mu")
    fIn_sig = ROOT.TFile("Tree_MC_Signal_mu.root")
    tree_sig = fIn_sig.Get("minitree_signal_mu")
    fOut = ROOT.TFile("outputs/Nominal_training_mu.root","RECREATE")
    #fOut = ROOT.TFile("outputs/Nominal_training_mu_withBjets.root","RECREATE")
    #fOut = ROOT.TFile("outputs/Nominal_training_mu_Wmass.root","RECREATE")

if not isMuon and not evaluate_BDT_systematic:
    fIn_bkg = ROOT.TFile("Tree_MC_Background_ele.root")
    tree_bkg = fIn_bkg.Get("minitree_background_ele")
    fIn_sig = ROOT.TFile("Tree_MC_Signal_ele.root")
    tree_sig = fIn_sig.Get("minitree_signal_ele")
    fOut = ROOT.TFile("outputs/Nominal_training_ele.root","RECREATE")
    #fOut = ROOT.TFile("outputs/Nominal_training_ele_withBjets.root","RECREATE")
    #fOut = ROOT.TFile("outputs/Nominal_training_ele_Wmass.root","RECREATE")

if isMuon and evaluate_BDT_systematic:
    fIn_bkg_training = ROOT.TFile("Tree_MC_Background_mu_shifted.root")
    tree_bkg_training = fIn_bkg_training.Get("minitree_background_mu")
    fIn_sig_training = ROOT.TFile("Tree_MC_Signal_mu_shifted.root")
    tree_sig_training = fIn_sig_training.Get("minitree_signal_mu")
    fIn_bkg_testing = ROOT.TFile("Tree_MC_Background_mu.root")
    tree_bkg_testing = fIn_bkg_testing.Get("minitree_background_mu")
    fIn_sig_testing = ROOT.TFile("Tree_MC_Signal_mu.root")
    tree_sig_testing = fIn_sig_testing.Get("minitree_signal_mu")
    fOut = ROOT.TFile("outputs/Nominal_training_mu_shifted.root","RECREATE")

if not isMuon and evaluate_BDT_systematic:
    fIn_bkg_training = ROOT.TFile("Tree_MC_Background_ele_shifted.root")
    tree_bkg_training = fIn_bkg_training.Get("minitree_background_ele")
    fIn_sig_training = ROOT.TFile("Tree_MC_Signal_ele_shifted.root")
    tree_sig_training = fIn_sig_training.Get("minitree_signal_ele")
    fIn_bkg_testing = ROOT.TFile("Tree_MC_Background_ele.root")
    tree_bkg_testing = fIn_bkg_testing.Get("minitree_background_ele")
    fIn_sig_testing = ROOT.TFile("Tree_MC_Signal_ele.root")
    tree_sig_testing = fIn_sig_testing.Get("minitree_signal_ele")
    fOut = ROOT.TFile("outputs/Nominal_training_ele_shifted.root","RECREATE")

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

if evaluate_BDT_systematic:
    dataloader.AddSignalTree(tree_sig_training, sig_weight, "Training")
    dataloader.AddSignalTree(tree_sig_testing, sig_weight, "Testing")
    dataloader.AddBackgroundTree(tree_bkg_training, bkg_weight, "Training")
    dataloader.AddBackgroundTree(tree_bkg_testing, bkg_weight, "Testing")

else:
    dataloader.AddSignalTree(tree_sig, sig_weight)
    dataloader.AddBackgroundTree(tree_bkg, bkg_weight)

dataloader.SetWeightExpression("weight")

mycutSig = ROOT.TCut("")
mycutBkg = ROOT.TCut("")

if isMuon:
    if not evaluate_BDT_systematic:
        dataloader.PrepareTrainingAndTestTree(mycutSig, mycutBkg, ":".join(["!V","nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0"]))
    method_btd  = factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "BDT", ":".join(["H","!V","NTrees=800", "MinNodeSize=2.5%","MaxDepth=3","BoostType=AdaBoost","AdaBoostBeta=0.25","nCuts=20","NegWeightTreatment=IgnoreNegWeightsInTraining"]))
else:
    if not evaluate_BDT_systematic:
        dataloader.PrepareTrainingAndTestTree(mycutSig, mycutBkg, ":".join(["!V","nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0"])) 
    method_btd  = factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "BDT", ":".join(["H","!V","NTrees=800", "MinNodeSize=2.5%","MaxDepth=3","BoostType=AdaBoost","AdaBoostBeta=0.25","nCuts=20","NegWeightTreatment=IgnoreNegWeightsInTraining"]))

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

fOut.Close()

weightfile_dir = "default/weights/TMVAClassification_BDT.weights.xml"

if evaluate_BDT_systematic:
    weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_shifted.xml"
    weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_shifted.xml"
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
