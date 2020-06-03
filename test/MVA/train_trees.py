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

evaluate_BDT_systematic     = False  # Train on shifted variables, test on nominal case
test_on_signal_shifted_up   = False
test_on_signal_shifted_down = False
test_on_signal_sin2         = False
test_on_signal_cos          = False

#---------------------------------#

if isMuon and not evaluate_BDT_systematic and not test_on_signal_shifted_up and not test_on_signal_shifted_down and not test_on_signal_sin2 and not test_on_signal_cos:
    fIn_bkg  = ROOT.TFile("Tree_MC_Background_mu.root")
    tree_bkg = fIn_bkg.Get("minitree_background_mu")
    fIn_sig  = ROOT.TFile("Tree_MC_Signal_mu.root")
    tree_sig = fIn_sig.Get("minitree_signal_mu")
    fOut     = ROOT.TFile("outputs/Nominal_training_mu.root","RECREATE")
    #fOut = ROOT.TFile("outputs/Nominal_training_mu_withBjets.root","RECREATE")
    #fOut = ROOT.TFile("outputs/Nominal_training_mu_Wmass.root","RECREATE")

if not isMuon and not evaluate_BDT_systematic and not test_on_signal_shifted_up and not test_on_signal_shifted_down and not test_on_signal_sin2 and not test_on_signal_cos:
    fIn_bkg  = ROOT.TFile("Tree_MC_Background_ele.root")
    tree_bkg = fIn_bkg.Get("minitree_background_ele")
    fIn_sig  = ROOT.TFile("Tree_MC_Signal_ele.root")
    tree_sig = fIn_sig.Get("minitree_signal_ele")
    fOut     = ROOT.TFile("outputs/Nominal_training_ele.root","RECREATE")
    #fOut = ROOT.TFile("outputs/Nominal_training_ele_withBjets.root","RECREATE")
    #fOut = ROOT.TFile("outputs/Nominal_training_ele_Wmass.root","RECREATE")

if isMuon and evaluate_BDT_systematic and not test_on_signal_shifted_up and not test_on_signal_shifted_down and not test_on_signal_sin2 and not test_on_signal_cos:
    fIn_bkg_training = ROOT.TFile("Tree_MC_Background_mu_shifted.root")
    tree_bkg_training = fIn_bkg_training.Get("minitree_background_mu")
    fIn_sig_training = ROOT.TFile("Tree_MC_Signal_mu_shifted.root")
    tree_sig_training = fIn_sig_training.Get("minitree_signal_mu")
    fIn_bkg_testing = ROOT.TFile("Tree_MC_Background_mu.root")
    tree_bkg_testing = fIn_bkg_testing.Get("minitree_background_mu")
    fIn_sig_testing = ROOT.TFile("Tree_MC_Signal_mu.root")
    tree_sig_testing = fIn_sig_testing.Get("minitree_signal_mu")
    fOut = ROOT.TFile("outputs/Nominal_training_mu_shifted.root","RECREATE")

if not isMuon and evaluate_BDT_systematic and not test_on_signal_shifted_up and not test_on_signal_shifted_down and not test_on_signal_sin2 and not test_on_signal_cos:
    fIn_bkg_training = ROOT.TFile("Tree_MC_Background_ele_shifted.root")
    tree_bkg_training = fIn_bkg_training.Get("minitree_background_ele")
    fIn_sig_training = ROOT.TFile("Tree_MC_Signal_ele_shifted.root")
    tree_sig_training = fIn_sig_training.Get("minitree_signal_ele")
    fIn_bkg_testing = ROOT.TFile("Tree_MC_Background_ele.root")
    tree_bkg_testing = fIn_bkg_testing.Get("minitree_background_ele")
    fIn_sig_testing = ROOT.TFile("Tree_MC_Signal_ele.root")
    tree_sig_testing = fIn_sig_testing.Get("minitree_signal_ele")
    fOut = ROOT.TFile("outputs/Nominal_training_ele_shifted.root","RECREATE")

if isMuon and test_on_signal_shifted_up:
    fIn_bkg = ROOT.TFile("Tree_MC_Background_mu.root")
    tree_bkg = fIn_bkg.Get("minitree_background_mu")
    fIn_sig_training = ROOT.TFile("Tree_MC_Signal_mu_Pythia_nominal.root")
    tree_sig_training = fIn_sig_training.Get("minitree_signal_mu")
    fIn_sig_testing = ROOT.TFile("Tree_MC_Signal_mu_Pythia_up.root")
    tree_sig_testing = fIn_sig_testing.Get("minitree_signal_mu")
    fOut = ROOT.TFile("outputs/Nominal_training_mu_Pythia_up.root","RECREATE")

if not isMuon and test_on_signal_shifted_up:
    fIn_bkg = ROOT.TFile("Tree_MC_Background_ele.root")
    tree_bkg = fIn_bkg.Get("minitree_background_ele")
    fIn_sig_training = ROOT.TFile("Tree_MC_Signal_ele_Pythia_nominal.root")
    tree_sig_training = fIn_sig_training.Get("minitree_signal_ele")
    fIn_sig_testing = ROOT.TFile("Tree_MC_Signal_ele_Pythia_up.root")
    tree_sig_testing = fIn_sig_testing.Get("minitree_signal_ele")
    fOut = ROOT.TFile("outputs/Nominal_training_ele_Pythia_up.root","RECREATE")

if isMuon and test_on_signal_shifted_down:
    fIn_bkg = ROOT.TFile("Tree_MC_Background_mu.root")
    tree_bkg = fIn_bkg.Get("minitree_background_mu")
    fIn_sig_training = ROOT.TFile("Tree_MC_Signal_mu_Pythia_nominal.root")
    tree_sig_training = fIn_sig_training.Get("minitree_signal_mu")
    fIn_sig_testing = ROOT.TFile("Tree_MC_Signal_mu_Pythia_down.root")
    tree_sig_testing = fIn_sig_testing.Get("minitree_signal_mu")
    fOut = ROOT.TFile("outputs/Nominal_training_mu_Pythia_down.root","RECREATE")

if not isMuon and test_on_signal_shifted_down:
    fIn_bkg = ROOT.TFile("Tree_MC_Background_ele.root")
    tree_bkg = fIn_bkg.Get("minitree_background_ele")
    fIn_sig_training = ROOT.TFile("Tree_MC_Signal_ele_Pythia_nominal.root")
    tree_sig_training = fIn_sig_training.Get("minitree_signal_ele")
    fIn_sig_testing = ROOT.TFile("Tree_MC_Signal_ele_Pythia_down.root")
    tree_sig_testing = fIn_sig_testing.Get("minitree_signal_ele")
    fOut = ROOT.TFile("outputs/Nominal_training_ele_Pythia_down.root","RECREATE")

if isMuon and test_on_signal_sin2:
    fIn_bkg = ROOT.TFile("Tree_MC_Background_mu.root")
    tree_bkg = fIn_bkg.Get("minitree_background_mu")
    fIn_sig_training = ROOT.TFile("Tree_MC_Signal_mu_no_polarization.root")
    tree_sig_training = fIn_sig_training.Get("minitree_signal_mu")
    fIn_sig_testing = ROOT.TFile("Tree_MC_Signal_mu_sin2.root")
    tree_sig_testing = fIn_sig_testing.Get("minitree_signal_mu")
    fOut = ROOT.TFile("outputs/Nominal_training_mu_sin2.root","RECREATE")

if not isMuon and test_on_signal_sin2:
    fIn_bkg = ROOT.TFile("Tree_MC_Background_ele.root")
    tree_bkg = fIn_bkg.Get("minitree_background_ele")
    fIn_sig_training = ROOT.TFile("Tree_MC_Signal_ele_no_polarization.root")
    tree_sig_training = fIn_sig_training.Get("minitree_signal_ele")
    fIn_sig_testing = ROOT.TFile("Tree_MC_Signal_ele_sin2.root")
    tree_sig_testing = fIn_sig_testing.Get("minitree_signal_ele")
    fOut = ROOT.TFile("outputs/Nominal_training_ele_sin2.root","RECREATE")

if isMuon and test_on_signal_cos:
    fIn_bkg = ROOT.TFile("Tree_MC_Background_mu.root")
    tree_bkg = fIn_bkg.Get("minitree_background_mu")
    fIn_sig_training = ROOT.TFile("Tree_MC_Signal_mu_no_polarization.root")
    tree_sig_training = fIn_sig_training.Get("minitree_signal_mu")
    fIn_sig_testing = ROOT.TFile("Tree_MC_Signal_mu_cos.root")
    tree_sig_testing = fIn_sig_testing.Get("minitree_signal_mu")
    fOut = ROOT.TFile("outputs/Nominal_training_mu_cos.root","RECREATE")

if not isMuon and test_on_signal_cos:
    fIn_bkg = ROOT.TFile("Tree_MC_Background_ele.root")
    tree_bkg = fIn_bkg.Get("minitree_background_ele")
    fIn_sig_training = ROOT.TFile("Tree_MC_Signal_ele_no_polarization.root")
    tree_sig_training = fIn_sig_training.Get("minitree_signal_ele")
    fIn_sig_testing = ROOT.TFile("Tree_MC_Signal_ele_cos.root")
    tree_sig_testing = fIn_sig_testing.Get("minitree_signal_ele")
    fOut = ROOT.TFile("outputs/Nominal_training_ele_cos.root","RECREATE")


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

sig_weight = 1.
bkg_weight = 1.

if evaluate_BDT_systematic:
    dataloader.AddSignalTree(tree_sig_training, sig_weight, "Training")
    dataloader.AddSignalTree(tree_sig_testing, sig_weight, "Testing")
    dataloader.AddBackgroundTree(tree_bkg_training, bkg_weight, "Training")
    dataloader.AddBackgroundTree(tree_bkg_testing, bkg_weight, "Testing")

elif (test_on_signal_shifted_up or test_on_signal_shifted_down):
    dataloader.AddSignalTree(tree_sig_training, sig_weight, "Training")
    dataloader.AddSignalTree(tree_sig_testing, sig_weight, "Testing")
    dataloader.AddBackgroundTree(tree_bkg, bkg_weight)

elif (test_on_signal_sin2 or test_on_signal_cos):
    dataloader.AddSignalTree(tree_sig_training, sig_weight, "Training")
    dataloader.AddSignalTree(tree_sig_testing, sig_weight, "Testing")
    dataloader.AddBackgroundTree(tree_bkg, bkg_weight)

else:
    dataloader.AddSignalTree(tree_sig, sig_weight)
    dataloader.AddBackgroundTree(tree_bkg, bkg_weight)

    # dataloader.AddSignalTree(tree_sig_training, sig_weight, "Training")
    # dataloader.AddSignalTree(tree_sig_testing, sig_weight, "Testing")
    # dataloader.AddBackgroundTree(tree_bkg_training, bkg_weight, "Training")
    # dataloader.AddBackgroundTree(tree_bkg_testing, bkg_weight, "Testing")

dataloader.SetWeightExpression("weight")

mycutSig = ROOT.TCut("")
mycutBkg = ROOT.TCut("")

if isMuon:

    if not evaluate_BDT_systematic and not test_on_signal_shifted_up and not test_on_signal_shifted_down and not test_on_signal_sin2 and not test_on_signal_cos:
        dataloader.PrepareTrainingAndTestTree(mycutSig, mycutBkg, ":".join(["!V","nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0"]))

    if test_on_signal_shifted_up or test_on_signal_shifted_down or test_on_signal_sin2 or test_on_signal_cos:
        dataloader.PrepareTrainingAndTestTree(mycutBkg, ":".join(["!V","nTrain_Background=0:nTest_Background=0"]))

    method_btd  = factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "BDT", ":".join(["H","!V","NTrees=800", "MinNodeSize=2.5%","MaxDepth=3","BoostType=AdaBoost","AdaBoostBeta=0.25","nCuts=20","NegWeightTreatment=IgnoreNegWeightsInTraining"]))

else:

    if not evaluate_BDT_systematic and not test_on_signal_shifted_up and not test_on_signal_shifted_down and not test_on_signal_sin2 and not test_on_signal_cos:
        dataloader.PrepareTrainingAndTestTree(mycutSig, mycutBkg, ":".join(["!V","nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0"])) 
    if test_on_signal_shifted_up or test_on_signal_shifted_down or test_on_signal_sin2 or test_on_signal_cos:
        dataloader.PrepareTrainingAndTestTree(mycutBkg, ":".join(["!V","nTrain_Background=0:nTest_Background=0"]))

    method_btd  = factory.BookMethod(dataloader, ROOT.TMVA.Types.kBDT, "BDT", ":".join(["H","!V","NTrees=800", "MinNodeSize=2.5%","MaxDepth=3","BoostType=AdaBoost","AdaBoostBeta=0.25","nCuts=20","NegWeightTreatment=IgnoreNegWeightsInTraining"]))

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

fOut.Close()

weightfile_dir = "default/weights/TMVAClassification_BDT.weights.xml"

if evaluate_BDT_systematic:
    weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_shifted.xml"
    weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_shifted.xml"
elif test_on_signal_shifted_up:
    weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_Pythia_up.xml"
    weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_Pythia_up.xml"
elif test_on_signal_shifted_down:
    weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_Pythia_down.xml"
    weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_Pythia_down.xml"
elif test_on_signal_sin2:
    weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_sin2_Top_Pi.xml"
    weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_sin2_Top_Pi.xml"
elif test_on_signal_cos:
    weightfile_mu  = "default/weights/TMVAClassification_BDT.weights_mu_cos.xml"
    weightfile_ele = "default/weights/TMVAClassification_BDT.weights_ele_cos.xml"
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
