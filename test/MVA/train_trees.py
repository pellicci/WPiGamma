import ROOT

fIn_bkg = ROOT.TFile("Tree_MC_Background.root")
tree_bkg = fIn_bkg.Get("minitree_background")

fIn_sig = ROOT.TFile("Tree_MC_Signal.root")
tree_sig = fIn_sig.Get("minitree_signal")

fOut = ROOT.TFile("outputs/Nominal_training.root","RECREATE")

ROOT.TMVA.Tools.Instance()

factory = ROOT.TMVA.Factory("TMVAClassification", fOut,":".join(["!V","Transformations=I;D;P;G,D","AnalysisType=Classification"]))

factory.AddVariable("pi_pT","F")
factory.AddVariable("gamma_eT","F")
factory.AddVariable("n_bjets","I")
factory.AddVariable("deltaphi_lep_pi","F")

#factory.AddSpectator("j4_btag","F")

factory.AddSignalTree(tree_sig)
factory.AddBackgroundTree(tree_bkg)

factory.SetWeightExpression("weight")

mycuts = ROOT.TCut("weight > 0.")
mycutb = ROOT.TCut("weight > 0.")

factory.PrepareTrainingAndTestTree(mycuts, mycutb, ":".join(["!V"]) )

method_cuts = factory.BookMethod(ROOT.TMVA.Types.kCuts,"Cuts",":".join(["!H","!V","FitMethod=MC","EffSel","SampleSize=200000","VarProp=FSmart"]))
method_btd  = factory.BookMethod(ROOT.TMVA.Types.kBDT, "BDT", ":".join(["!H","!V","NTrees=800", "MinNodeSize=2.5%","MaxDepth=3","BoostType=AdaBoost","AdaBoostBeta=0.5","nCuts=20"]))

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

fOut.Close()
