from ROOT import *

fileIn  = TFile("totalMC_medium_2018.root")
fileOut = TFile("bTagEff_medium_2018.root","RECREATE")

h2_Num = fileIn.Get("WPiGammaAnalysis/h2_BTaggingEff_Num_b")
h2_Den = fileIn.Get("WPiGammaAnalysis/h2_BTaggingEff_Denom_b")

h2_Num.Divide(h2_Den)

canvas = TCanvas()
gStyle.SetOptStat(0)
h2_Num.Draw("COLZ")

fileOut.cd()
h2_Num.Write("bTagEfficiency")

raw_input()
