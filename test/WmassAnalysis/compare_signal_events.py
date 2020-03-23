import ROOT

ROOT.gROOT.ProcessLineSync(".L dCB/RooDoubleCBFast.cc+")
#12665
#9627
#9740

#Retrive the sample
fInput_2016 = ROOT.TFile("Signal_model_0.root")
fInput_2016.cd()

massplot_2016 = fInput_2016.Get("massplot")

fInput_2017 = ROOT.TFile("Signal_model_1.root")
fInput_2017.cd()

massplot_2017 = fInput_2017.Get("massplot")

fInput_2018 = ROOT.TFile("Signal_model_2.root")
fInput_2018.cd()

massplot_2018 = fInput_2018.Get("massplot")

canvas = ROOT.TCanvas()
massplot_2016.Draw()
massplot_2017.Draw("SAME")
massplot_2018.Draw("SAME")
canvas.SaveAs("plots/comparison_signal_different_years.pdf")

raw_input()
