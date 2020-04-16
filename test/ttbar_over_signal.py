from ROOT import *

fIn_ttbar  = TFile("ttbar_pT.root")
fIn_signal = TFile("signal_pT.root")

h_ttbar  = fIn_ttbar.Get("W_ttbar")
h_signal = fIn_signal.Get("W_signal")

h_ttbar.Scale(1./h_ttbar.Integral())
h_signal.Scale(1./h_signal.Integral())

h_ttbar.Divide(h_signal)

canvas = TCanvas()
gStyle.SetOptStat(0)
h_ttbar.SetTitle("")
h_ttbar.SetMarkerStyle(21)
h_ttbar.SetLineColor(46)
h_ttbar.GetXaxis().SetTitle("p_{T}^{W} (GeV)")
h_ttbar.Draw()
canvas.SaveAs("ttbar_signal_ratio_2018.pdf")

fOut = TFile.Open("ttbar_signal_ratio_2018.root","RECREATE")
fOut.cd()
h_ttbar.Write("ttbar_signal_ratio")
fOut.Close()
