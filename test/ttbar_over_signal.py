from ROOT import *

fIn_ttbar  = TFile("ttbar_pT.root")
fIn_signal = TFile("signal_pT.root")

h_ttbar  = fIn_ttbar.Get("W_ttbar")
h_signal = fIn_signal.Get("W_signal")

h_ttbar.Scale(1./h_ttbar.Integral())
h_signal.Scale(1./h_signal.Integral())

h_ttbar.Divide(h_signal)

fOut = TFile.Open("ttbar_signal_ratio.root","RECREATE")
fOut.cd()
h_ttbar.Write("ttbar_signal_ratio")
fOut.Close()
