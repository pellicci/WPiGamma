import ROOT
import os
import math
import numpy as np
import copy

isMuon = True

if isMuon:
    Nsig_passed = 1
    Nbkg_passed = 112774
else:
    Nsig_passed = 1
    Nbkg_passed = 172296

if isMuon:
    BDT_file = ROOT.TFile("outputs/Nominal_training_mu_met.root")
else:
    BDT_file = ROOT.TFile("outputs/Nominal_training_ele_met.root")

h_BDT_effB_effS = BDT_file.Get("Method_BDT/BDT/MVA_BDT_effBvsS")

#h_sign = ROOT.TF1("h_sign","signif vs effS", 100,0,1)

canvas1 = ROOT.TCanvas()
sig_eff = []
bkg_eff = []
signif = []
#print "nbins: ", BDT_hist.GetNbinsX()
for jbin in range(1,h_BDT_effB_effS.GetNbinsX()):
    if h_BDT_effB_effS.GetBinCenter(jbin) > 0.1:
        sig_eff.append(h_BDT_effB_effS.GetBinCenter(jbin))
        bkg_eff.append(h_BDT_effB_effS.GetBinContent(jbin))
    #if h_BDT_effB_effS.GetBinContent(jbin) == 0:
    #    signif.append(0) # In fact it would be infinite
    #else:
        signif.append(Nsig_passed*h_BDT_effB_effS.GetBinCenter(jbin)/(math.sqrt(Nbkg_passed*h_BDT_effB_effS.GetBinContent(jbin))))

    #if not h_BDT_effB_effS.GetBinContent(jbin) == 0:
    #    print "sig_eff: ", h_BDT_effB_effS.GetBinCenter(jbin), "bkg_eff: ", h_BDT_effB_effS.GetBinContent(jbin), "signif: ", Nsig_passed*h_BDT_effB_effS.GetBinCenter(jbin)/(math.sqrt(Nbkg_passed*h_BDT_effB_effS.GetBinContent(jbin)))


sig_eff_array = np.array(sig_eff)
signif_array = np.array(signif)
sign = ROOT.TGraph(89,sig_eff_array,signif_array)
sign.GetXaxis().SetTitle("Signal efficiency")
sign.GetYaxis().SetTitle("Significance")
sign.SetMaximum(0.15)
sign.SetMarkerStyle(8)
sign.SetMarkerColor(4)
sign.Draw("AP")

if isMuon:
    canvas1.SaveAs("plots/signif_vs_effS_mu_met.pdf")
else:
    canvas1.SaveAs("plots/signif_vs_effS_ele_met.pdf")

#----Now find the BDT output corresponding to the highest significance

h_BDT_effS = BDT_file.Get("Method_BDT/BDT/MVA_BDT_effS")
signif_maximizing_eff = sig_eff_array[np.argmax(signif_array)]
BDT_output = 0.

for entry in xrange(h_BDT_effS.GetNbinsX()):

    effS = h_BDT_effS.GetBinContent(entry)
    effS = float(format(effS, '.3f'))
    signif_maximizing_eff = float(format(signif_maximizing_eff, '.3f'))
    #print "effS: ", effS#, "signif_max_eff: ", signif_maximizing_eff
    #if effS == signif_maximizing_eff:
    if effS == 0.550:
        BDT_output =  h_BDT_effS.GetBinCenter(entry)

print "The BDT output to maximize the significance is :", BDT_output
