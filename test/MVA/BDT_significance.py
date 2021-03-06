import ROOT
import math
import numpy as np
import argparse

#---------------------------------#
p = argparse.ArgumentParser(description='Select whether the plots will regard the muon or the electron channel')
p.add_argument('isMuon_option', help='Type <<muon>> or <<electron>>')
args = p.parse_args()

# Switch from muon to electron channel
if args.isMuon_option == "muon":
    isMuon = True
if args.isMuon_option == "electron":
    isMuon = False

#---------------------------------#

if isMuon:
    Nsig_passed = 1.24 + 1.34 + 2.06# Number of signal and background events from the sum of the weights (before applying BDT cuts)
    Nbkg_passed = 58602.88 + 68965.63 + 110551.74
    #Nsig_passed = 2.05938
    #Nbkg_passed = 40448.7
else:
    Nsig_passed = 0.73 + 0.78 + 1.22
    Nbkg_passed = 34051.46 + 35861.46 + 76162.44
    #Nsig_passed = 1.17629
    #Nbkg_passed = 15148.9

if isMuon:
    BDT_file = ROOT.TFile("outputs/Nominal_training_mu.root")
    #BDT_file = ROOT.TFile("outputs/Nominal_training_mu_Wmass.root")
else:
    BDT_file = ROOT.TFile("outputs/Nominal_training_ele.root")
    #BDT_file = ROOT.TFile("outputs/Nominal_training_ele_Wmass.root")

h_BDT_effB_effS = BDT_file.Get("default/Method_BDT/BDT/MVA_BDT_effBvsS")

canvas1 = ROOT.TCanvas()
sig_eff = []
bkg_eff = []
signif = []
_effS = 0

for jbin in range(1,h_BDT_effB_effS.GetNbinsX()+1):
    if h_BDT_effB_effS.GetBinCenter(jbin) > 0.3:
        sig_eff.append(h_BDT_effB_effS.GetBinCenter(jbin))
        if h_BDT_effB_effS.GetBinContent(jbin) < 0.:
            bkg_eff.append(0.)
            signif.append(0)
        else:
            bkg_eff.append(h_BDT_effB_effS.GetBinContent(jbin))
            signif.append(Nsig_passed*h_BDT_effB_effS.GetBinCenter(jbin)/math.sqrt(Nbkg_passed*h_BDT_effB_effS.GetBinContent(jbin)))


sig_eff_array = np.array(sig_eff)
bkg_eff_array = np.array(bkg_eff)
signif_array = np.array(signif)
#print "signif_len: ", len(signif_array)
sign = ROOT.TGraph(70,sig_eff_array,signif_array)
sign.SetTitle("")
sign.GetXaxis().SetTitle("#varepsilon_{S}^{BDT}")
sign.GetYaxis().SetTitle("Z")
sign.SetMaximum(0.20)
sign.SetMarkerStyle(8)
sign.SetMarkerColor(4)
sign.Draw("AP")

canvas2 = ROOT.TCanvas()
sign_vs_bkg = ROOT.TGraph(70,bkg_eff_array,signif_array)
sign_vs_bkg.SetTitle("")
sign_vs_bkg.GetXaxis().SetTitle("#varepsilon_{B}^{BDT}")
sign_vs_bkg.GetYaxis().SetTitle("Z")
sign_vs_bkg.GetXaxis().SetRangeUser(0.,0.02)
sign_vs_bkg.SetMaximum(0.20)
sign_vs_bkg.SetMarkerStyle(8)
sign_vs_bkg.SetMarkerColor(4)
sign_vs_bkg.Draw("AP")

if isMuon:
    canvas1.SaveAs("plots/signif_vs_effS_mu.pdf")
    canvas2.SaveAs("plots/signif_vs_effB_mu.pdf")
else:
    canvas1.SaveAs("plots/signif_vs_effS_ele.pdf")
    canvas2.SaveAs("plots/signif_vs_effB_ele.pdf")

#----Now find the BDT output corresponding to the highest significance

h_BDT_effS = BDT_file.Get("default/Method_BDT/BDT/MVA_BDT_effS")
signif_maximizing_eff = sig_eff_array[np.argmax(signif_array)]
BDT_output = 0.

for entry in xrange(h_BDT_effS.GetNbinsX()):

    effS = h_BDT_effS.GetBinContent(entry)
    effS = float(format(effS, '.3f'))
    signif_maximizing_eff = float(format(signif_maximizing_eff, '.3f'))
    #print "effS: ", effS#, "signif_max_eff: ", signif_maximizing_eff
    #if effS == signif_maximizing_eff:
    if effS == 0.55:
        BDT_output =  h_BDT_effS.GetBinCenter(entry)
        _effS = effS

print "For a signal efficiency of ", _effS, "the BDT output is :", BDT_output

raw_input()
