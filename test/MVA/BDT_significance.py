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
    Nsig_passed = 1.26 # Number of signal and background events from the sum of the weights (before applying BDT cuts)
    Nbkg_passed = 86898.70
else:
    Nsig_passed = 0.83
    Nbkg_passed = 44931.54

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

#print "nbins: ",h_BDT_effB_effS.GetNbinsX()
for jbin in range(1,h_BDT_effB_effS.GetNbinsX()+1):
    if h_BDT_effB_effS.GetBinCenter(jbin) > 0.3:
        sig_eff.append(h_BDT_effB_effS.GetBinCenter(jbin))
        if h_BDT_effB_effS.GetBinContent(jbin) < 0.:
            bkg_eff.append(0.)
        else:
            bkg_eff.append(h_BDT_effB_effS.GetBinContent(jbin))
        #print h_BDT_effB_effS.GetBinContent(jbin)

    #if h_BDT_effB_effS.GetBinContent(jbin) == 0:
    #    signif.append(0) # In fact it would be infinite
    #else:
        if h_BDT_effB_effS.GetBinContent(jbin) < 0.:
            signif.append(0)
        else:
            signif.append(Nsig_passed*h_BDT_effB_effS.GetBinCenter(jbin)/math.sqrt(Nbkg_passed*h_BDT_effB_effS.GetBinContent(jbin)))

    #if not h_BDT_effB_effS.GetBinContent(jbin) == 0:
    #    print "sig_eff: ", h_BDT_effB_effS.GetBinCenter(jbin), "bkg_eff: ", h_BDT_effB_effS.GetBinContent(jbin), "signif: ", Nsig_passed*h_BDT_effB_effS.GetBinCenter(jbin)/(math.sqrt(Nbkg_passed*h_BDT_effB_effS.GetBinContent(jbin)))


sig_eff_array = np.array(sig_eff)
signif_array = np.array(signif)
#print "signif_len: ", len(signif_array)
sign = ROOT.TGraph(70,sig_eff_array,signif_array)
sign.SetTitle("")
sign.GetXaxis().SetTitle("Signal efficiency")
sign.GetYaxis().SetTitle("Significance")
sign.SetMaximum(0.10)
sign.SetMarkerStyle(8)
sign.SetMarkerColor(4)
sign.Draw("AP")

if isMuon:
    canvas1.SaveAs("plots/signif_vs_effS_mu.pdf")
else:
    canvas1.SaveAs("plots/signif_vs_effS_ele.pdf")

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
    if effS == 0.685:
        BDT_output =  h_BDT_effS.GetBinCenter(entry)
        _effS = effS

print "For a signal efficiency of ", _effS, "the BDT output is :", BDT_output

raw_input()
