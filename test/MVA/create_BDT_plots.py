from DataFormats.FWLite import Events, Handle
import os
import math
import sys
import numpy as np
from ROOT import TH1F, TCanvas, TFile, TLegend, gStyle

isMuon = True

def BDT_output():

    if isMuon:
        f = TFile("outputs/Nominal_training_mu_sidebands.root")
    else:
        f = TFile("outputs/Nominal_training_ele_sidebands.root")

    h_BDT_sig = f.Get("Method_BDT/BDT/MVA_BDT_S")
    h_BDT_bkg = f.Get("Method_BDT/BDT/MVA_BDT_B")

    leg1 = TLegend(0.65,0.7,0.8,0.95)
    leg1.SetHeader(" ")
    leg1.SetFillColor(0)
    leg1.SetBorderSize(0)
    leg1.SetLineColor(1)
    leg1.SetLineStyle(1)
    leg1.SetLineWidth(1)
    leg1.SetFillStyle(0)
    leg1.AddEntry(h_BDT_sig,"Signal","f")
    leg1.AddEntry(h_BDT_bkg,"Background","f")
    
    gStyle.SetOptStat(0)
    canvas1 = TCanvas()
    h_BDT_sig.SetTitle("BDT Output")
    h_BDT_sig.SetFillColor(38)
    #h_BDT_sig.SetLineColor(1)
    h_BDT_bkg.SetFillColor(2)
    h_BDT_bkg.SetLineColor(2)
    h_BDT_bkg.SetFillStyle(3002)
    h_BDT_sig.Draw("hist")
    h_BDT_bkg.Draw("SAME, hist")
    h_BDT_sig.SetMaximum(6)
    leg1.Draw("SAME")
    #canvas1.Print("~rselvati/www/WPiGamma/MVA/06_12_2017/double_histo.pdf")

    if isMuon:
        canvas1.Print("plots/double_histo_mu_sidebands.pdf")
    else:
        canvas1.Print("plots/double_histo_ele_sidebands.pdf")

def rejB_vs_S():

    if isMuon:
        f1 = TFile("outputs/Nominal_training_mu_sidebands.root")
        f2 = TFile("outputs/Nominal_training_mu_MC.root")
    else:
        f1 = TFile("outputs/Nominal_training_ele_sidebands.root")
        f2 = TFile("outputs/Nominal_training_ele_MC.root")

    h_rejB_vs_S_1 = f1.Get("Method_BDT/BDT/MVA_BDT_rejBvsS")
    h_rejB_vs_S_2 = f2.Get("Method_BDT/BDT/MVA_BDT_rejBvsS")

    leg2 = TLegend(0.2,0.6,0.6,0.7)
    #leg1.SetHeader(" ")
    leg2.SetFillColor(0)
    leg2.SetBorderSize(0)
    leg2.SetLineColor(1)
    leg2.SetLineStyle(1)
    leg2.SetLineWidth(1)
    leg2.SetFillStyle(0)
    leg2.AddEntry(h_rejB_vs_S_1,"with data sidebands","l")
    leg2.AddEntry(h_rejB_vs_S_2,"with MC","l")

    gStyle.SetOptStat(0)
    canvas2 = TCanvas()
    h_rejB_vs_S_1.SetTitle("Background rejection vs Signal efficiency")
    h_rejB_vs_S_1.SetLineColor(2)
    h_rejB_vs_S_1.SetLineWidth(3)
    h_rejB_vs_S_2.SetLineColor(1)    
    h_rejB_vs_S_2.SetLineWidth(3)
    h_rejB_vs_S_1.Draw("hist")
    h_rejB_vs_S_2.Draw("SAME, hist")
    leg2.Draw("SAME")

    if isMuon:
        canvas2.Print("plots/rejBvsS_mu.pdf")
    else:
        canvas2.Print("plots/rejBvsS_ele.pdf")

if __name__ == "__main__":
    rejB_vs_S()
