from DataFormats.FWLite import Events, Handle
import os
import math
import sys
import numpy as np
from ROOT import TH1F, TCanvas, TFile, TLegend, gStyle
import argparse

#---------------------------------#
p = argparse.ArgumentParser(description='Select whether the plots will regard the muon or the electron channels')
p.add_argument('isMuon_option', help='Type <<muon>> or <<electron>>')
#p.add_argument('year_option', help='Type <<2016>> or <<2017>> or <<2018>>')
args = p.parse_args()

# Switch from muon to electron channel
if args.isMuon_option == "muon":
    isMuon = True
if args.isMuon_option == "electron":
    isMuon = False

#year = args.year_option

#---------------------------------#

def BDT_output():

    if isMuon:
        f = TFile("outputs/Nominal_training_mu.root")
    else:
        f = TFile("outputs/Nominal_training_ele.root")

    h_BDT_sig = f.Get("default/Method_BDT/BDT/MVA_BDT_S")
    h_BDT_bkg = f.Get("default/Method_BDT/BDT/MVA_BDT_B")

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
    h_BDT_sig.SetTitle("")
    h_BDT_bkg.SetTitle("")
    h_BDT_sig.SetFillColor(38)
    #h_BDT_sig.SetLineColor(1)
    h_BDT_bkg.SetFillColor(2)
    h_BDT_bkg.SetLineColor(2)
    h_BDT_bkg.SetFillStyle(3002)
    h_BDT_sig.Draw("hist")
    h_BDT_bkg.Draw("SAME, hist")
    h_BDT_sig.SetMaximum(4.5)
    leg1.Draw("SAME")


    if isMuon:
        canvas1.Print("plots/BDT_output_mu.pdf")
    else:
        canvas1.Print("plots/BDT_output_ele.pdf")

    raw_input()

def rejB_vs_S():

    if isMuon:
        #f1 = TFile("outputs/Nominal_training_mu_Pythia_down.root")
        #f1 = TFile("outputs/Nominal_training_mu_Wmass.root")
        f1 = TFile("outputs/Nominal_training_mu_shifted.root")
        f2 = TFile("outputs/Nominal_training_mu.root")
    else:
        #f1 = TFile("outputs/Nominal_training_ele_Pythia_down.root")
        #f1 = TFile("outputs/Nominal_training_ele_Wmass.root")
        f1 = TFile("outputs/Nominal_training_ele_shifted.root")
        f2 = TFile("outputs/Nominal_training_ele.root")

    h_rejB_vs_S_1 = f1.Get("default/Method_BDT/BDT/MVA_BDT_rejBvsS")
    h_rejB_vs_S_2 = f2.Get("default/Method_BDT/BDT/MVA_BDT_rejBvsS")

    leg2 = TLegend(0.2,0.6,0.6,0.7)
    #leg1.SetHeader(" ")
    leg2.SetFillColor(0)
    leg2.SetBorderSize(0)
    leg2.SetLineColor(1)
    leg2.SetLineStyle(1)
    leg2.SetLineWidth(1)
    leg2.SetFillStyle(0)
    #leg2.AddEntry(h_rejB_vs_S_1,"with data sidebands","l")
    #leg2.AddEntry(h_rejB_vs_S_1,"with #Delta#varphi_{l,#gamma}","l")
    #leg2.AddEntry(h_rejB_vs_S_1,"with m_{#pi#gamma}","l")
    leg2.AddEntry(h_rejB_vs_S_1,"shifted","l")
    #leg2.AddEntry(h_rejB_vs_S_1,"with signal scaled down","l")
    #leg2.AddEntry(h_rejB_vs_S_2,"with MC","l")
    #leg2.AddEntry(h_rejB_vs_S_2,"without #Delta#varphi_{l,#gamma}","l")
    #leg2.AddEntry(h_rejB_vs_S_2,"without m_{#pi#gamma}","l")
    leg2.AddEntry(h_rejB_vs_S_2,"nominal","l")
    #leg2.AddEntry(h_rejB_vs_S_2,"with signal non scaled","l")

    gStyle.SetOptStat(0)
    canvas2 = TCanvas()
    h_rejB_vs_S_1.SetTitle(" ")
    h_rejB_vs_S_1.SetLineColor(1)
    h_rejB_vs_S_1.SetLineWidth(3)
    h_rejB_vs_S_2.SetLineColor(2)    
    h_rejB_vs_S_2.SetLineWidth(3)
    h_rejB_vs_S_1.Draw("hist")
    h_rejB_vs_S_2.Draw("SAME, hist")
    leg2.Draw("SAME")

    if isMuon:
        #canvas2.Print("plots/rejBvsS_mu_Wmass.pdf")
        canvas2.Print("plots/rejBvsS_mu_shifted.pdf")
        #canvas2.Print("plots/rejBvsS_mu_Pythia_down.pdf")
    else:
        #canvas2.Print("plots/rejBvsS_ele_Wmass.pdf")
        canvas2.Print("plots/rejBvsS_ele_shifted.pdf")
        #canvas2.Print("plots/rejBvsS_ele_Pythia_down.pdf")

    raw_input()

if __name__ == "__main__":

    rejB_vs_S()
    #BDT_output()
