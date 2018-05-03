from DataFormats.FWLite import Events, Handle
import os
import math
import sys
import numpy as np
import ROOT
from ROOT import TH1F, TCanvas, TFile, gStyle, gPad

mult_if_signal_mu_1 = TH1F("Mult1 mu","Lepton multiplicity with tag muon",5,-0.5,4.5)
mult_if_signal_mu_2 = TH1F("Muons are blu","Lepton multiplicity with tag muon",5,-0.5,4.5)
mult_if_signal_mu_3 = TH1F("Mupn ch. with loose electrons","Lepton multiplicity with tag muon",5,-0.5,4.5)
mult_if_signal_el_1 = TH1F("Electrons are red","Lepton multiplicity with tag electron",5,-0.5,4.5)
mult_if_signal_el_2 = TH1F("Mult2 el","Lepton multiplicity with tag electron",5,-0.5,4.5)
mult_if_signal_el_3 = TH1F("Ele ch. with loose electrons","Lepton multiplicity with tag muon",5,-0.5,4.5)

f = TFile.Open("LeptonMultiplicity_output_v2.root")
mytree = f.Get("LeptonMultiplicity/mytree")


nevent = 0

for entry in mytree:
    nevent += 1
    if entry.isMuonSignal == True:
        mult_if_signal_mu_1.Fill(mytree.nMuons)
        #mult_if_signal_mu_2.SetLineColor(2)
        mult_if_signal_mu_2.Fill(mytree.nElectrons)
        mult_if_signal_mu_3.Fill(mytree.nElectronsLoose)

nevent = 0

for entry in mytree:
    nevent += 1
    if entry.isEleSignal == True:
        mult_if_signal_el_1.Fill(mytree.nMuons)
        #mult_if_signal_el_2.SetLineColor(2)
        mult_if_signal_el_2.Fill(mytree.nElectrons)
        mult_if_signal_el_3.Fill(mytree.nElectronsLoose)

f.Close()

#newfile = TFile("proof_mult.root","recreate")
legend = ROOT.TLegend(0.6868687,0.6120093,0.9511784,0.9491917)
legend.SetHeader(" ")
legend.SetFillColor(0)
legend.SetBorderSize(0)
legend.SetLineColor(1)
legend.SetLineStyle(1)
legend.SetLineWidth(1)
legend.SetFillStyle(0)
legend.SetTextFont(12)
legend.SetTextSize(0.04)
legend.AddEntry(mult_if_signal_mu_1,"#mu","l")
legend.AddEntry(mult_if_signal_mu_2,"e","l")
legend.AddEntry(mult_if_signal_mu_3,"e (loose)","l")

gStyle.SetOptStat(0)
canvas1 = TCanvas()
mult_if_signal_mu_1.SetTitle("")
mult_if_signal_mu_1.SetLineWidth(2)
mult_if_signal_mu_2.SetLineWidth(2)
mult_if_signal_mu_2.SetLineColor(2)
mult_if_signal_mu_3.SetLineWidth(2)
mult_if_signal_mu_3.SetLineColor(3)
mult_if_signal_mu_3.SetLineStyle(10)
mult_if_signal_mu_1.GetXaxis().SetTitle("Number of leptons")
mult_if_signal_mu_1.GetYaxis().SetTitle("Events")
mult_if_signal_mu_1.GetYaxis().SetTitleOffset(1.5)
mult_if_signal_mu_1.Draw()
mult_if_signal_mu_2.Draw("SAME")
mult_if_signal_mu_3.Draw("SAME")
gPad.SetLogy()
legend.Draw("SAME")
canvas1.Print("Muon_signal_logscale.pdf")

canvas2 = TCanvas()
mult_if_signal_el_1.SetTitle("")
mult_if_signal_el_1.SetLineWidth(2)
mult_if_signal_el_2.SetLineWidth(2)
mult_if_signal_el_2.SetLineColor(2)
mult_if_signal_el_2.SetLineWidth(2)
mult_if_signal_el_3.SetLineColor(3)
mult_if_signal_el_3.SetLineStyle(10)
mult_if_signal_el_1.GetXaxis().SetTitle("Number of leptons")
mult_if_signal_el_1.GetYaxis().SetTitle("Events")
mult_if_signal_el_1.GetYaxis().SetTitleOffset(1.5)
mult_if_signal_el_1.Draw()
mult_if_signal_el_2.Draw("SAME")
mult_if_signal_el_3.Draw("SAME")
gPad.SetLogy()
legend.Draw("SAME")
canvas2.Print("Electron_signal_logscale.pdf")

