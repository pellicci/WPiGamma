from ROOT import *

file_MC = TFile("WmassAnalysis/Tree_input_massfit_MC_3.root")
mytree  = file_MC.Get("minitree")

file_Data   = TFile("WmassAnalysis/Tree_input_massfit_Data_3.root")
mytree_Data = file_Data.Get("minitree")

h_MW_vs_BDT_sig_mu  = TH2F("sig_mu", "",100,-0.8,0.8,90,45.,105.)
h_MW_vs_BDT_sig_ele = TH2F("sig_ele", "",100,-0.8,0.8,90,45.,105.)
h_MW_vs_BDT_bkg_mu  = TH2F("bkg_mu", "",100,-0.8,0.8,90,45.,105.)
h_MW_vs_BDT_bkg_ele = TH2F("bkg_ele", "",100,-0.8,0.8,90,45.,105.)

h_MW_vs_BDT_data_mu  = TH2F("data_mu", "",100,-0.8,0.8,90,45.,105.)
h_MW_vs_BDT_data_ele = TH2F("data_ele", "",100,-0.8,0.8,90,45.,105.)

for event in mytree:

    if event.isMuon and event.isSignal:
        h_MW_vs_BDT_sig_mu.Fill(event.BDT_out,event.Wmass)
    if (not event.isMuon) and event.isSignal:
        h_MW_vs_BDT_sig_ele.Fill(event.BDT_out,event.Wmass)
    if event.isMuon and not event.isSignal:
        h_MW_vs_BDT_bkg_mu.Fill(event.BDT_out,event.Wmass)
    if not event.isMuon and not event.isSignal:
        h_MW_vs_BDT_bkg_ele.Fill(event.BDT_out,event.Wmass)

for event in mytree_Data:

    if event.isMuon:
        h_MW_vs_BDT_data_mu.Fill(event.BDT_out,event.Wmass)
    else:
        h_MW_vs_BDT_data_ele.Fill(event.BDT_out,event.Wmass)

gStyle.SetOptStat(0)
canvas_sig_mu  = TCanvas()
canvas_sig_mu.cd()
h_MW_vs_BDT_sig_mu.SetMarkerColor(4)
h_MW_vs_BDT_sig_mu.GetXaxis().SetTitle("BDT output")
h_MW_vs_BDT_sig_mu.GetYaxis().SetTitle("M_{#pi#gamma} (GeV)")
h_MW_vs_BDT_sig_mu.Draw("box")
canvas_sig_mu.SaveAs("plots/latest_production/2016_2017_2018/Wmass_vs_BDT_mu_sig.pdf")

canvas_sig_ele = TCanvas()
canvas_sig_ele.cd()
h_MW_vs_BDT_sig_ele.SetMarkerColor(4)
h_MW_vs_BDT_sig_ele.GetXaxis().SetTitle("BDT output")
h_MW_vs_BDT_sig_ele.GetYaxis().SetTitle("M_{#pi#gamma} (GeV)")
h_MW_vs_BDT_sig_ele.Draw("box")
canvas_sig_ele.SaveAs("plots/latest_production/2016_2017_2018/Wmass_vs_BDT_ele_sig.pdf")

canvas_bkg_mu  = TCanvas()
canvas_bkg_mu.cd()
h_MW_vs_BDT_bkg_mu.SetMarkerColor(4)
h_MW_vs_BDT_bkg_mu.GetXaxis().SetTitle("BDT output")
h_MW_vs_BDT_bkg_mu.GetYaxis().SetTitle("M_{#pi#gamma} (GeV)")
h_MW_vs_BDT_bkg_mu.Draw("box")
canvas_bkg_mu.SaveAs("plots/latest_production/2016_2017_2018/Wmass_vs_BDT_mu_bkg.pdf")

canvas_bkg_ele = TCanvas()
canvas_bkg_ele.cd()
h_MW_vs_BDT_bkg_ele.SetMarkerColor(4)
h_MW_vs_BDT_sig_ele.SetMarkerColor(4)
h_MW_vs_BDT_bkg_ele.GetXaxis().SetTitle("BDT output")
h_MW_vs_BDT_bkg_ele.GetYaxis().SetTitle("M_{#pi#gamma} (GeV)")
h_MW_vs_BDT_bkg_ele.Draw("box")
canvas_bkg_ele.SaveAs("plots/latest_production/2016_2017_2018/Wmass_vs_BDT_ele_bkg.pdf")

canvas_data_mu  = TCanvas()
canvas_data_mu.cd()
h_MW_vs_BDT_data_mu.SetMarkerColor(4)
h_MW_vs_BDT_data_mu.GetXaxis().SetTitle("BDT discriminant")
h_MW_vs_BDT_data_mu.GetYaxis().SetTitle("m_{#pi#gamma} (GeV)")
h_MW_vs_BDT_data_mu.Draw("box")
canvas_data_mu.SaveAs("plots/latest_production/2016_2017_2018/Wmass_vs_BDT_mu_data.pdf")

canvas_data_ele = TCanvas()
canvas_data_ele.cd()
h_MW_vs_BDT_data_ele.SetMarkerColor(4)
h_MW_vs_BDT_sig_ele.SetMarkerColor(4)
h_MW_vs_BDT_data_ele.GetXaxis().SetTitle("BDT discriminant")
h_MW_vs_BDT_data_ele.GetYaxis().SetTitle("m_{#pi#gamma} (GeV)")
h_MW_vs_BDT_data_ele.Draw("box")
canvas_data_ele.SaveAs("plots/latest_production/2016_2017_2018/Wmass_vs_BDT_ele_data.pdf")


raw_input()
