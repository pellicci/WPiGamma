import ROOT
import copy
import sys

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

runningEra = int(sys.argv[1])
signal_magnify = int(sys.argv[2])

list_inputfiles = []
for filename in sys.argv[3:]:
    list_inputfiles.append(filename)

if runningEra == 0:
    output_dir = "plots/latest_production/2016/"
elif runningEra == 1:
    output_dir = "plots/latest_production/2017/"

hstack = dict()
hsignal = dict()
hdata = dict()
canvas = dict()
histo_container = [] #just for memory management

list_histos = ["h_mupt", "h_elept", "h_pipt", "h_gammaet", "h_mueta", "h_eleeta","h_pieta","h_gammaeta", "h_nBjets_25","h_deltaphi_mu_pi","h_deltaphi_ele_pi","h_deltaphi_mu_W","h_deltaphi_ele_W","h_deltaeta_mu_pi","h_deltaeta_ele_pi","h_Wmass","h_Wmass_flag_mu","h_Wmass_flag_ele","h_mu_gamma_InvMass","h_ele_gamma_InvMass","h_piRelIso_05_mu_ch","h_piRelIso_05_mu","h_piRelIso_05_ele_ch","h_piRelIso_05_ele","h_met_mu","h_met_ele","h_met_puppi","h_nPV_mu","h_nPV_ele","h_deltaphi_mu_gamma","h_deltaphi_ele_gamma"]

for hname in list_histos:
    hstack[hname] = ROOT.THStack("hstack_" + hname,"")

# Color mask must have the same number of entries as non-QCD backgrounds
colors_mask = dict()
colors_mask["QCD"]                 = 12
colors_mask["ttbar"]               = ROOT.kYellow-8
colors_mask["ttbarlnu"]            = ROOT.kAzure+7
colors_mask["ttbarWQQ"]            = ROOT.kOrange-3
colors_mask["ttbarZlnu"]           = ROOT.kTeal-5
colors_mask["DY"]                  = ROOT.kViolet-6
colors_mask["ttbarZQQ"]            = ROOT.kSpring-1
colors_mask["GammaJets"]           = ROOT.kOrange+7
colors_mask["STtW"]                = ROOT.kMagenta+1
colors_mask["WJetsToLNu"]          = ROOT.kGreen+2
colors_mask["WZ"]                  = ROOT.kBlue-7
colors_mask["WGToLNuG"]            = ROOT.kRed-7
colors_mask["TTGJets"]             = ROOT.kOrange-2
colors_mask["ZGTo2LG"]             = ROOT.kYellow+3
colors_mask["WW"]                  = ROOT.kPink+1
colors_mask["ttbarWlnu"]           = ROOT.kOrange+8

# leg1 = ROOT.TLegend(0.15,0.6120093,0.34,0.9491917) #left positioning
leg1 = ROOT.TLegend(0.6868687,0.6120093,0.9511784,0.9491917) #right positioning
leg1.SetHeader(" ")
leg1.SetFillColor(0)
leg1.SetBorderSize(0)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)

for filename in list_inputfiles:
    fileIn = ROOT.TFile(filename)

    sample_name = filename.split("_")[2]

    for histo_name in list_histos:
        histo = fileIn.Get(histo_name)

        histo_container.append(copy.copy(histo))

        if "Signal" in sample_name:
            histo_container[-1].SetLineStyle(2)   #dashed
            histo_container[-1].SetLineColor(2)   #red
            histo_container[-1].SetLineWidth(4)   #kind of thick
            hsignal[histo_name] = histo_container[-1]
        elif "Data" in sample_name:
            histo_container[-1].SetMarkerStyle(20)   #point
            hdata[histo_name] = histo_container[-1]
        else:
            histo_container[-1].SetFillColor(colors_mask[sample_name])
            hstack[histo_name].Add(histo_container[-1])

        if histo_name == "h_gammaet":
            leg1.AddEntry(histo_container[-1],sample_name,"f")

    fileIn.Close()

for histo_name in list_histos:

    canvas[histo_name] = ROOT.TCanvas("Cavas_" + histo_name,"",200,106,600,600)
    canvas[histo_name].cd()

    hstack[histo_name].SetTitle("")
    hstack[histo_name].Draw("histo")

    if histo_name == "h_Wmass" or histo_name == "h_Wmass_flag_mu" or histo_name == "h_Wmass_flag_ele":
        hstack[histo_name].GetXaxis().SetTitle("m_{#pi#gamma} (GeV)")

    if histo_name == "h_mupt":
        hstack[histo_name].GetXaxis().SetTitle("p_{T}^{#mu} (GeV)")

    if histo_name == "h_elept":
        hstack[histo_name].GetXaxis().SetTitle("p_{T}^{e} (GeV)")
    
    if histo_name == "h_pipt":
        hstack[histo_name].GetXaxis().SetTitle("p_{T}^{#pi} (GeV)")

    if histo_name == "h_gammaet":
        hstack[histo_name].GetXaxis().SetTitle("E_{T}^{#gamma} (GeV)")

    if histo_name == "h_gammaeta":
        hstack[histo_name].GetXaxis().SetTitle("#eta^{#gamma}")

    if histo_name == "h_mueta":
        hstack[histo_name].GetXaxis().SetTitle("#eta^{#mu}")

    if histo_name == "h_eleeta":
        hstack[histo_name].GetXaxis().SetTitle("#eta^{e}")

    if histo_name == "h_pieta":
        hstack[histo_name].GetXaxis().SetTitle("#eta^{#pi}")

    if histo_name == "h_deltaphi_mu_pi":
        hstack[histo_name].GetXaxis().SetTitle("#Delta#varphi_{#mu-#pi}")

    if histo_name == "h_deltaphi_ele_pi":
        hstack[histo_name].GetXaxis().SetTitle("#Delta#varphi_{e-#pi}")

    if histo_name == "h_ele_gamma_InvMass":
        hstack[histo_name].GetXaxis().SetTitle("m_{e#gamma} (GeV/c^{2})")

    if histo_name == "h_nBjets":
        hstack[histo_name].GetXaxis().SetTitle("Number of b-jets")

    if histo_name == "h_nBjets_25":
        hstack[histo_name].GetXaxis().SetTitle("Number of b-jets (p_{T}>25 GeV/c)")

    if "h_piRelIso" in histo_name:
        hstack[histo_name].GetXaxis().SetTitle("#pi_{Iso}/p_{T}^{#pi}")
    
    if signal_magnify != 1:
        hsignal[histo_name].Scale(signal_magnify)      

    hstack[histo_name].Draw("histo")
    hsignal[histo_name].Draw("SAME,hist")
    hdata[histo_name].Draw("SAME,E1")

    leg1.Draw()

    canvas[histo_name].SaveAs(output_dir + histo_name + ".pdf")

