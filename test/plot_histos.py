import ROOT
import math
import copy
import sys
import tdrstyle, CMS_lumi

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   
runningEra = int(sys.argv[1])
signal_magnify = int(sys.argv[2])

list_inputfiles = []
for filename in sys.argv[3:]:
    list_inputfiles.append(filename)

#print list_inputfiles
#sample_order = [1,5,0,2,3,4,6,7,8,9,10,11,12,13,14,15,16]
sample_order = [1,4,0,3,5,6,7,8,2]
list_inputfiles = [list_inputfiles[i] for i in sample_order]

#print len(list_inputfiles)

#CMS-style plotting 
tdrstyle.setTDRStyle()
iPeriod = 4
iPos = 11
CMS_lumi.lumiTextSize = 0.9
CMS_lumi.cmsTextSize = 1.5

if runningEra == 0:
    output_dir = "plots/latest_production/2016/"
    CMS_lumi.lumi_13TeV = "35.9 fb^{-1}"
elif runningEra == 1:
    output_dir = "plots/latest_production/2017/"
    CMS_lumi.lumi_13TeV = "41.5 fb^{-1}"
elif runningEra == 2:
    output_dir = "plots/latest_production/2018/"
    CMS_lumi.lumi_13TeV = "59.7 fb^{-1}"
elif runningEra == 3:
    output_dir = "plots/latest_production/2016_2017_2018/"
    CMS_lumi.lumi_13TeV = "137 fb^{-1}"

hstack  = dict()
hsignal = dict()
hdata   = dict()
canvas  = dict()
histo_container = [] #just for memory management

h_nBjet_ratio   = ROOT.TH1F("nBjets_ratio", "nBjets ratio", 6,0,6)
h_met_mu_ratio  = ROOT.TH1F("met_mu_ratio", "met mu ratio",20,0.,200.)
h_met_ele_ratio = ROOT.TH1F("met_ele_ratio", "met ele ratio",20,0.,200.)

list_histos = ["h_mupt", "h_elept", "h_pipt", "h_gammaet", "h_mueta", "h_eleeta","h_pieta","h_gammaeta", "h_nBjets_25","h_deltaphi_mu_pi","h_deltaphi_ele_pi","h_deltaphi_mu_W","h_deltaphi_ele_W","h_deltaeta_mu_pi","h_deltaeta_ele_pi","h_Wmass","h_Wmass_flag_mu","h_Wmass_flag_ele","h_mu_gamma_InvMass","h_ele_gamma_InvMass","h_piRelIso_05_mu_ch","h_piRelIso_05_mu","h_piRelIso_05_ele_ch","h_piRelIso_05_ele","h_met_mu","h_met_ele","h_met_puppi","h_Wmass_alternative_mu","h_Wmass_alternative_ele","h_nPV_mu","h_nPV_ele","h_deltaphi_mu_gamma","h_deltaphi_ele_gamma","h_deltaR_mu_gamma","h_deltaR_ele_gamma","h_lepton_eta","h_lepton_pt","h_piRelIso_05_ch","h_deltaR_mu_pi","h_deltaR_ele_pi","h_nBjets_scaled","h_met_mu_scaled","h_met_ele_scaled","h_njets","h_nBjets_vs_njets","MCT_deltaR_lep_gamma","h_mu_met_mT","h_ele_met_mT","h_met","h_pi_ph_met_InvMass"]#,"h_piRelIso_03","h_piIso_05_mu","h_piIso_05_ele"]

for hname in list_histos:
    hstack[hname] = ROOT.THStack("hstack_" + hname,"")

# Color mask must have the same number of entries as non-QCD backgrounds + 1 (that is the cumulative QCD background)
colors_mask = dict()
colors_mask["QCD"]                 = 12
colors_mask["QCDEM"]               = ROOT.kSpring-3
#colors_mask["ttbar"]               = ROOT.kAzure+7
colors_mask["ttbar"]               = ROOT.kGreen-7
colors_mask["ttbarlnu"]            = ROOT.kYellow-8
colors_mask["ttbarWQQ"]            = ROOT.kOrange-3
colors_mask["ttbarZlnu"]           = ROOT.kTeal-5
colors_mask["DY"]                  = ROOT.kViolet-6
colors_mask["ttbarZQQ"]            = ROOT.kSpring-1
colors_mask["GammaJets"]           = ROOT.kOrange+7
colors_mask["STtW"]                = ROOT.kMagenta+1
colors_mask["WJetsToLNu"]          = ROOT.kGreen+2
colors_mask["WZ"]                  = ROOT.kBlue-7
colors_mask["WGToLNuG"]            = ROOT.kOrange-3
colors_mask["TTGJets"]             = ROOT.kRed-7
colors_mask["ZGTo2LG"]             = ROOT.kAzure+1
colors_mask["WW"]                  = ROOT.kPink+1
colors_mask["ttbarWlnu"]           = ROOT.kOrange+8
colors_mask["Other"]               = ROOT.kMagenta-9

names_mask = dict()
names_mask["QCD"]                 = "QCD multijet"
names_mask["ttbar"]               = "t#bar{t}"
names_mask["DY"]                  = "Drell-Yan"
names_mask["GammaJets"]           = "#bf{#gamma}+jets"
names_mask["STtW"]                = "tW"
names_mask["WJetsToLNu"]          = "W+jets"
names_mask["WZ"]                  = "WZ"
names_mask["WGToLNuG"]            = "W#bf{#gamma}#rightarrowl#bf{#nu#gamma}"
names_mask["TTGJets"]             = "t#bar{t}#gamma+jets"
names_mask["ZGTo2LG"]             = "Z#bf{#gamma}#rightarrowl^{+}l^{-}#bf{#gamma}"
names_mask["WW"]                  = "WW"
names_mask["Other"]               = "Other"

if signal_magnify == 10000:
    leg1 = ROOT.TLegend(0.35,0.62,0.98,0.95) #right positioning
else:
    leg1 = ROOT.TLegend(0.321,0.58,0.981,0.95) #right positioning
leg1.SetHeader(" ")
leg1.SetFillColorAlpha(0,0.)
leg1.SetBorderSize(0)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)
leg1.SetNColumns(2)

for filename in list_inputfiles:
    fileIn = ROOT.TFile(filename)

    #sample_name = filename.split("_")[5]
    sample_name = filename.split("_")[2]
    for histo_name in list_histos:
        histo = fileIn.Get(histo_name)

        # Set to 0 the bins containing negative values, due to negative weights
        hsize = histo.GetSize() - 2 # GetSize() returns the number of bins +2 (that is + overflow + underflow) 
        for bin in range(1,hsize+1): # The +1 is in order to get the last bin
            bincontent = histo.GetBinContent(bin)
            if bincontent < 0.:
                histo.SetBinContent(bin,0.)

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
            histo_container[-1].SetLineColor(colors_mask[sample_name])
            hstack[histo_name].Add(histo_container[-1])

        if histo_name == "h_gammaet": #Add the legend only once (gammaet is just a random variable)

            if histo.Integral() > float(signal_magnify)/12. or sample_name == "Signal": #Only plot in the legend those samples which have some contribution
                if not sample_name == "Data" and not sample_name == "Signal":
                    leg1.AddEntry(histo_container[-1],names_mask[sample_name],"f")
                elif sample_name == "Data":
                    leg1.AddEntry(histo_container[-1],sample_name,"ep")
                elif sample_name == "Signal":
                    #leg1.AddEntry(histo_container[-1],sample_name + " x" + str(signal_magnify),"l")
                    if signal_magnify == 10000:
                        leg1.AddEntry(histo_container[-1],"W#rightarrow#pi#gamma (B = 1%)","l")
                    else:
                        leg1.AddEntry(histo_container[-1],"W#rightarrow#pi#gamma (B = 10^{-4})","l")

    fileIn.Close()

for histo_name in list_histos:

    canvas[histo_name] = ROOT.TCanvas("Canvas_" + histo_name,"",200,106,600,600)
    canvas[histo_name].cd()
 
    ##########################################
    pad1 = ROOT.TPad("pad_" + histo_name,"",0,0.28,1,1.)
    pad2 = ROOT.TPad("pad_" + histo_name,'',0,0.01,1,0.27)
    pad1.SetTopMargin(0.047)
    pad1.SetBottomMargin(0.02)
    pad1.SetBorderMode(0)
    pad1.SetBorderSize(0)
    pad1.SetFrameBorderSize(0)
    pad2.SetBorderSize(0)
    pad2.SetFrameBorderSize(0)
    pad2.SetBottomMargin(0.3)
    pad2.SetBorderMode(0)
    pad1.Draw()
    pad2.Draw()
    ##########################################
    pad1.cd()
    hstack[histo_name].SetTitle("")
    hstack[histo_name].Draw("histo")
    hstack[histo_name].GetYaxis().SetTitleSize(0.07)
    #hstack[histo_name].GetYaxis().SetTitleOffset(1.15)
    hstack[histo_name].GetYaxis().SetTitleOffset(1.)
    hstack[histo_name].GetYaxis().SetTitle("Events/0.2")
    #hstack[histo_name].GetYaxis().SetTitle("Events/2.5 GeV")

    ##########################################
    #hstack[histo_name].GetXaxis().SetTickLength(0)
    hstack[histo_name].GetXaxis().SetLabelOffset(999)
    ##########################################

    if histo_name == "h_Wmass" or histo_name == "h_Wmass_flag_mu" or histo_name == "h_Wmass_flag_ele":
        hstack[histo_name].GetXaxis().SetTitle("#it{m}_{#pi#gamma} (GeV)")
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),30000))
        #hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),200.))

    if histo_name == "h_mupt":
        hstack[histo_name].GetXaxis().SetTitle("p_{T}^{#mu} (GeV)")
        #hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),12000))
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),60.))

    if histo_name == "h_elept":
        hstack[histo_name].GetXaxis().SetTitle("p_{T}^{e} (GeV)")
        #hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),4000.))
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),60.))

    if histo_name == "h_lepton_eta":
        hstack[histo_name].GetXaxis().SetTitle("\eta^{#ell}")
        #hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),40000.))
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),200.))

    if histo_name == "h_lepton_pt":
        hstack[histo_name].GetXaxis().SetTitle("p_{T}^{#ell} (GeV)")
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),65000))
        #hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),180.))
    
    if histo_name == "h_pipt":
        hstack[histo_name].GetXaxis().SetTitle("p_{T}^{#pi} (GeV)")
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),120000.))
        #hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),160.))

    if histo_name == "h_gammaet":
        hstack[histo_name].GetXaxis().SetTitle("p_{T}^{#gamma} (GeV)")
        #hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),210.))
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),120000.))

    if histo_name == "h_gammaeta":
        hstack[histo_name].GetXaxis().SetTitle("#eta^{#gamma}")
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),200.))

    if histo_name == "h_mueta":
        hstack[histo_name].GetXaxis().SetTitle("#eta^{#mu}")
        #hstack[histo_name].SetMaximum(60)
        #hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),60))

    if histo_name == "h_eleeta":
        hstack[histo_name].GetXaxis().SetTitle("#eta^{e}")
        #hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),50))
        #hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),60))

    if histo_name == "h_pieta":
        hstack[histo_name].GetXaxis().SetTitle("#eta^{#pi}")
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),180.))

    if histo_name == "h_deltaphi_mu_pi":
        hstack[histo_name].GetXaxis().SetTitle("#Delta#varphi(#mu,#pi)")
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),35000.))

    if histo_name == "h_deltaphi_ele_pi":
        hstack[histo_name].GetXaxis().SetTitle("#Delta#varphi(e,#pi)")
        #hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),45000))
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),11000.))

    if histo_name == "h_deltaphi_ele_gamma":
        hstack[histo_name].GetXaxis().SetTitle("#Delta#varphi(e,#gamma)")
        #canvas[histo_name].SetLogy()
        #hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),15000))

    if histo_name == "h_deltaphi_mu_gamma":
        hstack[histo_name].GetXaxis().SetTitle("#Delta#varphi(#mu,#gamma)")
        #canvas[histo_name].SetLogy()
        #hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),7000))

    if histo_name == "h_deltaR_mu_pi":
        hstack[histo_name].GetXaxis().SetTitle("#Delta R(#mu,#pi)")

    if histo_name == "h_deltaR_ele_pi":
        hstack[histo_name].GetXaxis().SetTitle("#Delta R(e,#pi)")

    if histo_name == "h_deltaR_mu_gamma":
        hstack[histo_name].GetXaxis().SetTitle("#Delta R(#mu,#gamma)")

    if histo_name == "h_deltaR_ele_gamma":
        hstack[histo_name].GetXaxis().SetTitle("#Delta R(e,#gamma)")
        #hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),60000.))

    if histo_name == "h_ele_gamma_InvMass":
        hstack[histo_name].GetXaxis().SetTitle("m_{e#gamma} (GeV/c^{2})")
        #hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),10000))
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),80))

    if histo_name == "h_mu_gamma_InvMass":
        hstack[histo_name].GetXaxis().SetTitle("m_{#mu#gamma} (GeV/c^{2})")
        #hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),7000))

    if histo_name == "h_nBjets":
        hstack[histo_name].GetXaxis().SetTitle("Number of b-jets")

    if histo_name == "h_nBjets_25":
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),200000))
        #hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),850))
        hstack[histo_name].GetXaxis().SetTitle("#it{n}_{b}")

    if "h_piRelIso" in histo_name:
        hstack[histo_name].GetXaxis().SetTitle("#pi_{Iso}/p_{T}^{#pi}")
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),800.))    

    if histo_name == "h_met":
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),70000.))    

    if histo_name == "h_mu_met_mT":
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),80.))    

    if histo_name == "h_ele_met_mT":
        hstack[histo_name].SetMaximum(max(hstack[histo_name].GetHistogram().GetMaximum(),60.))    

    if signal_magnify != 1:
        hsignal[histo_name].Scale(signal_magnify)      

    hstack[histo_name].Draw("SAME,histo")
    hsignal[histo_name].Draw("SAME,hist")

    if signal_magnify == 100 and histo_name == "h_Wmass":
        fIn_Wmass_frame = ROOT.TFile("WmassAnalysis/Wmass_frame.root")
        Wmass_frame = fIn_Wmass_frame.Get("massplot")
        Wmass_frame.Draw("SAME")
        #fIn_Wmass_frame.Close()
        fIn_dummy_plot = ROOT.TFile("dummy_plot.root")
        dummy_plot = fIn_dummy_plot.Get("f1")
        #leg1.AddEntry(dummy_plot,"Fit","l")
    else:
        hdata[histo_name].Draw("SAME,E1,X0")

    hMCErr = copy.deepcopy(hstack[histo_name].GetStack().Last())
    #Add systematic uncertainties to the background
    hMCErr_size = hMCErr.GetSize() - 2
    if signal_magnify == 100:
        for bin in range(1,hMCErr_size+1):
            new_BinError = hMCErr.GetBinError(bin) + (0.01+0.014+0.015)*hMCErr.GetBinContent(bin) 
            hMCErr.SetBinError(bin,new_BinError)

    hMCErr.SetFillStyle(3005)
    hMCErr.SetMarkerStyle(0)
    hMCErr.SetFillColor(ROOT.kBlack)
    hMCErr.SetLineColor(0)
    hMCErr.Draw("sameE2")

    if histo_name == "h_gammaet":# or histo_name == "h_pipt": #Add the legend only once (gammaet is just a random variable)
        leg1.AddEntry(hMCErr,"MC unc (stat + syst)","f")
    leg1.Draw()
    CMS_lumi.CMS_lumi(pad1, iPeriod, iPos) #Print integrated lumi and energy information

    ################################################
    pad2.cd()
    pad2.SetTopMargin(0.03)
    pad2.SetFillColor(0)
    pad2.SetFillStyle(0)
    ROOT.gStyle.SetOptStat(0)
    totalMC = copy.deepcopy(hMCErr)
    totalData = copy.deepcopy(hdata[histo_name])
    totalData_forErrors = copy.deepcopy(hdata[histo_name])
    totalData.Divide(totalMC)

    for bin in range(1,hMCErr_size+1):
        
        #Set MC error band to MC relative uncertainty
        if not totalMC.GetBinContent(bin) == 0:
            new_MC_BinError = totalMC.GetBinError(bin)/totalMC.GetBinContent(bin)
        else:
            new_MC_binError = 0.

        #Set data/MC ratio points error bar to data relative uncertainty
        if not totalData_forErrors.GetBinContent(bin) == 0:
            new_Data_BinError = totalData_forErrors.GetBinError(bin)/totalData_forErrors.GetBinContent(bin)
        else:
            new_Data_BinError = 0.

        totalMC.SetBinError(bin,new_MC_BinError)
        totalMC.SetBinContent(bin,1.)
        totalData.SetBinError(bin,new_Data_BinError)
    
    # if histo_name == "h_piRelIso_05_ch":
    #     deviation = 0.
    #     totalData_size = totalData.GetSize() - 2
    #     for bin in range(1,totalData_size+1): # The +1 is in order to get the last bin
    #         totalData_content = totalData.GetBinContent(bin)
    #         deviation += math.fabs(1. - totalData_content)
    #     deviation = deviation/totalData_size
    #     print "this is the deviation: ", deviation

    totalData.SetTitle("")
    totalData.SetMarkerStyle(8)
    totalData.SetMarkerColor(1)
    totalData.SetLineColor(1)
    totalData.GetYaxis().SetLabelSize(0.1)
    totalData.GetYaxis().SetTitle("Data/MC")
    totalData.GetYaxis().SetTitleSize(0.16)
    #totalData.GetYaxis().SetTitleOffset(0.48)
    totalData.GetYaxis().SetTitleOffset(0.42)
    totalData.GetYaxis().SetRangeUser(0.,2.)
    totalData.GetYaxis().SetNdivisions(502,ROOT.kFALSE)
    totalData.GetXaxis().SetLabelSize(0.10)
    totalData.GetXaxis().SetTitleSize(0.12)
    totalData.GetXaxis().SetTitleOffset(1.0)

    totalMC.SetTitle("")
    totalMC.SetFillStyle(3002)

    if histo_name == "h_Wmass" or histo_name == "h_Wmass_flag_mu" or histo_name == "h_Wmass_flag_ele":
        totalData.GetXaxis().SetTitle("#it{m}_{#pi#gamma} (GeV)")
        totalData.GetXaxis().SetTitleSize(0.17)
        totalData.GetXaxis().SetTitleOffset(0.78)

    if histo_name == "h_mupt":
        totalData.GetXaxis().SetTitle("p_{T}^{#mu} (GeV)")

    if histo_name == "h_elept":
        totalData.GetXaxis().SetTitle("p_{T}^{e} (GeV)")
    
    if histo_name == "h_pipt":
        totalData.GetXaxis().SetTitle("#it{p}_{T}^{#pi} (GeV)")

    if histo_name == "h_gammaet":
        totalData.GetXaxis().SetTitle("#it{p}_{T}^{#gamma} (GeV)")

    if histo_name == "h_gammaeta":
        totalData.GetXaxis().SetTitle("#eta^{#gamma}")

    if histo_name == "h_mueta":
        totalData.GetXaxis().SetTitle("#eta^{#mu}")

    if histo_name == "h_eleeta":
        totalData.GetXaxis().SetTitle("#eta^{e}")

    if histo_name == "h_lepton_eta":
        totalData.GetXaxis().SetTitle("\eta^{l}")

    if histo_name == "h_lepton_pt":
        totalData.GetXaxis().SetTitle("#it{p}_{T}^{l} (GeV)")

    if histo_name == "h_pieta":
        totalData.GetXaxis().SetTitle("#eta^{#pi}")

    if histo_name == "h_deltaphi_mu_pi":
        totalData.GetXaxis().SetTitle("#Delta#varphi(#mu,#pi)")

    if histo_name == "h_deltaphi_ele_pi":
        totalData.GetXaxis().SetTitle("#Delta#varphi(e,#pi)")

    if histo_name == "h_deltaphi_ele_gamma":
        totalData.GetXaxis().SetTitle("#Delta#varphi(e,#gamma)")

    if histo_name == "h_deltaphi_mu_gamma":
        totalData.GetXaxis().SetTitle("#Delta#varphi(#mu,#gamma)")

    if histo_name == "h_deltaR_mu_pi":
        totalData.GetXaxis().SetTitle("#DeltaR(#mu,#pi)")

    if histo_name == "h_deltaR_ele_pi":
        totalData.GetXaxis().SetTitle("#DeltaR(e,#pi)")

    if histo_name == "h_deltaR_ele_gamma":
        totalData.GetXaxis().SetTitle("#DeltaR(e,#gamma)")

    if histo_name == "h_deltaR_mu_gamma":
        totalData.GetXaxis().SetTitle("#DeltaR(#mu,#gamma)")

    if histo_name == "h_ele_gamma_InvMass":
        totalData.GetXaxis().SetTitle("m_{e#gamma} (GeV/c^{2})")

    if histo_name == "h_mu_gamma_InvMass":
        totalData.GetXaxis().SetTitle("m_{#mu#gamma} (GeV/c^{2})")

    if histo_name == "h_nBjets":
        totalData.GetXaxis().SetTitle("Number of b-jets")

    if histo_name == "h_nBjets_25":
        totalData.GetXaxis().SetTitle("#it{n}_{b}")

    if "h_piRelIso" in histo_name:
        totalData.GetXaxis().SetTitle("#Sigma#it{p}_{T}/#it{p}_{T}^{#pi}")
        totalData.GetXaxis().SetTitleSize(0.17)
        totalData.GetXaxis().SetTitleOffset(0.78)

    if histo_name == "h_mu_met_mT":
        totalData.GetXaxis().SetTitle("m_{T}(#mu,MET) (GeV)")

    if histo_name == "h_ele_met_mT":
        totalData.GetXaxis().SetTitle("m_{T}(e,MET) (GeV)")

    line_on_one = ROOT.TLine(totalData.GetXaxis().GetXmin(),1.,totalData.GetXaxis().GetXmax(),1.)
    line_on_one.SetLineColor(4)
    line_on_one.SetLineStyle(2)

    totalData.Draw("E1,X0")
    totalMC.Draw("sameE2")
    line_on_one.Draw("SAME")
    ################################################

    canvas[histo_name].SaveAs(output_dir + histo_name + ".pdf")
 
# Wmass_nominal_cut_mu = copy.deepcopy(hstack["h_Wmass_flag_mu"].GetStack().Last())
# Wmass_alternative_cut_mu = copy.deepcopy(hstack["h_Wmass_alternative_mu"].GetStack().Last())

# Wmass_nominal_cut_ele = copy.deepcopy(hstack["h_Wmass_flag_ele"].GetStack().Last())
# Wmass_alternative_cut_ele = copy.deepcopy(hstack["h_Wmass_alternative_ele"].GetStack().Last())

# Wmass_nominal_cut_mu.Scale(1./Wmass_nominal_cut_mu.Integral())
# Wmass_alternative_cut_mu.Scale(1./Wmass_alternative_cut_mu.Integral())

# Wmass_nominal_cut_ele.Scale(1./Wmass_nominal_cut_ele.Integral())
# Wmass_alternative_cut_ele.Scale(1./Wmass_alternative_cut_ele.Integral())

# Wmass_nominal_cut_mu.Divide(Wmass_nominal_cut_mu,Wmass_alternative_cut_mu,1.,1.,"B")
# Wmass_nominal_cut_ele.Divide(Wmass_nominal_cut_ele,Wmass_alternative_cut_ele,1.,1.,"B")

# canvas_Wmass_1 = ROOT.TCanvas()
# Wmass_nominal_cut_mu.SetMarkerStyle(21)
# Wmass_nominal_cut_mu.GetYaxis().SetRangeUser(0.,2.)
# Wmass_nominal_cut_mu.SetTitle("")
# Wmass_nominal_cut_mu.GetXaxis().SetTitle("m_{#pi#gamma} (GeV)")
# Wmass_nominal_cut_mu.Draw("Pe")
# canvas_Wmass_1.SaveAs(output_dir + "Wmass_ratio_mu_Wmass.pdf")

# canvas_Wmass_2 = ROOT.TCanvas()
# Wmass_nominal_cut_ele.SetMarkerStyle(21)
# Wmass_nominal_cut_ele.GetYaxis().SetRangeUser(0.,2.)
# Wmass_nominal_cut_ele.SetTitle("")
# Wmass_nominal_cut_ele.GetXaxis().SetTitle("m_{(#pi#gamma#)} (GeV)")
# Wmass_nominal_cut_ele.Draw("Pe")
# canvas_Wmass_2.SaveAs(output_dir + "Wmass_ratio_ele_Wmass.pdf")

# h_nBjets_ratio = copy.deepcopy(hstack["h_nBjets_25"].GetStack().Last())
# h_nBjets_ratio.Divide(h_nBjets_ratio,copy.deepcopy(hstack["h_nBjets_scaled"].GetStack().Last()),1.,1.,"B")
# canvas_Bjets = ROOT.TCanvas()
# h_nBjets_ratio.GetXaxis().SetTitle("Number of b-jets (p_{T} > 25 GeV)")
# h_nBjets_ratio.Draw()
# canvas_Bjets.SaveAs("plots/latest_production/2018/h_nBjets_ratio.pdf")

# h_met_mu_ratio = copy.deepcopy(hstack["h_met_mu"].GetStack().Last())
# h_met_mu_ratio.Divide(h_met_mu_ratio,copy.deepcopy(hstack["h_met_mu_scaled"].GetStack().Last()),1.,1.,"B")
# canvas_met_mu = ROOT.TCanvas()
# h_met_mu_ratio.Draw("hist")
# canvas_met_mu.SaveAs("plots/latest_production/2018/h_met_mu_ratio.pdf")

# h_met_ele_ratio = copy.deepcopy(hstack["h_met_ele"].GetStack().Last())
# h_met_ele_ratio.Divide(h_met_ele_ratio,copy.deepcopy(hstack["h_met_ele_scaled"].GetStack().Last()),1.,1.,"B")
# canvas_met_ele = ROOT.TCanvas()
# h_met_ele_ratio.Draw()
# canvas_met_ele.SaveAs("plots/latest_production/2018/h_met_ele_ratio.pdf")
