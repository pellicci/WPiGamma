from ROOT import TFile, TH1F, TCanvas, gStyle, TLegend, TPad, TLine, kFALSE
import copy

test_different_years = True

if test_different_years:

    fInput_2016 = TFile("Tree_input_massfit_MC_0.root")
    fInput_2017 = TFile("Tree_input_massfit_MC_1.root")
    fInput_2018 = TFile("Tree_input_massfit_MC_2.root")
    
    mytree_2016 = fInput_2016.Get("minitree")
    mytree_2017 = fInput_2017.Get("minitree")
    mytree_2018 = fInput_2018.Get("minitree")

    h_signal_2016 = TH1F("h_signal_2016","",100,50.,100.)
    h_signal_2017 = TH1F("h_signal_2017","",100,50.,100.)
    h_signal_2018 = TH1F("h_signal_2018","",100,50.,100.)

    h_ratio_hists = dict()
    h_ratio_list = ["h_ratio_2016_2017","h_ratio_2016_2018","h_ratio_2017_2018"]
    for hname in h_ratio_list:
        h_ratio_hists[hname] = TH1F(hname,"",100,50.,100.)
    
    for entry in mytree_2016:
        
        if not mytree_2016.isSignal:
            continue

        if mytree_2016.weight < 0.:
            Weight = 0.
        else:
            Weight = mytree_2016.weight

        if mytree_2016.Categorization == 1 or mytree_2016.Categorization == 3:
            h_signal_2016.Fill(mytree_2016.Wmass,Weight)

    for entry in mytree_2017:
        
        if not mytree_2017.isSignal:
            continue

        if mytree_2017.weight < 0.:
            Weight = 0.
        else:
            Weight = mytree_2017.weight

        if mytree_2017.Categorization == 1 or mytree_2017.Categorization == 3:
            h_signal_2017.Fill(mytree_2017.Wmass,Weight)

    for entry in mytree_2018:
        
        if not mytree_2018.isSignal:
            continue

        if mytree_2018.weight < 0.:
            Weight = 0.
        else:
            Weight = mytree_2018.weight

        if mytree_2018.Categorization == 1 or mytree_2018.Categorization == 3:
            h_signal_2018.Fill(mytree_2018.Wmass,Weight)


    #Chi2 test
    print "2016 with 2017: ", h_signal_2016.Chi2Test(h_signal_2017,"WW")
    print "2016 with 2018: ", h_signal_2016.Chi2Test(h_signal_2018,"WW")
    print "2017 with 2018: ", h_signal_2017.Chi2Test(h_signal_2018,"WW")

    h_ratio_hists["h_ratio_2016_2017"] = copy.deepcopy(h_signal_2016)
    h_ratio_hists["h_ratio_2016_2018"] = copy.deepcopy(h_signal_2016)
    h_ratio_hists["h_ratio_2017_2018"] = copy.deepcopy(h_signal_2017)

    gStyle.SetOptStat(0)
    canvas_years = TCanvas()
    pad1 = TPad("pad_main_plot","",0,0.28,1,1.)
    pad2 = TPad("pad_ratio_plot","",0,0.01,1,0.27)
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

    pad1.cd()
    legend_years = TLegend(0.15,0.7,0.3,0.95)
    legend_years.SetHeader(" ")
    legend_years.SetFillColor(0)
    legend_years.SetBorderSize(0)
    legend_years.SetLineColor(1)
    legend_years.SetLineStyle(1)
    legend_years.SetLineWidth(1)
    legend_years.SetFillStyle(0)
    legend_years.SetTextSize(0.03)
    legend_years.AddEntry(h_signal_2016,"2016","lep")
    legend_years.AddEntry(h_signal_2017,"2017","lep")
    legend_years.AddEntry(h_signal_2018,"2018","lep")

    scaling_factor_2016_2018 = (h_signal_2018.Integral())/(h_signal_2016.Integral())
    scaling_factor_2017_2018 = (h_signal_2018.Integral())/(h_signal_2017.Integral())
    h_signal_2016.Scale(scaling_factor_2016_2018)
    h_signal_2016.SetMarkerStyle(20)
    h_signal_2016.SetMarkerColor(3)
    h_signal_2017.Scale(scaling_factor_2017_2018)
    h_signal_2017.SetMarkerStyle(20)
    h_signal_2017.SetMarkerColor(2)
    h_signal_2018.SetMarkerStyle(20)
    h_signal_2018.SetMarkerColor(1)
    h_signal_2016.GetXaxis().SetLabelOffset(999)
    h_signal_2016.GetXaxis().SetTitle("#it{m}_{#pi#gamma} (GeV)")
    h_signal_2016.GetYaxis().SetTitle("Events")
    h_signal_2016.Draw()
    h_signal_2017.Draw("SAME")
    h_signal_2018.Draw("SAME")
    legend_years.Draw("SAME")

    pad2.cd()
    pad2.SetTopMargin(0.03)
    pad2.SetFillColor(0)
    pad2.SetFillStyle(0)
    
    for entry in h_ratio_list:
        h_ratio_hists[entry].SetTitle("")
        h_ratio_hists[entry].SetMarkerStyle(8)
        h_ratio_hists[entry].SetLineColor(1)
        h_ratio_hists[entry].GetYaxis().SetLabelSize(0.1)
        h_ratio_hists[entry].GetYaxis().SetTitle("Ratio")
        h_ratio_hists[entry].GetYaxis().SetTitleSize(0.12)
        #h_ratio_hists[entry].GetYaxis().SetTitleOffset(0.98)
        h_ratio_hists[entry].GetYaxis().SetTitleOffset(0.25)
        h_ratio_hists[entry].GetYaxis().SetRangeUser(0.5,1.5)
        h_ratio_hists[entry].GetYaxis().SetNdivisions(502,kFALSE)
        h_ratio_hists[entry].GetXaxis().SetLabelSize(0.10)
        h_ratio_hists[entry].GetXaxis().SetLabelOffset(0.03)
        h_ratio_hists[entry].GetXaxis().SetTitleSize(0.12)
        h_ratio_hists[entry].GetXaxis().SetTitleOffset(1.0)
        h_ratio_hists[entry].GetYaxis().SetRangeUser(0.,2.)

    line_on_one = TLine(h_ratio_hists["h_ratio_2016_2017"].GetXaxis().GetXmin(),1.,h_ratio_hists["h_ratio_2016_2017"].GetXaxis().GetXmax(),1.)
    line_on_one.SetLineColor(4)
    line_on_one.SetLineStyle(2)
    h_ratio_hists["h_ratio_2016_2017"].Divide(h_ratio_hists["h_ratio_2016_2017"],h_signal_2017,1.,1.,"B")
    h_ratio_hists["h_ratio_2016_2018"].Divide(h_ratio_hists["h_ratio_2016_2018"],h_signal_2018,1.,1.,"B")
    h_ratio_hists["h_ratio_2017_2018"].Divide(h_ratio_hists["h_ratio_2017_2018"],h_signal_2018,1.,1.,"B")
    h_ratio_hists["h_ratio_2016_2017"].SetMarkerColor(1)
    h_ratio_hists["h_ratio_2016_2018"].SetMarkerColor(1)
    h_ratio_hists["h_ratio_2017_2018"].SetMarkerColor(1)
    h_ratio_hists["h_ratio_2016_2017"].GetYaxis().SetRangeUser(0.,2.)
    h_ratio_hists["h_ratio_2016_2018"].GetYaxis().SetRangeUser(0.,2.)
    h_ratio_hists["h_ratio_2017_2018"].GetYaxis().SetRangeUser(0.,2.)
    #h_ratio_hists["h_ratio_2016_2017"].Draw()
    #h_ratio_hists["h_ratio_2016_2018"].Draw()
    h_ratio_hists["h_ratio_2017_2018"].Draw()
    line_on_one.Draw("SAME")
    canvas_years.SaveAs("plots/comparison_signal_different_years.pdf")


else:

    fInput = TFile("Tree_input_massfit_MC_3.root")
    
    h_signal_mu  = TH1F("h_signal_mu","",100,50.,100.)
    h_signal_ele = TH1F("h_signal_ele","",100,50.,100.)
    h_channels_ratio = TH1F("h_channels_ratio","",100,50.,100.)
    
    mytree = fInput.Get("minitree")
    
    for entry in mytree:
        
        if not mytree.isSignal:
            continue

        if mytree.weight < 0.:
            Weight = 0.
        else:
            Weight = mytree.weight

        if mytree.Categorization == 1:
            h_signal_mu.Fill(mytree.Wmass,Weight)
        if mytree.Categorization == 3:
            h_signal_ele.Fill(mytree.Wmass,Weight)

    #Chi2 test
    print h_signal_mu.Chi2Test(h_signal_ele,"WW")
    
    gStyle.SetOptStat(0)
    canvas_mu_ele = TCanvas()

    pad1 = TPad("pad_main_plot","",0,0.28,1,1.)
    pad2 = TPad("pad_ratio_plot","",0,0.01,1,0.27)
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

    pad1.cd()
    legend_mu_ele = TLegend(0.15,0.7,0.3,0.95)
    legend_mu_ele.SetHeader(" ")
    legend_mu_ele.SetFillColor(0)
    legend_mu_ele.SetBorderSize(0)
    legend_mu_ele.SetLineColor(1)
    legend_mu_ele.SetLineStyle(1)
    legend_mu_ele.SetLineWidth(1)
    legend_mu_ele.SetFillStyle(0)
    legend_mu_ele.SetTextSize(0.03)
    legend_mu_ele.AddEntry(h_signal_mu,"Muon channel","lep")
    legend_mu_ele.AddEntry(h_signal_ele,"Electron channel","lep")

    scaling_factor = (h_signal_ele.Integral())/(h_signal_mu.Integral())
    h_signal_mu.Scale(scaling_factor)
    h_signal_mu.SetMarkerStyle(20)
    h_signal_mu.SetMarkerColor(2)
    h_signal_ele.SetMarkerStyle(20)
    h_signal_ele.SetMarkerColor(4)
    h_signal_mu.GetXaxis().SetTitle("#it{m}_{#pi#gamma} (GeV)")
    h_signal_mu.GetXaxis().SetLabelOffset(999)
    h_signal_mu.GetYaxis().SetTitle("Events")
    h_signal_mu.GetYaxis().SetTitleSize(0.05)
    h_signal_mu.Draw()
    h_signal_ele.Draw("SAME")
    legend_mu_ele.Draw("SAME")

    pad2.cd()
    pad2.SetTopMargin(0.03)
    pad2.SetFillColor(0)
    pad2.SetFillStyle(0)
    h_channels_ratio = copy.deepcopy(h_signal_mu)
    h_channels_ratio.SetTitle("")
    h_channels_ratio.SetMarkerStyle(8)
    h_channels_ratio.SetMarkerColor(1)
    h_channels_ratio.SetLineColor(1)
    h_channels_ratio.GetYaxis().SetLabelSize(0.1)
    h_channels_ratio.GetYaxis().SetTitle("Ratio")
    h_channels_ratio.GetYaxis().SetTitleSize(0.12)
    #h_channels_ratio.GetYaxis().SetTitleOffset(0.98)
    h_channels_ratio.GetYaxis().SetTitleOffset(0.25)
    h_channels_ratio.GetYaxis().SetRangeUser(0.5,1.5)
    h_channels_ratio.GetYaxis().SetNdivisions(502,kFALSE)
    h_channels_ratio.GetXaxis().SetLabelSize(0.10)
    h_channels_ratio.GetXaxis().SetLabelOffset(0.03)
    h_channels_ratio.GetXaxis().SetTitleSize(0.12)
    h_channels_ratio.GetXaxis().SetTitleOffset(1.0)
    line_on_one = TLine(h_channels_ratio.GetXaxis().GetXmin(),1.,h_channels_ratio.GetXaxis().GetXmax(),1.)
    line_on_one.SetLineColor(4)
    line_on_one.SetLineStyle(2)
    h_channels_ratio.Divide(h_channels_ratio,h_signal_ele,1.,1.,"B")
    h_channels_ratio.SetMarkerColor(1)
    h_channels_ratio.GetYaxis().SetRangeUser(0.,2.)
    h_channels_ratio.Draw()
    line_on_one.Draw("SAME")
    canvas_mu_ele.SaveAs("plots/comparison_signal_different_channels.pdf")

raw_input()

