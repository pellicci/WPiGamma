from ROOT import TFile, TH1F, TCanvas, gStyle, TLegend

test_different_years = False

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

    # canvas_2016 = TCanvas()
    # canvas_2016.cd()
    # h_signal_2016.Draw()
    # canvas_2017 = TCanvas()
    # canvas_2017.cd()
    # h_signal_2017.Draw()
    # canvas_2018 = TCanvas()
    # canvas_2018.cd()
    # h_signal_2018.Draw()

    gStyle.SetOptStat(0)
    canvas_years = TCanvas()
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
    h_signal_2016.GetXaxis().SetTitle("m_{#pi#gamma} (GeV)")
    h_signal_2016.GetYaxis().SetTitle("Events")
    h_signal_2016.Draw()
    h_signal_2017.Draw("SAME")
    h_signal_2018.Draw("SAME")
    legend_years.Draw("SAME")
    canvas_years.SaveAs("plots/comparison_signal_different_years.pdf")


else:

    fInput = TFile("Tree_input_massfit_MC_3.root")
    
    h_signal_mu  = TH1F("h_signal_mu","",100,50.,100.)
    h_signal_ele = TH1F("h_signal_ele","",100,50.,100.)
    
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

    # canvas_mu = TCanvas()
    # canvas_mu.cd()
    # h_signal_mu.Draw()
    # canvas_ele = TCanvas()
    # canvas_ele.cd()
    # h_signal_ele.Draw()
    
    gStyle.SetOptStat(0)
    canvas_mu_ele = TCanvas()
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
    h_signal_mu.GetXaxis().SetTitle("m_{#pi#gamma} (GeV)")
    h_signal_mu.GetYaxis().SetTitle("Events")
    h_signal_mu.Draw()
    h_signal_ele.Draw("SAME")
    legend_mu_ele.Draw("SAME")
    canvas_mu_ele.SaveAs("plots/comparison_signal_different_channels.pdf")

raw_input()
