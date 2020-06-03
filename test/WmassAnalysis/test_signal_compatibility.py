from ROOT import TFile, TH1F, TCanvas

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

    canvas_2016 = TCanvas()
    canvas_2016.cd()
    h_signal_2016.Draw()
    canvas_2017 = TCanvas()
    canvas_2017.cd()
    h_signal_2017.Draw()
    canvas_2018 = TCanvas()
    canvas_2018.cd()
    h_signal_2018.Draw()


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

    canvas_mu = TCanvas()
    canvas_mu.cd()
    h_signal_mu.Draw()
    canvas_ele = TCanvas()
    canvas_ele.cd()
    h_signal_ele.Draw()

raw_input()
