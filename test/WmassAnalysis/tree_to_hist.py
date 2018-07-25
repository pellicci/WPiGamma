from ROOT import *

isData = False

if isData:
    fOpen = TFile("Tree_input_massfit_Data_prova.root")
else:
    fOpen = TFile("Tree_input_massfit_MC_prova.root")

mytree = fOpen.Get("minitree")

Wmass_mu = TH1F("Wmass_mu","Wmass mu",15,50,100)
Wmass_ele = TH1F("Wmass_ele","Wmass ele",15,50,100)

for event in mytree:
    if isData:
        if mytree.Categorization==1:
            Wmass_mu.Fill(mytree.Wmass)
        if mytree.Categorization==3:
            Wmass_ele.Fill(mytree.Wmass)
    else:
        if mytree.Categorization==1 and mytree.isSignal==0:
            Wmass_mu.Fill(mytree.Wmass,mytree.weight)
        if mytree.Categorization==3 and mytree.isSignal==0:
            Wmass_ele.Fill(mytree.Wmass,mytree.weight)

canvas1 = TCanvas("canvas1","canvas1",200,106,600,600)
gStyle.SetOptStat(0)
Wmass_mu.SetAxisRange(0.,65.,"Y")
Wmass_mu.Draw("hist")
if isData:
    canvas1.SaveAs("Wmass_mu_tree_to_hist_DATA.pdf")
else:
    canvas1.SaveAs("Wmass_mu_tree_to_hist.pdf")

canvas2 = TCanvas("canvas2","canvas2",200,106,600,600)
gStyle.SetOptStat(0)
Wmass_ele.SetAxisRange(0.,65.,"Y")
Wmass_ele.Draw("hist")
if isData:
    canvas2.SaveAs("Wmass_ele_tree_to_hist_DATA.pdf")
else:
    canvas2.SaveAs("Wmass_ele_tree_to_hist.pdf")
    
raw_input()
