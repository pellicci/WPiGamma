import ROOT
import os
import math
import numpy as np

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

from Workflow_Handler import Workflow_Handler
myWF = Workflow_Handler("Signal")

##Global constants
MU_MIN_PT    = 24.
ELE_MIN_PT   = 24.
PI_MIN_PT    = 24.
GAMMA_MIN_ET = 24.
#WMASS_MIN    = 10.
#WMASS_MAX    = 100.

#Normalize to this luminsity, in fb-1
luminosity_norm = 36.46

#Make signal histos larger
signal_magnify = 1500

output_dir = "plots"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#Here's the list of histos to plot
list_histos = ["h_mupt", "h_elept", "h_pipt", "h_gammaet", "h_Wmass", "h_nBjets", "h_mueta", "h_eleeta"]

#Color mask must have the same number of entries as non-QCD backgrounds
colors_mask = [1,400,840,616,860,432,880,416,800,900,820,920,405]   

def select_all_but_one(cutstring):

    selection_bools = dict()
    selection_bools["h_mupt"]    = lep_pt > MU_MIN_PT
    selection_bools["h_elept"]   = lep_pt > ELE_MIN_PT
    selection_bools["h_pipt"]    = pi_pt > PI_MIN_PT
    selection_bools["h_gammaet"] = gamma_et > GAMMA_MIN_ET
    #selection_bools["h_Wmass"]   = Wmass > WMASS_MIN and Wmass < WMASS_MAX

    result = True

    for hname in selection_bools:
        if cutstring == hname:
            continue
        else:
            result = result and selection_bools[hname]

    return result

##Here starts the program
Norm_Map = myWF.get_normalizations_map()

##Get the files and the names of the samples
samplename_list = myWF.get_samples_names()
root_file = myWF.get_root_files()

##Get the handlers for all the histos and graphics
hs      = dict()
h_base  = dict()
canvas  = dict()

for hname in list_histos:
    hs[hname] = ROOT.THStack("hs_" + hname,"")
    canvas[hname] = ROOT.TCanvas(hname,hname,200,106,600,600)

##Define the histos to be created
isQCDfirst = True
for sample_name in samplename_list:

    if "QCD" in sample_name:
        theSampleName = "QCD_"
        if not isQCDfirst:
            continue
        isQCDfirst = False
    else:
        theSampleName = sample_name

    h_base[theSampleName+list_histos[0]]  = ROOT.TH1F(theSampleName+list_histos[0], "p_{T} of the muon", 25, MU_MIN_PT, 300.)
    h_base[theSampleName+list_histos[1]]  = ROOT.TH1F(theSampleName+list_histos[1], "p_{T} of the electron", 25, ELE_MIN_PT, 300.)
    h_base[theSampleName+list_histos[2]]  = ROOT.TH1F(theSampleName+list_histos[2], "p_{T} of the pion", 25, PI_MIN_PT, 300.)
    h_base[theSampleName+list_histos[3]]  = ROOT.TH1F(theSampleName+list_histos[3], "e_{T} of the gamma", 25, GAMMA_MIN_ET, 300.)
    h_base[theSampleName+list_histos[4]]  = ROOT.TH1F(theSampleName+list_histos[4], "W mass", 40, 10, 125)
    h_base[theSampleName+list_histos[5]]  = ROOT.TH1F(theSampleName+list_histos[5], "n Bjets", 6, 0, 6.)
    h_base[theSampleName+list_histos[6]]  = ROOT.TH1F(theSampleName+list_histos[6], "eta of the muon", 50, -3, 3)
    h_base[theSampleName+list_histos[7]]  = ROOT.TH1F(theSampleName+list_histos[7], "eta of the electron", 50, -3, 3)

leg1 = ROOT.TLegend(0.6868687,0.6120093,0.9511784,0.9491917)
leg1.SetHeader(" ")
leg1.SetFillColor(0)
leg1.SetBorderSize(0)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)

Nsig_passed = 0.
Nbkg_passed = 0.

##Loop on samples, and then on events, and merge QCD stuff
idx_sample = 0
isFirstQCDlegend = True

for name_sample in samplename_list:

    theSampleName = name_sample

    if "QCD" in name_sample:
        QCDflag = True
        theSampleName = "QCD_"
    else:
        QCDflag = False

    if "DY_50" in name_sample:
        DYflag = True
        #theSampleName = "DY_"
    else:
        DYflag = False

    norm_factor = Norm_Map[name_sample]*luminosity_norm

    mytree = root_file[name_sample].Get("WPiGammaAnalysis/mytree")
 
    print "Processing Sample ", name_sample
    for jentry in xrange(mytree.GetEntriesFast()):
        ientry = mytree.LoadTree( jentry )
        if ientry < 0:
            break
        nb = mytree.GetEntry(jentry )
        if nb <= 0:
            continue

        Event_Weight = norm_factor   #Add other event weights here if necessary

        #This is how you access the tree variables
        lep_pt  = mytree.lepton_pT
        lep_eta = mytree.lepton_eta
        lep_phi = mytree.lepton_phi

        isMuon = mytree.is_muon

        pi_pt = mytree.pi_pT
        gamma_et = mytree.photon_eT

        Wmass = mytree.Wmass

        nBjets = mytree.nBjets

        if select_all_but_one("h_mupt") and isMuon:
            h_base[theSampleName+"h_mupt"].Fill(lep_pt,Event_Weight)

        if select_all_but_one("h_elept") and not isMuon:
            h_base[theSampleName+"h_elept"].Fill(lep_pt,Event_Weight)

        if isMuon:
            h_base[theSampleName+"h_mueta"].Fill(lep_eta,Event_Weight)

        if not isMuon:
            h_base[theSampleName+"h_eleeta"].Fill(lep_eta,Event_Weight)

        #if select_all_but_one("h_pipt"):
        h_base[theSampleName+"h_pipt"].Fill(pi_pt,Event_Weight)

        #if select_all_but_one("h_gammapt"):
        h_base[theSampleName+"h_gammaet"].Fill(gamma_et,Event_Weight)

        #if select_all_but_one("h_Wmass"):
        h_base[theSampleName+"h_Wmass"].Fill(Wmass,Event_Weight)

        #if select_all_but_one("h_nBjets"):
        h_base[theSampleName+"h_nBjets"].Fill(nBjets,Event_Weight)

        #Count the events
        if select_all_but_one("all cuts"):
            if name_sample == myWF.sig_samplename:
                Nsig_passed += Event_Weight
            else:
                Nbkg_passed += Event_Weight

    for idx_histo,hname in enumerate(list_histos):

        if QCDflag:
            h_base[theSampleName+hname].SetFillColor(12)
        elif name_sample == myWF.sig_samplename:
            h_base[theSampleName+hname].SetLineStyle(2)   #dashed
            h_base[theSampleName+hname].SetLineColor(4)   #blue
            h_base[theSampleName+hname].SetLineWidth(4)   #kind of thick
        else:
            h_base[theSampleName+hname].SetFillColor(colors_mask[idx_sample])
            #h_base[theSampleName+hname].SetLineStyle(1)
        #if DYflag:
        #    h_base[theSampleName+hname].SetLineStyle(1)
        #    h_base[theSampleName+hname].SetLineColor(2)
        #    h_base[theSampleName+hname].SetLineWidth(4)

        if idx_histo == 0:
            if QCDflag and isFirstQCDlegend:
                 leg1.AddEntry(h_base[theSampleName+hname],"QCD","f")
                 isFirstQCDlegend = False
            elif name_sample == myWF.sig_samplename:
                 sample_legend_name = "1500 x " + name_sample
                 leg1.AddEntry(h_base[name_sample+hname], sample_legend_name,"f")  #To comment when signal is has to be excluded.
            elif not QCDflag:
                 leg1.AddEntry(h_base[theSampleName+hname],theSampleName,"f")

        if not QCDflag:
            hs[hname].Add(h_base[theSampleName+hname])

    if not QCDflag and not name_sample == myWF.sig_samplename:
        idx_sample += 1

print "Finished runnning over samples!"

for idx_histo,hname in enumerate(list_histos):
    hs[hname].Add(h_base["QCD_"+hname])

x_axis_list = ["p_{T} of the muon", "p_{T} of the electron", "p_{T} of the pion", "p_{T} of the photon", "m_{W}", "eta of the muon", "eta of the electron"]

for hname in list_histos: 
    canvas[hname].cd()
    
    hs[hname].Draw("histo")
    
    #Graphic names
    
    hs[hname].SetTitle(" ")
    hs[hname].GetYaxis().SetTitle("Events")
    hs[hname].GetXaxis().SetTitle(hname)
    
    if signal_magnify != 1:
        h_base[myWF.sig_samplename+hname].Scale(signal_magnify)
        h_base[myWF.sig_samplename+hname].Draw("SAME,hist")
        
        leg1.Draw()
        
        canvas[hname].SaveAs("plots/tree_" + hname + ".gif")

print "Number of expected events for ", luminosity_norm, " in fb-1"
print "Number of signal events passed = ", Nsig_passed
print "Number of background events passed = ", Nbkg_passed
print "Significance S/sqrt(B) = ", Nsig_passed/math.sqrt(Nbkg_passed)
print "Significance S/sqrt(B + deltaB^2) = ", Nsig_passed/(math.sqrt(Nbkg_passed) + 0.2*Nbkg_passed)
print "Significance S/sqrt(S+B) = ", Nsig_passed/math.sqrt(Nsig_passed + Nbkg_passed)
print "\nAll the intresting plots have been produced..!"
