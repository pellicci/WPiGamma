import ROOT
import math
import numpy as np
import sys

##normalize to the 2016 lumi
luminosity_norm = 36.46

from Workflow_Handler import Workflow_Handler
myWF = Workflow_Handler("Signal","Data",isMedium = True)

#-----Some selection bools------#
is_pi_gamma        = True
is_ele_gamma       = False
is_mu_pT           = False
is_ele_pT          = False
is_deltaphi_mu_pi  = False
is_deltaphi_ele_pi = False

def is_Event_selected(Wmass,nBjets,lepton_iso,ele_gamma_InvMass,lep_pT):#,pi_pT,gamma_eT):
    """Save events according to some basic selection criteria"""
    bjet_cut = nBjets >= 2.

    mass_cut_down = Wmass >= 50.

    mass_cut_up = Wmass <= 100.

    pi_pT_cut = pi_pT >= 80.

    gamma_eT_cut = gamma_eT >= 60.

    if not isMuon:
        ele_iso_cut = lepton_iso <= 0.35
        ele_gamma_InvMass_cut = (ele_gamma_InvMass < 88.5 or ele_gamma_InvMass > 91.5)
        ele_pT_cut = lep_pT >= 29.
    else:
        mu_pT_cut = lep_pT >= 27.

    if isMuon:
        return mass_cut_down and mass_cut_up and bjet_cut and mu_pT_cut #and gamma_eT_cut and pi_pT_cut
    else:
        return mass_cut_down and mass_cut_up and bjet_cut and ele_iso_cut and ele_gamma_InvMass_cut and ele_pT_cut #and gamma_eT_cut and pi_pT_cut

##Here starts the program
Norm_Map = myWF.get_normalizations_map()

#----pi pT and gamma ET cuts----- 
if is_pi_gamma:
    steps_cut1 = 10
    cut1_init = 20.
    cut1_stepsize = 10.

    steps_cut2 = 10
    cut2_init = 20.
    cut2_stepsize = 10.

#----The next one is for lepton pT-----
if is_ele_pT:
    steps_cut1 = 26
    cut1_init = 26.
    cut1_stepsize = 1.

if is_mu_pT:
    steps_cut1 = 25
    cut1_init = 25.
    cut1_stepsize = 1.

#----The next one is for deltaphi------
if is_deltaphi_ele_pi or is_deltaphi_mu_pi:
    steps_cut1 = 30
    cut1_init = 0.
    cut1_stepsize = 0.1

#----The next one is for Z peak--------
if is_ele_gamma:
    steps_cut1 = 30
    cut1_init = 0.
    cut1_stepsize = 0.5

if not is_pi_gamma:
    steps_cut2 = 1
    cut2_init = 20.
    cut2_stepsize = 10.

cut_Nbkg = [[0 for x in range(steps_cut2)] for x in range(steps_cut1)]
cut_Nsig = [[0 for x in range(steps_cut2)] for x in range(steps_cut1)]

##collect all the root files
samplename_list = myWF.get_samples_names()
root_file = myWF.get_root_files()

##loop on the samples and on the cuts and calculate N_bkg and N_sig
for name_sample in samplename_list:

    mytree = root_file[name_sample].Get("WPiGammaAnalysis/mytree")

    if "Signal" in name_sample and not name_sample == myWF.sig_samplename:
        continue

    print "Working on sample ", name_sample

    for jentry in xrange(mytree.GetEntriesFast()):
        ientry = mytree.LoadTree( jentry )
        if ientry < 0:
            break

        nb = mytree.GetEntry( jentry )
        if nb <= 0:
            continue

        if "Data" in name_sample: 
            continue

        if name_sample == "ttbar" and mytree.isttbarlnu:
            continue

        isMuon = mytree.is_muon
        
        if is_ele_gamma and isMuon: 
            continue

        if is_ele_pT and isMuon:
            continue

        if is_mu_pT and not isMuon:
            continue

        if is_deltaphi_ele_pi and isMuon:
            continue

        if is_deltaphi_mu_pi and not isMuon:
            continue

        norm_factor = Norm_Map[name_sample]*luminosity_norm
        PU_Weight = mytree.PU_Weight
        Event_Weight = norm_factor*PU_Weight
        
        lep_pT  = mytree.lepton_pT
        lep_eta = mytree.lepton_eta
        lep_phi = mytree.lepton_phi
        lep_iso = mytree.lepton_iso
        ele_FourMomentum = ROOT.TLorentzVector()
        ele_FourMomentum.SetPtEtaPhiM(lep_pT,lep_eta,lep_phi,0.)
        pi_pT = mytree.pi_pT
        gamma_eta = mytree.photon_eta
        gamma_phi = mytree.photon_phi
        gamma_eT = mytree.photon_eT
        gamma_E = mytree.photon_energy
        gamma_FourMomentum = ROOT.TLorentzVector()
        gamma_FourMomentum.SetPtEtaPhiE(gamma_eT,gamma_eta,gamma_phi,gamma_E)
        ele_gamma_InvMass = (ele_FourMomentum + gamma_FourMomentum).M()

        deltaphi = math.fabs(mytree.lepton_phi-mytree.pi_phi)
        if deltaphi > 3.14:
            deltaphi = 6.28-deltaphi

        if not is_Event_selected(mytree.Wmass,mytree.nBjets_25,mytree.lepton_iso,ele_gamma_InvMass,lep_pT):#,pi_pT,gamma_eT):
            continue

        for icut1 in xrange(steps_cut1):
            cut1_value = cut1_init + cut1_stepsize*icut1

            if is_ele_gamma:
                if ele_gamma_InvMass > 90 - cut1_value and ele_gamma_InvMass < 90 + cut1_value:
                    continue
            
            if is_pi_gamma:
                if pi_pT < cut1_value:
                    continue

            if is_mu_pT or is_ele_pT:
                if lep_pT < cut1_value:
                    continue
            
            if is_deltaphi_mu_pi or is_deltaphi_ele_pi:
                if deltaphi < cut1_value:
                    continue

            for icut2 in xrange(steps_cut2):
                cut2_value = cut2_init + cut2_stepsize*icut2
                
                if is_pi_gamma:
                    if gamma_eT < cut2_value:
                        continue
                else:
                    if gamma_eT < 0: #-----a possible condition to use to optimize a single variable (the precedent)-----
                        continue

                if name_sample == myWF.sig_samplename:
                    cut_Nsig[icut1][icut2] += Event_Weight
                else:
                    cut_Nbkg[icut1][icut2] += Event_Weight

print "Done looping over the events"

##Calculate the significance
signif_list = [[0 for x in range(steps_cut2)] for x in range(steps_cut1)]

cut1_x_list = [cut1_init + cut1_stepsize*icut for icut in range(steps_cut1)]
cut2_x_list = [cut2_init + cut2_stepsize*icut for icut in range(steps_cut2)]

signif_max = -1.
cut1_max = -1
cut2_max = -1

for icut1 in xrange(steps_cut1):
    for icut2 in xrange(steps_cut2):

        if cut_Nbkg[icut1][icut2] != 0:
            signif_list[icut1][icut2] = cut_Nsig[icut1][icut2]/math.sqrt(cut_Nbkg[icut1][icut2])
        else:
            signif_list[icut1][icut2] = 0.

        if signif_list[icut1][icut2] > signif_max:
            signif_max = signif_list[icut1][icut2]
            cut1_max = icut1
            cut2_max = icut2

print "The cut1 value is ", cut1_init + cut1_stepsize*cut1_max
print "The cut2 value is ", cut2_init + cut2_stepsize*cut2_max
print "Number of signal events is ", cut_Nsig[cut1_max][cut2_max]
print "Number of background events is ", cut_Nbkg[cut1_max][cut2_max]
print "Max significance is ", signif_max

##Plot the significance as function of the cuts
cut1_x = np.array(cut1_x_list)
cut2_x = np.array(cut2_x_list)

cut1_y_list = [signif_list[icut][cut2_max] for icut in xrange(steps_cut1)]
cut2_y_list = [signif_list[cut1_max][icut] for icut in xrange(steps_cut2)]

cut1_y = np.array(cut1_y_list)
cut2_y = np.array(cut2_y_list)

graph_cut1 = ROOT.TGraph(steps_cut1,cut1_x,cut1_y)
graph_cut2 = ROOT.TGraph(steps_cut2,cut2_x,cut2_y)

c1 = ROOT.TCanvas("c1","c1",800,500)
c1.cd()
graph_cut1.Draw("A*")
graph_cut1.GetYaxis().SetTitleOffset(1.6)
#graph_cut1.SetTitle("Significance vs p_{T}^{#pi}; p_{T}^{#pi}; Significance")
#graph_cut1.SetTitle("Significance vs 90 #pm #varepsilon; 90 #pm #varepsilon; Significance")
#graph_cut1.SetTitle("Significance vs p_{T}^{e}; p_{T}^{e}; Significance")
#graph_cut1.SetTitle("; p_{T}^{e}; Significance")
#graph_cut1.SetTitle("Significance vs #Delta#varphi(#mu,#pi); #Delta#varphi(#mu,#pi); Significance") #for muons
#graph_cut1.SetTitle("Significance vs #Delta#varphi(e,#pi);  #Delta#varphi(e,#pi); Significance") #for electrons

if is_ele_gamma:
    graph_cut1.SetTitle("; 90 #pm #varepsilon (GeV/c^{2}); Significance")
    c1.SaveAs("~rselvati/www/WPiGamma/Significance/19_09_2017/Z_peak_thesis.pdf")
if is_ele_pT:
    graph_cut1.SetTitle("; p_{T}^{e} (GeV/c); Significance")
    c1.SaveAs("~rselvati/www/WPiGamma/Significance/19_09_2017/e_pT_thesis.pdf")
if is_mu_pT:
    graph_cut1.SetTitle("; p_{T}^{#mu} (GeV/c); Significance")
    c1.SaveAs("~rselvati/www/WPiGamma/Significance/19_09_2017/mu_pT_thesis.pdf")
if is_deltaphi_mu_pi:
    graph_cut1.SetTitle("; #Delta#varphi(#mu,#pi); Significance")
    c1.SaveAs("~rselvati/www/WPiGamma/Significance/19_09_2017/deltaphi_mu_pi_thesis.pdf")
if is_deltaphi_ele_pi:
    graph_cut1.SetTitle("; #Delta#varphi(e,#pi); Significance")
    c1.SaveAs("~rselvati/www/WPiGamma/Significance/19_09_2017/deltaphi_ele_pi_thesis.pdf")
if is_pi_gamma:
    graph_cut1.SetTitle("; p_{T}^{#pi} (GeV/c); Significance")
    c1.SaveAs("~rselvati/www/WPiGamma/Significance/19_09_2017/pi_pT_thesis.pdf")
    c2 = ROOT.TCanvas("c2","c2",800,500)
    c2.cd()
    graph_cut2.Draw("A*")
    graph_cut2.GetYaxis().SetTitleOffset(1.6)
    graph_cut2.SetTitle("; E_{T}^{#gamma} (GeV); Significance")
    c2.SaveAs("~rselvati/www/WPiGamma/Significance/19_09_2017/gamma_eT_thesis.pdf")
    
#c1.SaveAs("plots/pi_pT_signif.png")
#c1.SaveAs("plots/ele_pT_signif.png")
#c1.SaveAs("plots/deltaphi_mu_pi_signif.png") #for muons
#c1.SaveAs("plots/deltaphi_ele_pi_signif.png") #for electrons
"""
c2 = ROOT.TCanvas("c2","c2")
c2.cd()
graph_cut2.Draw("A*")
graph_cut2.GetYaxis().SetTitleOffset(1.5)
graph_cut2.SetTitle("Significance vs E_{T}^{#gamma}; E_{T}^{#gamma}; Significance")
c2.SaveAs("plots/gamma_eT_signif.png")
#c2.SaveAs("plots/nothing.png")
"""
