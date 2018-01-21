import ROOT
import math
import numpy as np
import sys
#sys.path.insert(0, '/afs/cern.ch/user/r/rselvati/wpigamma/CMSSW_8_0_28/src/StandardModel/WPiGamma/test/scripts')
from scripts.create_normalization_table import get_xsec_fromsample

##normalize to the 2016 lumi
luminosity_norm = 36.46

from Workflow_Handler import Workflow_Handler
myWF = Workflow_Handler("Signal","Data",isMedium = True)

#-----Some selection bools------#
is_pi_gamma        = False
is_ele_gamma       = True
is_mu_pT           = False
is_ele_pT          = False
is_deltaphi_mu_pi  = False
is_deltaphi_ele_pi = False
is_piIso_03        = False
is_piIso_05        = False


def is_Event_selected(Wmass,nBjets,lepton_iso,ele_gamma_InvMass):#,lep_pT):#,pi_pT,gamma_eT):
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
        return mass_cut_down and mass_cut_up #and bjet_cut #and mu_pT_cut #and gamma_eT_cut and pi_pT_cut
    else:
        return mass_cut_down and mass_cut_up and ele_iso_cut #and ele_gamma_InvMass_cut #and nbjet_cut #and ele_pT_cut #and gamma_eT_cut and pi_pT_cut

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
    cut1_stepsize = 0.1

if not is_pi_gamma:
    steps_cut2 = 1
    cut2_init = 20.
    cut2_stepsize = 10.

#----The next one is for pion isolation--------
if is_piIso_03 or is_piIso_05:
    steps_cut1 = 40
    cut1_init = 0.
    cut1_stepsize = 0.05

cut_Nbkg = [[0 for x in range(steps_cut2)] for x in range(steps_cut1)]
cut_Nsig = [[0 for x in range(steps_cut2)] for x in range(steps_cut1)]
cut_Nbkg_err = [[0 for x in range(steps_cut2)] for x in range(steps_cut1)]
cut_Nsig_err = [[0 for x in range(steps_cut2)] for x in range(steps_cut1)]
sigma_Nexp_sig =  [[0 for x in range(steps_cut2)] for x in range(steps_cut1)]
sigma_Nexp_bkg =  [[0 for x in range(steps_cut2)] for x in range(steps_cut1)]

##collect all the root files
samplename_list = myWF.get_samples_names()
root_file = myWF.get_root_files()

##loop on the samples and on the cuts and calculate N_bkg and N_sig
for name_sample in samplename_list:
   
    mytree = root_file[name_sample].Get("WPiGammaAnalysis/mytree")

    if "Signal" in name_sample and not name_sample == myWF.sig_samplename:
        continue

    print "Working on sample ", name_sample

    Nsig_passed = [[0 for x in range(steps_cut2)] for x in range(steps_cut1)]
    Nbkg_passed = [[0 for x in range(steps_cut2)] for x in range(steps_cut1)]

    if "Data" in name_sample: continue

    norm_factor = Norm_Map[name_sample]*luminosity_norm
    xsec = float(get_xsec_fromsample(name_sample))*1000
    N_evts_processed = xsec/Norm_Map[name_sample] # Obtain the number of events processed for each sample

    for jentry in xrange(mytree.GetEntriesFast()):
        ientry = mytree.LoadTree( jentry )
        if ientry < 0:
            break

        nb = mytree.GetEntry( jentry )
        if nb <= 0:
            continue

        if name_sample == "ttbar" and mytree.isttbarlnu:
            continue

        #if name_sample == "DY_50": continue

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

        #if isMuon:
        #    continue

        lep_pT  = mytree.lepton_pT
        lep_eta = mytree.lepton_eta
        lep_phi = mytree.lepton_phi
        lep_iso = mytree.lepton_iso
        ele_FourMomentum = ROOT.TLorentzVector()
        ele_FourMomentum.SetPtEtaPhiM(lep_pT,lep_eta,lep_phi,0.)

        pi_pT = mytree.pi_pT
        piIso_03 = mytree.sum_pT_03/pi_pT
        piIso_05 = mytree.sum_pT_05/pi_pT

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

        #---------Determine the total event weight---------#
        if(not isMuon):
            ele_Weight = myWF.get_ele_scale(lep_pT,lep_eta)

        ph_Weight = myWF.get_photon_scale(gamma_eT,gamma_eta)
        PU_Weight = mytree.PU_Weight
        Event_Weight = norm_factor*PU_Weight*ph_Weight

        if not isMuon:
            Event_Weight = Event_Weight*ele_Weight

        Event_Weight_err = Event_Weight/math.sqrt(N_evts_processed)
        #--------------------------------------------------#

        if not is_Event_selected(mytree.Wmass,mytree.nBjets_25,mytree.lepton_iso,ele_gamma_InvMass):#,lep_pT):#,pi_pT,gamma_eT):
            continue

        for icut1 in xrange(steps_cut1):
            cut1_value = cut1_init + cut1_stepsize*icut1

            if is_ele_gamma:
                if ele_gamma_InvMass > 91. - cut1_value and ele_gamma_InvMass < 91. + cut1_value:
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

            if is_piIso_03:
                if piIso_03 > cut1_value:
                    continue

            if is_piIso_05:
                if piIso_05 > cut1_value:
                    continue

            for icut2 in xrange(steps_cut2):
                cut2_value = cut2_init + cut2_stepsize*icut2
                
                if is_pi_gamma:
                    if gamma_eT < cut2_value:
                        continue
                else:
                    if gamma_eT < 0: #-----a possible condition to use to optimize a single variable-----#
                        continue

                if name_sample == myWF.sig_samplename:
                    cut_Nsig[icut1][icut2] += Event_Weight
                    Nsig_passed[icut1][icut2] += 1
                    #----sum in quadrature of the errors----#
                    cut_Nsig_err[icut1][icut2] += Event_Weight_err*Event_Weight_err
 
                else:
                    cut_Nbkg[icut1][icut2] += Event_Weight
                    Nbkg_passed[icut1][icut2] += 1
                    #----sum in quadrature of the errors----#
                    cut_Nbkg_err[icut1][icut2] += Event_Weight_err*Event_Weight_err 



    for icut1 in xrange(steps_cut1):
        cut1_value = cut1_init + cut1_stepsize*icut1

        for icut2 in xrange(steps_cut2):
            cut2_value = cut2_init + cut2_stepsize*icut2


            sigma_Nexp_sig[icut1][icut2] += Nsig_passed[icut1][icut2]*norm_factor*norm_factor

            sigma_Nexp_bkg[icut1][icut2] += Nbkg_passed[icut1][icut2]*norm_factor*norm_factor

            #print "Nsig passed: ", Nsig_passed[icut1][icut2]
            #print "Nbkg passed: ", Nbkg_passed[icut1][icut2] 
            #print "sigma_Nexp_sig:  ",  sigma_Nexp_sig[icut1][icut2]
            #print "sigma_Nexp_bkg:  ",  sigma_Nexp_bkg[icut1][icut2]

print "Done looping over the events"

##Calculate the significance
signif_list = [[0 for x in range(steps_cut2)] for x in range(steps_cut1)]
err_list = [[0 for x in range(steps_cut2)] for x in range(steps_cut1)]

cut1_x_list = [cut1_init + cut1_stepsize*icut for icut in range(steps_cut1)]
cut2_x_list = [cut2_init + cut2_stepsize*icut for icut in range(steps_cut2)]

signif_max = -1.
cut1_max = -1
cut2_max = -1

for icut1 in xrange(steps_cut1):
    for icut2 in xrange(steps_cut2):

        if cut_Nbkg[icut1][icut2] != 0:
            signif_list[icut1][icut2] = cut_Nsig[icut1][icut2]/math.sqrt(cut_Nbkg[icut1][icut2])
            err_list[icut1][icut2] = math.sqrt((1/cut_Nbkg[icut1][icut2])*sigma_Nexp_sig[icut1][icut2]+math.pow(cut_Nsig[icut1][icut2],2)/(4*math.pow(cut_Nbkg[icut1][icut2],3))* sigma_Nexp_bkg[icut1][icut2])
        else:
            signif_list[icut1][icut2] = 0.

        if signif_list[icut1][icut2] > signif_max:
            signif_max = signif_list[icut1][icut2]
            cut1_max = icut1
            cut2_max = icut2

        print "signif: ", signif_list[icut1][icut2], "pm: ",  err_list[icut1][icut2]

print "The cut1 value is ", cut1_init + cut1_stepsize*cut1_max
print "The cut2 value is ", cut2_init + cut2_stepsize*cut2_max
print "Number of signal events is ", cut_Nsig[cut1_max][cut2_max], " pm ", math.sqrt(sigma_Nexp_sig[cut1_max][cut2_max])
print "Number of background events is ", cut_Nbkg[cut1_max][cut2_max], " pm ", math.sqrt(sigma_Nexp_bkg[cut1_max][cut2_max])
print "Max significance is ", signif_max

##Plot the significance as function of the cuts
cut1_x = np.array(cut1_x_list)
cut2_x = np.array(cut2_x_list)

cut1_y_list = [signif_list[icut][cut2_max] for icut in xrange(steps_cut1)]
cut2_y_list = [signif_list[cut1_max][icut] for icut in xrange(steps_cut2)]

cut1_y_list_err = [err_list[icut][cut2_max] for icut in xrange(steps_cut1)]
cut2_y_list_err = [err_list[cut1_max][icut] for icut in xrange(steps_cut2)]

cut1_y = np.array(cut1_y_list)
cut2_y = np.array(cut2_y_list)

cut1_x_err = np.zeros(steps_cut1)
cut2_x_err = np.zeros(steps_cut2)

cut1_y_err = np.array(cut1_y_list_err)
cut2_y_err = np.array(cut2_y_list_err)

graph_cut1 = ROOT.TGraphErrors(steps_cut1,cut1_x,cut1_y,cut1_x_err,cut1_y_err)
graph_cut2 = ROOT.TGraphErrors(steps_cut2,cut2_x,cut2_y,cut2_x_err,cut2_y_err)

c1 = ROOT.TCanvas("c1","c1",800,500)
c1.cd()
graph_cut1.SetMarkerStyle(20)
graph_cut1.Draw("AP")
graph_cut1.GetYaxis().SetTitleOffset(1.6)

if is_ele_gamma:
    graph_cut1.SetTitle("; #varepsilon (GeV/c^{2}); Significance")
    c1.SaveAs("plots/Z_peak.pdf")
if is_ele_pT:
    graph_cut1.SetTitle("; p_{T}^{e} (GeV/c); Significance")
    c1.SaveAs("plots/e_pT.pdf")
if is_mu_pT:
    graph_cut1.SetTitle("; p_{T}^{#mu} (GeV/c); Significance")
    c1.SaveAs("plots/mu_pT.pdf")
if is_deltaphi_mu_pi:
    graph_cut1.SetTitle("; #Delta#varphi(#mu,#pi); Significance")
    c1.SaveAs("plots/deltaphi_mu_pi.pdf")
if is_deltaphi_ele_pi:
    graph_cut1.SetTitle("; #Delta#varphi(e,#pi); Significance")
    c1.SaveAs("plots/deltaphi_ele_pi.pdf")
if is_pi_gamma:
    graph_cut1.SetTitle("; p_{T}^{#pi} (GeV/c); Significance")
    c1.SaveAs("plots/pi_pT.pdf")
    c2 = ROOT.TCanvas("c2","c2",800,500)
    c2.cd()
    graph_cut2.SetMarkerStyle(20)
    graph_cut2.Draw("AP")
    graph_cut2.GetYaxis().SetTitleOffset(1.6)
    graph_cut2.SetTitle("; E_{T}^{#gamma} (GeV); Significance")
    c2.SaveAs("plots/gamma_eT.pdf")
if is_piIso_03:
    graph_cut1.SetTitle("; Rel.Iso_03; Significance")
    c1.SaveAs("plots/piIso_03.pdf")
if is_piIso_05:
    graph_cut1.SetTitle("; Rel.Iso_05; Significance")
    c1.SaveAs("plots/piIso_05.pdf")
