import ROOT
import math
import numpy as np
import sys

##normalize to the 2016 lumi
luminosity_norm = 36.46

from Workflow_Handler import Workflow_Handler
myWF = Workflow_Handler("Signal")

def is_Event_selected(nBjets,Wmass,pi_pt,gamma_et):
    """Save events according to some basic selection criteria"""
    bjet_cut = nBjets > 0.

    mass_cut = Wmass > 20.

    pi_pt_cut = pi_pt > 60.

    gamma_et_cut = gamma_et > 40.

    return bjet_cut and mass_cut and pi_pt_cut and gamma_et_cut

##Here starts the program
Norm_Map = myWF.get_normalizations_map()

steps_cut1 = 10
cut1_init = 0.
cut1_stepsize = 0.3

#steps_cut2 = 10
#cut2_init = 20.
#cut2_stepsize = 10.

cut_Nbkg = [0 for x in range(steps_cut1)]
cut_Nsig = [0 for x in range(steps_cut1)]

##collect all the root files
samplename_list = myWF.get_samples_names()
root_file = myWF.get_root_files()

##loop on the samples and on the cuts and calculate N_bkg and N_sig
for name_sample in samplename_list:

    norm_factor = Norm_Map[name_sample]*luminosity_norm
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

        if not is_Event_selected(mytree.nBjets, mytree.Wmass, mytree.pi_pT, mytree.photon_eT):
            continue
        
        if not mytree.is_muon: #to be set to 'if mytree.is_muon' if you want to find the cut on electrons
            continue

        deltaphi = math.fabs(mytree.lepton_phi-mytree.pi_phi)
        if deltaphi > 3.14:
            deltaphi = 6.28-deltaphi

        for icut1 in xrange(steps_cut1):
            cut1_value = cut1_init + cut1_stepsize*icut1

            if deltaphi < cut1_value:
                continue

            #for icut2 in xrange(steps_cut2):
                #cut2_value = cut2_init + cut2_stepsize*icut2

                #if mytree.photon_eT < cut2_value:
                    #continue

            if name_sample == myWF.sig_samplename:
                cut_Nsig[icut1] += norm_factor
            else:
                cut_Nbkg[icut1] += norm_factor

print "Done looping over the events"

##Calculate the significance
signif_list = [0 for x in range(steps_cut1)]

cut1_x_list = [cut1_init + cut1_stepsize*icut for icut in range(steps_cut1)]
#cut2_x_list = [cut2_init + cut2_stepsize*icut for icut in range(steps_cut2)]

signif_max = -1.
cut1_max = -1
#cut2_max = -1

for icut1 in xrange(steps_cut1):
    #for icut2 in xrange(steps_cut2):

    if cut_Nbkg[icut1] != 0:
        signif_list[icut1] = cut_Nsig[icut1]/math.sqrt(cut_Nbkg[icut1])
    else:
        signif_list[icut1] = 0.

    if signif_list[icut1] > signif_max:
        signif_max = signif_list[icut1]
        cut1_max = icut1
        #cut2_max = icut2

print "The cut1 value is ", cut1_init + cut1_stepsize*cut1_max
#print "The cut2 value is ", cut2_init + cut2_stepsize*cut2_max
print "Number of signal events is ", cut_Nsig[cut1_max]
print "Number of background events is ", cut_Nbkg[cut1_max]
print "Max significance is ", signif_max

##Plot the significance as function of the cuts
cut1_x = np.array(cut1_x_list)
#cut2_x = np.array(cut2_x_list)

cut1_y_list = [signif_list[icut] for icut in xrange(steps_cut1)] #There was cut2_max before
#cut2_y_list = [signif_list[cut1_max][icut] for icut in xrange(steps_cut2)]

cut1_y = np.array(cut1_y_list)
#cut2_y = np.array(cut2_y_list)

graph_cut1 = ROOT.TGraph(steps_cut1,cut1_x,cut1_y)
#graph_cut2 = ROOT.TGraph(steps_cut2,cut2_x,cut2_y)

c1 = ROOT.TCanvas("c1","c1")
c1.cd()
graph_cut1.Draw("A*")
graph_cut1.SetTitle("Significance vs deltaphi_mu_pi; deltaphi mu-pi; Significance")

c1.SaveAs("plots/deltaphi_mu_pi_signif.png")

#c2 = ROOT.TCanvas("c2","c2")
#c2.cd()
#graph_cut2.Draw("A*")
#graph_cut2.SetTitle("Significance vs e_{T} of gamma; e_{T} of gamma; Significance")
#c2.SaveAs("plots/cut2_signif.png")