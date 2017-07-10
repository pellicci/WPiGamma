import ROOT
import os
import math
import numpy as np
from ROOT import TFile, TTree, TBranch

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

from Workflow_Handler import Workflow_Handler
myWF = Workflow_Handler("Signal")

Wmass = np.zeros(1, dtype=float)
isMuon = np.zeros(1, dtype=int)
isSignal = np.zeros(1, dtype=int)
weight = np.zeros(1, dtype=float)
f = TFile('WmassAnalysis/Tree_MC.root','recreate')
t = TTree('minitree','tree with branches')
t.Branch('Wmass',Wmass,'Wmass/D')
t.Branch('isMuon',isMuon,'isMuon/I')
t.Branch('isSignal',isSignal,'isSignal/I')
t.Branch('weight',weight,'weight/D')

##Global constants
MU_MIN_PT    = 25.
ELE_MIN_PT   = 26.
PI_MIN_PT    = 50.
GAMMA_MIN_ET = 50.
N_BJETS_MIN  = 0.
WMASS_MIN    = 50.
WMASS_MAX    = 100.
DELTAPHI_MU_PI_MIN = 1.2
DELTAPHI_ELE_PI_MIN = 1.8

#Normalize to this luminsity, in fb-1
luminosity_norm = 36.46

def select_all_but_one(cutstring):

    selection_bools = dict()
    if ismuon:
        selection_bools["mupt"]      = lep_pt > MU_MIN_PT
        selection_bools["deltaphi_mu_pi"]  = deltaphi_lep_pi > DELTAPHI_MU_PI_MIN
    if not ismuon:
        selection_bools["elept"]     = lep_pt > ELE_MIN_PT
        selection_bools["deltaphi_ele_pi"]  = deltaphi_lep_pi > DELTAPHI_ELE_PI_MIN
    selection_bools["pipt"]      = pi_pt > PI_MIN_PT
    selection_bools["gammaet"]   = gamma_et > GAMMA_MIN_ET
    selection_bools["nBjets"]    = nBjets > N_BJETS_MIN
    selection_bools["Wmass_bottom"]     = wmass > WMASS_MIN
    selection_bools["Wmass_top"]     = wmass < WMASS_MAX
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

for name_sample in samplename_list:

    theSampleName = name_sample

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
        
        if name_sample == "ttbar" and mytree.isttbarlnu:
            continue

        PU_Weight = mytree.PU_Weight
        Event_Weight = norm_factor*PU_Weight
        
        lep_pt  = mytree.lepton_pT
        lep_eta = mytree.lepton_eta
        lep_phi = mytree.lepton_phi
        
        ismuon = mytree.is_muon

        pi_pt = mytree.pi_pT
        pi_eta = mytree.pi_eta
        pi_phi = mytree.pi_phi
        
        gamma_et = mytree.photon_eT

        wmass = mytree.Wmass

        nBjets = mytree.nBjets
        
        deltaeta_lep_pi = lep_eta-pi_eta
        
        deltaphi_lep_pi = math.fabs(lep_phi-pi_phi)
        if deltaphi_lep_pi > 3.14:
            deltaphi_lep_pi = 6.28 - deltaphi_lep_pi

#-------- Filling tree -------------
        if select_all_but_one("all cuts"):
            isMuon[0] = mytree.is_muon
            Wmass[0] = mytree.Wmass
            weight[0] = Event_Weight
            if name_sample == myWF.sig_samplename:
                isSignal[0] = 1
            else:
                isSignal[0] = 0
            t.Fill()

print "Finished runnning over samples!"

f.Write()
print "file written"
f.Close()
