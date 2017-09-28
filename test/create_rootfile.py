import ROOT
import os
import math
import numpy as np
from ROOT import TFile, TTree, TBranch

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

from Workflow_Handler import Workflow_Handler
myWF = Workflow_Handler("Signal","Data",isMedium = True)

isData = False ##---------switch from DATA to MC and vice versa---------##

Wmass = np.zeros(1, dtype=float)
isMuon = np.zeros(1, dtype=int)
isSignal = np.zeros(1, dtype=int)
weight = np.zeros(1, dtype=float)

if not isData:
    f = TFile('WmassAnalysis/Tree_MC.root','recreate')
    t = TTree('minitree','tree with branches')
    t.Branch('Wmass',Wmass,'Wmass/D')
    t.Branch('isMuon',isMuon,'isMuon/I')
    t.Branch('isSignal',isSignal,'isSignal/I')
    t.Branch('weight',weight,'weight/D')

if isData:
    f = TFile('WmassAnalysis/Tree_Data.root','recreate')
    t = TTree('minitree','tree with branches')
    t.Branch('Wmass',Wmass,'Wmass/D')
    t.Branch('isMuon',isMuon,'isMuon/I')

##Global constants
MU_MIN_PT = 27.
ELE_MIN_PT = 29.
PI_MIN_PT = 50.
GAMMA_MIN_ET = 60.
N_BJETS_MIN = 2.
WMASS_MIN = 50.
WMASS_MAX  = 100.
DELTAPHI_MU_PI_MIN = 0.
DELTAPHI_ELE_PI_MIN = 0.
ELE_ISO_MAX = 0.35
ELE_GAMMA_INVMASS_MIN = 88.5
ELE_GAMMA_INVMASS_MAX = 91.5

#Normalize to this luminsity, in fb-1
luminosity_norm = 36.46

def select_all_but_one(cutstring):

    selection_bools = dict()
    if ismuon:
        selection_bools["mupt"]                = lep_pT >= MU_MIN_PT
        selection_bools["deltaphi_mu_pi"]      = deltaphi_lep_pi >= DELTAPHI_MU_PI_MIN
    if not ismuon:
        selection_bools["elept"]               = lep_pT >= ELE_MIN_PT
        selection_bools["deltaphi_ele_pi"]     = deltaphi_lep_pi >= DELTAPHI_ELE_PI_MIN
        selection_bools["h_ele_iso"]           = lep_iso <= ELE_ISO_MAX
        selection_bools["h_ele_gamma_InvMass"] = (ele_gamma_InvMass < ELE_GAMMA_INVMASS_MIN or ele_gamma_InvMass > ELE_GAMMA_INVMASS_MAX)
    selection_bools["pipt"]                    = pi_pT >= PI_MIN_PT
    selection_bools["gammaet"]                 = gamma_eT >= GAMMA_MIN_ET
    selection_bools["nBjets"]                  = nBjets_25 >= N_BJETS_MIN
    selection_bools["Wmass"]                   = (wmass >= WMASS_MIN and wmass <= WMASS_MAX)
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

    mytree = root_file[name_sample].Get("WPiGammaAnalysis/mytree")
 
    print "Processing Sample ", name_sample

    for jentry in xrange(mytree.GetEntriesFast()):
        ientry = mytree.LoadTree( jentry )
        if ientry < 0:
            break
        nb = mytree.GetEntry(jentry )
        if nb <= 0:
            continue

        if isData and not "Data" in name_sample: 
            continue
        
        if not isData and "Data" in name_sample:
            continue
        
        if name_sample == "ttbar" and mytree.isttbarlnu:
            continue

        Event_Weight = 1.
        if not isData:
            norm_factor = Norm_Map[name_sample]*luminosity_norm
            PU_Weight = mytree.PU_Weight
            Event_Weight = norm_factor*PU_Weight

        ismuon = mytree.is_muon

        lep_pT  = mytree.lepton_pT
        lep_eta = mytree.lepton_eta
        lep_phi = mytree.lepton_phi
        lep_iso = mytree.lepton_iso

        pi_pT = mytree.pi_pT
        pi_eta = mytree.pi_eta
        pi_phi = mytree.pi_phi
        pi_E = mytree.pi_energy
        pi_FourMomentum = ROOT.TLorentzVector()
        pi_FourMomentum.SetPtEtaPhiE(pi_pT,pi_eta,pi_phi,pi_E)
            
        gamma_eT = mytree.photon_eT
        gamma_eta = mytree.photon_eta
        gamma_phi = mytree.photon_phi
        gamma_E = mytree.photon_energy
        gamma_FourMomentum = ROOT.TLorentzVector()
        gamma_FourMomentum.SetPtEtaPhiE(gamma_eT,gamma_eta,gamma_phi,gamma_E)
        gamma_iso_ChHad = mytree.photon_iso_ChargedHadron
        gamma_iso_NeuHad = mytree.photon_iso_NeutralHadron
        gamma_iso_Ph = mytree.photon_iso_Photon
        gamma_iso_eArho = mytree.photon_iso_eArho

        W_phi = (pi_FourMomentum + gamma_FourMomentum).Phi()

        wmass = mytree.Wmass

        #if (wmass > 65. and wmass < 90.): continue  #-------Blind window----------

        lep_FourMomentum = ROOT.TLorentzVector()
        lep_FourMomentum.SetPtEtaPhiM(lep_pT,lep_eta,lep_phi,0.)

        if not ismuon:
            ele_gamma_InvMass = (lep_FourMomentum + gamma_FourMomentum).M()
        else:
            mu_gamma_InvMass = (lep_FourMomentum + gamma_FourMomentum).M()
        
        nBjets = mytree.nBjets
        nBjets_25 = mytree.nBjets_25

        deltaeta_lep_pi = lep_eta-pi_eta

        deltaphi_lep_pi = math.fabs(lep_phi-pi_phi)
        if deltaphi_lep_pi > 3.14:
            deltaphi_lep_pi = 6.28 - deltaphi_lep_pi

        deltaphi_lep_W = math.fabs(lep_phi-W_phi)
        if deltaphi_lep_W > 3.14:
            deltaphi_lep_W = 6.28 - deltaphi_lep_W        

#-------- Filling tree -------------
        if select_all_but_one("all cuts"):
            isMuon[0] = mytree.is_muon
            Wmass[0] = mytree.Wmass
            if not isData:
                weight[0] = Event_Weight
                if name_sample == myWF.sig_samplename:
                    isSignal[0] = 1
                else:
                    isSignal[0] = 0
            t.Fill()

print "Finished runnning over samples!"

f.Write()
print "File written"
f.Close()
