import ROOT
import os
import math
import numpy as np
from ROOT import TFile, TTree, TBranch

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)

from Workflow_Handler import Workflow_Handler
myWF = Workflow_Handler("Signal","Data",isMedium = True)

#------------------------------------------------------------------------#

isData = False ##---------switch from DATA to MC and vice versa---------##
split_MC = False #-------if True, MC signal sample is split in two for the training/testing of the BDT

#------------------------------------------------------------------------#

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
ELE_GAMMA_INVMASS_MIN = 90.4
ELE_GAMMA_INVMASS_MAX = 91.6

#Normalize to this luminsity, in fb-1
#luminosity_norm = 36.46
luminosity_norm = 35.86
luminosity_BtoF = 19.72
luminosity_GH   = 16.14

#Variables for the trees
Wmass    = np.zeros(1, dtype=float)
isMuon   = np.zeros(1, dtype=int)
isSignal = np.zeros(1, dtype=int)
weight   = np.zeros(1, dtype=float)

_gamma_eT          = np.zeros(1, dtype=float)
_pi_pT             = np.zeros(1, dtype=float)
_lep_pT            = np.zeros(1, dtype=float)
_lep_iso           = np.zeros(1, dtype=float)
_nBjets            = np.zeros(1, dtype=int)
_nBjets_25         = np.zeros(1, dtype=int)
_deltaphi_lep_pi   = np.zeros(1, dtype=float)
_piRelIso_03       = np.zeros(1, dtype=float)
_piRelIso_05       = np.zeros(1, dtype=float)
_ele_gamma_InvMass = np.zeros(1, dtype=float)
_met               = np.zeros(1, dtype=float)

if not isData:
    f = TFile('WmassAnalysis/Tree_MC.root','recreate')
    t = TTree('minitree','tree with branches')
    t.Branch('Wmass',Wmass,'Wmass/D')
    t.Branch('isMuon',isMuon,'isMuon/I')
    t.Branch('isSignal',isSignal,'isSignal/I')
    t.Branch('weight',weight,'weight/D')

    if split_MC:
        fMVA_signal_mu_training = TFile('MVA/Tree_MC_Signal_mu_training.root','recreate')
        tMVA_signal_mu_training = TTree('minitree_signal_mu_training','tree with branches')
        tMVA_signal_mu_training.Branch('isMuon',isMuon,'isMuon/I')
        tMVA_signal_mu_training.Branch('weight',weight,'weight/D')
        tMVA_signal_mu_training.Branch('gamma_eT',_gamma_eT,'gamma_eT/D')
        tMVA_signal_mu_training.Branch('pi_pT',_pi_pT,'pi_pT/D')
        tMVA_signal_mu_training.Branch('lep_pT',_lep_pT,'lep_pT/D')
        tMVA_signal_mu_training.Branch('lep_iso',_lep_iso,'lep_iso/D')
        tMVA_signal_mu_training.Branch('nBjets',_nBjets,'nBjets/I')
        tMVA_signal_mu_training.Branch('nBjets_25',_nBjets_25,'nBjets_25/I')
        tMVA_signal_mu_training.Branch('deltaphi_lep_pi',_deltaphi_lep_pi,'deltaphi_lep_pi/D')
        tMVA_signal_mu_training.Branch('piRelIso_03',_piRelIso_03,'piRelIso_03/D')
        tMVA_signal_mu_training.Branch('piRelIso_05',_piRelIso_05,'piRelIso_05/D')
        tMVA_signal_mu_training.Branch('MET',_met,'met/D')
        
        fMVA_signal_ele_training = TFile('MVA/Tree_MC_Signal_ele_training.root','recreate')
        tMVA_signal_ele_training = TTree('minitree_signal_ele_training','tree with branches')
        tMVA_signal_ele_training.Branch('isMuon',isMuon,'isMuon/I')
        tMVA_signal_ele_training.Branch('weight',weight,'weight/D')
        tMVA_signal_ele_training.Branch('gamma_eT',_gamma_eT,'gamma_eT/D')
        tMVA_signal_ele_training.Branch('pi_pT',_pi_pT,'pi_pT/D')
        tMVA_signal_ele_training.Branch('lep_pT',_lep_pT,'lep_pT/D')
        tMVA_signal_ele_training.Branch('lep_iso',_lep_iso,'lep_iso/D')
        tMVA_signal_ele_training.Branch('nBjets',_nBjets,'nBjets/I')
        tMVA_signal_ele_training.Branch('nBjets_25',_nBjets_25,'nBjets_25/I')
        tMVA_signal_ele_training.Branch('deltaphi_lep_pi',_deltaphi_lep_pi,'deltaphi_lep_pi/D')
        tMVA_signal_ele_training.Branch('piRelIso_03',_piRelIso_03,'piRelIso_03/D')
        tMVA_signal_ele_training.Branch('piRelIso_05',_piRelIso_05,'piRelIso_05/D')
        tMVA_signal_ele_training.Branch('MET',_met,'met/D')
        
        fMVA_signal_mu_test = TFile('MVA/Tree_MC_Signal_mu_test.root','recreate')
        tMVA_signal_mu_test = TTree('minitree_signal_mu_test','tree with branches')
        tMVA_signal_mu_test.Branch('isMuon',isMuon,'isMuon/I')
        tMVA_signal_mu_test.Branch('weight',weight,'weight/D')
        tMVA_signal_mu_test.Branch('gamma_eT',_gamma_eT,'gamma_eT/D')
        tMVA_signal_mu_test.Branch('pi_pT',_pi_pT,'pi_pT/D')
        tMVA_signal_mu_test.Branch('lep_pT',_lep_pT,'lep_pT/D')
        tMVA_signal_mu_test.Branch('lep_iso',_lep_iso,'lep_iso/D')
        tMVA_signal_mu_test.Branch('nBjets',_nBjets,'nBjets/I')
        tMVA_signal_mu_test.Branch('nBjets_25',_nBjets_25,'nBjets_25/I')
        tMVA_signal_mu_test.Branch('deltaphi_lep_pi',_deltaphi_lep_pi,'deltaphi_lep_pi/D')
        tMVA_signal_mu_test.Branch('piRelIso_03',_piRelIso_03,'piRelIso_03/D')
        tMVA_signal_mu_test.Branch('piRelIso_05',_piRelIso_05,'piRelIso_05/D')
        tMVA_signal_mu_test.Branch('MET',_met,'met/D')
        
        fMVA_signal_ele_test = TFile('MVA/Tree_MC_Signal_ele_test.root','recreate')
        tMVA_signal_ele_test = TTree('minitree_signal_ele_test','tree with branches')
        tMVA_signal_ele_test.Branch('isMuon',isMuon,'isMuon/I')
        tMVA_signal_ele_test.Branch('weight',weight,'weight/D')
        tMVA_signal_ele_test.Branch('gamma_eT',_gamma_eT,'gamma_eT/D')
        tMVA_signal_ele_test.Branch('pi_pT',_pi_pT,'pi_pT/D')
        tMVA_signal_ele_test.Branch('lep_pT',_lep_pT,'lep_pT/D')
        tMVA_signal_ele_test.Branch('lep_iso',_lep_iso,'lep_iso/D')
        tMVA_signal_ele_test.Branch('nBjets',_nBjets,'nBjets/I')
        tMVA_signal_ele_test.Branch('nBjets_25',_nBjets_25,'nBjets_25/I')
        tMVA_signal_ele_test.Branch('deltaphi_lep_pi',_deltaphi_lep_pi,'deltaphi_lep_pi/D')
        tMVA_signal_ele_test.Branch('piRelIso_03',_piRelIso_03,'piRelIso_03/D')
        tMVA_signal_ele_test.Branch('piRelIso_05',_piRelIso_05,'piRelIso_05/D')
        tMVA_signal_ele_test.Branch('MET',_met,'met/D')

    else:
        fMVA_signal_mu = TFile('MVA/Tree_MC_Signal_mu.root','recreate')
        tMVA_signal_mu = TTree('minitree_signal_mu','tree with branches')
        tMVA_signal_mu.Branch('isMuon',isMuon,'isMuon/I')
        tMVA_signal_mu.Branch('weight',weight,'weight/D')
        tMVA_signal_mu.Branch('gamma_eT',_gamma_eT,'gamma_eT/D')
        tMVA_signal_mu.Branch('pi_pT',_pi_pT,'pi_pT/D')
        tMVA_signal_mu.Branch('lep_pT',_lep_pT,'lep_pT/D')
        tMVA_signal_mu.Branch('lep_iso',_lep_iso,'lep_iso/D')
        tMVA_signal_mu.Branch('nBjets',_nBjets,'nBjets/I')
        tMVA_signal_mu.Branch('nBjets_25',_nBjets_25,'nBjets_25/I')
        tMVA_signal_mu.Branch('deltaphi_lep_pi',_deltaphi_lep_pi,'deltaphi_lep_pi/D')
        tMVA_signal_mu.Branch('piRelIso_03',_piRelIso_03,'piRelIso_03/D')
        tMVA_signal_mu.Branch('piRelIso_05',_piRelIso_05,'piRelIso_05/D')
        tMVA_signal_mu.Branch('MET',_met,'met/D')
        
        fMVA_signal_ele = TFile('MVA/Tree_MC_Signal_ele.root','recreate')
        tMVA_signal_ele = TTree('minitree_signal_ele','tree with branches')
        tMVA_signal_ele.Branch('isMuon',isMuon,'isMuon/I')
        tMVA_signal_ele.Branch('weight',weight,'weight/D')
        tMVA_signal_ele.Branch('gamma_eT',_gamma_eT,'gamma_eT/D')
        tMVA_signal_ele.Branch('pi_pT',_pi_pT,'pi_pT/D')
        tMVA_signal_ele.Branch('lep_pT',_lep_pT,'lep_pT/D')
        tMVA_signal_ele.Branch('lep_iso',_lep_iso,'lep_iso/D')
        tMVA_signal_ele.Branch('nBjets',_nBjets,'nBjets/I')
        tMVA_signal_ele.Branch('nBjets_25',_nBjets_25,'nBjets_25/I')
        tMVA_signal_ele.Branch('deltaphi_lep_pi',_deltaphi_lep_pi,'deltaphi_lep_pi/D')
        tMVA_signal_ele.Branch('piRelIso_03',_piRelIso_03,'piRelIso_03/D')
        tMVA_signal_ele.Branch('piRelIso_05',_piRelIso_05,'piRelIso_05/D')
        tMVA_signal_ele.Branch('MET',_met,'met/D')
    
    fMVA_background_mu = TFile('MVA/Tree_MC_Background_mu.root','recreate')
    tMVA_background_mu = TTree('minitree_background_mu','tree with branches')
    tMVA_background_mu.Branch('isMuon',isMuon,'isMuon/I')
    tMVA_background_mu.Branch('weight',weight,'weight/D')
    tMVA_background_mu.Branch('gamma_eT',_gamma_eT,'gamma_eT/D')
    tMVA_background_mu.Branch('pi_pT',_pi_pT,'pi_pT/D')
    tMVA_background_mu.Branch('lep_pT',_lep_pT,'lep_pT/D')
    tMVA_background_mu.Branch('lep_iso',_lep_iso,'lep_iso/D')
    tMVA_background_mu.Branch('nBjets',_nBjets,'nBjets/I')
    tMVA_background_mu.Branch('nBjets_25',_nBjets_25,'nBjets_25/I')
    tMVA_background_mu.Branch('deltaphi_lep_pi',_deltaphi_lep_pi,'deltaphi_lep_pi/D')
    tMVA_background_mu.Branch('piRelIso_03',_piRelIso_03,'piRelIso_03/D')
    tMVA_background_mu.Branch('piRelIso_05',_piRelIso_05,'piRelIso_05/D')
    tMVA_background_mu.Branch('MET',_met,'met/D')

    fMVA_background_ele = TFile('MVA/Tree_MC_Background_ele.root','recreate')
    tMVA_background_ele = TTree('minitree_background_ele','tree with branches')
    tMVA_background_ele.Branch('isMuon',isMuon,'isMuon/I')
    tMVA_background_ele.Branch('weight',weight,'weight/D')
    tMVA_background_ele.Branch('gamma_eT',_gamma_eT,'gamma_eT/D')
    tMVA_background_ele.Branch('pi_pT',_pi_pT,'pi_pT/D')
    tMVA_background_ele.Branch('lep_pT',_lep_pT,'lep_pT/D')
    tMVA_background_ele.Branch('lep_iso',_lep_iso,'lep_iso/D')
    tMVA_background_ele.Branch('nBjets',_nBjets,'nBjets/I')
    tMVA_background_ele.Branch('nBjets_25',_nBjets_25,'nBjets_25/I')
    tMVA_background_ele.Branch('deltaphi_lep_pi',_deltaphi_lep_pi,'deltaphi_lep_pi/D')
    tMVA_background_ele.Branch('piRelIso_03',_piRelIso_03,'piRelIso_03/D')
    tMVA_background_ele.Branch('piRelIso_05',_piRelIso_05,'piRelIso_05/D')
    tMVA_background_ele.Branch('MET',_met,'met/D')


if isData:

    f = TFile('WmassAnalysis/Tree_Data.root','recreate')
    t = TTree('minitree','tree with branches')
    t.Branch('Wmass',Wmass,'Wmass/D')
    t.Branch('isMuon',isMuon,'isMuon/I')

    fMVA_background_mu_DATA = TFile('MVA/Tree_MC_Background_mu_DATA.root','recreate')
    tMVA_background_mu_DATA = TTree('minitree_background_mu_DATA','tree with branches')
    tMVA_background_mu_DATA.Branch('isMuon',isMuon,'isMuon/I')
    tMVA_background_mu_DATA.Branch('weight',weight,'weight/D')
    tMVA_background_mu_DATA.Branch('gamma_eT',_gamma_eT,'gamma_eT/D')
    tMVA_background_mu_DATA.Branch('pi_pT',_pi_pT,'pi_pT/D')
    tMVA_background_mu_DATA.Branch('lep_pT',_lep_pT,'lep_pT/D')
    tMVA_background_mu_DATA.Branch('lep_iso',_lep_iso,'lep_iso/D')
    tMVA_background_mu_DATA.Branch('nBjets',_nBjets,'nBjets/I')
    tMVA_background_mu_DATA.Branch('nBjets_25',_nBjets_25,'nBjets_25/I')
    tMVA_background_mu_DATA.Branch('deltaphi_lep_pi',_deltaphi_lep_pi,'deltaphi_lep_pi/D')
    tMVA_background_mu_DATA.Branch('piRelIso_03',_piRelIso_03,'piRelIso_03/D')
    tMVA_background_mu_DATA.Branch('piRelIso_05',_piRelIso_05,'piRelIso_05/D')
    tMVA_background_mu_DATA.Branch('MET',_met,'met/D')

    fMVA_background_ele_DATA = TFile('MVA/Tree_MC_Background_ele_DATA.root','recreate')
    tMVA_background_ele_DATA = TTree('minitree_background_ele_DATA','tree with branches')
    tMVA_background_ele_DATA.Branch('isMuon',isMuon,'isMuon/I')
    tMVA_background_ele_DATA.Branch('weight',weight,'weight/D')
    tMVA_background_ele_DATA.Branch('gamma_eT',_gamma_eT,'gamma_eT/D')
    tMVA_background_ele_DATA.Branch('pi_pT',_pi_pT,'pi_pT/D')
    tMVA_background_ele_DATA.Branch('lep_pT',_lep_pT,'lep_pT/D')
    tMVA_background_ele_DATA.Branch('lep_iso',_lep_iso,'lep_iso/D')
    tMVA_background_ele_DATA.Branch('nBjets',_nBjets,'nBjets/I')
    tMVA_background_ele_DATA.Branch('nBjets_25',_nBjets_25,'nBjets_25/I')
    tMVA_background_ele_DATA.Branch('deltaphi_lep_pi',_deltaphi_lep_pi,'deltaphi_lep_pi/D')
    tMVA_background_ele_DATA.Branch('piRelIso_03',_piRelIso_03,'piRelIso_03/D')
    tMVA_background_ele_DATA.Branch('piRelIso_05',_piRelIso_05,'piRelIso_05/D')
    tMVA_background_ele_DATA.Branch('MET',_met,'met/D')


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

SignalTree_entries = 0
entry_index = 0

for name_sample in samplename_list:

    theSampleName = name_sample

    #i_event = 0
    is_ttbarlnu = 0

    mytree = root_file[name_sample].Get("WPiGammaAnalysis/mytree")
 
    print "Processing Sample ", name_sample

    #nEvts = mytree.GetEntriesFast() #Get the number of events per sample

    for jentry in xrange(mytree.GetEntriesFast()):
        #i_event += 1
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
            is_ttbarlnu += 1
            continue

        if name_sample == "Signal":
            SignalTree_entries = mytree.GetEntriesFast()
            entry_index += 1

        ismuon = mytree.is_muon

        #run_number = mytree.run_number
        isSingleMuTrigger_24 = mytree.isSingleMuTrigger_24
        isSingleMuTrigger_50 = mytree.isSingleMuTrigger_50

        lep_pT  = mytree.lepton_pT
        lep_eta = mytree.lepton_eta
        lep_phi = mytree.lepton_phi
        lep_iso = mytree.lepton_iso
        lep_FourMomentum = ROOT.TLorentzVector()
        lep_FourMomentum.SetPtEtaPhiM(lep_pT,lep_eta,lep_phi,0.)

        pi_pT = mytree.pi_pT
        pi_eta = mytree.pi_eta
        pi_phi = mytree.pi_phi
        pi_E = mytree.pi_energy
        pi_FourMomentum = ROOT.TLorentzVector()
        pi_FourMomentum.SetPtEtaPhiE(pi_pT,pi_eta,pi_phi,pi_E)
        piRelIso_03 = mytree.sum_pT_03/pi_pT
        piRelIso_05 = mytree.sum_pT_05/pi_pT
            
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

        met = mytree.met_pT

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

        #---------Retrieve the BDT output----------#

        BDT_out = myWF.get_BDT_output(pi_pT,gamma_eT,nBjets_25,deltaphi_lep_pi,lep_pT,piRelIso_05,isMuon)    

        #--------Determining the event weight--------#

        if ismuon:
            mu_weight_BtoF = myWF.get_muon_scale_BtoF(lep_pT,lep_eta,isSingleMuTrigger_24,isSingleMuTrigger_50)
            mu_weight_GH   = myWF.get_muon_scale_GH(lep_pT,lep_eta,isSingleMuTrigger_24,isSingleMuTrigger_50)
            mu_weight_tot  = mu_weight_BtoF*luminosity_BtoF/luminosity_norm + mu_weight_GH*luminosity_GH/luminosity_norm
        else:
            ele_weight = myWF.get_ele_scale(lep_pT,lep_eta)

        ph_weight = myWF.get_photon_scale(gamma_eT,gamma_eta)

        Event_Weight = 1.

        if not isData:
            norm_factor = Norm_Map[name_sample]*luminosity_norm
            PU_Weight = mytree.PU_Weight
            Event_Weight = norm_factor*PU_Weight*ph_weight
            if ismuon:
                Event_Weight = Event_Weight*mu_weight_tot
            else:
                Event_Weight = Event_Weight*ele_weight
  

        #-------- Filling mass tree -------------
        #if select_all_but_one("all cuts"):
        if (ismuon and BDT_out >= 0.14) or (not ismuon and BDT_out >= 0.15):
            isMuon[0] = ismuon
            Wmass[0] = wmass
            if not isData:
                weight[0] = Event_Weight
                if name_sample == myWF.sig_samplename :
                    isSignal[0] = 1
                else :
                    isSignal[0] = 0
            t.Fill()

        #------- Filling MVA tree ------------

        if not isData:
            
            _gamma_eT[0]        = gamma_eT
            _pi_pT[0]           = pi_pT
            _lep_pT[0]          = lep_pT
            _lep_iso[0]         = lep_iso
            _nBjets[0]          = nBjets
            _nBjets_25[0]       = nBjets_25
            weight[0]           = Event_Weight
            _deltaphi_lep_pi[0] = deltaphi_lep_pi
            isMuon[0]           = ismuon
            _piRelIso_03[0]     = piRelIso_03
            _piRelIso_05[0]     = piRelIso_05
            _met[0]             = met
            #if not ismuon:
            #    _ele_gamma_InvMass[0] = math.fabs(91.-ele_gamma_InvMass)

            if name_sample == myWF.sig_samplename and ismuon and split_MC and entry_index <= SignalTree_entries/2:
                tMVA_signal_mu_training.Fill()
            if name_sample == myWF.sig_samplename and not ismuon and split_MC and entry_index <= SignalTree_entries/2:
                tMVA_signal_ele_training.Fill()
            if name_sample == myWF.sig_samplename and ismuon and split_MC and entry_index > SignalTree_entries/2:
                tMVA_signal_mu_test.Fill()
            if name_sample == myWF.sig_samplename and not ismuon and split_MC and entry_index > SignalTree_entries/2:
                tMVA_signal_ele_test.Fill()
            if name_sample == myWF.sig_samplename and ismuon and not split_MC:
                tMVA_signal_mu.Fill()
            if name_sample == myWF.sig_samplename and not ismuon and not split_MC:
                tMVA_signal_ele.Fill()
            if not name_sample == myWF.sig_samplename and ismuon:
                tMVA_background_mu.Fill()
            if not name_sample == myWF.sig_samplename and not ismuon:
                tMVA_background_ele.Fill()

        if isData and (wmass < 65. or wmass > 90.):

            _gamma_eT[0]        = gamma_eT
            _pi_pT[0]           = pi_pT
            _lep_pT[0]          = lep_pT
            _lep_iso[0]         = lep_iso
            _nBjets[0]          = nBjets
            _nBjets_25[0]       = nBjets_25
            weight[0]           = 1
            _deltaphi_lep_pi[0] = deltaphi_lep_pi
            isMuon[0]           = ismuon
            _piRelIso_03[0]     = piRelIso_03
            _piRelIso_05[0]     = piRelIso_05
            _met[0]             = met
            
        if name_sample == "Data" and ismuon:
            tMVA_background_mu_DATA.Fill()
        if name_sample == "Data" and not ismuon:
            tMVA_background_ele_DATA.Fill()


print "Finished runnning over samples!"

f.cd()
t.Write()
f.Close()

if not isData :

    if split_MC:
        fMVA_signal_mu_training.cd()
        tMVA_signal_mu_training.Write()
        fMVA_signal_mu_training.Close()
        
        fMVA_signal_ele_training.cd()
        tMVA_signal_ele_training.Write()
        fMVA_signal_ele_training.Close()
        
        fMVA_signal_mu_test.cd()
        tMVA_signal_mu_test.Write()
        fMVA_signal_mu_test.Close()
        
        fMVA_signal_ele_test.cd()
        tMVA_signal_ele_test.Write()
        fMVA_signal_ele_test.Close()

    else:
        fMVA_signal_mu.cd()
        tMVA_signal_mu.Write()
        fMVA_signal_mu.Close()
        
        fMVA_signal_ele.cd()
        tMVA_signal_ele.Write()
        fMVA_signal_ele.Close()

    fMVA_background_mu.cd()
    tMVA_background_mu.Write()
    fMVA_background_mu.Close()

    fMVA_background_ele.cd()
    tMVA_background_ele.Write()
    fMVA_background_ele.Close()

if isData:
    
    fMVA_background_mu_DATA.cd()
    tMVA_background_mu_DATA.Write()
    fMVA_background_mu_DATA.Close()

    fMVA_background_ele_DATA.cd()
    tMVA_background_ele_DATA.Write()
    fMVA_background_ele_DATA.Close()

print "File written"

