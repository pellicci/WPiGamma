import ROOT
import os
import math
import numpy as np
import copy
from array import array

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

from Workflow_Handler import Workflow_Handler
myWF = Workflow_Handler("Signal","Data",isBDT_with_Wmass = False)

##Global constants
MU_MIN_PT = 27.
ELE_MIN_PT = 29.
PI_MIN_PT = 50.
GAMMA_MIN_ET = 60.
N_BJETS_MIN = 2.
WMASS_MIN = 50.
WMASS_MAX  = 100.
WMASS_MIN_1 = 65.
WMASS_MAX_1 = 90.
DELTAPHI_MU_PI_MIN = 0.
DELTAPHI_ELE_PI_MIN = 0.
ELE_ISO_MAX = 0.35
ELE_GAMMA_INVMASS_MIN = 88.5
ELE_GAMMA_INVMASS_MAX = 91.5

random_mu_SF  = False #------if True, muon scale factors are sampled from a Gaussian
random_ele_SF = False #------if True, electron scale factors are sampled from a Gaussian
random_ph_SF  = False #------if True, photon scale factors are sampled from a Gaussian

#Normalize to this luminsity, in fb-1
#luminosity_norm = 36.46
luminosity_norm = 35.86
luminosity_BtoF = 19.72
luminosity_GH   = 16.14

#Make signal histos larger
signal_magnify = 10000

output_dir = "plots"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#Here's the list of histos to plot
list_histos = ["h_mupt", "h_elept", "h_pipt", "h_gammaet", "h_Wmass", "h_nBjets", "h_mueta", "h_eleeta","h_deltaphi_mu_pi","h_deltaphi_ele_pi","h_deltaphi_mu_W","h_deltaphi_ele_W","h_deltaeta_mu_pi","h_deltaeta_ele_pi","h_Wmass_flag_mu","h_Wmass_flag_ele","h_mu_iso","h_ele_iso","h_gamma_iso_ChHad","h_gamma_iso_NeuHad","h_gamma_iso_Ph","h_gamma_iso_eArho","h_ele_gamma_InvMass","h_mu_gamma_InvMass","h_nBjets_25","h_evts_Bjets","h_piIso_03_mu","h_piIso_05_mu","h_piRelIso_05_mu_ch","h_piRelIso_05_mu","h_piIso_03_ele","h_piIso_05_ele","h_piRelIso_05_ele_ch","h_piRelIso_05_ele","h_mu_pi_InvMass","h_met","h_Wmass_mu_plus","h_Wmass_ele_plus","h_Wmass_mu_minus","h_Wmass_ele_minus","h_pieta","h_gammaeta","h_piRelIso_05_mu_ch_AfterCut","h_piRelIso_05_ele_ch_AfterCut","h_mupt_sig","h_elept_sig","h_pipt_sig","h_gammaet_sig","h_mueta_sig","h_eleeta_sig","h_pieta_sig","h_gammaeta_sig"]

h_events_sig_mu  = ROOT.TH1F("h_events_sig_mu","Progressive event loss (#mu)",5,0,5)
h_events_sig_ele = ROOT.TH1F("h_events_sig_ele","Progressive event loss (e)",7,0,7)
h_events_bkg_mu  = ROOT.TH1F("h_events_bkg_mu","Progressive event loss (#mu)",5,0,5)
h_events_bkg_ele = ROOT.TH1F("h_events_bkg_ele","Progressive event loss (e)",7,0,7)

Wmass_mu        = ROOT.TH1F("Wmass_mu","Wmass mu",10,40,100)
Wmass_ele       = ROOT.TH1F("Wmass_ele","Wmass ele",10,40,100)
Wmass_mu_plus   = ROOT.TH1F("Wmass_mu_plus","Wmass mu plus",10,40,100)
Wmass_mu_minus  = ROOT.TH1F("Wmass_mu_minus","Wmass mu minus",10,40,100)
Wmass_ele_plus  = ROOT.TH1F("Wmass_ele_plus","Wmass ele plus",10,40,100)
Wmass_ele_minus = ROOT.TH1F("Wmass_ele_minus","Wmass ele minus",10,40,100)

Wmass_mu.Sumw2()
Wmass_ele.Sumw2()
Wmass_mu_plus.Sumw2()
Wmass_mu_minus.Sumw2()
Wmass_ele_plus.Sumw2()
Wmass_ele_minus.Sumw2()

#Color mask must have the same number of entries as non-QCD backgrounds
colors_mask = [26,400,840,616,860,432,880,900,800,416,885,910,200,630,420,608,960,ROOT.kGreen+3,ROOT.kOrange+8]

def select_all_but_one(cutstring):

    selection_bools = dict()
    if isMuon:
        selection_bools["h_mupt"]             = lep_pT >= MU_MIN_PT
        selection_bools["h_deltaphi_mu_pi"]   = deltaphi_lep_pi >= DELTAPHI_MU_PI_MIN
    if not isMuon:
        selection_bools["h_elept"]            = lep_pT >= ELE_MIN_PT
        selection_bools["h_deltaphi_ele_pi"]  = deltaphi_lep_pi >= DELTAPHI_ELE_PI_MIN
        selection_bools["h_ele_iso"]          = lep_iso <= ELE_ISO_MAX
        selection_bools["h_ele_gamma_InvMass"]= (ele_gamma_InvMass < ELE_GAMMA_INVMASS_MIN or ele_gamma_InvMass > ELE_GAMMA_INVMASS_MAX)
    selection_bools["h_pipt"]                 = pi_pT >= PI_MIN_PT
    selection_bools["h_gammaet"]              = gamma_eT >= GAMMA_MIN_ET
    selection_bools["h_nBjets"]               = nBjets_25 >= N_BJETS_MIN
    selection_bools["h_Wmass"]                = (Wmass >= WMASS_MIN and Wmass <= WMASS_MAX)
    #selection_bools["h_Wmass"]                = (Wmass >= WMASS_MIN and Wmass <= WMASS_MIN_1) or (Wmass >= WMASS_MAX_1 and Wmass <= WMASS_MAX)
    #if not "Data" in name_sample:
        #selection_bools["h_Wmass"]            = Wmass >= WMASS_MAX or Wmass <= WMASS_MIN
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

    h_base[theSampleName+list_histos[0]]  = ROOT.TH1F(theSampleName+list_histos[0], "p_{T} of the muon", 15, 20, 100.)
    h_base[theSampleName+list_histos[1]]  = ROOT.TH1F(theSampleName+list_histos[1], "p_{T} of the electron", 15, 20, 100.)
    h_base[theSampleName+list_histos[2]]  = ROOT.TH1F(theSampleName+list_histos[2], "p_{T} of the pion", 15, 20, 100.)
    h_base[theSampleName+list_histos[3]]  = ROOT.TH1F(theSampleName+list_histos[3], "E_{T} of the gamma", 15, 20, 100.)
    h_base[theSampleName+list_histos[4]]  = ROOT.TH1F(theSampleName+list_histos[4], "W mass", 15, 50, 100)
    h_base[theSampleName+list_histos[5]]  = ROOT.TH1F(theSampleName+list_histos[5], "n Bjets", 6, 0, 6.)
    h_base[theSampleName+list_histos[6]]  = ROOT.TH1F(theSampleName+list_histos[6], "eta of the muon", 20, -3, 3)
    h_base[theSampleName+list_histos[7]]  = ROOT.TH1F(theSampleName+list_histos[7], "eta of the electron", 20, -3, 3)
    h_base[theSampleName+list_histos[8]]  = ROOT.TH1F(theSampleName+list_histos[8], "deltaphi mu-pi", 10, 0, 3.14)
    h_base[theSampleName+list_histos[9]]  = ROOT.TH1F(theSampleName+list_histos[9], "deltaphi ele-pi", 10, 0, 3.14)
    h_base[theSampleName+list_histos[10]] = ROOT.TH1F(theSampleName+list_histos[10], "deltaphi mu-W", 10, 0, 3.14)
    h_base[theSampleName+list_histos[11]] = ROOT.TH1F(theSampleName+list_histos[11], "deltaphi ele-W", 10, 0, 3.14)
    h_base[theSampleName+list_histos[12]] = ROOT.TH1F(theSampleName+list_histos[12], "deltaeta mu-pi", 20, -5, 5)
    h_base[theSampleName+list_histos[13]] = ROOT.TH1F(theSampleName+list_histos[13], "deltaeta ele-pi", 20, -5, 5)
    h_base[theSampleName+list_histos[14]] = ROOT.TH1F(theSampleName+list_histos[14], "W mass if flag mu", 15, 50, 100)
    h_base[theSampleName+list_histos[15]] = ROOT.TH1F(theSampleName+list_histos[15], "W mass if flag ele", 15, 50, 100)
    h_base[theSampleName+list_histos[16]] = ROOT.TH1F(theSampleName+list_histos[16], "muon isolation", 20, 0, 0.3)
    h_base[theSampleName+list_histos[17]] = ROOT.TH1F(theSampleName+list_histos[17], "electron isolation", 20, 0, 0.3)
    h_base[theSampleName+list_histos[18]] = ROOT.TH1F(theSampleName+list_histos[18], "Photon isolation - ChargedHadron", 50, 0, 5)
    h_base[theSampleName+list_histos[19]] = ROOT.TH1F(theSampleName+list_histos[19], "Photon isolation - NeutralHadron", 50, 0, 5)
    h_base[theSampleName+list_histos[20]] = ROOT.TH1F(theSampleName+list_histos[20], "Photon isolation - Photon", 50, 0, 5)
    h_base[theSampleName+list_histos[21]] = ROOT.TH1F(theSampleName+list_histos[21], "Photon isolation - eArho", 50, 0, 5)
    h_base[theSampleName+list_histos[22]] = ROOT.TH1F(theSampleName+list_histos[22], "ele-gamma InvMass", 50, 0, 300)
    h_base[theSampleName+list_histos[23]] = ROOT.TH1F(theSampleName+list_histos[23], "mu-gamma InvMass", 50, 0, 300)
    h_base[theSampleName+list_histos[24]] = ROOT.TH1F(theSampleName+list_histos[24], "n Bjets 25", 6, 0, 6)
    h_base[theSampleName+list_histos[25]] = ROOT.TH1F(theSampleName+list_histos[25], "events with Bjets 25", 4, 0, 4)
    h_base[theSampleName+list_histos[26]] = ROOT.TH1F(theSampleName+list_histos[26], "Pion isolation 03 - mu", 75, 0, 150)
    h_base[theSampleName+list_histos[27]] = ROOT.TH1F(theSampleName+list_histos[27], "Pion isolation 05 - mu", 75, 0, 150)
    h_base[theSampleName+list_histos[28]] = ROOT.TH1F(theSampleName+list_histos[28], "Pion rel. isolation 05 - mu - ch", 50, 0, 10)
    h_base[theSampleName+list_histos[29]] = ROOT.TH1F(theSampleName+list_histos[29], "Pion rel. isolation 05 - mu", 50, 0, 10)
    h_base[theSampleName+list_histos[30]] = ROOT.TH1F(theSampleName+list_histos[30], "Pion isolation 03 - ele", 75, 0, 150)
    h_base[theSampleName+list_histos[31]] = ROOT.TH1F(theSampleName+list_histos[31], "Pion isolation 05 - ele", 75, 0, 150)
    h_base[theSampleName+list_histos[32]] = ROOT.TH1F(theSampleName+list_histos[32], "Pion rel. isolation 05 - ele - ch", 50, 0, 10)
    h_base[theSampleName+list_histos[33]] = ROOT.TH1F(theSampleName+list_histos[33], "Pion rel. isolation 05 - ele", 50, 0, 10)
    h_base[theSampleName+list_histos[34]] = ROOT.TH1F(theSampleName+list_histos[34], "mu-pi InvMass", 50, 0, 200)
    h_base[theSampleName+list_histos[35]] = ROOT.TH1F(theSampleName+list_histos[35], "met", 150, 0, 200)
    h_base[theSampleName+list_histos[36]] = ROOT.TH1F(theSampleName+list_histos[36], "Wmass mu plus", 10, 40, 100)
    h_base[theSampleName+list_histos[37]] = ROOT.TH1F(theSampleName+list_histos[37], "Wmass ele plus", 10, 40, 100)
    h_base[theSampleName+list_histos[38]] = ROOT.TH1F(theSampleName+list_histos[38], "Wmass mu minus", 10, 40, 100)
    h_base[theSampleName+list_histos[39]] = ROOT.TH1F(theSampleName+list_histos[39], "Wmass ele minus", 10, 40, 100)
    h_base[theSampleName+list_histos[40]] = ROOT.TH1F(theSampleName+list_histos[40], "eta of the pion", 20, -3, 3)
    h_base[theSampleName+list_histos[41]] = ROOT.TH1F(theSampleName+list_histos[41], "eta of the photon", 20, -3, 3)
    h_base[theSampleName+list_histos[42]] = ROOT.TH1F(theSampleName+list_histos[42], "Pion rel. isolation 05 - mu - ch -AfterCut", 10, 0, 1)
    h_base[theSampleName+list_histos[43]] = ROOT.TH1F(theSampleName+list_histos[43], "Pion rel. isolation 05 - ele - ch -AfterCut", 10, 0, 1)
    h_base[theSampleName+list_histos[44]] = ROOT.TH1F(theSampleName+list_histos[44], "p_{T} of the muon - sig", 15, 20, 100.)
    h_base[theSampleName+list_histos[45]] = ROOT.TH1F(theSampleName+list_histos[45], "p_{T} of the electron -sig", 15, 20, 100.)
    h_base[theSampleName+list_histos[46]] = ROOT.TH1F(theSampleName+list_histos[46], "p_{T} of the pion -sig", 15, 20, 100.)
    h_base[theSampleName+list_histos[47]] = ROOT.TH1F(theSampleName+list_histos[47], "E_{T} of the gamma -sig", 15, 20, 100.)
    h_base[theSampleName+list_histos[48]] = ROOT.TH1F(theSampleName+list_histos[48], "eta of the muon -sig", 20, -3, 3)
    h_base[theSampleName+list_histos[49]] = ROOT.TH1F(theSampleName+list_histos[49], "eta of the electron -sig", 20, -3, 3)
    h_base[theSampleName+list_histos[50]] = ROOT.TH1F(theSampleName+list_histos[50], "eta of the pion -sig", 20, -3, 3)
    h_base[theSampleName+list_histos[51]] = ROOT.TH1F(theSampleName+list_histos[51], "eta of the gamma -sig", 20, -3, 3)


#leg1 = ROOT.TLegend(0.15,0.6120093,0.34,0.9491917) #left positioning
leg1 = ROOT.TLegend(0.6868687,0.6120093,0.9511784,0.9491917) #right positioning
leg1.SetHeader(" ")
leg1.SetFillColor(0)
leg1.SetBorderSize(0)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)

leg2 = ROOT.TLegend(0.6,0.8,0.85,0.86)
leg2.SetHeader(" ")
leg2.SetFillColor(0)
leg2.SetBorderSize(0)
leg2.SetLineColor(1)
leg2.SetLineStyle(1)
leg2.SetLineWidth(1)
leg2.SetFillStyle(1001)

Nsig_passed = 0.
Ndata_passed = 0.
Nbkg_passed = 0.
mu_sig_events = 0
ele_sig_events = 0
mu_sig_events_deltaphi_mu_pi = 0
mu_sig_events_pi_pT = 0
mu_sig_events_gamma_eT = 0
mu_sig_events_nBjets = 0
mu_sig_events_Wmass = 0
ele_sig_events_deltaphi_ele_pi = 0
ele_sig_events_lep_iso = 0
ele_sig_events_elegamma = 0
ele_sig_events_pi_pT = 0
ele_sig_events_gamma_eT = 0
ele_sig_events_nBjets = 0
ele_sig_events_Wmass = 0
mu_bkg_events = 0
ele_bkg_events = 0
mu_bkg_events_deltaphi_mu_pi = 0
mu_bkg_events_pi_pT = 0
mu_bkg_events_gamma_eT = 0
mu_bkg_events_nBjets = 0
mu_bkg_events_Wmass = 0
ele_bkg_events_deltaphi_ele_pi = 0
ele_bkg_events_lep_iso = 0
ele_bkg_events_elegamma = 0
ele_bkg_events_pi_pT = 0
ele_bkg_events_gamma_eT = 0
ele_bkg_events_nBjets = 0
ele_bkg_events_Wmass = 0
piIso_03 = 0
piIso_05 = 0
ttbar_mu = 0
ttbar_ele = 0
#isMuon_evt = 0
#isNotMuon_evt = 0
Sevts_mu_SFvariation  = 0 # Counters for the number of signal events (weighted) when variating scale factors
Sevts_ele_SFvariation = 0
Sevts_tot = 0
Bevts_tot = 0
Sevts_weighted_mu = 0
Bevts_weighted_mu = 0
Sevts_weighted_ele = 0
Bevts_weighted_ele = 0
_Nrandom_for_SF = ROOT.TRandom3(44317)
_Nrandom_for_Gaus_SF = ROOT.TRandom3(44329)
N_WGToLNuG_mu = 0.

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


    if not "Data" in name_sample:
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

        if name_sample == "ttbar" and mytree.isttbarlnu: # Avoid double-counting of the ttbarlnu background
            continue

        #if "Data" in name_sample: continue  #-------------Excluding data-------------#
        #if not name_sample == "Signal":
        #    continue


        if "Signal" in name_sample:
            Sevts_tot += 1
        else:
            Bevts_tot += 1

        #---------- Access the tree variables ----------#

        isMuon = mytree.is_muon

        #run_number = mytree.run_number
        isSingleMuTrigger_24 = mytree.isSingleMuTrigger_24
        isSingleMuTrigger_50 = mytree.isSingleMuTrigger_50

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
        piRelIso_03 = mytree.sum_pT_03/pi_pT
        piRelIso_05 = mytree.sum_pT_05/pi_pT
        piRelIso_05_ch = mytree.sum_pT_05_ch/pi_pT
            
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

        met = mytree.met_pT
       
        Wmass = mytree.Wmass

        lep_FourMomentum = ROOT.TLorentzVector()
        lep_FourMomentum.SetPtEtaPhiM(lep_pT,lep_eta,lep_phi,0.)

        if not isMuon:
            ele_gamma_InvMass = (lep_FourMomentum + gamma_FourMomentum).M()
        else:
            mu_gamma_InvMass = (lep_FourMomentum + gamma_FourMomentum).M()
            mu_pi_InvMass = (lep_FourMomentum + pi_FourMomentum).M()
        
        nBjets = mytree.nBjets
        nBjets_25 = mytree.nBjets_25

        deltaeta_lep_pi = lep_eta-pi_eta
        
        deltaphi_lep_pi = math.fabs(lep_phi-pi_phi)
        if deltaphi_lep_pi > 3.14:
            deltaphi_lep_pi = 6.28 - deltaphi_lep_pi

        deltaphi_lep_W = math.fabs(lep_phi-W_phi)
        if deltaphi_lep_W > 3.14:
            deltaphi_lep_W = 6.28 - deltaphi_lep_W

        #---------Determine the total event weight---------#
        
        if isMuon: # Get muon scale factors, which are different for two groups of datasets, and weight them for the respective integrated lumi 
            mu_weight_BtoF, mu_weight_BtoF_err = myWF.get_muon_scale_BtoF(lep_pT,lep_eta,isSingleMuTrigger_24)
            mu_weight_GH, mu_weight_GH_err     = myWF.get_muon_scale_GH(lep_pT,lep_eta,isSingleMuTrigger_24)

            # Use a random number to select which muon scale factor to use, depending on the associated lumi fraction
            Nrandom_for_SF = _Nrandom_for_SF.Rndm()

            if Nrandom_for_SF <= (luminosity_BtoF/luminosity_norm):  # Accessing muon SF, B to F

                if random_mu_SF:
                    mu_weight = _Nrandom_for_Gaus_SF.Gaus(mu_weight_BtoF,mu_weight_BtoF_err)

                else:
                    mu_weight = mu_weight_BtoF

            else: #Accessing muon SF, G and H
                
                if random_mu_SF:
                    mu_weight = _Nrandom_for_Gaus_SF.Gaus(mu_weight_GH,mu_weight_GH_err)

                else:
                    mu_weight = mu_weight_GH

        else:
            ele_weight, ele_weight_err = myWF.get_ele_scale(lep_pT,lep_eta)

            if random_ele_SF:
                ele_weight = _Nrandom_for_Gaus_SF.Gaus(ele_weight,ele_weight_err) 

        ph_weight, ph_weight_err = myWF.get_photon_scale(gamma_eT,gamma_eta)
        
        if random_ph_SF:
            ph_weight = _Nrandom_for_Gaus_SF.Gaus(ph_weight,ph_weight_err)
        

        if not "Data" in name_sample:
            PU_Weight = mytree.PU_Weight        
            Event_Weight = norm_factor*PU_Weight*ph_weight

            if isMuon:
                Event_Weight = Event_Weight*mu_weight
            else:
                Event_Weight = Event_Weight*ele_weight




            # Obtaining the number of sig and bkg events (weighted)
            
            if "Signal" in name_sample and isMuon:
                Sevts_weighted_mu += Event_Weight
            if not "Signal" in name_sample and isMuon:
                Bevts_weighted_mu += Event_Weight
            if "Signal" in name_sample and not isMuon:
                Sevts_weighted_ele += Event_Weight
            if not "Signal" in name_sample and not isMuon:
                Bevts_weighted_ele += Event_Weight
            

        else:
            Event_Weight = 1.
 

        #---------Retrieve the BDT output----------#

        BDT_out = myWF.get_BDT_output(pi_pT,gamma_eT,nBjets_25,lep_pT,piRelIso_05_ch,met,isMuon)    
        # BDT_out = myWF.get_BDT_output(pi_pT/Wmass,gamma_eT/Wmass,nBjets_25,lep_pT,piRelIso_05_ch,met,isMuon)    

        #---------- filling histos ------------
        
        # if select_all_but_one("h_mupt") and isMuon:
        #     h_base[theSampleName+"h_mupt"].Fill(lep_pT,Event_Weight)
                    
        # if select_all_but_one("h_elept") and not isMuon:
        #     h_base[theSampleName+"h_elept"].Fill(lep_pT,Event_Weight)
                        
        # if select_all_but_one("h_pipt"):
        #     h_base[theSampleName+"h_pipt"].Fill(pi_pT,Event_Weight)
                            
        if select_all_but_one("h_nBjets"):
            h_base[theSampleName+"h_nBjets"].Fill(nBjets,Event_Weight)
        #if not "Data" in name_sample:
        h_base[theSampleName+"h_nBjets_25"].Fill(nBjets_25,Event_Weight)
                                
        # if select_all_but_one("h_gammaet"):
        #     h_base[theSampleName+"h_gammaet"].Fill(gamma_eT,Event_Weight)
                                  
        # if select_all_but_one("h_Wmass"):
        
        #     if "Data" in name_sample and (Wmass < 65. or Wmass > 90):
        #         h_base[theSampleName+"h_Wmass"].Fill(Wmass,Event_Weight)
        #     if not "Data" in name_sample:
        #         h_base[theSampleName+"h_Wmass"].Fill(Wmass,Event_Weight)

        #if select_all_but_one("h_deltaphi_mu_pi") and isMuon:
        #    h_base[theSampleName+"h_deltaphi_mu_pi"].Fill(deltaphi_lep_pi,Event_Weight)

        if select_all_but_one("all cuts") and isMuon:
            h_base[theSampleName+"h_deltaeta_mu_pi"].Fill(deltaeta_lep_pi,Event_Weight)
            #h_base[theSampleName+"h_mueta"].Fill(lep_eta,Event_Weight)
                    
        #if select_all_but_one("h_deltaphi_ele_pi") and not isMuon:
        #    h_base[theSampleName+"h_deltaphi_ele_pi"].Fill(deltaphi_lep_pi,Event_Weight)

        if select_all_but_one("all cuts") and not isMuon:
            h_base[theSampleName+"h_deltaeta_ele_pi"].Fill(deltaeta_lep_pi,Event_Weight)
            #h_base[theSampleName+"h_eleeta"].Fill(lep_eta,Event_Weight)

        #if select_all_but_one("h_ele_gamma_InvMass") and not isMuon:
        if not isMuon:
            h_base[theSampleName+"h_ele_gamma_InvMass"].Fill(ele_gamma_InvMass,Event_Weight)
        
        
        # if select_all_but_one("h_Wmass") and isMuon:
        
        #     if "Data" in name_sample and (Wmass < 65. or Wmass > 90):
        #         h_base[theSampleName+"h_Wmass_flag_mu"].Fill(Wmass,Event_Weight)
        #     if not "Data" in name_sample:
        #         h_base[theSampleName+"h_Wmass_flag_mu"].Fill(Wmass,Event_Weight)

        # if select_all_but_one("h_Wmass") and not isMuon:
        
        #     if "Data" in name_sample and (Wmass < 65. or Wmass > 90):
        #         h_base[theSampleName+"h_Wmass_flag_ele"].Fill(Wmass,Event_Weight)
        #     if not "Data" in name_sample:
        #         h_base[theSampleName+"h_Wmass_flag_ele"].Fill(Wmass,Event_Weight)
        
        if select_all_but_one("all cuts") and isMuon:
            h_base[theSampleName+"h_mu_iso"].Fill(lep_iso,Event_Weight)
            h_base[theSampleName+"h_deltaphi_mu_W"].Fill(deltaphi_lep_W,Event_Weight)
            h_base[theSampleName+"h_mu_gamma_InvMass"].Fill(mu_gamma_InvMass,Event_Weight)   

        if select_all_but_one("h_ele_iso") and not isMuon:
            h_base[theSampleName+"h_ele_iso"].Fill(lep_iso,Event_Weight)
            h_base[theSampleName+"h_deltaphi_ele_W"].Fill(deltaphi_lep_W,Event_Weight)

        if select_all_but_one("all cuts"):
            h_base[theSampleName+"h_gamma_iso_ChHad"].Fill(gamma_iso_ChHad,Event_Weight)
            h_base[theSampleName+"h_gamma_iso_NeuHad"].Fill(gamma_iso_NeuHad,Event_Weight)
            h_base[theSampleName+"h_gamma_iso_Ph"].Fill(gamma_iso_Ph,Event_Weight)
            h_base[theSampleName+"h_gamma_iso_eArho"].Fill(gamma_iso_eArho,Event_Weight)
        
        if select_all_but_one("h_nBjets"):
            if nBjets_25 == 0:
                h_base[theSampleName+"h_evts_Bjets"].Fill(0.5,1)
            if nBjets > 0:
                h_base[theSampleName+"h_evts_Bjets"].Fill(1.5,1)
            if nBjets > 2:
                h_base[theSampleName+"h_evts_Bjets"].Fill(2.5,1)
            if nBjets > 3:
                h_base[theSampleName+"h_evts_Bjets"].Fill(3.5,1)

        if isMuon:
            h_base[theSampleName+"h_piRelIso_05_mu"].Fill(piRelIso_05,Event_Weight)
            h_base[theSampleName+"h_mu_pi_InvMass"].Fill(mu_pi_InvMass,Event_Weight)
            h_base[theSampleName+"h_piRelIso_05_mu_ch"].Fill(piRelIso_05_ch,Event_Weight)

        h_base[theSampleName+"h_met"].Fill(met,Event_Weight)


        if not isMuon:
            h_base[theSampleName+"h_piRelIso_05_ele"].Fill(piRelIso_05,Event_Weight)
            h_base[theSampleName+"h_piRelIso_05_ele_ch"].Fill(piRelIso_05_ch,Event_Weight)

        if name_sample == "Signal":
            h_base[theSampleName+"h_pipt_sig"].Fill(pi_pT,Event_Weight)
            h_base[theSampleName+"h_pieta_sig"].Fill(pi_eta,Event_Weight)
            h_base[theSampleName+"h_gammaet_sig"].Fill(gamma_eT,Event_Weight)
            h_base[theSampleName+"h_gammaeta_sig"].Fill(gamma_eta,Event_Weight)
            if isMuon:
                h_base[theSampleName+"h_mueta_sig"].Fill(lep_eta,Event_Weight)
                h_base[theSampleName+"h_mupt_sig"].Fill(lep_pT,Event_Weight)
            else:
                h_base[theSampleName+"h_eleeta_sig"].Fill(lep_eta,Event_Weight)
                h_base[theSampleName+"h_elept_sig"].Fill(lep_pT,Event_Weight)

        if isMuon:
            h_base[theSampleName+"h_Wmass_flag_mu"].Fill(Wmass,Event_Weight)
            h_base[theSampleName+"h_mupt"].Fill(lep_pT,Event_Weight)
            h_base[theSampleName+"h_mueta"].Fill(lep_eta,Event_Weight)

        if not isMuon:
            h_base[theSampleName+"h_Wmass_flag_ele"].Fill(Wmass,Event_Weight)
            h_base[theSampleName+"h_elept"].Fill(lep_pT,Event_Weight)
            h_base[theSampleName+"h_eleeta"].Fill(lep_eta,Event_Weight)

        h_base[theSampleName+"h_Wmass"].Fill(Wmass,Event_Weight)
        h_base[theSampleName+"h_pipt"].Fill(pi_pT,Event_Weight)
        h_base[theSampleName+"h_pieta"].Fill(pi_eta,Event_Weight)
        h_base[theSampleName+"h_gammaet"].Fill(gamma_eT,Event_Weight)
        h_base[theSampleName+"h_gammaeta"].Fill(gamma_eta,Event_Weight)


        #---------------------Here's where the BDT selection starts---------------------#
      
        # if (isMuon and BDT_out >= 0.255) or (not isMuon and BDT_out >= 0.250):
        # # if (isMuon and BDT_out >= 0.223) or (not isMuon and BDT_out >= 0.238): #Wmass
        #     if (Wmass >= 50. and Wmass <= 100.):
        #         if "Data" in name_sample and (Wmass < 65. or Wmass > 90):
        #             h_base[theSampleName+"h_Wmass"].Fill(Wmass,Event_Weight)
        #         if not "Data" in name_sample:
        #             h_base[theSampleName+"h_Wmass"].Fill(Wmass,Event_Weight)

        #         h_base[theSampleName+"h_pipt"].Fill(pi_pT,Event_Weight)
        #         h_base[theSampleName+"h_gammaet"].Fill(gamma_eT,Event_Weight)
        #         h_base[theSampleName+"h_pieta"].Fill(pi_eta,Event_Weight)
        #         h_base[theSampleName+"h_gammaeta"].Fill(gamma_eta,Event_Weight)

        #         if isMuon:
        #             if "Data" in name_sample and (Wmass < 65. or Wmass > 90):
        #                 h_base[theSampleName+"h_Wmass_flag_mu"].Fill(Wmass,Event_Weight)
        #             if not "Data" in name_sample:
        #                 h_base[theSampleName+"h_Wmass_flag_mu"].Fill(Wmass,Event_Weight)
                        
        #                 if "Signal" in name_sample:
        #                     Sevts_mu_SFvariation += Event_Weight

        #             h_base[theSampleName+"h_mupt"].Fill(lep_pT,Event_Weight)
        #             h_base[theSampleName+"h_mueta"].Fill(lep_eta,Event_Weight)
        #             h_base[theSampleName+"h_piRelIso_05_mu_ch_AfterCut"].Fill(piRelIso_05_ch,Event_Weight)

        #             if name_sample == "WGToLNuG":
        #                 N_WGToLNuG_mu += Event_Weight

 
        #         if (not isMuon) and lep_iso <= 0.35:
        #             if "Data" in name_sample and (Wmass < 65. or Wmass > 90):
        #                 h_base[theSampleName+"h_Wmass_flag_ele"].Fill(Wmass,Event_Weight)
        #             if not "Data" in name_sample:
        #                 h_base[theSampleName+"h_Wmass_flag_ele"].Fill(Wmass,Event_Weight)

        #                 if "Signal" in name_sample:
        #                     Sevts_ele_SFvariation += Event_Weight
                
        #             h_base[theSampleName+"h_elept"].Fill(lep_pT,Event_Weight)
        #             h_base[theSampleName+"h_eleeta"].Fill(lep_eta,Event_Weight)
        #             h_base[theSampleName+"h_piRelIso_05_ele_ch_AfterCut"].Fill(piRelIso_05_ch,Event_Weight)


        #-------BDT cut variation------#

        if (isMuon and BDT_out >= -0.1 and BDT_out < 0.255) or (not isMuon and BDT_out >= -0.1 and BDT_out < 0.250):
        # if (isMuon and BDT_out >= -0.1 and BDT_out < 0.223) or (not isMuon and BDT_out >= -0.1 and BDT_out < 0.238): #Wmass
            if (Wmass >= 50. and Wmass <= 100.) and isMuon:
                h_base[theSampleName+"h_Wmass_mu_minus"].Fill(Wmass,Event_Weight)
            if (Wmass >= 50. and Wmass <= 100.) and not isMuon and lep_iso <= 0.35:
                h_base[theSampleName+"h_Wmass_ele_minus"].Fill(Wmass,Event_Weight)

        # if (isMuon and BDT_out >= 0.295) or (not isMuon and BDT_out >= 0.270):
        # if (isMuon and BDT_out >= 0.265) or (not isMuon and BDT_out >= 0.295): #Wmass
            # if (Wmass >= 50. and Wmass <= 100.) and isMuon:
            #     h_base[theSampleName+"h_Wmass_mu_plus"].Fill(Wmass,Event_Weight)
            # if (Wmass >= 50. and Wmass <= 100.) and not isMuon and lep_iso <= 0.35:
            #     h_base[theSampleName+"h_Wmass_ele_plus"].Fill(Wmass,Event_Weight)

        
        #Count the events
        if select_all_but_one("all cuts"):
            if name_sample == myWF.sig_samplename:
                Nsig_passed += Event_Weight
            elif name_sample == "Data":
                Ndata_passed += Event_Weight
            else:
                Nbkg_passed += Event_Weight

  #--------- progressive histos ----------
        if name_sample == myWF.sig_samplename:
            if isMuon:
                mu_sig_events += 1
                if deltaphi_lep_pi >= DELTAPHI_MU_PI_MIN: mu_sig_events_deltaphi_mu_pi += 1
                else: continue
                if pi_pT >= PI_MIN_PT: mu_sig_events_pi_pT += 1
                else: continue
                if gamma_eT >= GAMMA_MIN_ET: mu_sig_events_gamma_eT += 1
                else: continue
                if nBjets_25 >= N_BJETS_MIN: mu_sig_events_nBjets += 1
                else: continue
                if (Wmass >= WMASS_MIN and Wmass <= WMASS_MAX): mu_sig_events_Wmass += 1
            if not isMuon:
                ele_sig_events += 1
                if deltaphi_lep_pi >= DELTAPHI_ELE_PI_MIN: ele_sig_events_deltaphi_ele_pi += 1
                else: continue
                if lep_iso <= ELE_ISO_MAX: ele_sig_events_lep_iso += 1
                else: continue
                if (ele_gamma_InvMass < ELE_GAMMA_INVMASS_MIN or ele_gamma_InvMass > ELE_GAMMA_INVMASS_MAX): ele_sig_events_elegamma += 1
                else: continue
                if pi_pT >= PI_MIN_PT: ele_sig_events_pi_pT += 1
                else: continue
                if gamma_eT >= GAMMA_MIN_ET: ele_sig_events_gamma_eT += 1
                else: continue
                if nBjets_25 >= N_BJETS_MIN: ele_sig_events_nBjets += 1
                else: continue
                if (Wmass >= WMASS_MIN and Wmass <= WMASS_MAX): ele_sig_events_Wmass += 1
        else:
            if not "Data" in name_sample:
                if isMuon:
                    mu_bkg_events += 1
                    if deltaphi_lep_pi >= DELTAPHI_MU_PI_MIN: mu_bkg_events_deltaphi_mu_pi += 1
                    else: continue
                    if pi_pT >= PI_MIN_PT: mu_bkg_events_pi_pT += 1
                    else: continue
                    if gamma_eT >= GAMMA_MIN_ET: mu_bkg_events_gamma_eT += 1
                    else: continue
                    if nBjets_25 >= N_BJETS_MIN: mu_bkg_events_nBjets += 1
                    else: continue
                    if (Wmass >= WMASS_MIN and Wmass <= WMASS_MAX): mu_bkg_events_Wmass += 1
                if not isMuon:
                    ele_bkg_events += 1
                    if deltaphi_lep_pi >= DELTAPHI_ELE_PI_MIN: ele_bkg_events_deltaphi_ele_pi += 1
                    else:continue
                    if lep_iso <= ELE_ISO_MAX: ele_bkg_events_lep_iso += 1
                    else: continue
                    if (ele_gamma_InvMass < ELE_GAMMA_INVMASS_MIN or ele_gamma_InvMass > ELE_GAMMA_INVMASS_MAX): ele_bkg_events_elegamma += 1
                    else: continue
                    if pi_pT >= PI_MIN_PT: ele_bkg_events_pi_pT += 1
                    else: continue
                    if gamma_eT >= GAMMA_MIN_ET: ele_bkg_events_gamma_eT += 1
                    else: continue
                    if nBjets_25 >= N_BJETS_MIN: ele_bkg_events_nBjets += 1
                    else: continue
                    if (Wmass >= WMASS_MIN and Wmass <= WMASS_MAX): ele_bkg_events_Wmass += 1
                

    for idx_histo,hname in enumerate(list_histos):

        if QCDflag:
            h_base[theSampleName+hname].SetFillColor(12)
        elif name_sample == myWF.sig_samplename:
            h_base[theSampleName+hname].SetLineStyle(2)   #dashed
            h_base[theSampleName+hname].SetLineColor(2)   #red
            h_base[theSampleName+hname].SetLineWidth(4)   #kind of thick
        elif name_sample == myWF.data_samplename:
            h_base[theSampleName+hname].SetMarkerStyle(20)   #dashed
        else:
            h_base[theSampleName+hname].SetFillColor(colors_mask[idx_sample])


        if idx_histo == 0:
            if QCDflag and isFirstQCDlegend:
                leg1.AddEntry(h_base[theSampleName+hname],"QCD","f")
                isFirstQCDlegend = False
            elif name_sample == myWF.sig_samplename:
                sample_legend_name = str(signal_magnify) + " x " + name_sample
                leg1.AddEntry(h_base[name_sample+hname], sample_legend_name,"f")  #To comment when signal has to be excluded.
                #leg2.AddEntry(h_base[name_sample+hname], sample_legend_name,"f")
            elif name_sample == myWF.data_samplename:
                leg1.AddEntry(h_base[name_sample+hname], name_sample,"lep") # lep shows on the TLegend a point with errors to indicate data
            elif not QCDflag and not name_sample == myWF.data_samplename:
                leg1.AddEntry(h_base[theSampleName+hname],theSampleName,"f")

        if not QCDflag and not name_sample == myWF.sig_samplename and not name_sample == myWF.data_samplename:
            hs[hname].Add(h_base[theSampleName+hname])

    if not QCDflag and not name_sample == myWF.sig_samplename and not name_sample == myWF.data_samplename:
        idx_sample += 1

print "Finished runnning over samples!"

for idx_histo,hname in enumerate(list_histos):
    hs[hname].Add(h_base["QCD_"+hname])

for hname in list_histos:

    canvas[hname].cd()

    hs[hname].Draw("histo")
    # if "h_Wmass_" in hname:
    #     hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),65.))
    # if hname == "h_Wmass":
    #     hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),100.))
    if hname == "h_piRelIso_05_ele_ch":
        hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),27000.))
    if hname == "h_piRelIso_05_mu_ch":
        hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),12000.))
    if hname == "h_piRelIso_05_mu_ch_AfterCut":
        hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),130.))
    if hname == "h_piRelIso_05_ele_ch_AfterCut":
        hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),130.))
    # if "eta" in hname:
    #     hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),90.))
    # if "pt" in hname:
    #     hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),90.))
    # if hname == "h_gammaet":
    #     hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),100.))
    if hname == "h_pieta_sig" or hname == "h_pipt_sig" or hname == "h_gammaeta_sig" or hname == "h_gammaet_sig":
        hs[hname].SetMaximum(30.)
    if hname == "h_mueta_sig" or hname == "h_mupt_sig" or hname == "h_eleeta_sig" or hname == "h_elept_sig":
        hs[hname].SetMaximum(20.)


    down = ROOT.gPad.GetUymin()
    up   = ROOT.gPad.GetUymax()

    #---------Histos names---------#
    
    hs[hname].SetTitle(" ")
    hs[hname].GetYaxis().SetTitleOffset(1.7)
    #hs[hname].GetYaxis().SetTitle("Events")
    
    if hname == "h_Wmass" or hname == "h_Wmass_flag_mu" or hname == "h_Wmass_flag_ele":
        hs[hname].GetXaxis().SetTitle("m_{#pi#gamma} (GeV)")
        line = ROOT.TLine(WMASS_MIN, 0, WMASS_MIN, 72)
        line.SetLineColor(8)
        line.SetLineStyle(9)
        line.SetLineWidth(4)
        #line.Draw()

    if hname == "h_mupt" or hname == "h_mupt_sig":
        hs[hname].GetXaxis().SetTitle("p_{T}^{#mu} (GeV)")
        line = ROOT.TLine(MU_MIN_PT, 0, MU_MIN_PT, 520)
        line.SetLineColor(8)
        line.SetLineStyle(9)
        line.SetLineWidth(4)
        # line.Draw()

    if hname == "h_elept" or hname == "h_elept_sig":
        hs[hname].GetXaxis().SetTitle("p_{T}^{e} (GeV)")
        line = ROOT.TLine(ELE_MIN_PT, 0, ELE_MIN_PT, 520)
        line.SetLineColor(8)
        line.SetLineStyle(9)
        line.SetLineWidth(4)
        # line.Draw()
    
    if hname == "h_pipt" or hname == "h_pipt_sig":
        hs[hname].GetXaxis().SetTitle("p_{T}^{#pi} (GeV)")
        line = ROOT.TLine(PI_MIN_PT, 0, PI_MIN_PT, 520)
        line.SetLineColor(8)
        line.SetLineStyle(9)
        line.SetLineWidth(4)
        # line.Draw()

    if hname == "h_gammaet" or hname == "h_gammaet_sig":
        hs[hname].GetXaxis().SetTitle("E_{T}^{#gamma} (GeV)")
        line = ROOT.TLine(GAMMA_MIN_ET, 0, GAMMA_MIN_ET, 520)
        line.SetLineColor(8)
        line.SetLineStyle(9)
        line.SetLineWidth(4)
        # line.Draw()

    if hname == "h_gammaeta" or hname == "h_gammaeta_sig":
        hs[hname].GetXaxis().SetTitle("#eta^{#gamma}")

    if hname == "h_mueta" or hname == "h_mueta_sig":
        hs[hname].GetXaxis().SetTitle("#eta^{#mu}")

    if hname == "h_eleeta" or hname == "h_eleeta_sig":
        hs[hname].GetXaxis().SetTitle("#eta^{e}")

    if hname == "h_pieta" or hname == "h_pieta_sig":
        hs[hname].GetXaxis().SetTitle("#eta^{#pi}")

    if hname == "h_deltaphi_mu_pi":
        hs[hname].GetXaxis().SetTitle("#Delta#varphi_{#mu-#pi})")

    if hname == "h_deltaphi_ele_pi":
        hs[hname].GetXaxis().SetTitle("#Delta#varphi_{e-#pi})")

    if hname == "h_ele_gamma_InvMass":
        hs[hname].GetXaxis().SetTitle("m_{e#gamma} (GeV/c^{2})")

    if hname == "h_nBjets":
        hs[hname].GetXaxis().SetTitle("Number of b-jets")

    if hname == "h_nBjets_25":
        hs[hname].GetXaxis().SetTitle("Number of b-jets (p_{T}>25 GeV/c)")

    if "h_piRelIso" in hname:
        hs[hname].GetXaxis().SetTitle("#pi_{Iso}/p_{T}^{#pi}")
    
    if signal_magnify != 1:
        h_base[myWF.sig_samplename+hname].Scale(signal_magnify)


    #---Wmass ratio plots---#
    if hname == "h_Wmass_flag_mu":
        Wmass_mu = hs[hname].GetStack().Last()

    if hname == "h_Wmass_flag_ele":
        Wmass_ele = hs[hname].GetStack().Last()

    # if hname == "h_Wmass_mu_plus":
    #     Wmass_mu_plus = hs[hname].GetStack().Last()

    if hname == "h_Wmass_mu_minus":
        Wmass_mu_minus = hs[hname].GetStack().Last()

    # if hname == "h_Wmass_ele_plus":
    #     Wmass_ele_plus = hs[hname].GetStack().Last()

    if hname == "h_Wmass_ele_minus":
        Wmass_ele_minus = hs[hname].GetStack().Last()


    h_base[myWF.sig_samplename+hname].Draw("SAME,hist")
    h_base[myWF.data_samplename+hname].Draw("SAME,E1")

    hMCErr = copy.deepcopy(hs[hname].GetStack().Last())
    hMCErr.SetFillStyle(3005)
    hMCErr.SetMarkerStyle(1)
    hMCErr.SetFillColor(ROOT.kBlack)
    hMCErr.Draw("sameE2")
    
    if not "_sig" in hname:
        leg1.Draw()
 
        
    canvas[hname].SaveAs("plots/" + hname + ".pdf")


#------draw progressive histo------
canvas1 = ROOT.TCanvas()
h_events_sig_mu.Fill(0.5,mu_sig_events)
#h_events_sig_mu.Fill(1.5,mu_sig_events_deltaphi_mu_pi)
h_events_sig_mu.Fill(1.5,mu_sig_events_pi_pT)
h_events_sig_mu.Fill(2.5,mu_sig_events_gamma_eT)
h_events_sig_mu.Fill(3.5,mu_sig_events_nBjets)
h_events_sig_mu.Fill(4.5,mu_sig_events_Wmass)

ROOT.gStyle.SetOptStat(0)
h_events_sig_mu.GetXaxis().SetBinLabel(1,"initial events")
#h_events_sig_mu.GetXaxis().SetBinLabel(2,"#Delta#varphi(#mu-#pi)")
h_events_sig_mu.GetXaxis().SetBinLabel(2,"p_{T}^{#pi}")
h_events_sig_mu.GetXaxis().SetBinLabel(3,"E_{T}^{#gamma}")
h_events_sig_mu.GetXaxis().SetBinLabel(4,"nBjets")
h_events_sig_mu.GetXaxis().SetBinLabel(5,"Wmass")
h_events_sig_mu.Draw("hist")
ROOT.gPad.SetLogy()
canvas1.SaveAs("plots/h_events_sig_mu.pdf")


canvas2 = ROOT.TCanvas()
h_events_sig_ele.Fill(0.5,ele_sig_events)
#h_events_sig_ele.Fill(1.5,ele_sig_events_deltaphi_ele_pi)
h_events_sig_ele.Fill(1.5,ele_sig_events_lep_iso)
h_events_sig_ele.Fill(2.5,ele_sig_events_elegamma)
h_events_sig_ele.Fill(3.5,ele_sig_events_pi_pT)
h_events_sig_ele.Fill(4.5,ele_sig_events_gamma_eT)
h_events_sig_ele.Fill(5.5,ele_sig_events_nBjets)
h_events_sig_ele.Fill(6.5,ele_sig_events_Wmass)

ROOT.gStyle.SetOptStat(0)
h_events_sig_ele.GetXaxis().SetBinLabel(1,"initial events")
#h_events_sig_ele.GetXaxis().SetBinLabel(2,"#Delta#varphi(e-#pi)")
h_events_sig_ele.GetXaxis().SetBinLabel(2,"ele iso")
h_events_sig_ele.GetXaxis().SetBinLabel(3,"e-#gamma inv mass")
h_events_sig_ele.GetXaxis().SetBinLabel(4,"p_{T}^{#pi}")
h_events_sig_ele.GetXaxis().SetBinLabel(5,"E_{T}^{#gamma}")
h_events_sig_ele.GetXaxis().SetBinLabel(6,"nBjets")
h_events_sig_ele.GetXaxis().SetBinLabel(7,"Wmass")
h_events_sig_ele.Draw("hist")
ROOT.gPad.SetLogy()
canvas2.SaveAs("plots/h_events_sig_ele.pdf")


canvas3 = ROOT.TCanvas()
h_events_bkg_mu.Fill(0.5,mu_bkg_events)
#h_events_bkg_mu.Fill(1.5,mu_bkg_events_deltaphi_mu_pi)
h_events_bkg_mu.Fill(1.5,mu_bkg_events_pi_pT)
h_events_bkg_mu.Fill(2.5,mu_bkg_events_gamma_eT)
h_events_bkg_mu.Fill(3.5,mu_bkg_events_nBjets)
h_events_bkg_mu.Fill(4.5,mu_bkg_events_Wmass)

ROOT.gStyle.SetOptStat(0)
h_events_bkg_mu.GetXaxis().SetBinLabel(1,"initial events")
#h_events_bkg_mu.GetXaxis().SetBinLabel(2,"#Delta#varphi(#mu-#pi)")
h_events_bkg_mu.GetXaxis().SetBinLabel(2,"p_{T}^{#pi}")
h_events_bkg_mu.GetXaxis().SetBinLabel(3,"E_{T}^{#gamma}")
h_events_bkg_mu.GetXaxis().SetBinLabel(4,"nBjets")
h_events_bkg_mu.GetXaxis().SetBinLabel(5,"Wmass")
h_events_bkg_mu.Draw("hist")
ROOT.gPad.SetLogy()
canvas3.SaveAs("plots/h_events_bkg_mu.pdf")


canvas4 = ROOT.TCanvas()
h_events_bkg_ele.Fill(0.5,ele_bkg_events)
#h_events_bkg_ele.Fill(1.5,ele_bkg_events_deltaphi_ele_pi)
h_events_bkg_ele.Fill(1.5,ele_bkg_events_lep_iso)
h_events_bkg_ele.Fill(2.5,ele_bkg_events_elegamma)
h_events_bkg_ele.Fill(3.5,ele_bkg_events_pi_pT)
h_events_bkg_ele.Fill(4.5,ele_bkg_events_gamma_eT)
h_events_bkg_ele.Fill(5.5,ele_bkg_events_nBjets)
h_events_bkg_ele.Fill(6.5,ele_bkg_events_Wmass)

ROOT.gStyle.SetOptStat(0)
h_events_bkg_ele.GetXaxis().SetBinLabel(1,"initial events")
#h_events_bkg_ele.GetXaxis().SetBinLabel(2,"#Delta#varphi(e-#pi)")
h_events_bkg_ele.GetXaxis().SetBinLabel(2,"ele iso")
h_events_bkg_ele.GetXaxis().SetBinLabel(3,"e-#gamma inv mass")
h_events_bkg_ele.GetXaxis().SetBinLabel(4,"p_{T}^{#pi}")
h_events_bkg_ele.GetXaxis().SetBinLabel(5,"E_{T}^{#gamma}")
h_events_bkg_ele.GetXaxis().SetBinLabel(6,"nBjets")
h_events_bkg_ele.GetXaxis().SetBinLabel(7,"Wmass")
h_events_bkg_ele.Draw("hist")
ROOT.gPad.SetLogy()
canvas4.SaveAs("plots/h_events_bkg_ele.pdf")

Wmass_mu.Scale(1/Wmass_mu.Integral())
#Wmass_mu_plus.Scale(1/Wmass_mu_plus.Integral())
Wmass_mu_minus.Scale(1/Wmass_mu_minus.Integral())
Wmass_ele.Scale(1/Wmass_ele.Integral())
# Wmass_ele_plus.Scale(1/Wmass_ele_plus.Integral())
Wmass_ele_minus.Scale(1/Wmass_ele_minus.Integral())

#Wmass_mu_plus.Divide(Wmass_mu_plus,Wmass_mu,1.0,1.0,"B")
Wmass_mu_minus.Divide(Wmass_mu_minus,Wmass_mu,1.0,1.0,"B")
#Wmass_ele_plus.Divide(Wmass_ele_plus,Wmass_ele,1.0,1.0,"B")
Wmass_ele_minus.Divide(Wmass_ele_minus,Wmass_ele,1.0,1.0,"B")

# canvas5 = ROOT.TCanvas()
# Wmass_mu_plus.SetMarkerStyle(21)
# Wmass_mu_plus.Draw("Pe")
# canvas5.SaveAs("plots/Wmass_mu_plus.pdf")

canvas6 = ROOT.TCanvas()
Wmass_mu_minus.SetMarkerStyle(21)
Wmass_mu_minus.SetTitle(" ")
Wmass_mu_minus.GetXaxis().SetTitle("m_{#pi#gamma} (GeV)")
Wmass_mu_minus.Draw("Pe")
canvas6.SaveAs("plots/Wmass_mu_minus.pdf")

# canvas7 = ROOT.TCanvas()
# Wmass_ele_plus.SetMarkerStyle(21)
# Wmass_ele_plus.Draw("Pe")
# canvas7.SaveAs("plots/Wmass_ele_plus.pdf")

canvas8 = ROOT.TCanvas()
Wmass_ele_minus.SetMarkerStyle(21)
Wmass_ele_minus.SetTitle(" ")
Wmass_ele_minus.GetXaxis().SetTitle("m_{#pi#gamma} (GeV)")
Wmass_ele_minus.Draw("Pe")
canvas8.SaveAs("plots/Wmass_ele_minus.pdf")

# canvas9 = ROOT.TCanvas()
# Wmass_mu.SetMarkerStyle(21)
# Wmass_mu.Draw("Pe")
# canvas9.SaveAs("plots/Wmass_mu.pdf")

print "Number of expected events for ", luminosity_norm, " in fb-1"
print "Number of signal events passed = ", Nsig_passed
print "Number of background events passed = ", Nbkg_passed
print "Number of data events passed = ", Ndata_passed
# print "Significance S/sqrt(B) = ", Nsig_passed/math.sqrt(Nbkg_passed)
# print "Significance S/sqrt(B + deltaB^2) = ", Nsig_passed/(math.sqrt(Nbkg_passed) + 0.2*Nbkg_passed)
# print "Significance S/sqrt(S+B) = ", Nsig_passed/math.sqrt(Nsig_passed + Nbkg_passed)
# print "\nAll the intresting plots have been produced..!"

print "number of Signal events in muon channel: ", Sevts_mu_SFvariation
print "number of Signal events in electron channel: ", Sevts_ele_SFvariation
#print "total number of S evts: ", Sevts_tot
#print "total number of B evts: ", Bevts_tot
print "total number of S evts weighted -mu: ", Sevts_weighted_mu
print "total number of B evts weighted -mu: ", Bevts_weighted_mu
print "total number of S evts weighted -ele: ", Sevts_weighted_ele
print "total number of B evts weighted -ele: ", Bevts_weighted_ele

print "N_WGToLNuG_mu", N_WGToLNuG_mu 
