import ROOT
import os
import math
import numpy as np
from ROOT import TFile, TTree, TBranch, TCanvas, TH1F
from Simplified_Workflow_Handler import Simplified_Workflow_Handler
import argparse

#---------------------------------#
p = argparse.ArgumentParser(description='Select whether to fill the mass tree to be used for fits, or the MVA trees')
p.add_argument('create_mass_tree_option', help='Type <<mass>> or <<MVA>>')
p.add_argument('runningEra_option', help='Type <<0>> for 2016, <<1>> for 2017 and <<2>> for 2018')
args = p.parse_args()

# Switch from muon to electron channel
if args.create_mass_tree_option == "mass":
    create_mass_tree = True
if args.create_mass_tree_option == "MVA":
    create_mass_tree = False
runningEra = int(args.runningEra_option)

#------------------------------------------------------------------------#

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)

isData   = False # Switch from DATA to MC and vice versa
isBDT_with_Wmass = False # If true, pT(pi) and ET(gamma) in the BDT are normalized to Wmass 
split_MC = False # If True, MC signal sample is split in two for the training/testing of the BDT
scale_signal_up = False # If true, scales the MC signal weight up (Pythia modeling) 
scale_signal_down = False # If true, scales the MC signal weight down (Pythia modeling) 
scale_signal_sin2 = False
scale_signal_cos = False
only_signal_relevant_weights = False

myWF = Simplified_Workflow_Handler("Signal","Data",create_mass_tree,isBDT_with_Wmass,runningEra)

############################################################################
#                                                                          #
#-------------------------- Integrated luminosity -------------------------#
#                                                                          #
############################################################################
#Normalize to this luminsity, in fb-1
luminosity_norm_list = [35.86,41.529,59.69] #2016,2017,2018

#############---------------- BDT score cut values ----------------#############

BDT_OUT_MU  = 0.281
BDT_OUT_ELE = 0.269

############################################################################
#                                                                          #
#------------------- Get rootfiles and name of samples --------------------#
#                                                                          #
############################################################################

samplename_list = myWF.get_samples_names()
root_file = myWF.get_root_files()

############################################################################
#                                                                          #
#------------------------------ Tree variables ----------------------------#
#                                                                          #
############################################################################

#Variables for the FIT trees
_Wmass_fit          = np.zeros(1, dtype=float)
_isMuon_fit         = np.zeros(1, dtype=int)
_isSignal_fit       = np.zeros(1, dtype=int)
_weight_fit         = np.zeros(1, dtype=float)
_BDT_out_fit        = np.zeros(1, dtype=float)
_isSignalRegion_fit = np.zeros(1, dtype=int)
_Categorization_fit = np.zeros(1, dtype=int)

#Variables for the MVA trees
_isMuon             = np.zeros(1, dtype=int)
_gamma_eT           = np.zeros(1, dtype=float)
_pi_pT              = np.zeros(1, dtype=float)
_lep_pT             = np.zeros(1, dtype=float)
_lep_iso            = np.zeros(1, dtype=float)
_nBjets             = np.zeros(1, dtype=int)
_nBjets_25          = np.zeros(1, dtype=int)
_deltaphi_lep_pi    = np.zeros(1, dtype=float)
_deltaphi_lep_gamma = np.zeros(1, dtype=float)
_piRelIso_05        = np.zeros(1, dtype=float)
_piRelIso_05_ch     = np.zeros(1, dtype=float)
_pi_dxy             = np.zeros(1, dtype=float)
_ele_gamma_InvMass  = np.zeros(1, dtype=float)
_weight             = np.zeros(1, dtype=float)
_met                = np.zeros(1, dtype=float)
_met_puppi          = np.zeros(1, dtype=float)
_Wmass              = np.zeros(1, dtype=float)

_Nrandom_for_BDT_systematic = ROOT.TRandom3(24593)

nSig_mu  = 0.
nBkg_mu  = 0.
nSig_ele = 0.
nBkg_ele = 0.
nBkg_processed_mu  = 0.
nBkg_processed_ele = 0.

############################################################################
#                                                                          #
#------------------------ Create rootfiles and trees ----------------------#
#                                                                          #
############################################################################

if not isData:

    if create_mass_tree:
        f = TFile('WmassAnalysis/Tree_input_massfit_MC_' + str(runningEra) + '.root','recreate')        
        t = TTree('minitree','tree with branches')
        t.Branch('Wmass',_Wmass_fit,'Wmass/D')
        t.Branch('isMuon',_isMuon_fit,'isMuon/I')
        t.Branch('isSignal',_isSignal_fit,'isSignal/I')
        t.Branch('weight',_weight_fit,'weight/D')
        t.Branch('BDT_out',_BDT_out_fit,'BDT_out/D')
        t.Branch('isSignalRegion',_isSignalRegion_fit,'isSignalRegion/I')
        t.Branch('Categorization',_Categorization_fit,'Categorization/I')

    else:
        fMVA_signal_mu = TFile('MVA/Tree_MC_Signal_mu_' + str(runningEra) + '.root','recreate')
        tMVA_signal_mu = TTree('minitree_signal_mu','tree with branches')
        tMVA_signal_mu.Branch('isMuon',_isMuon,'isMuon/I')
        tMVA_signal_mu.Branch('weight',_weight,'weight/D')
        tMVA_signal_mu.Branch('gamma_eT',_gamma_eT,'gamma_eT/D')
        tMVA_signal_mu.Branch('pi_pT',_pi_pT,'pi_pT/D')
        tMVA_signal_mu.Branch('lep_pT',_lep_pT,'lep_pT/D')
        tMVA_signal_mu.Branch('lep_iso',_lep_iso,'lep_iso/D')
        tMVA_signal_mu.Branch('nBjets',_nBjets,'nBjets/I')
        tMVA_signal_mu.Branch('nBjets_25',_nBjets_25,'nBjets_25/I')
        tMVA_signal_mu.Branch('deltaphi_lep_pi',_deltaphi_lep_pi,'deltaphi_lep_pi/D')
        tMVA_signal_mu.Branch('deltaphi_lep_gamma',_deltaphi_lep_gamma,'deltaphi_lep_gamma/D')
        tMVA_signal_mu.Branch('piRelIso_05_ch',_piRelIso_05_ch,'piRelIso_05_ch/D')
        tMVA_signal_mu.Branch('piRelIso_05',_piRelIso_05,'piRelIso_05/D')
        tMVA_signal_mu.Branch('pi_dxy',_pi_dxy,'pi_dxy/D')
        tMVA_signal_mu.Branch('MET',_met,'MET/D')
        tMVA_signal_mu.Branch('MET_puppi',_met_puppi,'MET_puppi/D')
        tMVA_signal_mu.Branch('Wmass',_Wmass,'Wmass/D')
        
        fMVA_signal_ele = TFile('MVA/Tree_MC_Signal_ele_' + str(runningEra) + '.root','recreate')
        tMVA_signal_ele = TTree('minitree_signal_ele','tree with branches')
        tMVA_signal_ele.Branch('isMuon',_isMuon,'isMuon/I')
        tMVA_signal_ele.Branch('weight',_weight,'weight/D')
        tMVA_signal_ele.Branch('gamma_eT',_gamma_eT,'gamma_eT/D')
        tMVA_signal_ele.Branch('pi_pT',_pi_pT,'pi_pT/D')
        tMVA_signal_ele.Branch('lep_pT',_lep_pT,'lep_pT/D')
        tMVA_signal_ele.Branch('lep_iso',_lep_iso,'lep_iso/D')
        tMVA_signal_ele.Branch('nBjets',_nBjets,'nBjets/I')
        tMVA_signal_ele.Branch('nBjets_25',_nBjets_25,'nBjets_25/I')
        tMVA_signal_ele.Branch('deltaphi_lep_pi',_deltaphi_lep_pi,'deltaphi_lep_pi/D')
        tMVA_signal_ele.Branch('deltaphi_lep_gamma',_deltaphi_lep_gamma,'deltaphi_lep_gamma/D')
        tMVA_signal_ele.Branch('piRelIso_05_ch',_piRelIso_05_ch,'piRelIso_05_ch/D')
        tMVA_signal_ele.Branch('piRelIso_05',_piRelIso_05,'piRelIso_05/D')
        tMVA_signal_ele.Branch('pi_dxy',_pi_dxy,'pi_dxy/D')
        tMVA_signal_ele.Branch('MET',_met,'MET/D')
        tMVA_signal_ele.Branch('MET_puppi',_met_puppi,'MET_puppi/D')
        tMVA_signal_ele.Branch('Wmass',_Wmass,'Wmass/D')

        fMVA_background_mu = TFile('MVA/Tree_MC_Background_mu_' + str(runningEra) + '.root','recreate')
        tMVA_background_mu = TTree('minitree_background_mu','tree with branches')
        tMVA_background_mu.Branch('isMuon',_isMuon,'isMuon/I')
        tMVA_background_mu.Branch('weight',_weight,'weight/D')
        tMVA_background_mu.Branch('gamma_eT',_gamma_eT,'gamma_eT/D')
        tMVA_background_mu.Branch('pi_pT',_pi_pT,'pi_pT/D')
        tMVA_background_mu.Branch('lep_pT',_lep_pT,'lep_pT/D')
        tMVA_background_mu.Branch('lep_iso',_lep_iso,'lep_iso/D')
        tMVA_background_mu.Branch('nBjets',_nBjets,'nBjets/I')
        tMVA_background_mu.Branch('nBjets_25',_nBjets_25,'nBjets_25/I')
        tMVA_background_mu.Branch('deltaphi_lep_pi',_deltaphi_lep_pi,'deltaphi_lep_pi/D')
        tMVA_background_mu.Branch('deltaphi_lep_gamma',_deltaphi_lep_gamma,'deltaphi_lep_gamma/D')
        tMVA_background_mu.Branch('piRelIso_05_ch',_piRelIso_05_ch,'piRelIso_05_ch/D')
        tMVA_background_mu.Branch('piRelIso_05',_piRelIso_05,'piRelIso_05/D')
        tMVA_background_mu.Branch('pi_dxy',_pi_dxy,'pi_dxy/D')
        tMVA_background_mu.Branch('MET',_met,'MET/D')
        tMVA_background_mu.Branch('MET_puppi',_met_puppi,'MET_puppi/D')
        tMVA_background_mu.Branch('Wmass',_Wmass,'Wmass/D')
        
        fMVA_background_ele = TFile('MVA/Tree_MC_Background_ele_' + str(runningEra) + '.root','recreate')
        tMVA_background_ele = TTree('minitree_background_ele','tree with branches')
        tMVA_background_ele.Branch('isMuon',_isMuon,'isMuon/I')
        tMVA_background_ele.Branch('weight',_weight,'weight/D')
        tMVA_background_ele.Branch('gamma_eT',_gamma_eT,'gamma_eT/D')
        tMVA_background_ele.Branch('pi_pT',_pi_pT,'pi_pT/D')
        tMVA_background_ele.Branch('lep_pT',_lep_pT,'lep_pT/D')
        tMVA_background_ele.Branch('lep_iso',_lep_iso,'lep_iso/D')
        tMVA_background_ele.Branch('nBjets',_nBjets,'nBjets/I')
        tMVA_background_ele.Branch('nBjets_25',_nBjets_25,'nBjets_25/I')
        tMVA_background_ele.Branch('deltaphi_lep_pi',_deltaphi_lep_pi,'deltaphi_lep_pi/D')
        tMVA_background_ele.Branch('deltaphi_lep_gamma',_deltaphi_lep_gamma,'deltaphi_lep_gamma/D')
        tMVA_background_ele.Branch('piRelIso_05_ch',_piRelIso_05_ch,'piRelIso_05_ch/D')
        tMVA_background_ele.Branch('piRelIso_05',_piRelIso_05,'piRelIso_05/D')
        tMVA_background_ele.Branch('pi_dxy',_pi_dxy,'pi_dxy/D')
        tMVA_background_ele.Branch('MET',_met,'MET/D')
        tMVA_background_ele.Branch('MET_puppi',_met_puppi,'MET_puppi/D')
        tMVA_background_ele.Branch('Wmass',_Wmass,'Wmass/D')


if isData:

    if create_mass_tree:
        f = TFile('WmassAnalysis/Tree_input_massfit_Data_' + str(runningEra) + '.root','recreate')
        t = TTree('minitree','tree with branches')
        t.Branch('Wmass',_Wmass_fit,'Wmass/D')
        t.Branch('isMuon',_isMuon_fit,'isMuon/I')
        t.Branch('BDT_out',_BDT_out_fit,'BDT_out/D')
        t.Branch('Categorization',_Categorization_fit,'Categorization/I')

    else:
        fMVA_background_mu_DATA = TFile('MVA/Tree_Background_mu_DATA.root','recreate')
        tMVA_background_mu_DATA = TTree('minitree_background_mu_DATA','tree with branches')
        tMVA_background_mu_DATA.Branch('isMuon',_isMuon,'isMuon/I')
        tMVA_background_mu_DATA.Branch('weight',_weight,'weight/D')
        tMVA_background_mu_DATA.Branch('gamma_eT',_gamma_eT,'gamma_eT/D')
        tMVA_background_mu_DATA.Branch('pi_pT',_pi_pT,'pi_pT/D')
        tMVA_background_mu_DATA.Branch('lep_pT',_lep_pT,'lep_pT/D')
        tMVA_background_mu_DATA.Branch('lep_iso',_lep_iso,'lep_iso/D')
        tMVA_background_mu_DATA.Branch('nBjets',_nBjets,'nBjets/I')
        tMVA_background_mu_DATA.Branch('nBjets_25',_nBjets_25,'nBjets_25/I')
        tMVA_background_mu_DATA.Branch('deltaphi_lep_pi',_deltaphi_lep_pi,'deltaphi_lep_pi/D')
        tMVA_background_mu_DATA.Branch('deltaphi_lep_gamma',_deltaphi_lep_gamma,'deltaphi_lep_gamma/D')
        tMVA_background_mu_DATA.Branch('piRelIso_05_ch',_piRelIso_05_ch,'piRelIso_05_ch/D')
        tMVA_background_mu_DATA.Branch('piRelIso_05',_piRelIso_05,'piRelIso_05/D')
        tMVA_background_mu_DATA.Branch('pi_dxy',_pi_dxy,'pi_dxy/D')
        tMVA_background_mu_DATA.Branch('MET',_met,'MET/D')
        tMVA_background_mu_DATA.Branch('MET_puppi',_met_puppi,'MET_puppi/D')
        tMVA_background_mu_DATA.Branch('Wmass',_Wmass,'Wmass/D')
        
        fMVA_background_ele_DATA = TFile('MVA/Tree_Background_ele_DATA.root','recreate')
        tMVA_background_ele_DATA = TTree('minitree_background_ele_DATA','tree with branches')
        tMVA_background_ele_DATA.Branch('isMuon',_isMuon,'isMuon/I')
        tMVA_background_ele_DATA.Branch('weight',_weight,'weight/D')
        tMVA_background_ele_DATA.Branch('gamma_eT',_gamma_eT,'gamma_eT/D')
        tMVA_background_ele_DATA.Branch('pi_pT',_pi_pT,'pi_pT/D')
        tMVA_background_ele_DATA.Branch('lep_pT',_lep_pT,'lep_pT/D')
        tMVA_background_ele_DATA.Branch('lep_iso',_lep_iso,'lep_iso/D')
        tMVA_background_ele_DATA.Branch('nBjets',_nBjets,'nBjets/I')
        tMVA_background_ele_DATA.Branch('nBjets_25',_nBjets_25,'nBjets_25/I')
        tMVA_background_ele_DATA.Branch('deltaphi_lep_pi',_deltaphi_lep_pi,'deltaphi_lep_pi/D')
        tMVA_background_ele_DATA.Branch('deltaphi_lep_gamma',_deltaphi_lep_gamma,'deltaphi_lep_gamma/D')
        tMVA_background_ele_DATA.Branch('piRelIso_05_ch',_piRelIso_05_ch,'piRelIso_05_ch/D')
        tMVA_background_ele_DATA.Branch('piRelIso_05',_piRelIso_05,'piRelIso_05/D')
        tMVA_background_ele_DATA.Branch('pi_dxy',_pi_dxy,'pi_dxy/D')
        tMVA_background_ele_DATA.Branch('MET',_met,'MET/D')
        tMVA_background_ele_DATA.Branch('MET_puppi',_met_puppi,'MET_puppi/D')
        tMVA_background_ele_DATA.Branch('Wmass',_Wmass,'Wmass/D')

## ttbar-signal calibration ##
ttbar_sig_calib = myWF.get_ttbar_signal_reweight()


for full_sample_name in samplename_list:

    sample_name = full_sample_name.split("_")[0]
    sample_era  = full_sample_name.split("_")[1]
    if sample_era == "2016":
        sample_era = 0
    if sample_era == "2017":
        sample_era = 1
    if sample_era == "2018":
        sample_era = 2

    print "full_sample_name: ", full_sample_name, "sample_era: ", sample_era

    if runningEra == 0 and not sample_era == 0:
        continue
    if runningEra == 1 and not sample_era == 1:
        continue
    if runningEra == 2 and not sample_era == 2:
        continue

    Norm_Map = myWF.get_normalizations_map(sample_era)
    luminosity_norm = luminosity_norm_list[sample_era]
    print "lumi norm: ", luminosity_norm

    if not "Data" in sample_name:
        norm_factor = Norm_Map[sample_name]*luminosity_norm

    mytree = root_file[full_sample_name].Get("WPiGammaAnalysis/mytree")
    #print "root_file[sample_name]: ", root_file[full_sample_name]

    sample_entries = mytree.GetEntriesFast() 
    print "Processing Sample ", sample_name, ", number of events: ", sample_entries
    SignalTree_entries = 0
    entry_index = 0.

    #nEvts = mytree.GetEntriesFast() #Get the number of events per sample

    for jentry in xrange(mytree.GetEntriesFast()):

        ientry = mytree.LoadTree( jentry )
        if ientry < 0:
            break
        nb = mytree.GetEntry(jentry )
        if nb <= 0:
            continue

        ############################################################################
        #                                                                          #
        #-------------------------- Samples to be excluded ------------------------#
        #                                                                          #
        ############################################################################

        if isData and not "Data" in sample_name: 
            continue
        
        if not isData and "Data" in sample_name:
            continue

        if mytree.is_muon and sample_name == "QCDDoubleEMEnriched30toInf":
            continue

        if not mytree.is_muon and sample_name == "QCDHT200to300":
            continue

        ############################################################################
        #                                                                          #
        #------------------------ Access the tree variables -----------------------#
        #                                                                          #
        ############################################################################

        entry_index += 1.
        
        #if entry_index > sample_entries/2.:
        #    continue

        # Print entry_index

        isMuon = mytree.is_muon

        isTriggerMatched = mytree.isTriggerMatched
        isSingleMuTrigger_24 = mytree.isSingleMuTrigger_24
        isSingleMuTrigger_50 = mytree.isSingleMuTrigger_50
        isSingleMuTrigger_27 = mytree.isSingleMuTrigger_27
        isSingleEleTrigger_32 = mytree.isSingleEleTrigger_32
        isSingleEleTrigger_32_DoubleEG = mytree.isSingleEleTrigger_32_DoubleEG
            
        lep_pT  = mytree.lepton_pT
        lep_eta = mytree.lepton_eta
        lep_phi = mytree.lepton_phi
        lep_iso = mytree.lepton_iso
        ele_etaSC = mytree.lepton_etaSC
        lep_FourMomentum = ROOT.TLorentzVector()
        lep_FourMomentum.SetPtEtaPhiM(lep_pT,lep_eta,lep_phi,0.)

        pi_pT = mytree.pi_pT
        pi_eta = mytree.pi_eta
        pi_phi = mytree.pi_phi
        pi_E = mytree.pi_energy
        pi_FourMomentum = ROOT.TLorentzVector()
        pi_FourMomentum.SetPtEtaPhiE(pi_pT,pi_eta,pi_phi,pi_E)
        piRelIso_05 = mytree.sum_pT_05/pi_pT
        piRelIso_05_ch = mytree.sum_pT_05_ch/pi_pT
        pi_dxy = mytree.pi_dxy
            
        gamma_eT = mytree.photon_eT
        gamma_eta = mytree.photon_eta
        gamma_etaSC = mytree.photon_etaSC
        gamma_phi = mytree.photon_phi
        gamma_E = mytree.photon_energy
        gamma_FourMomentum = ROOT.TLorentzVector()
        gamma_FourMomentum.SetPtEtaPhiE(gamma_eT,gamma_eta,gamma_phi,gamma_E)

        W_phi = (pi_FourMomentum + gamma_FourMomentum).Phi()
        W_pT  = (pi_FourMomentum + gamma_FourMomentum).Pt()

        if not "Data" in sample_name:
            Wplus_pT = mytree.Wplus_pT
            Wminus_pT = mytree.Wminus_pT

        Wmass = mytree.Wmass

        met = mytree.met_pT
        met_puppi = mytree.metpuppi_pT

        if not isMuon:
            ele_gamma_InvMass = (lep_FourMomentum + gamma_FourMomentum).M()
        else:
            mu_gamma_InvMass = (lep_FourMomentum + gamma_FourMomentum).M()
        
        nBjets_25 = mytree.nBjets_25

        deltaeta_lep_pi = math.fabs(lep_eta-pi_eta)

        deltaphi_lep_pi = math.fabs(lep_phi - pi_phi)
        if deltaphi_lep_pi > 3.14:
            deltaphi_lep_pi = 6.28 - deltaphi_lep_pi

        deltaphi_lep_W = math.fabs(lep_phi - W_phi)
        if deltaphi_lep_W > 3.14:
            deltaphi_lep_W = 6.28 - deltaphi_lep_W  
            
        deltaphi_lep_W = math.fabs(lep_phi - W_phi)
        if deltaphi_lep_W > 3.14:
            deltaphi_lep_W = 6.28 - deltaphi_lep_W

        deltaphi_lep_gamma = math.fabs(lep_phi - gamma_phi)
        if deltaphi_lep_gamma > 3.14:
            deltaphi_lep_gamma = 6.28 - deltaphi_lep_gamma

        deltaeta_lep_gamma = math.fabs(lep_eta - gamma_eta)
        deltaR_lep_gamma = math.sqrt(deltaphi_lep_gamma*deltaphi_lep_gamma + deltaeta_lep_gamma*deltaeta_lep_gamma)

        if "WJetsToLNu" in sample_name:
            
            if isMuon:
                MCT_lep_eta = mytree.MCT_HpT_mu_eta
                MCT_lep_phi = mytree.MCT_HpT_mu_phi
                MCT_lep_pT  = mytree.MCT_HpT_mu_pT
            else:
                MCT_lep_eta = mytree.MCT_HpT_ele_eta
                MCT_lep_phi = mytree.MCT_HpT_ele_phi
                MCT_lep_pT  = mytree.MCT_HpT_ele_pT
                
            MCT_deltaeta_lep_gamma = math.fabs(MCT_lep_eta - mytree.MCT_HeT_ph_eta)
            MCT_deltaphi_lep_gamma = math.fabs(MCT_lep_phi - mytree.MCT_HeT_ph_phi)

            if MCT_deltaphi_lep_gamma > 3.14:
                MCT_deltaphi_lep_gamma = 6.28 - MCT_deltaphi_lep_gamma

            MCT_deltaR_lep_gamma = math.sqrt(MCT_deltaeta_lep_gamma*MCT_deltaeta_lep_gamma + MCT_deltaphi_lep_gamma*MCT_deltaphi_lep_gamma)

            if MCT_deltaR_lep_gamma > 0.2 and not MCT_lep_pT < 0. and not mytree.MCT_HeT_ph_eT < 0.:
                continue

        if sample_name == "Signal":
            is_signal_Wplus = mytree.is_signal_Wplus
            is_signal_Wminus = mytree.is_signal_Wminus
            genPi_pT  = mytree.gen_pi_pT 
            genPi_eta = mytree.gen_pi_eta 
            genPi_phi = mytree.gen_pi_phi
            genPi_energy = mytree.gen_pi_energy
            genPh_pT  = mytree.gen_ph_pT 
            genPh_eta = mytree.gen_ph_eta
            genPh_phi = mytree.gen_ph_phi
            genPh_energy = mytree.gen_ph_energy
            genTop_pT  = mytree.gen_pi_grandmother_pT 
            genTop_eta = mytree.gen_pi_grandmother_eta 
            genTop_phi = mytree.gen_pi_grandmother_phi
            genTop_energy = mytree.gen_pi_grandmother_energy
        
            genTop_FourMomentum = ROOT.TLorentzVector()
            genTop_FourMomentum.SetPtEtaPhiE(genTop_pT,genTop_eta,genTop_phi,genTop_energy)
            genPi_FourMomentum = ROOT.TLorentzVector()
            genPi_FourMomentum.SetPtEtaPhiE(genPi_pT,genPi_eta,genPi_phi,genPi_energy)
            genPh_FourMomentum = ROOT.TLorentzVector()
            genPh_FourMomentum.SetPtEtaPhiE(genPh_pT,genPh_eta,genPh_phi,genPh_energy)
        
            genW_pT      = (genPi_FourMomentum + genPh_FourMomentum).Pt()
            genW_eta     = (genPi_FourMomentum + genPh_FourMomentum).Eta()
            genW_phi     = (genPi_FourMomentum + genPh_FourMomentum).Phi()
            genW_energy  = (genPi_FourMomentum + genPh_FourMomentum).E()
            genW_FourMomentum = ROOT.TLorentzVector()
            genW_FourMomentum.SetPtEtaPhiE(genW_pT,genW_eta,genW_phi,genW_energy)

            genW_ThreeMomentum = ROOT.TVector3()
            genW_ThreeMomentum = genW_FourMomentum.Vect()

            boosted_genTop_FourMomentum = ROOT.TLorentzVector()
            boosted_genTop_FourMomentum.SetPtEtaPhiE(genTop_pT,genTop_eta,genTop_phi,genTop_energy)
            boosted_genTop_FourMomentum.Boost(-genW_FourMomentum.BoostVector()) #The sign - is required
            boosted_genPi_FourMomentum = ROOT.TLorentzVector()
            boosted_genPi_FourMomentum.SetPtEtaPhiE(genPi_pT,genPi_eta,genPi_phi,genPi_energy)
            boosted_genPi_FourMomentum.Boost(-genW_FourMomentum.BoostVector()) #The sign - is required

            angle_Top_Pi = boosted_genTop_FourMomentum.Angle(boosted_genPi_FourMomentum.Vect())
            angle_Pi_W = boosted_genPi_FourMomentum.Angle(genW_FourMomentum.Vect())

        ############################################################################
        #                                                                          #
        #----------------------- Some post-preselection cuts ----------------------#
        #                                                                          #
        ############################################################################

        if myWF.post_preselection_cuts(lep_eta,lep_pT,ele_etaSC,gamma_etaSC,isMuon,deltaphi_lep_pi,deltaphi_lep_gamma,isTriggerMatched,sample_era):
            continue

        ############################################################################
        #                                                                          #
        #--------------------- Determine the total event weight -------------------#
        #                                                                          #
        ############################################################################

        ################ MUON SFs ################

        if isMuon: # Get muon scale factors, which are different for two groups of datasets, and weight them for the respective integrated lumi 
            if sample_era == 0 or sample_era == 2:
                isSingleMuTrigger_LOW = isSingleMuTrigger_24
            if sample_era == 1:
                isSingleMuTrigger_LOW = isSingleMuTrigger_27

            lep_weight, lep_weight_err = myWF.get_muon_scale(lep_pT,lep_eta,isSingleMuTrigger_LOW,sample_era)
            
        ############## ELECTRON SFs ##############

        else:
            lep_weight, lep_weight_err = myWF.get_ele_scale(lep_pT,ele_etaSC,sample_era)
            
        ############### PHOTON SFs ###############

        ph_weight, ph_weight_err = myWF.get_photon_scale(gamma_eT,gamma_etaSC,sample_era)

        ############### Multiply weights and SFs for MC. Set weight to 1 for data ###############

        if not "Data" in sample_name:
            MC_Weight = mytree.MC_Weight # Add MC weight        
            PU_Weight = mytree.PU_Weight # Add Pile Up weight

            if not runningEra == 2: # Prefiring weight not to be applied to 2018 MC
                Prefiring_Weight = mytree.Prefiring_Weight # Add prefiring weight

                if only_signal_relevant_weights: #Suppress non pt/eta-dependent weights
                    Event_Weight = ph_weight*PU_Weight*Prefiring_Weight
                else:
                    Event_Weight = norm_factor*ph_weight*MC_Weight*PU_Weight*Prefiring_Weight/math.fabs(MC_Weight) # Just take the sign of the gen weight

            else:
                if only_signal_relevant_weights: #Suppress non pt/eta-dependent weights
                    Event_Weight = ph_weight*PU_Weight
                else:
                    Event_Weight = norm_factor*ph_weight*MC_Weight*PU_Weight/math.fabs(MC_Weight) # Just take the sign of the gen weight                

            Event_Weight = Event_Weight*lep_weight

            ################ Zvtx inefficiency weight for 2017 MC ################
            #https://twiki.cern.ch/twiki/bin/view/CMS/Egamma2017DataRecommendations

            if sample_era == 1:
                if isSingleEleTrigger_32_DoubleEG and not isSingleMuTrigger_27 and not isSingleMuTrigger_50:
                    Event_Weight = Event_Weight*0.991 #uniform penalty for all the 2017 eras
                    
            ################ Correct for the difference in pT of the generated W in Pythia and Madgraph samples ################

            if sample_name == "Signal":
                if is_signal_Wplus:
                    local_W_pT = Wplus_pT
                else:
                    local_W_pT = Wminus_pT

                if local_W_pT > 300.:
                    local_W_pT = 300.

                Event_Weight = Event_Weight*ttbar_sig_calib.GetBinContent(ttbar_sig_calib.GetXaxis().FindBin(local_W_pT))

                ################ Create up- and down-scalded event weights for Pythia signal modeling ################

                Event_Weight_up   = Event_Weight*myWF.get_Pythia_pT_modeling_reweight(True,genW_pT)
                Event_Weight_down = Event_Weight*myWF.get_Pythia_pT_modeling_reweight(False,genW_pT)
                Event_Weight_sin2 = Event_Weight*myWF.get_Pythia_polarization_modeling_reweight(True,angle_Pi_W)
                Event_Weight_cos  = Event_Weight*myWF.get_Pythia_polarization_modeling_reweight(False,angle_Pi_W)

            #####################################################################################################

            if "Signal" in sample_name and isMuon:
                nSig_mu += Event_Weight
            if not "Signal" in sample_name and isMuon:
                nBkg_mu += Event_Weight
                nBkg_processed_mu  += 1.
            if "Signal" in sample_name and not isMuon:
                nSig_ele += Event_Weight
            if not "Signal" in sample_name and not isMuon:
                nBkg_ele += Event_Weight
                nBkg_processed_ele  += 1.

        else:
            Event_Weight = 1.

        # Remove QCD events with abnormal weight
        if not create_mass_tree and "QCD" in sample_name and Event_Weight >= 1600.:
            continue
        if create_mass_tree and "QCD" in sample_name and Event_Weight >= 30.:
            continue


        ############################################################################
        #                                                                          #
        #------------------------------ Fill mass tree ----------------------------#
        #                                                                          #
        ############################################################################


        #---------Retrieve the BDT output----------#

        if create_mass_tree:
            BDT_out = myWF.get_BDT_output(pi_pT,gamma_eT,nBjets_25,lep_pT,piRelIso_05_ch,met,isMuon)

        #---------- Fill the tree ----------#

        if create_mass_tree:

            #if ((Wmass >=50. and Wmass<=65.) or (Wmass>=90. and Wmass<=100.)):
            if (Wmass >= 50. and Wmass <= 100.): 
                _isMuon_fit[0]  = isMuon
                _Wmass_fit[0]   = Wmass
                _BDT_out_fit[0] = BDT_out
            

                if (isMuon and BDT_out >= BDT_OUT_MU) or ( (not isMuon) and BDT_out >= BDT_OUT_ELE):
                    _isSignalRegion_fit[0] = 1
                else:
                    _isSignalRegion_fit[0] = 0

                #------------------------ Fill Tree ------------------------#
                if isMuon and BDT_out < BDT_OUT_MU:
                    _Categorization_fit[0] = 0 # Control Region muon channel = 0
                if isMuon and BDT_out >= BDT_OUT_MU:
                    _Categorization_fit[0] = 1 # Signal Region muon channel = 1
                    if not sample_name=="Signal":
                        if not Wmass == _Wmass_fit[0]:
                            print "WARNING!!! mu channel"

                if not isMuon and BDT_out < BDT_OUT_ELE:
                    _Categorization_fit[0] = 2 # Control Region electron channel = 2
                if not isMuon and BDT_out >= BDT_OUT_ELE:
                    _Categorization_fit[0] = 3 # Signal Region electron channel = 3
                    if not sample_name=="Signal":
                        if not Wmass == _Wmass_fit[0]:
                            print "WARNING!!! ele channel"

                if not isData:
                    _weight_fit[0] = Event_Weight
                    if sample_name == "Signal":
                        _isSignal_fit[0] = 1
                    else:
                        _isSignal_fit[0] = 0
  
                t.Fill() #Filling the tree for FIT



            ############################################################################
            #                                                                          #
            #------------------------------- Fill MVA tree ----------------------------#
            #                                                                          #
            ############################################################################

        else:

            if not isData:
                
                _gamma_eT[0]           = gamma_eT
                _pi_pT[0]              = pi_pT
                _lep_pT[0]             = lep_pT
                #_lep_iso[0]            = lep_iso
                _nBjets_25[0]          = nBjets_25
                #_deltaphi_lep_pi[0]    = deltaphi_lep_pi
                #_deltaphi_lep_gamma[0] = deltaphi_lep_gamma
                _isMuon[0]             = isMuon
                #_piRelIso_05[0]        = piRelIso_05
                _piRelIso_05_ch[0]     = piRelIso_05_ch
                _pi_dxy[0]             = pi_dxy
                _met[0]                = met
                _met_puppi[0]          = met_puppi
                _Wmass[0]              = Wmass
                if sample_name == "Signal" and scale_signal_up:
                    _weight[0]         = Event_Weight_up
                elif sample_name == "Signal" and scale_signal_down:
                    _weight[0]         = Event_Weight_down
                elif sample_name == "Signal" and scale_signal_sin2:
                    _weight[0]         = Event_Weight_sin2
                elif sample_name == "Signal" and scale_signal_cos:
                    _weight[0]         = Event_Weight_cos
                else:
                    _weight[0]         = Event_Weight
                #-----SCALED VARIABLES-----#
                # _gamma_eT[0]           = _Nrandom_for_BDT_systematic.Gaus(gamma_eT,gamma_eT*0.05)
                # _pi_pT[0]              = _Nrandom_for_BDT_systematic.Gaus(pi_pT,pi_pT*0.05)
                # _lep_pT[0]             = _Nrandom_for_BDT_systematic.Gaus(lep_pT,lep_pT*0.05)
                # _piRelIso_05_ch[0]     = _Nrandom_for_BDT_systematic.Gaus(piRelIso_05_ch,piRelIso_05_ch*0.1)
                # _met[0]                = _Nrandom_for_BDT_systematic.Gaus(met,met*0.05)
                # _Nrandom_for_BDT_systematic_nBjets_25 = _Nrandom_for_BDT_systematic.Rndm()
                # if _Nrandom_for_BDT_systematic_nBjets_25 <= 0.90:
                #     _nBjets_25[0]          = nBjets_25
                # elif (_Nrandom_for_BDT_systematic_nBjets_25 > 0.90 and _Nrandom_for_BDT_systematic_nBjets_25 <= 0.95):
                #     _nBjets_25[0]          = nBjets_25-1
                # else:
                #     _nBjets_25[0]          = nBjets_25+1
                
                if sample_name == "Signal" and isMuon and (not split_MC):
                    tMVA_signal_mu.Fill()
                if sample_name == "Signal" and (not isMuon) and (not split_MC):
                    tMVA_signal_ele.Fill()
                if (not sample_name == "Signal") and isMuon:
                    tMVA_background_mu.Fill()
                if (not sample_name == "Signal") and (not isMuon):
                    tMVA_background_ele.Fill()

            if isData:

                _gamma_eT[0]        = gamma_eT
                _pi_pT[0]           = pi_pT
                _lep_pT[0]          = lep_pT
                _lep_iso[0]         = lep_iso
                _nBjets_25[0]       = nBjets_25
                _weight[0]          = 1
                _deltaphi_lep_pi[0] = deltaphi_lep_pi
                _isMuon[0]          = isMuon
                _piRelIso_05[0]     = piRelIso_05
                _piRelIso_05_ch[0]  = piRelIso_05_ch
                _pi_dxy[0]          = pi_dxy
                _met[0]             = met
                _met_puppi[0]       = met_puppi
                _Wmass[0]           = Wmass
                
                if sample_name == "Data" and isMuon:
                    tMVA_background_mu_DATA.Fill()
                if sample_name == "Data" and not isMuon:
                    tMVA_background_ele_DATA.Fill()

print "Finished runnning over samples!"

############################################################################
#                                                                          #
#------------------------------- Write trees ------------------------------#
#                                                                          #
############################################################################

if create_mass_tree:
    f.cd()
    t.Write()
    f.Close()

else:
    if not isData :
                    
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

print "nSig_mu: ", nSig_mu, "nBkg_mu: ", nBkg_mu 
print "nSig_ele: ", nSig_ele, "nBkg_ele: ", nBkg_ele 
print "nBkg_processed_mu: ", nBkg_processed_mu, "nBkg_processed_ele: ", nBkg_processed_ele
