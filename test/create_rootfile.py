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

isBDT_with_Wmass = False # If true, pT(pi) and ET(gamma) in the BDT are normalized to Wmass 
isData   = False # Switch from DATA to MC and vice versa
split_MC = False # If True, MC signal sample is split in two for the training/testing of the BDT

myWF = Simplified_Workflow_Handler("Signal","Data",create_mass_tree,isBDT_with_Wmass,runningEra)

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

############################################################################
#                                                                          #
#-------------------------- Integrated luminosity -------------------------#
#                                                                          #
############################################################################
#Normalize to this luminsity, in fb-1
luminosity_norm_list = [35.86,41.529,59.69] #2016,2017,2018

luminosity_norm_2017_Ele32_WPTight = 27.13
_Nrandom_for_Ele_32_WPTight_exclusion = ROOT.TRandom3(64524)

#############---------------- BDT score cut values ----------------#############

BDT_OUT_MU  = 0.220
BDT_OUT_ELE = 0.180

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

#_Nrandom_for_SF = ROOT.TRandom3(44317)

Wmass_mu   = ROOT.TH1F("Wmass_mu","Wmass mu",15,50,100)
Wmass_ele  = ROOT.TH1F("Wmass_ele","Wmass ele",15,50,100)
lep_pt_mu  = ROOT.TH1F("lep_pt_mu","lep pt mu",15,25,100)
lep_pt_ele = ROOT.TH1F("lep_pt_ele","lep pt ele",15,28,100)

Wmass_mu_WJets  = ROOT.TH1F("Wmass_mu_WJets","Wmass mu WJets",15,50,100)
Wmass_ele_WJets = ROOT.TH1F("Wmass_ele_WJets","Wmass ele WJets",15,50,100)

N_WGToLNuG_mu = 0.

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


def select_all_but_one(cutstring):

    selection_bools = dict()
    if isMuon:
        selection_bools["mupt"]                = lep_pT >= MU_MIN_PT
        selection_bools["deltaphi_mu_pi"]      = deltaphi_lep_pi >= DELTAPHI_MU_PI_MIN
    if not isMuon:
        selection_bools["elept"]               = lep_pT >= ELE_MIN_PT
        selection_bools["deltaphi_ele_pi"]     = deltaphi_lep_pi >= DELTAPHI_ELE_PI_MIN
        selection_bools["h_ele_iso"]           = lep_iso <= ELE_ISO_MAX
        selection_bools["h_ele_gamma_InvMass"] = (ele_gamma_InvMass < ELE_GAMMA_INVMASS_MIN or ele_gamma_InvMass > ELE_GAMMA_INVMASS_MAX)
    selection_bools["pipt"]                    = pi_pT >= PI_MIN_PT
    selection_bools["gammaet"]                 = gamma_eT >= GAMMA_MIN_ET
    selection_bools["nBjets"]                  = nBjets_25 >= N_BJETS_MIN
    selection_bools["Wmass"]                   = (Wmass >= WMASS_MIN and Wmass <= WMASS_MAX)
    result = True

    for hname in selection_bools:
        if cutstring == hname:
            continue
        else:
            result = result and selection_bools[hname]

    return result

## ttbar-signal calibration ##
ttbar_sig_calib = myWF.get_ttbar_signal_reweight()


SignalTree_entries = 0
entry_index = 0

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
 
    print "Processing Sample ", sample_name

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
        
        #if runningEra == 0 and sample_name == "ttbar" and mytree.isttbarlnu: # Avoid double-counting of the ttbarlnu background
        #    continue

        if mytree.is_muon and sample_name == "QCDDoubleEMEnriched30toInf":
            continue

        if not mytree.is_muon and sample_name == "QCDHT200to300":
            continue

        ############################################################################
        #                                                                          #
        #------------------------ Access the tree variables -----------------------#
        #                                                                          #
        ############################################################################

        if sample_name == "Signal":
            SignalTree_entries = mytree.GetEntriesFast()
            entry_index += 1

        isMuon = mytree.is_muon
        LepPiOppositeCharge = mytree.LepPiOppositeCharge

        isSingleMuTrigger_24 = mytree.isSingleMuTrigger_24
        isSingleMuTrigger_50 = mytree.isSingleMuTrigger_50

        if sample_era == 1:
            isSingleMuTrigger_27 = mytree.isSingleMuTrigger_27
            isSingleEleTrigger_32 = mytree.isSingleEleTrigger_32
            isSingleEleTrigger_32_DoubleEG = mytree.isSingleEleTrigger_32_DoubleEG
            
            if sample_name == "Data":
                run_number = mytree.run_number
                #Use only Ele32_WPTight trigger for the period it is on
                if run_number > 302026 and not isSingleMuTrigger_27 and not isSingleMuTrigger_50 and not isSingleEleTrigger_32:
                    continue
            else: #Use only Ele32_WPTight trigger for the fraction of luminosity it is on
                if not isSingleMuTrigger_27 and not isSingleMuTrigger_50:
                    if _Nrandom_for_Ele_32_WPTight_exclusion.Rndm() <= (luminosity_norm_2017_Ele32_WPTight/luminosity_norm):
                        if not isSingleEleTrigger_32:
                            continue

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

        gamma_iso_ChHad = mytree.photon_iso_ChargedHadron
        gamma_iso_NeuHad = mytree.photon_iso_NeutralHadron
        gamma_iso_Ph = mytree.photon_iso_Photon
        gamma_iso_eArho = mytree.photon_iso_eArho

        W_phi = (pi_FourMomentum + gamma_FourMomentum).Phi()

        if not "Data" in sample_name:
            Wplus_pT = mytree.Wplus_pT
            Wminus_pT = mytree.Wminus_pT

        if sample_name == "Signal":
            is_signal_Wplus = mytree.is_signal_Wplus
            is_signal_Wminus = mytree.is_signal_Wminus

        Wmass = mytree.Wmass

        met = mytree.met_pT
        met_puppi = mytree.metpuppi_pT

        if not isMuon:
            ele_etaSC = mytree.lepton_etaSC
            ele_gamma_InvMass = (lep_FourMomentum + gamma_FourMomentum).M()
        else:
            mu_gamma_InvMass = (lep_FourMomentum + gamma_FourMomentum).M()
        
        nBjets = mytree.nBjets
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

            if MCT_deltaR_lep_gamma > 0.5 and not MCT_lep_pT < 0. and not mytree.MCT_HeT_ph_eT < 0.:
                continue


        ############################################################################
        #                                                                          #
        #----------------------- Some post-preselection cuts ----------------------#
        #                                                                          #
        ############################################################################

        if myWF.post_preselection_cuts(lep_eta,lep_pT,isMuon,LepPiOppositeCharge,deltaphi_lep_gamma,sample_era):
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
                Event_Weight = norm_factor*ph_weight*MC_Weight*PU_Weight/math.fabs(MC_Weight)*Prefiring_Weight # Just take the sign of the gen weight
            else:
                Event_Weight = norm_factor*ph_weight*MC_Weight*PU_Weight/math.fabs(MC_Weight) # Just take the sign of the gen weight


            Event_Weight = Event_Weight*lep_weight

            ################ Zvtx inefficiency weight for 2017 MC ################
            #https://twiki.cern.ch/twiki/bin/view/CMS/Egamma2017DataRecommendations

            if sample_era == 1:
                if (isSingleEleTrigger_32_DoubleEG or isSingleEleTrigger_32) and not isSingleMuTrigger_27 and not isSingleMuTrigger_50:
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

            if "Signal" in sample_name and isMuon:
                nSig_mu += Event_Weight
            if not "Signal" in sample_name and isMuon:
                nBkg_mu += Event_Weight
                nBkg_processed_mu  += 1.
                # if "WJetsToLNu" in sample_name:
                #     Wmass_mu_WJets.Fill(Wmass,Event_Weight)
            if "Signal" in sample_name and not isMuon:
                nSig_ele += Event_Weight
            if not "Signal" in sample_name and not isMuon:
                nBkg_ele += Event_Weight
                nBkg_processed_ele  += 1.
                # if "WJetsToLNu" in sample_name:
                #     Wmass_ele_WJets.Fill(Wmass,Event_Weight)

        else:
            Event_Weight = 1.


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

                #------------------------ 2016 ------------------------#
                if sample_era == 0 and isMuon and BDT_out < BDT_OUT_MU:
                    _Categorization_fit[0] = 0 # Control Region muon channel = 0
                if sample_era == 0 and isMuon and BDT_out >= BDT_OUT_MU:
                    _Categorization_fit[0] = 1 # Signal Region muon channel = 1
                    if not sample_name=="Signal":
                        Wmass_mu.Fill(Wmass,Event_Weight)
                        #print "Wmass: ", Wmass, "  Wmass_fit[0]", _Wmass_fit[0]
                        if not Wmass == _Wmass_fit[0]:
                            print "WARNING!!! mu channel"

                if sample_era == 0 and (not isMuon) and BDT_out < BDT_OUT_ELE:
                    _Categorization_fit[0] = 2 # Control Region electron channel = 2
                if sample_era == 0 and (not isMuon) and BDT_out >= BDT_OUT_ELE:
                    _Categorization_fit[0] = 3 # Signal Region electron channel = 3
                    if not sample_name=="Signal":
                        Wmass_ele.Fill(Wmass,Event_Weight)
                        #print "Wmass: ", Wmass, "  Wmass_fit[0]", _Wmass_fit[0]
                        if not Wmass == _Wmass_fit[0]:
                            print "WARNING!!! ele channel"

                #------------------------ 2017 ------------------------#
                if sample_era == 1 and isMuon and BDT_out < BDT_OUT_MU:
                    _Categorization_fit[0] = 4 # Control Region muon channel = 0
                if sample_era == 1 and isMuon and BDT_out >= BDT_OUT_MU:
                    _Categorization_fit[0] = 5 # Signal Region muon channel = 1
                    if not sample_name=="Signal":
                        Wmass_mu.Fill(Wmass,Event_Weight)
                        #print "Wmass: ", Wmass, "  Wmass_fit[0]", _Wmass_fit[0]
                        if not Wmass == _Wmass_fit[0]:
                            print "WARNING!!! mu channel"

                if sample_era == 1 and (not isMuon) and BDT_out < BDT_OUT_ELE:
                    _Categorization_fit[0] = 6 # Control Region electron channel = 2
                if sample_era == 1 and (not isMuon) and BDT_out >= BDT_OUT_ELE:
                    _Categorization_fit[0] = 7 # Signal Region electron channel = 3
                    if not sample_name=="Signal":
                        Wmass_ele.Fill(Wmass,Event_Weight)
                        #print "Wmass: ", Wmass, "  Wmass_fit[0]", _Wmass_fit[0]
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
                _lep_iso[0]            = lep_iso
                _nBjets[0]             = nBjets
                _nBjets_25[0]          = nBjets_25
                _weight[0]             = Event_Weight
                _deltaphi_lep_pi[0]    = deltaphi_lep_pi
                _deltaphi_lep_gamma[0] = deltaphi_lep_gamma
                _isMuon[0]             = isMuon
                _piRelIso_05[0]        = piRelIso_05
                _piRelIso_05_ch[0]     = piRelIso_05_ch
                _pi_dxy[0]             = pi_dxy
                _met[0]                = met
                _met_puppi[0]          = met_puppi
                _Wmass[0]              = Wmass
                
                
                if sample_name == "Signal" and isMuon and split_MC and entry_index <= SignalTree_entries/2:
                    tMVA_signal_mu_training.Fill()
                if sample_name == "Signal" and (not isMuon) and split_MC and entry_index <= SignalTree_entries/2:
                    tMVA_signal_ele_training.Fill()
                if sample_name == "Signal" and isMuon and split_MC and entry_index > SignalTree_entries/2:
                    tMVA_signal_mu_test.Fill()
                if sample_name == "Signal" and (not isMuon) and split_MC and entry_index > SignalTree_entries/2:
                    tMVA_signal_ele_test.Fill()
                if sample_name == "Signal" and isMuon and (not split_MC):
                    tMVA_signal_mu.Fill()
                if sample_name == "Signal" and (not isMuon) and (not split_MC):
                    tMVA_signal_ele.Fill()
                if (not sample_name == "Signal") and isMuon:
                    lep_pt_mu.Fill(lep_pT,Event_Weight)
                    tMVA_background_mu.Fill()
                if (not sample_name == "Signal") and (not isMuon):
                    lep_pt_ele.Fill(lep_pT,Event_Weight)
                    tMVA_background_ele.Fill()

            if isData:# and (Wmass < 65. or Wmass > 90.):

                _gamma_eT[0]        = gamma_eT
                _pi_pT[0]           = pi_pT
                _lep_pT[0]          = lep_pT
                _lep_iso[0]         = lep_iso
                _nBjets[0]          = nBjets
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

#print "N_WGToLNuG_mu", N_WGToLNuG_mu 

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

print "nSig_mu: ", nSig_mu, "nBkg_mu: ", nBkg_mu 
print "nSig_ele: ", nSig_ele, "nBkg_ele: ", nBkg_ele 
print "nBkg_processed_mu: ", nBkg_processed_mu, "nBkg_processed_ele: ", nBkg_processed_ele

canvas1 = TCanvas("canvas1","canvas1",200,106,600,600)
ROOT.gStyle.SetOptStat(0)
Wmass_mu.SetAxisRange(0.,65.,"Y")
Wmass_mu.Draw("hist")
canvas1.SaveAs("Wmass_mu_create_rootfile.pdf")

canvas2 = TCanvas("canvas2","canvas2",200,106,600,600)
ROOT.gStyle.SetOptStat(0)
Wmass_ele.SetAxisRange(0.,65.,"Y")
Wmass_ele.Draw("hist")
canvas2.SaveAs("Wmass_ele_create_rootfile.pdf")

canvas3 = TCanvas("canvas3","canvas3",200,106,600,600)
ROOT.gStyle.SetOptStat(0)
lep_pt_mu.Draw("hist")
canvas3.SaveAs("lep_pt_mu_create_rootfile.pdf")

canvas4 = TCanvas("canvas4","canvas4",200,106,600,600)
ROOT.gStyle.SetOptStat(0)
lep_pt_ele.Draw("hist")
canvas4.SaveAs("lep_pt_ele_create_rootfile.pdf")

# canvas5 = TCanvas("canvas5","canvas5",200,106,600,600)
# ROOT.gStyle.SetOptStat(0)
# Wmass_mu_WJets.SetAxisRange(0.,7000.,"Y")
# Wmass_mu_WJets.Draw("hist")
# print "mu integral: ", Wmass_mu_WJets.Integral()
# canvas5.SaveAs("Wmass_mu_WJets_create_rootfile.pdf")

# canvas6= TCanvas("canvas6","canvas6",200,106,600,600)
# ROOT.gStyle.SetOptStat(0)
# Wmass_ele_WJets.SetAxisRange(0.,5000.,"Y")
# Wmass_ele_WJets.Draw("hist")
# print "ele integral: ", Wmass_ele_WJets.Integral()
# canvas6.SaveAs("Wmass_ele_WJets_create_rootfile.pdf")
