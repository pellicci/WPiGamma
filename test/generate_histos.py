import ROOT
import math
import argparse
import numpy as np
import array as arr
from Simplified_Workflow_Handler import Simplified_Workflow_Handler

############################################################################
#                                                                          #
#----------------------- Some bools to be initialized ---------------------#
#                                                                          #
############################################################################
p = argparse.ArgumentParser(description='Select whether to fill the histograms after pre-selection or after BDT')
p.add_argument('isBDT_option', help='Type <<preselection>> or <<BDT>>')
p.add_argument('runningEra_option', help='Type <<0>> for 2016, <<1>> for 2017, <<2>> for 2018, <<3>> for combination 2016+2017')
p.add_argument('inputfile_option', help='Provide input file name')
p.add_argument('outputfile_option', help='Provide output file name')
args = p.parse_args()

isBDT = False
if args.isBDT_option == "BDT":
    isBDT = True

runningEra = int(args.runningEra_option)
input_filename = args.inputfile_option
output_filename = args.outputfile_option

#Configuration bools---------------------#
isBDT_with_Wmass = False # If true, pT(pi) and ET(gamma) in the BDT are normalized to Wmass 
random_mu_SF     = False #------if True, muon scale factors are sampled from a Gaussian
random_ele_SF    = False #------if True, electron scale factors are sampled from a Gaussian
random_ph_SF     = False #------if True, photon scale factors are sampled from a Gaussian

#-------------------------#
myWF = Simplified_Workflow_Handler("Signal","Data",isBDT,isBDT_with_Wmass,runningEra)

############################################################################
#                                                                          #
#-------------------------- Integrated luminosity -------------------------#
#                                                                          #
############################################################################
#Normalize to this luminsity, in fb-1
if runningEra == 0:
    luminosity_norm = 35.86
if runningEra == 1:
    luminosity_norm = 41.529
    #luminosity_norm = 27.13 #lumi for Ele32_WPTight only
if runningEra == 2:
    luminosity_norm = 59.69

#############---------------- BDT score cut values ----------------#############
BDT_OUT_MU  = 0.281
BDT_OUT_ELE = 0.269
#BDT_OUT_MU  = 0.273 
#BDT_OUT_ELE = 0.258 #Wmass

############################################################################
#                                                                          #
#-------------------- Get files and normalization map  --------------------#
#                                                                          #
############################################################################

# Get the files and the names of the samples
sample_name = input_filename.split("_")[2]
#sample_name = input_filename.split("_")[6]

# Get the normalization
Norm_Map = myWF.get_normalizations_map(runningEra)

############################################################################
#                                                                          #
#------------------------------ Create histos -----------------------------#
#                                                                          #
############################################################################

#Here's the list of histos to plot
W_signal_hist = ROOT.TH1F("W_signal","W signal",8,0,300)
W_ttbar_hist  = ROOT.TH1F("W_ttbar","W ttbar",8,0,300)

h_PUdistrib = ROOT.TH1F("pile up", "pile up",75,0,75)

W_signal_hist.Sumw2()
W_ttbar_hist.Sumw2()
h_PUdistrib.Sumw2()

ttbar_sig_calib = myWF.get_ttbar_signal_reweight()

h_genW_pT       = ROOT.TH1F("h_genW_pT", "", 15, 0., 300.)
h_genW_eta      = ROOT.TH1F("h_genW_eta", "", 20, -2.6, 2.6)
h_genW_phi      = ROOT.TH1F("h_genW_phi", "", 15, -3.14, 3.14)
h_genTop_pT     = ROOT.TH1F("h_genTop_pT", "", 15, 0., 300.)
h_genTop_eta    = ROOT.TH1F("h_genTop_eta", "", 20, -2.6, 2.6)
h_genTop_phi    = ROOT.TH1F("h_genTop_phi", "", 15, -3.14, 3.14)
h_angle_Top_Pi  = ROOT.TH1F("h_angle_Top_Pi", "", 25, 0., 3.14)
#h_angle_Top_Pi  = ROOT.TH1F("h_angle_Top_Pi", "", 20, -1., 1.)
h_angle_Pi_Ph   = ROOT.TH1F("h_angle_Pi_Ph", "", 35, -3.14, 6.28)
h_angle_Pi_W    = ROOT.TH1F("h_angle_Pi_W", "", 20, -1., 1.)
#h_angle_Pi_W    = ROOT.TH1F("h_angle_Pi_W", "", 35, 0., 3.14)

h_genW_pT_up    = ROOT.TH1F("h_genW_pT_up", "", 15, 0., 300.)
h_genW_pT_down  = ROOT.TH1F("h_genW_pT_down", "", 15, 0., 300.)
h_genW_pT_up_over_nominal   = ROOT.TH1F("h_genW_pT_up_over_nominal", "", 15, 0., 300.)
h_genW_pT_down_over_nominal = ROOT.TH1F("h_genW_pT_down_over_nominal", "", 15, 0., 300.)


pT_ele_binning = arr.array('f',[33.,40.,50.,100.,200.,500.])

##Get the handlers for all the histos and graphics
h_base  = dict()

list_histos = ["h_mupt", "h_elept", "h_pipt", "h_gammaet", "h_mueta", "h_eleeta","h_pieta","h_gammaeta", "h_nBjets_25","h_deltaphi_mu_pi","h_deltaphi_ele_pi","h_deltaphi_mu_W","h_deltaphi_ele_W","h_deltaeta_mu_pi","h_deltaeta_ele_pi","h_Wmass","h_Wmass_flag_mu","h_Wmass_flag_ele","h_mu_gamma_InvMass","h_ele_gamma_InvMass","h_piIso_05_mu","h_piRelIso_05_mu_ch","h_piRelIso_05_mu","h_piIso_05_ele","h_piRelIso_05_ele_ch","h_piRelIso_05_ele","h_met_mu","h_met_ele","h_met_puppi","h_Wmass_alternative_mu","h_Wmass_alternative_ele","h_nPV_mu","h_nPV_ele","h_deltaphi_mu_gamma","h_deltaphi_ele_gamma","h_deltaR_mu_gamma","h_deltaR_ele_gamma","h_lepton_eta","h_lepton_pt","h_piRelIso_05_ch","h_deltaR_mu_pi","h_deltaR_ele_pi","h_nBjets_scaled","h_met_mu_scaled","h_met_ele_scaled","h_njets","h_nBjets_vs_njets","MCT_deltaR_lep_gamma","h_BDT_out_mu","h_BDT_out_ele","h_mu_met_mT","h_ele_met_mT","h_met","h_pi_ph_met_InvMass","h_piRelIso_03"]

h_base[list_histos[0]]  = ROOT.TH1F(list_histos[0], "p_{T} of the muon", 15, 25, 100.)
h_base[list_histos[1]]  = ROOT.TH1F(list_histos[1], "p_{T} of the electron", 15, 28, 100.)
#h_base[list_histos[1]]  = ROOT.TH1F(list_histos[1], "p_{T} of the electron", 5, pT_ele_binning)
h_base[list_histos[2]]  = ROOT.TH1F(list_histos[2], "p_{T} of the pion", 15, 20., 95.)
h_base[list_histos[3]]  = ROOT.TH1F(list_histos[3], "E_{T} of the gamma", 15, 20., 95.)
h_base[list_histos[4]]  = ROOT.TH1F(list_histos[4], "eta of the muon", 20, -3, 3)
h_base[list_histos[5]]  = ROOT.TH1F(list_histos[5], "eta of the electron", 20, -3, 3)
h_base[list_histos[6]]  = ROOT.TH1F(list_histos[6], "eta of the pion", 20, -3, 3)
h_base[list_histos[7]]  = ROOT.TH1F(list_histos[7], "eta of the photon", 20, -3, 3)
h_base[list_histos[8]]  = ROOT.TH1F(list_histos[8], "n Bjets 25", 6, -0.5, 5.5)
h_base[list_histos[9]]  = ROOT.TH1F(list_histos[9], "deltaphi mu-pi", 50, 0, 3.14)
h_base[list_histos[10]] = ROOT.TH1F(list_histos[10], "deltaphi ele-pi", 50, 0, 3.14)
h_base[list_histos[11]] = ROOT.TH1F(list_histos[11], "deltaphi mu-W", 10, 0, 3.14)
h_base[list_histos[12]] = ROOT.TH1F(list_histos[12], "deltaphi ele-W", 10, 0, 3.14)
h_base[list_histos[13]] = ROOT.TH1F(list_histos[13], "deltaeta mu-pi", 20, -5, 5)
h_base[list_histos[14]] = ROOT.TH1F(list_histos[14], "deltaeta ele-pi", 20, -5, 5)
h_base[list_histos[15]] = ROOT.TH1F(list_histos[15], "W mass", 20, 50, 100)
h_base[list_histos[16]] = ROOT.TH1F(list_histos[16], "W mass if flag mu", 20, 50, 100)
h_base[list_histos[17]] = ROOT.TH1F(list_histos[17], "W mass if flag ele", 20, 50, 100)
h_base[list_histos[18]] = ROOT.TH1F(list_histos[18], "mu-gamma InvMass", 20, 0, 300)
h_base[list_histos[19]] = ROOT.TH1F(list_histos[19], "ele-gamma InvMass", 20, 0, 300)
h_base[list_histos[20]] = ROOT.TH1F(list_histos[20], "Pion isolation 05 - mu", 75, 0, 150)
h_base[list_histos[21]] = ROOT.TH1F(list_histos[21], "Pion rel. isolation 05 - mu - ch", 50, 0, 10)
h_base[list_histos[22]] = ROOT.TH1F(list_histos[22], "Pion rel. isolation 05 - mu", 50, 0, 10)
h_base[list_histos[23]] = ROOT.TH1F(list_histos[23], "Pion isolation 05 - ele", 75, 0, 150)
h_base[list_histos[24]] = ROOT.TH1F(list_histos[24], "Pion rel. isolation 05 - ele - ch", 50, 0, 10)
h_base[list_histos[25]] = ROOT.TH1F(list_histos[25], "Pion rel. isolation 05 - ele", 50, 0, 10)
h_base[list_histos[26]] = ROOT.TH1F(list_histos[26], "met mu", 20, 0, 200)
h_base[list_histos[27]] = ROOT.TH1F(list_histos[27], "met ele", 20, 0, 200)
h_base[list_histos[28]] = ROOT.TH1F(list_histos[28], "met puppi", 45, 0, 300)
h_base[list_histos[29]] = ROOT.TH1F(list_histos[29], "Wmass ratio mu", 10, 50, 100)
h_base[list_histos[30]] = ROOT.TH1F(list_histos[30], "Wmass ratio ele", 10, 50, 100)
h_base[list_histos[31]] = ROOT.TH1F(list_histos[31], "nPV - mu", 15, 0, 50)
h_base[list_histos[32]] = ROOT.TH1F(list_histos[32], "nPV - ele", 15, 0, 50)
# h_base[list_histos[33]] = ROOT.TH1F(list_histos[33], "deltaphi mu-gamma", 30, 0., 3.14)
# h_base[list_histos[34]] = ROOT.TH1F(list_histos[34], "deltaphi ele-gamma", 30, 0., 3.14)
h_base[list_histos[33]] = ROOT.TH1F(list_histos[33], "deltaphi mu-gamma", 50, 0., 3.14)
h_base[list_histos[34]] = ROOT.TH1F(list_histos[34], "deltaphi ele-gamma", 50, 0., 3.14)
h_base[list_histos[35]] = ROOT.TH1F(list_histos[35], "deltaR mu-gamma", 60, 0., 6.)
h_base[list_histos[36]] = ROOT.TH1F(list_histos[36], "deltaR ele-gamma", 60, 0., 6.)
h_base[list_histos[37]] = ROOT.TH1F(list_histos[37], "lepton eta", 20, -3., 3.)
h_base[list_histos[38]] = ROOT.TH1F(list_histos[38], "p_{T} of the lepton", 15, 25., 100.)
#h_base[list_histos[39]] = ROOT.TH1F(list_histos[39], "Pion rel. isolation 05 - ch", 30, 0, 6)
h_base[list_histos[39]] = ROOT.TH1F(list_histos[39], "Pion rel. isolation 05 - ch", 16, 0, 0.8)
h_base[list_histos[40]] = ROOT.TH1F(list_histos[40], "deltaR mu-pi", 60, 0., 6.)
h_base[list_histos[41]] = ROOT.TH1F(list_histos[41], "deltaR ele-pi", 60, 0., 6.)
h_base[list_histos[42]] = ROOT.TH1F(list_histos[42], "n Bjets scaled", 6, -0.5, 5.5)
h_base[list_histos[43]] = ROOT.TH1F(list_histos[43], "met mu scaled", 20, 0, 200)
h_base[list_histos[44]] = ROOT.TH1F(list_histos[44], "met ele scaled", 20, 0, 200)
h_base[list_histos[45]] = ROOT.TH1F(list_histos[45], "n jets", 10, -0.5, 9.5)
h_base[list_histos[46]] = ROOT.TH2F(list_histos[46], "nBjets vs njets", 20, -0.5, 19.5, 10, -0.5, 9.5,)
h_base[list_histos[47]] = ROOT.TH1F(list_histos[47], "MCT deltaR lep gamma", 20, 0., 2.)
h_base[list_histos[48]] = ROOT.TH1F(list_histos[48], "BDT out_mu", 40, -0.7, 0.7)
h_base[list_histos[49]] = ROOT.TH1F(list_histos[49], "BDT out_ele", 40, -0.7, 0.7)
h_base[list_histos[50]] = ROOT.TH1F(list_histos[50], "mu-met mT", 40, 0., 200.)
h_base[list_histos[51]] = ROOT.TH1F(list_histos[51], "ele-met mT", 40, 0., 200.)
h_base[list_histos[52]] = ROOT.TH1F(list_histos[52], "met", 20, 0., 200.)
h_base[list_histos[53]] = ROOT.TH1F(list_histos[53], "Pion - photon - MET invariant mass", 45, 0., 300.)
h_base[list_histos[54]] = ROOT.TH1F(list_histos[54], "Pion rel. isolation 03", 30, 0, 6)

_Nrandom_for_Gaus_SF = ROOT.TRandom3(44329)

##Loop on events

if not "Data" in sample_name:
    norm_factor = Norm_Map[sample_name]*luminosity_norm
    print "Norm_Map[", sample_name, "]: ", Norm_Map[sample_name]
    
root_file = ROOT.TFile(input_filename)
mytree = root_file.Get("WPiGammaAnalysis/mytree")
print "Processing Sample ", sample_name

Nevts_per_sample   = 0. # Count the number of events in input per each sample processed
Nevts_selected_mu  = 0. # Count the number of events survived per each sample processed (muon channel)
Nevts_selected_ele = 0. # Count the number of events survived per each sample processed (electron channel)
Nevts_expected     = 0. # Number of expected events from weights
Nevts_expected_mu  = 0. # Number of expected events from weights
Nevts_expected_ele = 0. # Number of expected events from weights

nPionsTrue = 0.
nPionsTot = 0.
nPileUp = 0.

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

    if mytree.is_muon and sample_name == "QCDDoubleEMEnriched30toInf":
        continue

    if not mytree.is_muon and sample_name == "QCDHT200to300":
        continue

    ############################################################################
    #                                                                          #
    #------------------------ Access the tree variables -----------------------#
    #                                                                          #
    ############################################################################

    Nevts_per_sample = Nevts_per_sample + 1

    isMuon = mytree.is_muon

    nPV = mytree.nPV

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
    gamma_etaSC = mytree.photon_etaSC
    gamma_phi = mytree.photon_phi
    gamma_E = mytree.photon_energy
    gamma_FourMomentum = ROOT.TLorentzVector()
    gamma_FourMomentum.SetPtEtaPhiE(gamma_eT,gamma_eta,gamma_phi,gamma_E)

    met = mytree.met_pT
    met_scaled = mytree.met_pT_scaled
    met_puppi = mytree.metpuppi_pT
       
    Wmass = mytree.Wmass
    W_phi = (pi_FourMomentum + gamma_FourMomentum).Phi()

    lep_FourMomentum = ROOT.TLorentzVector()
    lep_FourMomentum.SetPtEtaPhiM(lep_pT,lep_eta,lep_phi,0.)

    met_FourMomentum = ROOT.TLorentzVector()
    met_FourMomentum.SetPtEtaPhiE(met,0.,mytree.met_phi,met)

    pi_ph_met_InvMass = (pi_FourMomentum + gamma_FourMomentum + met_FourMomentum).M()

    #lep_met_TransverseMass = (lep_FourMomentum + met_FourMomentum).Mt()

    deltaphi_lep_met = math.fabs(lep_phi - mytree.met_phi)
    if deltaphi_lep_met > 3.14:
       deltaphi_lep_met = 6.28 - deltaphi_lep_met

    lep_met_TransverseMass = math.sqrt(2*lep_pT*met*(1. - math.cos(deltaphi_lep_met)))

    nBjets_30 = mytree.nBjets_30
    nBjets_25 = mytree.nBjets_25
    nBjets_scaled = mytree.nBjets_scaled

    deltaeta_lep_pi = math.fabs(lep_eta - pi_eta)
        
    deltaphi_lep_pi = math.fabs(lep_phi - pi_phi)
    if deltaphi_lep_pi > 3.14:
        deltaphi_lep_pi = 6.28 - deltaphi_lep_pi

    deltaphi_lep_W = math.fabs(lep_phi - W_phi)
    if deltaphi_lep_W > 3.14:
        deltaphi_lep_W = 6.28 - deltaphi_lep_W

    deltaphi_lep_gamma = math.fabs(lep_phi - gamma_phi)
    if deltaphi_lep_gamma > 3.14:
        deltaphi_lep_gamma = 6.28 - deltaphi_lep_gamma

    deltaeta_lep_gamma = math.fabs(lep_eta - gamma_eta)

    deltaR_lep_pi    = math.sqrt(deltaphi_lep_pi*deltaphi_lep_pi + deltaeta_lep_pi*deltaeta_lep_pi)
    deltaR_lep_gamma = math.sqrt(deltaphi_lep_gamma*deltaphi_lep_gamma + deltaeta_lep_gamma*deltaeta_lep_gamma)

    #Avoid double counting between WJetsToLNu and WGToLNuG, based on the DeltaR(lepton,photon)
    if "WJetsToLNu" in sample_name:# or "WGToLNuG" in sample_name:

        if isMuon:
            MCT_lep_eta = mytree.MCT_HpT_mu_eta
            MCT_lep_phi = mytree.MCT_HpT_mu_phi
            MCT_lep_pT  = mytree.MCT_HpT_mu_pT
            # if MCT_lep_eta == -1000. or MCT_lep_phi == -1000.:
            #     print "NO MCT MUON IN THIS EVENT"
        else:
            MCT_lep_eta = mytree.MCT_HpT_ele_eta
            MCT_lep_phi = mytree.MCT_HpT_ele_phi
            MCT_lep_pT  = mytree.MCT_HpT_ele_pT
            # if MCT_lep_eta == -1000. or MCT_lep_phi == -1000.:
            #     print "NO MCT ELECTRON IN THIS EVENT"
        
        MCT_deltaeta_lep_gamma = math.fabs(MCT_lep_eta - mytree.MCT_HeT_ph_eta)
        MCT_deltaphi_lep_gamma = math.fabs(MCT_lep_phi - mytree.MCT_HeT_ph_phi)
        # if not (mytree.MCT_HeT_ph_eta == -1000. or mytree.MCT_HeT_ph_phi == -1000. or mytree.MCT_HeT_ph_eT == -1000.):
        #     print "NO MCT PHOTON IN THIS EVENT"

        if MCT_deltaphi_lep_gamma > 3.14:
            MCT_deltaphi_lep_gamma = 6.28 - MCT_deltaphi_lep_gamma

        MCT_deltaR_lep_gamma = math.sqrt(MCT_deltaeta_lep_gamma*MCT_deltaeta_lep_gamma + MCT_deltaphi_lep_gamma*MCT_deltaphi_lep_gamma)
        
        if MCT_deltaR_lep_gamma > 0.2 and not MCT_lep_pT < 0. and not mytree.MCT_HeT_ph_eT < 0.:
            continue


    if not "Data" in sample_name:
        Wplus_pT = mytree.Wplus_pT
        Wminus_pT = mytree.Wminus_pT

    if sample_name == "Signal":
        is_signal_Wplus = mytree.is_signal_Wplus
        is_signal_Wminus = mytree.is_signal_Wminus
        isMuonSignal_fromTau = mytree.isMuonSignal_fromTau
        isEleSignal_fromTau = mytree.isEleSignal_fromTau

        #if isMuonSignal_fromTau or isEleSignal_fromTau:
        #    continue

    if not isMuon:
        ele_gamma_InvMass = (lep_FourMomentum + gamma_FourMomentum).M()
    else:
        mu_gamma_InvMass = (lep_FourMomentum + gamma_FourMomentum).M()
        mu_pi_InvMass    = (lep_FourMomentum + pi_FourMomentum).M()

    ############## For Signal Pythia Modeling Only ##############
    if sample_name == "Signal":
        genPi_pT  = mytree.gen_pi_pT 
        genPi_eta = mytree.gen_pi_eta 
        genPi_phi = mytree.gen_pi_phi
        genPi_energy = mytree.gen_pi_energy
        genTop_pT  = mytree.gen_pi_grandmother_pT 
        genTop_eta = mytree.gen_pi_grandmother_eta 
        genTop_phi = mytree.gen_pi_grandmother_phi
        genTop_energy = mytree.gen_pi_grandmother_energy
        genPh_pT  = mytree.gen_ph_pT 
        genPh_eta = mytree.gen_ph_eta
        genPh_phi = mytree.gen_ph_phi
        genPh_energy = mytree.gen_ph_energy
        
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

        boosted_genW_FourMomentum = ROOT.TLorentzVector()
        boosted_genW_FourMomentum.SetPtEtaPhiE(genW_pT,genW_eta,genW_phi,genW_energy)
        boosted_genW_FourMomentum.Boost(-genW_FourMomentum.BoostVector()) #The sign - is required  
        
        boosted_genTop_FourMomentum = ROOT.TLorentzVector()
        boosted_genTop_FourMomentum.SetPtEtaPhiE(genTop_pT,genTop_eta,genTop_phi,genTop_energy)
        boosted_genTop_FourMomentum.Boost(-genW_FourMomentum.BoostVector()) #The sign - is required
        boosted_genPi_FourMomentum = ROOT.TLorentzVector()
        boosted_genPi_FourMomentum.SetPtEtaPhiE(genPi_pT,genPi_eta,genPi_phi,genPi_energy)
        boosted_genPi_FourMomentum.Boost(-genW_FourMomentum.BoostVector()) #The sign - is required
        boosted_genPh_FourMomentum = ROOT.TLorentzVector()
        boosted_genPh_FourMomentum.SetPtEtaPhiE(genPh_pT,genPh_eta,genPh_phi,genPh_energy)
        boosted_genPh_FourMomentum.Boost(-genW_FourMomentum.BoostVector()) #The sign - is required 
  
        #print boosted_genW_FourMomentum.Px(), boosted_genW_FourMomentum.Py(), boosted_genW_FourMomentum.Pz()
        angle_Top_Pi = boosted_genTop_FourMomentum.Angle(boosted_genPi_FourMomentum.Vect())
        angle_Pi_Ph  = boosted_genPi_FourMomentum.Angle(boosted_genPh_FourMomentum.Vect())
        angle_Pi_W   = boosted_genPi_FourMomentum.Angle(genW_FourMomentum.Vect())
        #print "Boosted eta: ", boosted_genPi_FourMomentum.Eta(), " energy: ", boosted_genPi_FourMomentum.E(), " pT: ", boosted_genPi_FourMomentum.Pt(), " pT da variabile: ", genPi_pT
        #print
    ############################################################

    ############################################################################
    #                                                                          #
    #----------------------- Some post-preselection cuts ----------------------#
    #                                                                          #
    ############################################################################

    if myWF.post_preselection_cuts(lep_eta,lep_pT,ele_etaSC,gamma_etaSC,isMuon,deltaphi_lep_pi,deltaphi_lep_gamma,isTriggerMatched,runningEra):
        continue

    ############################################################################
    #                                                                          #
    #--------------------- Determine the total event weight -------------------#
    #                                                                          #
    ############################################################################

    ################ MUON SFs ################

    if isMuon: # Get muon scale factors, which are different for two groups of datasets, and weight them for the respective integrated lumi 

        isSingleMuTrigger_LOW = isSingleMuTrigger_24
        if runningEra == 1:
            isSingleMuTrigger_LOW = isSingleMuTrigger_27

        lep_weight, lep_weight_err = myWF.get_muon_scale(lep_pT,lep_eta,isSingleMuTrigger_LOW,runningEra)

        if random_mu_SF:
            lep_weight = _Nrandom_for_Gaus_SF.Gaus(lep_weight,lep_weight_err)

    ############## ELECTRON SFs ##############

    else:
        lep_weight, lep_weight_err = myWF.get_ele_scale(lep_pT,ele_etaSC,runningEra)

        if random_ele_SF:
            lep_weight = _Nrandom_for_Gaus_SF.Gaus(lep_weight,lep_weight_err) 

    ############### PHOTON SFs ###############

    ph_weight, ph_weight_err = myWF.get_photon_scale(gamma_eT,gamma_etaSC,runningEra)

    if random_ph_SF:
        ph_weight = _Nrandom_for_Gaus_SF.Gaus(ph_weight,ph_weight_err)
        

    ############### Multiply weights and SFs for MC. Set weight to 1 for data ###############

    if not "Data" in sample_name:
        MC_Weight   = mytree.MC_Weight # Add MC weight        
        PU_Weight   = mytree.PU_Weight # Add Pile Up weight
        # bTag_Weight = mytree.bTag_Weight # Add b-tagging weight
        # if math.isnan(bTag_Weight) or bTag_Weight == np.inf or bTag_Weight == -np.inf: #FIXME
        #     bTag_Weight = 1.
        
        Event_Weight = norm_factor*MC_Weight*lep_weight*ph_weight*PU_Weight/math.fabs(MC_Weight) # Just take the sign of the gen weight

        if not runningEra == 2: # Prefiring weight NOT to be applied to 2018 MC
            Prefiring_Weight = mytree.Prefiring_Weight # Add prefiring weight 
            Event_Weight *= Prefiring_Weight


        ################ Zvtx inefficiency weight for 2017 MC ################
        #https://twiki.cern.ch/twiki/bin/view/CMS/Egamma2017DataRecommendations

        if runningEra == 1:
            if isSingleEleTrigger_32_DoubleEG and not isSingleMuTrigger_27 and not isSingleMuTrigger_50:
                Event_Weight *= 0.991 #uniform penalty for all the 2017 eras

        ################ Correct for the difference in pT of the generated W in Pythia and PowHeg samples ################

        if sample_name == "Signal":
            if is_signal_Wplus:
                local_W_pT = Wplus_pT
            else:
                local_W_pT = Wminus_pT

            if local_W_pT > 300.:
                local_W_pT = 300.

            Event_Weight = Event_Weight*ttbar_sig_calib.GetBinContent(ttbar_sig_calib.GetXaxis().FindBin(local_W_pT))

    else:
        Event_Weight = 1.

    # Remove QCD events with abnormal weight
    if not isBDT and "QCD" in sample_name and Event_Weight >= 1600.:
        continue
    if isBDT and "QCD" in sample_name and Event_Weight >= 30.:
        continue

    ############################################################################
    #                                                                          #
    #-------------------------- Retrieve BDT output  --------------------------#
    #                                                                          #
    ############################################################################
    if isBDT_with_Wmass:
        BDT_out = myWF.get_BDT_output(pi_pT/Wmass,gamma_eT/Wmass,nBjets_25,lep_pT,piRelIso_05_ch,met,isMuon)
    elif isBDT:
        BDT_out = myWF.get_BDT_output(pi_pT,gamma_eT,nBjets_25,lep_pT,piRelIso_05_ch,met,isMuon)
        if isMuon:
            h_base["h_BDT_out_mu"].Fill(BDT_out)
        if not isMuon:
            h_base["h_BDT_out_ele"].Fill(BDT_out)
    ############################################################################
    #                                                                          #
    #------------------------------- Fill histos ------------------------------#
    #                                                                          #
    ############################################################################

    if isBDT and ((isMuon and BDT_out > -0.1 and BDT_out < BDT_OUT_MU) or (not isMuon and BDT_out > -0.1 and BDT_out < BDT_OUT_ELE)): # Alternative cut on BDT score. Two histos will be filled. Eventually, they will be the ratio of Wmass distributions: one with pi_pT and gamma_ET normalized to Wmass, one not
        if (Wmass >= 50. and Wmass <= 100.) and isMuon:
            h_base["h_Wmass_alternative_mu"].Fill(Wmass,Event_Weight)
        if (Wmass >= 50. and Wmass <= 100.) and not isMuon:
            h_base["h_Wmass_alternative_ele"].Fill(Wmass,Event_Weight)

    if (isBDT and isMuon and BDT_out < BDT_OUT_MU) or (isBDT and not isMuon and BDT_out < BDT_OUT_ELE): # Cut on BDT output
        continue
    if isBDT and (Wmass < 50. or Wmass > 100.): # General Wmass condition
        continue
    # if isBDT and sample_name == "Data" and (Wmass >= 65. and Wmass <= 90.): # Exclude data in the Blind Window
    #     continue

    Nevts_expected += Event_Weight # Increment the number of events survived in the analyzed sample


    #if piRelIso_05_ch < 0.2:
    # nPionsTot += Event_Weight
    # if mytree.isPionTrue:
    #     nPionsTrue += Event_Weight
    #     if math.fabs(mytree.pi_gen_mother_ID == 15) and mytree.pi_gen_nDaughters > 2: 
    #         print "trackID", mytree.pi_gen_ID, "  track mother: ", mytree.pi_gen_mother_ID, "  nDaughters", mytree.pi_gen_nDaughters, "  pT: ", pi_pT
    #         nPileUp += Event_Weight
        #else:
            #if math.fabs(mytree.pi_gen_ID == 1) or math.fabs(mytree.pi_gen_ID == 2) or math.fabs(mytree.pi_gen_ID == 3) or math.fabs(mytree.pi_gen_ID == 4):
            #    nPileUp += Event_Weight
            #print "trackID", mytree.pi_gen_ID, "  track mother: ", mytree.pi_gen_mother_ID, "  nDaughters", mytree.pi_gen_nDaughters, "  pT: ", pi_pT

    if isMuon:
        Nevts_selected_mu = Nevts_selected_mu + 1
        #if Wmass >= 65. and Wmass <= 90.:
        Nevts_expected_mu += Event_Weight
    else:
        Nevts_selected_ele = Nevts_selected_ele + 1
        #if Wmass >= 65. and Wmass <= 90.:
        Nevts_expected_ele += Event_Weight


    # if isBDT and sample_name == "Data" and (Wmass < 65. or Wmass > 90.): #Only for blind analysis
    #     h_base["h_Wmass"].Fill(Wmass,Event_Weight)
    # if isBDT and not sample_name == "Data":
    #     h_base["h_Wmass"].Fill(Wmass,Event_Weight)
    # if not isBDT:
    #     h_base["h_Wmass"].Fill(Wmass,Event_Weight)

    #if "WJetsToLNu" in sample_name or "WGToLNuG" in sample_name:
    #    h_base["MCT_deltaR_lep_gamma"].Fill(MCT_deltaR_lep_gamma,Event_Weight)
    h_base["h_Wmass"].Fill(Wmass,Event_Weight)
    h_base["h_nBjets_25"].Fill(nBjets_25,Event_Weight)
    h_base["h_njets"].Fill(mytree.njets,Event_Weight)
    h_base["h_nBjets_vs_njets"].Fill(mytree.njets,nBjets_25,Event_Weight)
    h_base["h_nBjets_scaled"].Fill(nBjets_scaled,Event_Weight)
    h_base["h_met_puppi"].Fill(met_puppi,Event_Weight)
    h_base["h_met"].Fill(met,Event_Weight)
    h_base["h_pipt"].Fill(pi_pT,Event_Weight)
    h_base["h_pieta"].Fill(pi_eta,Event_Weight)
    h_base["h_gammaet"].Fill(gamma_eT,Event_Weight)
    h_base["h_gammaeta"].Fill(gamma_eta,Event_Weight)
    h_base["h_lepton_eta"].Fill(lep_eta,Event_Weight)
    h_base["h_lepton_pt"].Fill(lep_pT,Event_Weight)
    h_base["h_piRelIso_05_ch"].Fill(piRelIso_05_ch,Event_Weight)
    h_base["h_piRelIso_03"].Fill(piRelIso_03,Event_Weight)
    h_base["h_pi_ph_met_InvMass"].Fill(pi_ph_met_InvMass,Event_Weight)

    if sample_name == "Data":
        h_PUdistrib.Fill(nPV,Event_Weight)      

    if isMuon:
        h_base["h_Wmass_flag_mu"].Fill(Wmass,Event_Weight)
        h_base["h_mupt"].Fill(lep_pT,Event_Weight)
        h_base["h_mueta"].Fill(lep_eta,Event_Weight)
        h_base["h_nPV_mu"].Fill(nPV,Event_Weight)
        h_base["h_piRelIso_05_mu"].Fill(piRelIso_05,Event_Weight)
        h_base["h_piRelIso_05_mu_ch"].Fill(piRelIso_05_ch,Event_Weight)
        h_base["h_met_mu"].Fill(met,Event_Weight)
        h_base["h_met_mu_scaled"].Fill(met_scaled,Event_Weight)
        h_base["h_deltaphi_mu_W"].Fill(deltaphi_lep_W,Event_Weight)
        h_base["h_mu_gamma_InvMass"].Fill(mu_gamma_InvMass,Event_Weight)  
        h_base["h_deltaeta_mu_pi"].Fill(deltaeta_lep_pi,Event_Weight)
        h_base["h_deltaphi_mu_pi"].Fill(deltaphi_lep_pi,Event_Weight)
        h_base["h_deltaphi_mu_gamma"].Fill(deltaphi_lep_gamma,Event_Weight)
        h_base["h_deltaR_mu_pi"].Fill(deltaR_lep_pi,Event_Weight)
        h_base["h_mu_met_mT"].Fill(lep_met_TransverseMass,Event_Weight)
        #if "WJetsToLNu" in sample_name:
        h_base["h_deltaR_mu_gamma"].Fill(deltaR_lep_gamma,Event_Weight)
        h_base["h_piIso_05_mu"].Fill(mytree.sum_pT_05_ch,Event_Weight)

    else:
        if math.fabs(ele_etaSC) < 0.8:
            h_base["h_elept"].Fill(lep_pT,Event_Weight)
        h_base["h_Wmass_flag_ele"].Fill(Wmass,Event_Weight)
        h_base["h_eleeta"].Fill(lep_eta,Event_Weight)
        h_base["h_nPV_ele"].Fill(nPV,Event_Weight)
        h_base["h_piRelIso_05_ele"].Fill(piRelIso_05,Event_Weight)
        h_base["h_piRelIso_05_ele_ch"].Fill(piRelIso_05_ch,Event_Weight)
        h_base["h_ele_gamma_InvMass"].Fill(ele_gamma_InvMass,Event_Weight)
        h_base["h_met_ele"].Fill(met,Event_Weight)
        h_base["h_met_ele_scaled"].Fill(met_scaled,Event_Weight)
        h_base["h_deltaeta_ele_pi"].Fill(deltaeta_lep_pi,Event_Weight)
        h_base["h_deltaphi_ele_pi"].Fill(deltaphi_lep_pi,Event_Weight)
        h_base["h_deltaphi_ele_gamma"].Fill(deltaphi_lep_gamma,Event_Weight)
        h_base["h_deltaphi_ele_W"].Fill(deltaphi_lep_W,Event_Weight)
        h_base["h_deltaR_ele_pi"].Fill(deltaR_lep_pi,Event_Weight)
        h_base["h_ele_met_mT"].Fill(lep_met_TransverseMass,Event_Weight)
        #if "WJetsToLNu" in sample_name:
        h_base["h_deltaR_ele_gamma"].Fill(deltaR_lep_gamma,Event_Weight)
        h_base["h_piIso_05_ele"].Fill(mytree.sum_pT_05_ch,Event_Weight)
    #########    #########    #########    #########    #########    #########    #########
    if sample_name == "Signal":
        h_genW_pT.Fill(genW_pT)
        h_genW_pT_up.Fill(genW_pT*1.05)
        h_genW_pT_down.Fill(genW_pT*0.95)
        h_genW_eta.Fill(genW_eta)
        h_genW_phi.Fill(genW_phi)
        h_genTop_pT.Fill(genTop_pT)
        h_genTop_eta.Fill(genTop_eta)
        h_genTop_phi.Fill(genTop_phi)
        h_angle_Top_Pi.Fill(angle_Top_Pi)
        if angle_Pi_W > 3.14:
            angle_Pi_W = 6.28 - angle_Pi_W
        h_angle_Pi_W.Fill(math.cos(angle_Pi_W))
        h_angle_Pi_Ph.Fill(angle_Pi_Ph)


    ############################################################################
    #                                                                          #
    #------------------- Fill Signal-ttbar comparison histos ------------------#
    #                                                                          #
    ############################################################################
    ########----------- To be filled in "preselection" mode! -----------########
#     if sample_name == "Signal":
#         if is_signal_Wplus:
#             local_pT = Wplus_pT
#         else:
#             local_pT = Wminus_pT

#         if local_pT > 300.:
#             W_signal_hist.Fill(300.,Event_Weight)
#         else:
#             W_signal_hist.Fill(local_pT,Event_Weight)


#     if sample_name == "ttbarToSemiLeptonic":
#         if Wplus_pT > 300.:
#             W_ttbar_hist.Fill(300.,Event_Weight)
#         else:
#             W_ttbar_hist.Fill(Wplus_pT,Event_Weight)

#         if Wminus_pT > 300.:
#             W_ttbar_hist.Fill(300.,Event_Weight)
#         else:
#             W_ttbar_hist.Fill(Wminus_pT,Event_Weight)

# if sample_name == "Signal":
#     fOut_signal = ROOT.TFile("signal_pT.root","RECREATE")
#     fOut_signal.cd()
#     W_signal_hist.Write()
#     fOut_signal.Close()

# if sample_name == "ttbarToSemiLeptonic":
#     fOut_ttbar = ROOT.TFile("ttbar_pT.root","RECREATE")
#     fOut_ttbar.cd()
#     W_ttbar_hist.Write()
#     fOut_ttbar.Close()


    ############################################################################                                                 
    #                                                                          #                                                 
    #------------------- Fill Pythia Signal Modeling histos -------------------#                                                 
    #                                                                          #                                                 
    ############################################################################  

# ROOT.gStyle.SetOptStat(0)
# canvas_W_pT = ROOT.TCanvas()
# h_genW_pT.GetXaxis().SetTitle("p_{T}^{W} (GeV)")
# h_genW_pT.Draw("hist")
# canvas_W_pT.SaveAs("plots/latest_production/2018/genW_pT.pdf")
# canvas_W_eta = ROOT.TCanvas()
# h_genW_eta.Draw("hist")
# canvas_W_eta.SaveAs("plots/latest_production/2018/genW_eta.pdf")
# canvas_W_phi = ROOT.TCanvas()
# h_genW_phi.Draw("hist")
# canvas_W_phi.SaveAs("plots/latest_production/2018/genW_phi.pdf")
# canvas_Top_pT = ROOT.TCanvas()
# h_genTop_pT.Draw("hist")
# canvas_Top_pT.SaveAs("plots/latest_production/2018/genTop_pT.pdf")
# canvas_Top_eta = ROOT.TCanvas()
# h_genTop_eta.Draw("hist")
# canvas_Top_eta.SaveAs("plots/latest_production/2018/genTop_eta.pdf")
# canvas_Top_phi = ROOT.TCanvas()
# h_genTop_phi.Draw("hist")
# canvas_Top_phi.SaveAs("plots/latest_production/2018/genTop_phi.pdf")
# canvas_angle_Top_Pi = ROOT.TCanvas()
# h_angle_Top_Pi.GetXaxis().SetTitle("#theta^{h}")
# h_angle_Top_Pi.Draw("hist")
# canvas_angle_Top_Pi.SaveAs("plots/latest_production/2018/angle_Top_Pi.pdf")
# fOut_angle_Top_Pi = ROOT.TFile("cosine_angle_Top_Pi_2016.root","RECREATE")
# fOut_angle_Top_Pi.cd()                                                                                                    
# h_angle_Top_Pi.Write()                                                                                                      
# fOut_angle_Top_Pi.Close()  
# canvas_angle_Pi_Ph = ROOT.TCanvas()
# h_angle_Pi_Ph.Draw("hist")
# canvas_angle_Pi_Ph.SaveAs("plots/latest_production/2018/angle_Pi_Ph.pdf")
# canvas_angle_Pi_W = ROOT.TCanvas()
# h_angle_Pi_W.GetXaxis().SetTitle("#theta^{*}")
# h_angle_Pi_W.Draw("hist")
# canvas_angle_Pi_W.SaveAs("plots/latest_production/2018/angle_Pi_W.pdf")
# fOut_angle_Pi_W = ROOT.TFile("angle_Pi_W_2018.root","RECREATE")
# fOut_angle_Pi_W.cd()                                                                                                    
# h_angle_Pi_W.Write()                                                                                                      
# fOut_angle_Pi_W.Close()  

# ROOT.gStyle.SetOptStat(0)
# canvas1 = ROOT.TCanvas()
# h_genW_pT.Draw("hist")
# h_genW_pT.GetXaxis().SetTitle("p_{T}^{W} (GeV)")
# h_genW_pT_up.SetLineColor(2)
# h_genW_pT_down.SetLineColor(3)
# h_genW_pT_up.Draw("SAME,hist")
# h_genW_pT_down.Draw("SAME,hist")


# leg = ROOT.TLegend(0.6,0.6,0.8,0.8)
# leg.SetHeader(" ")
# leg.SetFillColor(0)
# leg.SetBorderSize(0)
# leg.SetLineColor(1)
# leg.SetLineStyle(1)
# leg.SetLineWidth(1)
# leg.SetFillStyle(1001)
# leg.AddEntry(h_genW_pT,"gen W","l")
# leg.AddEntry(h_genW_pT_up,"gen W up","l")
# leg.AddEntry(h_genW_pT_down,"gen W down","l")
# leg.Draw("SAME")

# canvas1.SaveAs("plots/latest_production/2018/h_genW_pT.pdf")

# h_genW_pT.Scale(1./h_genW_pT.Integral())
# h_genW_pT_up.Scale(1./h_genW_pT_up.Integral())
# h_genW_pT_down.Scale(1./h_genW_pT_down.Integral())
# h_genW_pT_up_over_nominal = h_genW_pT_up
# h_genW_pT_down_over_nominal = h_genW_pT_down
# h_genW_pT_up_over_nominal.Divide(h_genW_pT)
# h_genW_pT_down_over_nominal.Divide(h_genW_pT)
# canvas2 = ROOT.TCanvas()
# h_genW_pT_up_over_nominal.Draw()
# canvas2.SaveAs("plots/latest_production/2018/h_genW_pT_up_over_nominal.pdf")
# canvas3 = ROOT.TCanvas()
# h_genW_pT_down_over_nominal.Draw()
# canvas3.SaveAs("plots/latest_production/2018/h_genW_pT_down_over_nominal.pdf")
# fOutPythia = ROOT.TFile("Pythia_pT_modeling_2018.root","RECREATE")
# fOutPythia.cd()
# h_genW_pT_up_over_nominal.Write()
# h_genW_pT_down_over_nominal.Write()
# fOutPythia.Close()

    #########    #########    #########    #########    #########    #########    #########

fOut = ROOT.TFile(output_filename,"RECREATE")
fOut.cd()
for hist_name in list_histos:
    h_base[hist_name].Write()
fOut.Close()


print "###################"
print "Number of expected events for ", luminosity_norm, " in fb-1, for sample " , sample_name
print "Number of events processed = ", Nevts_per_sample
print "Number of events selected in muon channel = ", Nevts_selected_mu
print "Number of events selected in electron channel = ", Nevts_selected_ele
print "Number of events expected = ", Nevts_expected
print "Number of events expected in muon channel = ", Nevts_expected_mu
print "Number of events expected in electron channel = ", Nevts_expected_ele
print "nPionsTot = ", nPionsTot
print "nPionsTrue = ", nPionsTrue
print "nPileUp = ", nPileUp

print "###################"

