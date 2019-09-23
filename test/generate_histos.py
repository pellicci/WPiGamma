import ROOT
import math
import argparse
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

luminosity_norm_2017 = 41.529
luminosity_norm_2017_Ele32_WPTight = 27.13

_Nrandom_for_Ele_32_WPTight_exclusion = ROOT.TRandom3(64524)

#############---------------- BDT score cut values ----------------#############
BDT_OUT_MU  = 0.240
BDT_OUT_ELE = 0.190

############################################################################
#                                                                          #
#-------------------- Get files and normalization map  --------------------#
#                                                                          #
############################################################################

# Get the files and the names of the samples
sample_name = input_filename.split("_")[2]

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

##Get the handlers for all the histos and graphics
h_base  = dict()

list_histos = ["h_mupt", "h_elept", "h_pipt", "h_gammaet", "h_mueta", "h_eleeta","h_pieta","h_gammaeta", "h_nBjets_25","h_deltaphi_mu_pi","h_deltaphi_ele_pi","h_deltaphi_mu_W","h_deltaphi_ele_W","h_deltaeta_mu_pi","h_deltaeta_ele_pi","h_Wmass","h_Wmass_flag_mu","h_Wmass_flag_ele","h_mu_gamma_InvMass","h_ele_gamma_InvMass","h_piIso_05_mu","h_piRelIso_05_mu_ch","h_piRelIso_05_mu","h_piIso_05_ele","h_piRelIso_05_ele_ch","h_piRelIso_05_ele","h_met_mu","h_met_ele","h_met_puppi","h_Wmass_ratio_mu","h_Wmass_ratio_ele","h_nPV_mu","h_nPV_ele","h_deltaphi_mu_gamma","h_deltaphi_ele_gamma","h_deltaR_mu_gamma","h_deltaR_ele_gamma"]

h_base[list_histos[0]]  = ROOT.TH1F(list_histos[0], "p_{T} of the muon", 15, 25, 100.)
h_base[list_histos[1]]  = ROOT.TH1F(list_histos[1], "p_{T} of the electron", 15, 28, 100.)
h_base[list_histos[2]]  = ROOT.TH1F(list_histos[2], "p_{T} of the pion", 15, 20, 100.)
h_base[list_histos[3]]  = ROOT.TH1F(list_histos[3], "E_{T} of the gamma", 15, 20, 100.)
h_base[list_histos[4]]  = ROOT.TH1F(list_histos[4], "eta of the muon", 20, -3, 3)
h_base[list_histos[5]]  = ROOT.TH1F(list_histos[5], "eta of the electron", 20, -3, 3)
h_base[list_histos[6]]  = ROOT.TH1F(list_histos[6], "eta of the pion", 20, -3, 3)
h_base[list_histos[7]]  = ROOT.TH1F(list_histos[7], "eta of the photon", 20, -3, 3)
h_base[list_histos[8]]  = ROOT.TH1F(list_histos[8], "n Bjets 25", 6, 0, 6.)
h_base[list_histos[9]]  = ROOT.TH1F(list_histos[9], "deltaphi mu-pi", 10, 0, 3.14)
h_base[list_histos[10]] = ROOT.TH1F(list_histos[10], "deltaphi ele-pi", 10, 0, 3.14)
h_base[list_histos[11]] = ROOT.TH1F(list_histos[11], "deltaphi mu-W", 10, 0, 3.14)
h_base[list_histos[12]] = ROOT.TH1F(list_histos[12], "deltaphi ele-W", 10, 0, 3.14)
h_base[list_histos[13]] = ROOT.TH1F(list_histos[13], "deltaeta mu-pi", 20, -5, 5)
h_base[list_histos[14]] = ROOT.TH1F(list_histos[14], "deltaeta ele-pi", 20, -5, 5)
h_base[list_histos[15]] = ROOT.TH1F(list_histos[15], "W mass", 10, 50, 100)
h_base[list_histos[16]] = ROOT.TH1F(list_histos[16], "W mass if flag mu", 10, 50, 100)
h_base[list_histos[17]] = ROOT.TH1F(list_histos[17], "W mass if flag ele", 10, 50, 100)
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
h_base[list_histos[33]] = ROOT.TH1F(list_histos[33], "deltaphi mu-gamma", 100, 0., 0.5)
h_base[list_histos[34]] = ROOT.TH1F(list_histos[34], "deltaphi ele-gamma", 100, 0., 0.5)
h_base[list_histos[35]] = ROOT.TH1F(list_histos[35], "deltaR mu-gamma", 50, 0., 5.)
h_base[list_histos[36]] = ROOT.TH1F(list_histos[36], "deltaR ele-gamma", 50, 0., 5.)

_Nrandom_for_Gaus_SF = ROOT.TRandom3(44329)

##Loop on events

if not "Data" in sample_name:
    norm_factor = Norm_Map[sample_name]*luminosity_norm
    print "Norm_Map[", sample_name, "]: ", Norm_Map[sample_name]
    
root_file = ROOT.TFile(input_filename)
mytree = root_file.Get("WPiGammaAnalysis/mytree")
print "Processing Sample ", sample_name

Nevts_per_sample   = 0. # Count the number of events in input per each sample processed
Nevts_selected     = 0. # Count the number of events survived per each sample processed
Nevts_expected     = 0. # Number of expected events from weights
Nevts_expected_mu  = 0. # Number of expected events from weights
Nevts_expected_ele = 0. # Number of expected events from weights

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

    if runningEra == 0 and sample_name == "ttbar" and mytree.isttbarlnu: # Avoid double-counting of the ttbarlnu background
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

    Nevts_per_sample = Nevts_per_sample + 1

    isMuon = mytree.is_muon
    LepPiOppositeCharge = mytree.LepPiOppositeCharge

    nPV = mytree.nPV

    isSingleMuTrigger_24 = mytree.isSingleMuTrigger_24
    isSingleMuTrigger_50 = mytree.isSingleMuTrigger_50
    if runningEra == 1:
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
                if _Nrandom_for_Ele_32_WPTight_exclusion.Rndm() <= (luminosity_norm_2017_Ele32_WPTight/luminosity_norm_2017):
                    if not isSingleEleTrigger_32:
                        continue


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
    gamma_etaSC = mytree.photon_etaSC
    gamma_phi = mytree.photon_phi
    gamma_E = mytree.photon_energy
    gamma_FourMomentum = ROOT.TLorentzVector()
    gamma_FourMomentum.SetPtEtaPhiE(gamma_eT,gamma_eta,gamma_phi,gamma_E)
    gamma_iso_ChHad = mytree.photon_iso_ChargedHadron
    gamma_iso_NeuHad = mytree.photon_iso_NeutralHadron
    gamma_iso_Ph = mytree.photon_iso_Photon
    gamma_iso_eArho = mytree.photon_iso_eArho

    met = mytree.met_pT
    met_puppi = mytree.metpuppi_pT
       
    Wmass = mytree.Wmass
    W_phi = (pi_FourMomentum + gamma_FourMomentum).Phi()

    lep_FourMomentum = ROOT.TLorentzVector()
    lep_FourMomentum.SetPtEtaPhiM(lep_pT,lep_eta,lep_phi,0.)

    nBjets = mytree.nBjets
    nBjets_25 = mytree.nBjets_25

    deltaeta_lep_pi = math.fabs(lep_eta-pi_eta)
        
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
    deltaR_lep_gamma = math.sqrt(deltaphi_lep_gamma*deltaphi_lep_gamma + deltaeta_lep_gamma*deltaeta_lep_gamma)

    #Avoid double counting between WJetsToLNu and WGToLNuG, based on the DeltaR(lepton,photon)
    if "WJetsToLNu" in sample_name:

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

        if MCT_deltaR_lep_gamma > 0.5 and not MCT_lep_pT < 0. and not mytree.MCT_HeT_ph_eT < 0.:
            continue
        # if not MCT_lep_pT < 0. and not mytree.MCT_HeT_ph_eT < 0.:
        #     continue


    if not "Data" in sample_name:
        Wplus_pT = mytree.Wplus_pT
        Wminus_pT = mytree.Wminus_pT

    if sample_name == "Signal":
        is_signal_Wplus = mytree.is_signal_Wplus
        is_signal_Wminus = mytree.is_signal_Wminus

    if not isMuon:
        ele_gamma_InvMass = (lep_FourMomentum + gamma_FourMomentum).M()
        ele_etaSC = mytree.lepton_etaSC
    else:
        mu_gamma_InvMass = (lep_FourMomentum + gamma_FourMomentum).M()
        mu_pi_InvMass    = (lep_FourMomentum + pi_FourMomentum).M()


    ############################################################################
    #                                                                          #
    #----------------------- Some post-preselection cuts ----------------------#
    #                                                                          #
    ############################################################################

    if myWF.post_preselection_cuts(lep_eta,lep_pT,isMuon,LepPiOppositeCharge,deltaphi_lep_gamma,runningEra):
        continue
        
    Nevts_selected = Nevts_selected + 1

    ############################################################################
    #                                                                          #
    #--------------------- Determine the total event weight -------------------#
    #                                                                          #
    ############################################################################

    ################ MUON SFs ################

    if isMuon: # Get muon scale factors, which are different for two groups of datasets, and weight them for the respective integrated lumi 
        if runningEra == 0:
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
        MC_Weight = mytree.MC_Weight # Add MC weight        
        PU_Weight = mytree.PU_Weight # Add Pile Up weight
        Prefiring_Weight = mytree.Prefiring_Weight # Add prefiring weight  
        Event_Weight = norm_factor*ph_weight*MC_Weight*PU_Weight/math.fabs(MC_Weight)*Prefiring_Weight # Just take the sign of the gen weight

        Event_Weight = Event_Weight*lep_weight

        ################ Zvtx inefficiency weight for 2017 MC ################
        #https://twiki.cern.ch/twiki/bin/view/CMS/Egamma2017DataRecommendations

        if runningEra == 1:
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

    else:
        Event_Weight = 1.

    if "QCD" in sample_name and Event_Weight >= 4:
        continue

    Nevts_expected += Event_Weight # Increment the number of events survived in the analyzed sample
    if isMuon:
        Nevts_expected_mu += Event_Weight
    else:
        Nevts_expected_ele += Event_Weight

    ############################################################################
    #                                                                          #
    #-------------------------- Retrieve BDT output  --------------------------#
    #                                                                          #
    ############################################################################
    if isBDT_with_Wmass:
        #BDT_out = myWF.get_BDT_output(pi_pT/Wmass,gamma_eT/Wmass,nBjets_25,lep_pT,piRelIso_05_ch,met,deltaphi_lep_gamma,isMuon)
        BDT_out = myWF.get_BDT_output(pi_pT/Wmass,gamma_eT/Wmass,nBjets_25,lep_pT,piRelIso_05_ch,met,isMuon)
    elif isBDT:
        #BDT_out = myWF.get_BDT_output(pi_pT,gamma_eT,nBjets_25,lep_pT,piRelIso_05_ch,met,deltaphi_lep_gamma,isMuon)
        BDT_out = myWF.get_BDT_output(pi_pT,gamma_eT,nBjets_25,lep_pT,piRelIso_05_ch,met,isMuon)
   
    ############################################################################
    #                                                                          #
    #------------------------------- Fill histos ------------------------------#
    #                                                                          #
    ############################################################################

    if isBDT and ((isMuon and BDT_out > -0.1 and BDT_out < BDT_OUT_MU) or (not isMuon and BDT_out > -0.1 and BDT_out < BDT_OUT_ELE)): # Alternative cut on BDT score. Two histos will be filled. Eventually, they will be the ratio of Wmass distributions: one with pi_pT and gamma_ET normalized to Wmass, one not
        if (Wmass >= 50. and Wmass <= 100.) and isMuon:
            h_base["h_Wmass_ratio_mu"].Fill(Wmass,Event_Weight)
        if (Wmass >= 50. and Wmass <= 100.) and not isMuon:
            h_base["h_Wmass_ratio_ele"].Fill(Wmass,Event_Weight)

    if (isBDT and isMuon and BDT_out < BDT_OUT_MU) or (isBDT and not isMuon and BDT_out < BDT_OUT_ELE): # Cut on BDT output
        continue
    if isBDT and (Wmass < 50. or Wmass > 100.): # General Wmass condition
        continue
    if isBDT and sample_name == "Data" and (Wmass >= 65. and Wmass <= 90.): # Exclude data in the Blind Window
        continue

    h_base["h_nBjets_25"].Fill(nBjets_25,Event_Weight)
    h_base["h_met_puppi"].Fill(met_puppi,Event_Weight)
    h_base["h_Wmass"].Fill(Wmass,Event_Weight)
    h_base["h_pipt"].Fill(pi_pT,Event_Weight)
    h_base["h_pieta"].Fill(pi_eta,Event_Weight)
    h_base["h_gammaet"].Fill(gamma_eT,Event_Weight)
    h_base["h_gammaeta"].Fill(gamma_eta,Event_Weight)

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
        h_base["h_deltaphi_mu_W"].Fill(deltaphi_lep_W,Event_Weight)
        h_base["h_mu_gamma_InvMass"].Fill(mu_gamma_InvMass,Event_Weight)  
        h_base["h_deltaeta_mu_pi"].Fill(deltaeta_lep_pi,Event_Weight)
        h_base["h_deltaphi_mu_pi"].Fill(deltaphi_lep_pi,Event_Weight)
        h_base["h_deltaphi_mu_gamma"].Fill(deltaphi_lep_gamma,Event_Weight)
        if "WJetsToLNu" in sample_name:
            h_base["h_deltaR_mu_gamma"].Fill(MCT_deltaR_lep_gamma,Event_Weight)

    else:
        h_base["h_Wmass_flag_ele"].Fill(Wmass,Event_Weight)
        h_base["h_elept"].Fill(lep_pT,Event_Weight)
        h_base["h_eleeta"].Fill(lep_eta,Event_Weight)
        h_base["h_nPV_ele"].Fill(nPV,Event_Weight)
        h_base["h_piRelIso_05_ele"].Fill(piRelIso_05,Event_Weight)
        h_base["h_piRelIso_05_ele_ch"].Fill(piRelIso_05_ch,Event_Weight)
        h_base["h_ele_gamma_InvMass"].Fill(ele_gamma_InvMass,Event_Weight)
        h_base["h_met_ele"].Fill(met,Event_Weight)
        h_base["h_deltaeta_ele_pi"].Fill(deltaeta_lep_pi,Event_Weight)
        h_base["h_deltaphi_ele_pi"].Fill(deltaphi_lep_pi,Event_Weight)
        h_base["h_deltaphi_ele_gamma"].Fill(deltaphi_lep_gamma,Event_Weight)
        h_base["h_deltaphi_ele_W"].Fill(deltaphi_lep_W,Event_Weight)
        if "WJetsToLNu" in sample_name:
            h_base["h_deltaR_ele_gamma"].Fill(MCT_deltaR_lep_gamma,Event_Weight)


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
#     fOut_ttbar = ROOT.TFile("ttbar_SemiLeptonic_pT.root","RECREATE")
#     fOut_ttbar.cd()
#     W_ttbar_hist.Write()
#     fOut_ttbar.Close()

fOut = ROOT.TFile(output_filename,"RECREATE")
fOut.cd()
for hist_name in list_histos:
    h_base[hist_name].Write()
fOut.Close()


print "###################"
print "Number of expected events for ", luminosity_norm, " in fb-1, for sample " , sample_name
print "Number of events processed = ", Nevts_per_sample
print "Number of events selected = ", Nevts_selected
print "Number of events expected = ", Nevts_expected
print "Number of events expected in muon channel = ", Nevts_expected_mu
print "Number of events expected in electron channel = ", Nevts_expected_ele
print "###################"
