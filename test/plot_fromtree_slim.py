import ROOT
import os
import math
import numpy as np
import copy
import argparse
from array import array
from Workflow_Handler import Workflow_Handler

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

############################################################################
#                                                                          #
#----------------------- Some bools to be initialized ---------------------#
#                                                                          #
############################################################################

#---------------------------------#
p = argparse.ArgumentParser(description='Select whether to fill the histograms after pre-selection or after BDT')
p.add_argument('isBDT_option', help='Type <<preselection>> or <<BDT>>')
args = p.parse_args()

# Switch from muon to electron channel
if args.isBDT_option == "preselection":
    isBDT = False
if args.isBDT_option == "BDT":
    isBDT = True
#---------------------------------#

isBDT_with_Wmass = False # If true, pT(pi) and ET(gamma) in the BDT are normalized to Wmass 
myWF = Workflow_Handler("Signal","Data",isBDT_with_Wmass)

##Bools for rundom SF variation
random_mu_SF  = False #------if True, muon scale factors are sampled from a Gaussian
random_ele_SF = False #------if True, electron scale factors are sampled from a Gaussian
random_ph_SF  = False #------if True, photon scale factors are sampled from a Gaussian


############################################################################
#                                                                          #
#-------------------------- Integrated luminosity -------------------------#
#                                                                          #
############################################################################

#Normalize to this luminsity, in fb-1
#luminosity_norm = 36.46
luminosity_norm = 35.86
luminosity_BtoF = 19.72
luminosity_GH   = 16.14


############################################################################
#                                                                          #
#------------------------------ Create histos -----------------------------#
#                                                                          #
############################################################################

#Make signal histos larger
if isBDT:
    signal_magnify = 100
else:
    signal_magnify = 10000

output_dir = "plots"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#Here's the list of histos to plot
list_histos = ["h_mupt", "h_elept", "h_pipt", "h_gammaet","h_mupt_sig","h_elept_sig","h_pipt_sig","h_gammaet_sig", "h_mueta", "h_eleeta","h_pieta","h_gammaeta","h_mueta_sig","h_eleeta_sig","h_pieta_sig","h_gammaeta_sig", "h_nBjets_25","h_deltaphi_mu_pi","h_deltaphi_ele_pi","h_deltaphi_mu_W","h_deltaphi_ele_W","h_deltaeta_mu_pi","h_deltaeta_ele_pi","h_Wmass","h_Wmass_flag_mu","h_Wmass_flag_ele","h_mu_gamma_InvMass","h_ele_gamma_InvMass","h_piIso_05_mu","h_piRelIso_05_mu_ch","h_piRelIso_05_mu","h_piIso_05_ele","h_piRelIso_05_ele_ch","h_piRelIso_05_ele","h_met_mu","h_met_ele","h_met_puppi","h_Wmass_ratio_mu","h_Wmass_ratio_ele","h_nPV_mu","h_nPV_ele","h_motherID_mu","h_motherID_ele"]

Wmass_mu        = ROOT.TH1F("Wmass_mu","Wmass mu",10,40,100)
Wmass_ele       = ROOT.TH1F("Wmass_ele","Wmass ele",10,40,100)
Wmass_ratio_mu  = ROOT.TH1F("Wmass_ratio_mu","Wmass ratio mu",10,40,100)
Wmass_ratio_ele = ROOT.TH1F("Wmass_ratio_ele","Wmass ratio ele",10,40,100)

W_signal_hist = ROOT.TH1F("W_signal"," W signal",8,0,300)
W_ttbar_hist  = ROOT.TH1F("W_ttbar"," W ttbar",8,0,300)

Wmass_mu.Sumw2()
Wmass_ele.Sumw2()
Wmass_ratio_mu.Sumw2()
Wmass_ratio_ele.Sumw2()

W_signal_hist.Sumw2()
W_ttbar_hist.Sumw2()

# Color mask must have the same number of entries as non-QCD backgrounds
colors_mask = [ROOT.kYellow-8,ROOT.kAzure+7,ROOT.kOrange-3,ROOT.kTeal-5,ROOT.kViolet-6,ROOT.kSpring-1,ROOT.kOrange+6+1,ROOT.kMagenta+1,ROOT.kYellow-2,ROOT.kGreen+2,ROOT.kBlue-7,ROOT.kRed-7,ROOT.kOrange-2,ROOT.kYellow+3,ROOT.kBlue-4,ROOT.kGreen+3,ROOT.kCyan-7,ROOT.kPink+10,ROOT.kOrange+8,ROOT.kRed-7]

# Get the files and the names of the samples
samplename_list = myWF.get_samples_names()
root_file = myWF.get_root_files()

# Get the normalization
Norm_Map = myWF.get_normalizations_map()
ttbar_sig_calib_file = ROOT.TFile.Open("ttbar_signal_ratio.root")
ttbar_sig_calib = ttbar_sig_calib_file.Get("ttbar_signal_ratio")

##Get the handlers for all the histos and graphics
hs      = dict()
h_base  = dict()
canvas  = dict()

#####################################################################
#Event counters
motherID = dict()
mother_list = ["1","2","3","4","5","13","15","21","22","24","111","213","221","223","313","323","331","333","423","2212","3214"]

for mother in mother_list:
    motherID[mother] = 0.

#####################################################################

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

    h_base[theSampleName+list_histos[0]]  = ROOT.TH1F(theSampleName+list_histos[0], "p_{T} of the muon", 15, 25, 100.)
    h_base[theSampleName+list_histos[1]]  = ROOT.TH1F(theSampleName+list_histos[1], "p_{T} of the electron", 15, 26, 100.)
    h_base[theSampleName+list_histos[2]]  = ROOT.TH1F(theSampleName+list_histos[2], "p_{T} of the pion", 15, 20, 100.)
    h_base[theSampleName+list_histos[3]]  = ROOT.TH1F(theSampleName+list_histos[3], "E_{T} of the gamma", 15, 20, 100.)
    h_base[theSampleName+list_histos[4]]  = ROOT.TH1F(theSampleName+list_histos[4], "p_{T} of the muon - sig", 15, 20, 100.)
    h_base[theSampleName+list_histos[5]]  = ROOT.TH1F(theSampleName+list_histos[5], "p_{T} of the electron -sig", 15, 20, 100.)
    h_base[theSampleName+list_histos[6]]  = ROOT.TH1F(theSampleName+list_histos[6], "p_{T} of the pion -sig", 15, 20, 100.)
    h_base[theSampleName+list_histos[7]]  = ROOT.TH1F(theSampleName+list_histos[7], "E_{T} of the gamma -sig", 15, 20, 100.)
    h_base[theSampleName+list_histos[8]]  = ROOT.TH1F(theSampleName+list_histos[8], "eta of the muon", 20, -3, 3)
    h_base[theSampleName+list_histos[9]]  = ROOT.TH1F(theSampleName+list_histos[9], "eta of the electron", 20, -3, 3)
    h_base[theSampleName+list_histos[10]] = ROOT.TH1F(theSampleName+list_histos[10], "eta of the pion", 20, -3, 3)
    h_base[theSampleName+list_histos[11]] = ROOT.TH1F(theSampleName+list_histos[11], "eta of the photon", 20, -3, 3)
    h_base[theSampleName+list_histos[12]] = ROOT.TH1F(theSampleName+list_histos[12], "eta of the muon -sig", 20, -3, 3)
    h_base[theSampleName+list_histos[13]] = ROOT.TH1F(theSampleName+list_histos[13], "eta of the electron -sig", 20, -3, 3)
    h_base[theSampleName+list_histos[14]] = ROOT.TH1F(theSampleName+list_histos[14], "eta of the pion -sig", 20, -3, 3)
    h_base[theSampleName+list_histos[15]] = ROOT.TH1F(theSampleName+list_histos[15], "eta of the gamma -sig", 20, -3, 3)
    h_base[theSampleName+list_histos[16]] = ROOT.TH1F(theSampleName+list_histos[16], "n Bjets 25", 6, 0, 6.)
    h_base[theSampleName+list_histos[17]] = ROOT.TH1F(theSampleName+list_histos[17], "deltaphi mu-pi", 10, 0, 3.14)
    h_base[theSampleName+list_histos[18]] = ROOT.TH1F(theSampleName+list_histos[18], "deltaphi ele-pi", 10, 0, 3.14)
    h_base[theSampleName+list_histos[19]] = ROOT.TH1F(theSampleName+list_histos[19], "deltaphi mu-W", 10, 0, 3.14)
    h_base[theSampleName+list_histos[20]] = ROOT.TH1F(theSampleName+list_histos[20], "deltaphi ele-W", 10, 0, 3.14)
    h_base[theSampleName+list_histos[21]] = ROOT.TH1F(theSampleName+list_histos[21], "deltaeta mu-pi", 20, -5, 5)
    h_base[theSampleName+list_histos[22]] = ROOT.TH1F(theSampleName+list_histos[22], "deltaeta ele-pi", 20, -5, 5)
    h_base[theSampleName+list_histos[23]] = ROOT.TH1F(theSampleName+list_histos[23], "W mass", 15, 50, 100)
    h_base[theSampleName+list_histos[24]] = ROOT.TH1F(theSampleName+list_histos[24], "W mass if flag mu", 15, 50, 100)
    h_base[theSampleName+list_histos[25]] = ROOT.TH1F(theSampleName+list_histos[25], "W mass if flag ele", 15, 50, 100)
    h_base[theSampleName+list_histos[26]] = ROOT.TH1F(theSampleName+list_histos[26], "mu-gamma InvMass", 50, 0, 300)
    h_base[theSampleName+list_histos[27]] = ROOT.TH1F(theSampleName+list_histos[27], "ele-gamma InvMass", 50, 0, 300)
    h_base[theSampleName+list_histos[28]] = ROOT.TH1F(theSampleName+list_histos[28], "Pion isolation 05 - mu", 75, 0, 150)
    h_base[theSampleName+list_histos[29]] = ROOT.TH1F(theSampleName+list_histos[29], "Pion rel. isolation 05 - mu - ch", 50, 0, 10)
    h_base[theSampleName+list_histos[30]] = ROOT.TH1F(theSampleName+list_histos[30], "Pion rel. isolation 05 - mu", 50, 0, 10)
    h_base[theSampleName+list_histos[31]] = ROOT.TH1F(theSampleName+list_histos[31], "Pion isolation 05 - ele", 75, 0, 150)
    h_base[theSampleName+list_histos[32]] = ROOT.TH1F(theSampleName+list_histos[32], "Pion rel. isolation 05 - ele - ch", 50, 0, 10)
    h_base[theSampleName+list_histos[33]] = ROOT.TH1F(theSampleName+list_histos[33], "Pion rel. isolation 05 - ele", 50, 0, 10)
    h_base[theSampleName+list_histos[34]] = ROOT.TH1F(theSampleName+list_histos[34], "met mu", 40, 0, 200)
    h_base[theSampleName+list_histos[35]] = ROOT.TH1F(theSampleName+list_histos[35], "met ele", 40, 0, 200)
    h_base[theSampleName+list_histos[36]] = ROOT.TH1F(theSampleName+list_histos[36], "met puppi", 45, 0, 300)
    h_base[theSampleName+list_histos[37]] = ROOT.TH1F(theSampleName+list_histos[37], "Wmass ratio mu", 10, 40, 100)
    h_base[theSampleName+list_histos[38]] = ROOT.TH1F(theSampleName+list_histos[38], "Wmass ratio ele", 10, 40, 100)
    h_base[theSampleName+list_histos[39]] = ROOT.TH1F(theSampleName+list_histos[39], "nPV - mu", 50, 0, 50)
    h_base[theSampleName+list_histos[40]] = ROOT.TH1F(theSampleName+list_histos[40], "nPV - ele", 50, 0, 50)
    h_base[theSampleName+list_histos[41]] = ROOT.TH1F(theSampleName+list_histos[41], "mother ID - mu", 1000, 0, 1000)
    h_base[theSampleName+list_histos[42]] = ROOT.TH1F(theSampleName+list_histos[42], "mother ID - ele", 1000, 0, 1000)



# leg1 = ROOT.TLegend(0.15,0.6120093,0.34,0.9491917) #left positioning
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
Sevts_mu_SFvariation  = 0 # Counters for the number of signal events (weighted) when variating scale factors
Sevts_ele_SFvariation = 0
Sevts_weighted_mu = 0
Bevts_weighted_mu = 0
Sevts_weighted_ele = 0
Bevts_weighted_ele = 0
_Nrandom_for_SF = ROOT.TRandom3(44317)
_Nrandom_for_Gaus_SF = ROOT.TRandom3(44329)
N_WGToLNuG_mu = 0.
Nevts_per_sample = 0.
N_DoubleEMEnriched = 0.

#############---------------- BDT score cut values ----------------#############

BDT_OUT_MU = 0.190
BDT_OUT_ELE = 0.210

################################################################################

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

    Nevts_per_sample = 0. # Count the number of events survived per each sample processed

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

        if name_sample == "ttbar" and mytree.isttbarlnu: # Avoid double-counting of the ttbarlnu background
            continue

        # if "Data" in name_sample: continue  #-------------Excluding data-------------#

        # if not "Data" in name_sample: continue

        # if not (name_sample == "ttbar" or name_sample == "Signal"):
        #     continue

        # if name_sample == "ZGTo2LG" or name_sample == "TTGJets" or name_sample == "WGToLNuG":
        #     continue


        ############################################################################
        #                                                                          #
        #------------------------ Access the tree variables -----------------------#
        #                                                                          #
        ############################################################################

        isMuon = mytree.is_muon
        LepPiOppositeCharge = mytree.LepPiOppositeCharge

        nPV = mytree.nPV

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

        deltaeta_lep_pi = lep_eta-pi_eta
        
        deltaphi_lep_pi = math.fabs(lep_phi-pi_phi)
        if deltaphi_lep_pi > 3.14:
            deltaphi_lep_pi = 6.28 - deltaphi_lep_pi

        deltaphi_lep_W = math.fabs(lep_phi-W_phi)
        if deltaphi_lep_W > 3.14:
            deltaphi_lep_W = 6.28 - deltaphi_lep_W

        if not "Data" in name_sample:
            Wplus_pT = mytree.Wplus_pT
            Wminus_pT = mytree.Wminus_pT
            is_gen_ph = mytree.is_gen_ph
            gen_ph_pT = mytree.gen_ph_pT
            gen_ph_mother = str(math.fabs(mytree.gen_ph_mother))
            gen_ph_mother = gen_ph_mother.replace('.0','')

        if name_sample == "Signal":
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

        if myWF.post_preselection_cuts(lep_eta,lep_pT,isMuon,LepPiOppositeCharge):
            continue

        #if ("WJetsToLNu" in name_sample or "DY" in name_sample) and is_gen_ph:# and not gen_ph_mother == "13":
        #    continue
            
        if "Signal" in name_sample:
            Nsig_passed += 1
        if not "Signal" in name_sample and not "Data" in name_sample:
        #if "WJetsToLNu" in name_sample:
            Nbkg_passed += 1

        ############################################################################
        #                                                                          #
        #--------------------- Determine the total event weight -------------------#
        #                                                                          #
        ############################################################################
        
        if isMuon: # Get muon scale factors, which are different for two groups of datasets, and weight them for the respective integrated lumi 
            mu_weight_BtoF, mu_weight_BtoF_err = myWF.get_muon_scale_BtoF(lep_pT,lep_eta,isSingleMuTrigger_24)
            mu_weight_GH, mu_weight_GH_err     = myWF.get_muon_scale_GH(lep_pT,lep_eta,isSingleMuTrigger_24)

            # Use a random number to select which muon scale factor to use, depending on the respective lumi fraction
            Nrandom_for_SF = _Nrandom_for_SF.Rndm()

            if Nrandom_for_SF <= (luminosity_BtoF/luminosity_norm):  # Access muon SF: B to F

                if random_mu_SF:
                    mu_weight = _Nrandom_for_Gaus_SF.Gaus(mu_weight_BtoF,mu_weight_BtoF_err)
                else:
                    mu_weight = mu_weight_BtoF

            else: # Access muon SF: G and H
                
                if random_mu_SF:
                    mu_weight = _Nrandom_for_Gaus_SF.Gaus(mu_weight_GH,mu_weight_GH_err)
                else:
                    mu_weight = mu_weight_GH


        else:
            ele_weight, ele_weight_err = myWF.get_ele_scale(lep_pT,ele_etaSC)

            if random_ele_SF:
                ele_weight = _Nrandom_for_Gaus_SF.Gaus(ele_weight,ele_weight_err) 

        ph_weight, ph_weight_err = myWF.get_photon_scale(gamma_eT,gamma_etaSC)

        
        if random_ph_SF:
            ph_weight = _Nrandom_for_Gaus_SF.Gaus(ph_weight,ph_weight_err)
        

        if not "Data" in name_sample:
            PU_Weight = mytree.PU_Weight # Add Pile Up weight        
            Event_Weight = norm_factor*ph_weight*PU_Weight

            if isMuon:
                Event_Weight = Event_Weight*mu_weight
            else:
                Event_Weight = Event_Weight*ele_weight

            # Correct for the difference in pT of the generated W in Pythia and Madgraph samples
            if name_sample == "Signal" and is_signal_Wplus:
                local_Wplus_pT = Wplus_pT
                if Wplus_pT > 300.:
                    local_Wplus_pT = 300.

                Event_Weight = Event_Weight#*ttbar_sig_calib.GetBinContent(ttbar_sig_calib.GetXaxis().FindBin(local_Wplus_pT))


            if name_sample == "Signal" and is_signal_Wminus:
                local_Wminus_pT = Wminus_pT
                if Wminus_pT > 300.:
                    local_Wminus_pT = 300.

                Event_Weight = Event_Weight#*ttbar_sig_calib.GetBinContent(ttbar_sig_calib.GetXaxis().FindBin(local_Wminus_pT))

        else:
            Event_Weight = 1.


        Nevts_per_sample += Event_Weight # Increment the number of events survived in the analyzed sample


        ############################################################################
        #                                                                          #
        #-------------------------- Retrieve BDT output  --------------------------#
        #                                                                          #
        ############################################################################

        if isBDT_with_Wmass:
            BDT_out = myWF.get_BDT_output(pi_pT/Wmass,gamma_eT/Wmass,nBjets_25,lep_pT,piRelIso_05_ch,met,isMuon)  
        else:
            BDT_out = myWF.get_BDT_output(pi_pT,gamma_eT,nBjets_25,lep_pT,piRelIso_05_ch,met,isMuon)

   
        ############################################################################
        #                                                                          #
        #------------------------------- Fill histos ------------------------------#
        #                                                                          #
        ############################################################################

        if (isMuon and BDT_out >= -0.1 and BDT_out < BDT_OUT_MU) or (not isMuon and BDT_out >= -0.1 and BDT_out < BDT_OUT_ELE): # Alternative cut on BDT score. Two histos will be filled. Eventually, they will be the ratio of Wmass distributions: one with pi_pT and gamma_ET normalized to Wmass, one not
            if (Wmass >= 50. and Wmass <= 100.) and isMuon:
                h_base[theSampleName+"h_Wmass_ratio_mu"].Fill(Wmass,Event_Weight)
            if (Wmass >= 50. and Wmass <= 100.) and not isMuon:
                h_base[theSampleName+"h_Wmass_ratio_ele"].Fill(Wmass,Event_Weight)

            
        if (isBDT and isMuon and BDT_out < BDT_OUT_MU) or (isBDT and not isMuon and BDT_out < BDT_OUT_ELE): # Cut on BDT output
            continue
        if isBDT and (Wmass < 50. or Wmass > 100.): # General Wmass condition
            continue
        if isBDT and name_sample=="Data" and (Wmass >= 65. and Wmass <= 90.): # Exclude data in the Blind Window
            continue

        h_base[theSampleName+"h_nBjets_25"].Fill(nBjets_25,Event_Weight)
        h_base[theSampleName+"h_met_puppi"].Fill(met_puppi,Event_Weight)
        h_base[theSampleName+"h_Wmass"].Fill(Wmass,Event_Weight)
        h_base[theSampleName+"h_pipt"].Fill(pi_pT,Event_Weight)
        h_base[theSampleName+"h_pieta"].Fill(pi_eta,Event_Weight)
        h_base[theSampleName+"h_gammaet"].Fill(gamma_eT,Event_Weight)
        h_base[theSampleName+"h_gammaeta"].Fill(gamma_eta,Event_Weight)
        

        if isMuon:
            h_base[theSampleName+"h_Wmass_flag_mu"].Fill(Wmass,Event_Weight)
            h_base[theSampleName+"h_mupt"].Fill(lep_pT,Event_Weight)
            h_base[theSampleName+"h_mueta"].Fill(lep_eta,Event_Weight)
            h_base[theSampleName+"h_nPV_mu"].Fill(nPV,Event_Weight)
            h_base[theSampleName+"h_piRelIso_05_mu"].Fill(piRelIso_05,Event_Weight)
            h_base[theSampleName+"h_piRelIso_05_mu_ch"].Fill(piRelIso_05_ch,Event_Weight)
            h_base[theSampleName+"h_met_mu"].Fill(met,Event_Weight)
            h_base[theSampleName+"h_deltaphi_mu_W"].Fill(deltaphi_lep_W,Event_Weight)
            h_base[theSampleName+"h_mu_gamma_InvMass"].Fill(mu_gamma_InvMass,Event_Weight)  
            h_base[theSampleName+"h_deltaeta_mu_pi"].Fill(deltaeta_lep_pi,Event_Weight)

            if isBDT and"Signal" in name_sample:
                Sevts_mu_SFvariation += Event_Weight

        if not isMuon:
            h_base[theSampleName+"h_Wmass_flag_ele"].Fill(Wmass,Event_Weight)
            h_base[theSampleName+"h_elept"].Fill(lep_pT,Event_Weight)
            h_base[theSampleName+"h_eleeta"].Fill(lep_eta,Event_Weight)
            h_base[theSampleName+"h_nPV_ele"].Fill(nPV,Event_Weight)
            h_base[theSampleName+"h_piRelIso_05_ele"].Fill(piRelIso_05,Event_Weight)
            h_base[theSampleName+"h_piRelIso_05_ele_ch"].Fill(piRelIso_05_ch,Event_Weight)
            h_base[theSampleName+"h_ele_gamma_InvMass"].Fill(ele_gamma_InvMass,Event_Weight)
            h_base[theSampleName+"h_met_ele"].Fill(met,Event_Weight)
            h_base[theSampleName+"h_deltaeta_ele_pi"].Fill(deltaeta_lep_pi,Event_Weight)
            h_base[theSampleName+"h_deltaphi_ele_W"].Fill(deltaphi_lep_W,Event_Weight)

            if isBDT and "Signal" in name_sample:
                Sevts_ele_SFvariation += Event_Weight


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


        # Obtain the number of sig and bkg events (weighted)
        if not "Data" in name_sample:
            if "Signal" in name_sample and isMuon:
                Sevts_weighted_mu += Event_Weight
            #if "WJetsToLNu" in name_sample and isMuon:
            if not "Signal" in name_sample and isMuon:
                Bevts_weighted_mu += Event_Weight
            if "Signal" in name_sample and not isMuon:
                Sevts_weighted_ele += Event_Weight
            #if "WJetsToLNu" in name_sample and not isMuon:
            if not "Signal" in name_sample and not isMuon:
                Bevts_weighted_ele += Event_Weight


        if "DoubleEMEnriched" in name_sample:
            N_DoubleEMEnriched += Event_Weight

        if "WGToLNuG" in name_sample and is_gen_ph and gen_ph_pT > 30.:# and isMuon:
            motherID[gen_ph_mother] += Event_Weight 


        ############################################################################
        #                                                                          #
        #------------------- Fill Signal-ttbar comparison histos ------------------#
        #                                                                          #
        ############################################################################

        ########----------- To be filled in "preselection" mode! -----------########
        if name_sample == "Signal" and is_signal_Wplus:
            if Wplus_pT > 300.:
                W_signal_hist.Fill(300.,Event_Weight)
            else:
                W_signal_hist.Fill(Wplus_pT,Event_Weight)

        if name_sample == "Signal" and is_signal_Wminus:
            if Wminus_pT > 300.:
                W_signal_hist.Fill(300.,Event_Weight)
            else:
                W_signal_hist.Fill(Wminus_pT,Event_Weight)

        if name_sample == "ttbar":
            if Wplus_pT > 300.:
                W_ttbar_hist.Fill(300.,Event_Weight)
            else:
                W_ttbar_hist.Fill(Wplus_pT,Event_Weight)

            if Wminus_pT > 300.:
                W_ttbar_hist.Fill(300.,Event_Weight)
            else:
                W_ttbar_hist.Fill(Wminus_pT,Event_Weight)


        ############################################################################
        #                                                                          #
        #--------------------- Define histos and legend features ------------------#
        #                                                                          #
        ############################################################################


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


        if idx_histo == 0 and (Nevts_per_sample > 800 or name_sample == myWF.sig_samplename):
            if QCDflag and isFirstQCDlegend:
                leg1.AddEntry(h_base[theSampleName+hname],"QCD","f")
                isFirstQCDlegend = False
            elif name_sample == myWF.sig_samplename:
                sample_legend_name = str(signal_magnify) + " x " + name_sample
                leg1.AddEntry(h_base[name_sample+hname],sample_legend_name,"f")  #To comment when signal has to be excluded.
                #leg2.AddEntry(h_base[name_sample+hname],sample_legend_name,"f")
            elif name_sample == myWF.data_samplename:
                leg1.AddEntry(h_base[name_sample+hname],name_sample,"lep") # lep shows on the TLegend a point with errors to indicate data
            elif not QCDflag and not name_sample == myWF.data_samplename:
                leg1.AddEntry(h_base[theSampleName+hname],theSampleName,"f")

        if not QCDflag and not name_sample == myWF.sig_samplename and not name_sample == myWF.data_samplename:
            hs[hname].Add(h_base[theSampleName+hname])

    if not QCDflag and not name_sample == myWF.sig_samplename and not name_sample == myWF.data_samplename:
        idx_sample += 1

print "Finished runnning over samples!"



##################### Plotting all the backgrounds in separate histos for a given distribution ######################

color_index = 0
eleeta_canvas = dict()

for sample_name in samplename_list:
    
    if "QCD" in sample_name or "Data" in sample_name or "Signal" in sample_name:
        continue
    else:
        eleeta_canvas[sample_name] = ROOT.TCanvas(sample_name,sample_name,200,106,600,600)
        h_base[sample_name+"h_eleeta"].SetTitle(sample_name)
        h_base[sample_name+"h_eleeta"].SetFillColor(colors_mask[color_index])
        h_base[sample_name+"h_eleeta"].Draw("hist")
        print sample_name, " integral: ", h_base[sample_name+"h_eleeta"].Integral()
        color_index += 1
        eleeta_canvas[sample_name].SaveAs("plots/eleeta/h_eleeta_" + sample_name + ".pdf")

 ####################################################################################################################
        
    

for idx_histo,hname in enumerate(list_histos):
    hs[hname].Add(h_base["QCD_"+hname])

for hname in list_histos:

    canvas[hname].cd()

    hs[hname].Draw("histo")

    if "h_Wmass_" in hname:
        hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),65.))
    if hname == "h_Wmass":
        hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),100.))
    if hname == "h_piRelIso_05_ele_ch":
        hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),10000.))
    if hname == "h_piRelIso_05_mu_ch":
        hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),10000.))
    if hname == "h_piRelIso_05_mu_ch_AfterCut":
        hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),130.))
    if hname == "h_piRelIso_05_ele_ch_AfterCut":
        hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),130.))
    if hname == "h_mueta":
        hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),90.))
    if hname == "h_eleeta":
        hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),90.))
    if hname == "h_pieta":
        hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),90.))
    if hname == "h_gammaeta":
        hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),90.))
    if hname == "h_pipt":
        hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),90.))
    if "pt" in hname:
        hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),90.))
    if hname == "h_gammaet":
        hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),100.))
    #if hname == "h_pieta_sig" or hname == "h_pipt_sig" or hname == "h_gammaeta_sig" or hname == "h_gammaet_sig":
    #    hs[hname].SetMaximum(30.)
    #if hname == "h_mueta_sig" or hname == "h_mupt_sig" or hname == "h_eleeta_sig" or hname == "h_elept_sig":
    #    hs[hname].SetMaximum(20.)


    #---------Histos names---------#
    
    hs[hname].SetTitle(" ")
    hs[hname].GetYaxis().SetTitleOffset(1.7)
    #hs[hname].GetYaxis().SetTitle("Events")
    
    if hname == "h_Wmass" or hname == "h_Wmass_flag_mu" or hname == "h_Wmass_flag_ele":
        hs[hname].GetXaxis().SetTitle("m_{#pi#gamma} (GeV)")

    if hname == "h_mupt" or hname == "h_mupt_sig":
        hs[hname].GetXaxis().SetTitle("p_{T}^{#mu} (GeV)")

    if hname == "h_elept" or hname == "h_elept_sig":
        hs[hname].GetXaxis().SetTitle("p_{T}^{e} (GeV)")
    
    if hname == "h_pipt" or hname == "h_pipt_sig":
        hs[hname].GetXaxis().SetTitle("p_{T}^{#pi} (GeV)")

    if hname == "h_gammaet" or hname == "h_gammaet_sig":
        hs[hname].GetXaxis().SetTitle("E_{T}^{#gamma} (GeV)")

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

    if hname == "h_Wmass_ratio_mu":
        Wmass_ratio_mu = hs[hname].GetStack().Last()

    if hname == "h_Wmass_ratio_ele":
        Wmass_ratio_ele = hs[hname].GetStack().Last()


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

#print "Wmass_mu: ", Wmass_mu.Integral()
#print "Wmass_ele: ", Wmass_ele.Integral()
Wmass_mu.Scale(1/Wmass_mu.Integral())
Wmass_ratio_mu.Scale(1/Wmass_ratio_mu.Integral())
Wmass_ele.Scale(1/Wmass_ele.Integral())
Wmass_ratio_ele.Scale(1/Wmass_ratio_ele.Integral())

Wmass_ratio_mu.Divide(Wmass_ratio_mu,Wmass_mu,1.0,1.0,"B")
Wmass_ratio_ele.Divide(Wmass_ratio_ele,Wmass_ele,1.0,1.0,"B")


canvas5 = ROOT.TCanvas()
Wmass_ratio_mu.SetMarkerStyle(21)
Wmass_ratio_mu.SetTitle(" ")
Wmass_ratio_mu.GetXaxis().SetTitle("m_{#pi#gamma} (GeV)")
Wmass_ratio_mu.Draw("Pe")
canvas5.SaveAs("plots/Wmass_ratio_mu.pdf")


canvas6 = ROOT.TCanvas()
Wmass_ratio_ele.SetMarkerStyle(21)
Wmass_ratio_ele.SetTitle(" ")
Wmass_ratio_ele.GetXaxis().SetTitle("m_{#pi#gamma} (GeV)")
Wmass_ratio_ele.Draw("Pe")
canvas6.SaveAs("plots/Wmass_ratio_ele.pdf")


#---------- Signal-ttbar comparison ----------#

leg_sig_ttbar = ROOT.TLegend(0.6868687,0.6120093,0.86,0.85)
leg_sig_ttbar.SetHeader(" ")
leg_sig_ttbar.SetFillColor(0)
leg_sig_ttbar.SetBorderSize(0)
leg_sig_ttbar.SetLineColor(1)
leg_sig_ttbar.SetLineStyle(1)
leg_sig_ttbar.SetLineWidth(1)
leg_sig_ttbar.SetFillStyle(1001)
leg_sig_ttbar.SetTextSize(0.04)
leg_sig_ttbar.AddEntry(W_signal_hist,"Signal","l")
leg_sig_ttbar.AddEntry(W_ttbar_hist,"ttbar","l")

W_signal_hist.Scale(1/W_signal_hist.Integral())
W_ttbar_hist.Scale(1/W_ttbar_hist.Integral())

canvas7 = ROOT.TCanvas()
ROOT.gStyle.SetOptStat(0)
W_ttbar_hist.SetLineColor(46)
W_signal_hist.SetLineColor(38)
W_ttbar_hist.SetTitle(" ")
W_signal_hist.SetTitle(" ")
W_ttbar_hist.GetXaxis().SetTitle("p_{T}^{W} (GeV)")
W_signal_hist.GetXaxis().SetTitle("p_{T}^{W} (GeV)")
W_ttbar_hist.GetYaxis().SetTitle("Probability density")
W_signal_hist.GetYaxis().SetTitle("Probability density")
W_signal_hist.Draw("hist")
W_ttbar_hist.Draw("hist,SAME")
ROOT.gPad.SetLogy()
leg_sig_ttbar.Draw("SAME")
canvas7.SaveAs("plots/ttbar_signal.pdf")

W_ttbar_hist.Divide(W_signal_hist)

canvas8 = ROOT.TCanvas()
# ttbar_sig_file = ROOT.TFile("ttbar_signal_ratio.root","RECREATE")
ROOT.gStyle.SetOptStat(0)
W_ttbar_hist.SetMarkerStyle(21)
W_ttbar_hist.SetTitle(" ")
W_ttbar_hist.GetXaxis().SetTitle("p_{T}^{W} (GeV)")
W_ttbar_hist.Draw("Pe")
canvas8.SaveAs("plots/ttbar_signal_ratio.pdf")
# W_ttbar_hist.Write("ttbar_signal_ratio")
# ttbar_sig_file.Close()

print "Number of expected events for ", luminosity_norm, " in fb-1"
print "Number of signal events passed = ", Nsig_passed
print "Number of background events passed = ", Nbkg_passed
print "number of Signal events in muon channel: ", Sevts_mu_SFvariation
print "number of Signal events in electron channel: ", Sevts_ele_SFvariation
print "total number of S evts weighted -mu: ", Sevts_weighted_mu
print "total number of B (WJetsToLNu) evts weighted -mu: ", Bevts_weighted_mu
print "total number of S evts weighted -ele: ", Sevts_weighted_ele
print "total number of B (WJetsToLNu) evts weighted -ele: ", Bevts_weighted_ele

print "N_WGToLNuG_mu", N_WGToLNuG_mu
print "N_DoubleEMEnriched", N_DoubleEMEnriched

#for mother in mother_list:
#    print mother, " --- ", motherID[mother] 
