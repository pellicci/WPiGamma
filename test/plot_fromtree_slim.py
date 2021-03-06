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
p.add_argument('runningEra_option', help='Type <<0>> for 2016, <<1>> for 2017, <<2>> for 2018, <<3>> for combination 2016+2017')
args = p.parse_args()

# Switch from muon to electron channel, and from 2016 to 2017
if args.isBDT_option == "preselection":
    isBDT = False
if args.isBDT_option == "BDT":
    isBDT = True
runningEra = int(args.runningEra_option)
#---------------------------------#

isBDT_with_Wmass = False # If true, pT(pi) and ET(gamma) in the BDT are normalized to Wmass 
myWF = Workflow_Handler("Signal","Data",isBDT_with_Wmass,runningEra)

#Bools for rundom SF variation
random_mu_SF  = False #------if True, muon scale factors are sampled from a Gaussian
random_ele_SF = False #------if True, electron scale factors are sampled from a Gaussian
random_ph_SF  = False #------if True, photon scale factors are sampled from a Gaussian


############################################################################
#                                                                          #
#-------------------------- Integrated luminosity -------------------------#
#                                                                          #
############################################################################

#Normalize to this luminsity, in fb-1

luminosity_norm_2016 = 35.86
luminosity_norm_2017 = 41.529

if runningEra == 0:
    luminosity_norm = luminosity_norm_2016
if runningEra == 1:
    luminosity_norm = luminosity_norm_2017
if runningEra == 3:
    luminosity_norm = luminosity_norm_2016 + luminosity_norm_2017

luminosity_BtoF = 19.72 #For 2016
luminosity_GH   = 16.14 #For 2016


#############---------------- BDT score cut values ----------------#############

BDT_OUT_MU  = 0.220
BDT_OUT_ELE = 0.170

################################################################################






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
list_histos = ["h_mupt", "h_elept", "h_pipt", "h_gammaet","h_mupt_sig","h_elept_sig","h_pipt_sig","h_gammaet_sig", "h_mueta", "h_eleeta","h_pieta","h_gammaeta","h_mueta_sig","h_eleeta_sig","h_pieta_sig","h_gammaeta_sig", "h_nBjets_25","h_deltaphi_mu_pi","h_deltaphi_ele_pi","h_deltaphi_mu_W","h_deltaphi_ele_W","h_deltaeta_mu_pi","h_deltaeta_ele_pi","h_Wmass","h_Wmass_flag_mu","h_Wmass_flag_ele","h_mu_gamma_InvMass","h_ele_gamma_InvMass","h_piIso_05_mu","h_piRelIso_05_mu_ch","h_piRelIso_05_mu","h_piIso_05_ele","h_piRelIso_05_ele_ch","h_piRelIso_05_ele","h_met_mu","h_met_ele","h_met_puppi","h_Wmass_ratio_mu","h_Wmass_ratio_ele","h_nPV_mu","h_nPV_ele","h_deltaphi_mu_gamma","h_deltaphi_ele_gamma"]

Wmass_mu        = ROOT.TH1F("Wmass_mu","Wmass mu",10,50,100)
Wmass_ele       = ROOT.TH1F("Wmass_ele","Wmass ele",10,50,100)
Wmass_ratio_mu  = ROOT.TH1F("Wmass_ratio_mu","Wmass ratio mu",10,50,100)
Wmass_ratio_ele = ROOT.TH1F("Wmass_ratio_ele","Wmass ratio ele",10,50,100)

W_signal_hist = ROOT.TH1F("W_signal"," W signal",8,0,300)
W_ttbar_hist  = ROOT.TH1F("W_ttbar"," W ttbar",8,0,300)

h_PUdistrib = ROOT.TH1F("pile up", "pile up",75,0,75)
h_PUdistrib.Sumw2()

Wmass_mu.Sumw2()
Wmass_ele.Sumw2()
Wmass_ratio_mu.Sumw2()
Wmass_ratio_ele.Sumw2()

W_signal_hist.Sumw2()
W_ttbar_hist.Sumw2()

# Color mask must have the same number of entries as non-QCD backgrounds
colors_mask = dict()
colors_mask["ttbar"]               = ROOT.kYellow-8
colors_mask["ttbarlnu"]            = ROOT.kAzure+7
colors_mask["ttbarWQQ"]            = ROOT.kOrange-3
colors_mask["ttbarZlnu"]           = ROOT.kTeal-5
colors_mask["DY50"]                = ROOT.kViolet-6
colors_mask["ttbarZQQ"]            = ROOT.kSpring-1
colors_mask["GammaJets20to40"]     = ROOT.kOrange+7
colors_mask["SingleToptW"]         = ROOT.kMagenta+1
colors_mask["WJetsToLNu"]          = ROOT.kGreen+2
colors_mask["SingleAntiToptW"]     = ROOT.kYellow-2
colors_mask["WZ"]                  = ROOT.kBlue-7
colors_mask["WGToLNuG"]            = ROOT.kRed-7
colors_mask["TTGJets"]             = ROOT.kOrange-2
colors_mask["ZGTo2LG"]             = ROOT.kYellow+3
colors_mask["DY10to50"]            = ROOT.kBlue-4
colors_mask["GammaJets20toInf"]    = ROOT.kGreen+3
colors_mask["GammaJets40toInf"]    = ROOT.kCyan-7
colors_mask["WW"]                  = ROOT.kPink+1
colors_mask["ttbarWlnu"]           = ROOT.kOrange+8
colors_mask["ttbarToHadronic"]     = ROOT.kRed-7
colors_mask["ttbarToSemiLeptonic"] = ROOT.kBlue+3


# Get the files and the names of the samples
samplename_list = myWF.get_samples_names()
print "samplename_list: ", samplename_list#, "  sampleEra_list: ", sampleEra_list
root_file = myWF.get_root_files()
#print "root_file: ",root_file

# Get the normalization
#Norm_Map = myWF.get_normalizations_map()
Norm_Map_2016, Norm_Map_2017 = myWF.get_normalizations_map()


ttbar_sig_calib_file = ROOT.TFile.Open("ttbar_signal_ratio.root")
ttbar_sig_calib = ttbar_sig_calib_file.Get("ttbar_signal_ratio")

##Get the handlers for all the histos and graphics
hs      = dict()
h_base  = dict()
canvas  = dict()

#####################################################################

for hname in list_histos:
    hs[hname] = ROOT.THStack("hs_" + hname,"")
    canvas[hname] = ROOT.TCanvas(hname,hname,200,106,600,600)

##Define the histos to be created
isQCDfirst = True

double_samplenames = dict() #Identify the samples which appear twice or more with the same name
for name_sample_withYear in samplename_list:
    name_sample = name_sample_withYear.split("_")[0]
    double_samplenames[name_sample] = 0

#for sample_name, sampleEra in zip(samplename_list,sampleEra_list):
for name_sample_withYear in samplename_list:

    # if samplename_list.count(name_sample) == 1:
    #     print "name of samples which appear only once in the list: ", name_sample

    ###################################################
    #                                                 #
    #------- Select which era will be plotted --------#
    #                                                 #
    ###################################################

    # if runningEra == 0 and not sampleEra == "2016":
    #     continue
    # if runningEra == 1 and not sampleEra == "2017":
    #     continue
    # if runningEra == 2 and not sampleEra == "2018":
    #     continue
    # if runningEra == 3 and not (sampleEra == "2016" or sampleEra == "2017"):
    #     continue

    if runningEra == 0 and not "2016" in name_sample_withYear:
        continue
    if runningEra == 1 and not "2017" in name_sample_withYear:
        continue
    if runningEra == 2 and not "2018" in name_sample_withYear:
        continue
    if runningEra == 3 and not ("2016" or "2017") in name_sample_withYear:
        continue

    name_sample = name_sample_withYear.split("_")[0]
    print "nome del sample: ", name_sample
    sampleEra   = name_sample_withYear.split("_")[1]
    print "era del sample: ", sampleEra

    #Prevent memory leaks from creation of two or more histos with same name. At the same time, account for samples which are unique among the eras
    double_samplenames[name_sample] += 1 
    if double_samplenames[name_sample] > 1:
        continue
    #######################################################################

    if "QCD" in name_sample:
        theSampleName = "QCD_"
        if not isQCDfirst:
            continue
        isQCDfirst = False
    else:
        theSampleName = name_sample

    h_base[theSampleName+list_histos[0]]  = ROOT.TH1F(theSampleName+list_histos[0], "p_{T} of the muon", 15, 25, 100.)
    h_base[theSampleName+list_histos[1]]  = ROOT.TH1F(theSampleName+list_histos[1], "p_{T} of the electron", 15, 28, 100.)
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
    h_base[theSampleName+list_histos[23]] = ROOT.TH1F(theSampleName+list_histos[23], "W mass", 10, 50, 100)
    h_base[theSampleName+list_histos[24]] = ROOT.TH1F(theSampleName+list_histos[24], "W mass if flag mu", 10, 50, 100)
    h_base[theSampleName+list_histos[25]] = ROOT.TH1F(theSampleName+list_histos[25], "W mass if flag ele", 10, 50, 100)
    h_base[theSampleName+list_histos[26]] = ROOT.TH1F(theSampleName+list_histos[26], "mu-gamma InvMass", 20, 0, 300)
    h_base[theSampleName+list_histos[27]] = ROOT.TH1F(theSampleName+list_histos[27], "ele-gamma InvMass", 20, 0, 300)
    h_base[theSampleName+list_histos[28]] = ROOT.TH1F(theSampleName+list_histos[28], "Pion isolation 05 - mu", 75, 0, 150)
    h_base[theSampleName+list_histos[29]] = ROOT.TH1F(theSampleName+list_histos[29], "Pion rel. isolation 05 - mu - ch", 50, 0, 10)
    h_base[theSampleName+list_histos[30]] = ROOT.TH1F(theSampleName+list_histos[30], "Pion rel. isolation 05 - mu", 50, 0, 10)
    h_base[theSampleName+list_histos[31]] = ROOT.TH1F(theSampleName+list_histos[31], "Pion isolation 05 - ele", 75, 0, 150)
    h_base[theSampleName+list_histos[32]] = ROOT.TH1F(theSampleName+list_histos[32], "Pion rel. isolation 05 - ele - ch", 50, 0, 10)
    h_base[theSampleName+list_histos[33]] = ROOT.TH1F(theSampleName+list_histos[33], "Pion rel. isolation 05 - ele", 50, 0, 10)
    h_base[theSampleName+list_histos[34]] = ROOT.TH1F(theSampleName+list_histos[34], "met mu", 20, 0, 200)
    h_base[theSampleName+list_histos[35]] = ROOT.TH1F(theSampleName+list_histos[35], "met ele", 20, 0, 200)
    h_base[theSampleName+list_histos[36]] = ROOT.TH1F(theSampleName+list_histos[36], "met puppi", 45, 0, 300)
    h_base[theSampleName+list_histos[37]] = ROOT.TH1F(theSampleName+list_histos[37], "Wmass ratio mu", 10, 50, 100)
    h_base[theSampleName+list_histos[38]] = ROOT.TH1F(theSampleName+list_histos[38], "Wmass ratio ele", 10, 50, 100)
    h_base[theSampleName+list_histos[39]] = ROOT.TH1F(theSampleName+list_histos[39], "nPV - mu", 15, 0, 50)
    h_base[theSampleName+list_histos[40]] = ROOT.TH1F(theSampleName+list_histos[40], "nPV - ele", 15, 0, 50)
    h_base[theSampleName+list_histos[41]] = ROOT.TH1F(theSampleName+list_histos[41], "deltaphi mu-gamma", 10, 0, 3.14)
    h_base[theSampleName+list_histos[42]] = ROOT.TH1F(theSampleName+list_histos[42], "deltaphi ele-gamma", 10, 0, 3.14)




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
#_Nrandom_for_SF = ROOT.TRandom3(44317)
_Nrandom_for_Gaus_SF = ROOT.TRandom3(44329)
N_WGToLNuG_mu = 0.
Nevts_per_sample = 0.
N_DoubleEMEnriched = 0.

##Loop on samples, and then on events, and merge QCD stuff
idx_sample = 0
isFirstQCDlegend = True


# for name_sample, sampleEra in zip(samplename_list,sampleEra_list):
for name_sample_withYear in samplename_list:

    ###################################################
    #                                                 #
    #------- Select which era will be plotted --------#
    #                                                 #
    ###################################################

    # if runningEra == 0 and not sampleEra == "2016":
    #     continue
    # if runningEra == 1 and not sampleEra == "2017":
    #     continue
    # if runningEra == 2 and not sampleEra == "2018":
    #     continue
    # if runningEra == 3 and not (sampleEra == "2016" or sampleEra == "2017"):
    #     continue

    if runningEra == 0 and not "2016" in name_sample_withYear:
        continue
    if runningEra == 1 and not "2017" in name_sample_withYear:
        continue
    if runningEra == 2 and not "2018" in name_sample_withYear:
        continue
    if runningEra == 3 and not ("2016" or "2017") in name_sample_withYear:
        continue

    name_sample = name_sample_withYear.split("_")[0]
    sampleEra   = name_sample_withYear.split("_")[1]

    theSampleName = name_sample

    if "QCD" in name_sample:
        QCDflag = True
        theSampleName = "QCD_"
    else:
        QCDflag = False

    #######################################################################

    if not "Data" in name_sample:
        
        #norm_factor = Norm_Map[name_sample]*luminosity_norm
        if sampleEra == "2016":
            norm_factor = Norm_Map_2016[name_sample]*luminosity_norm
        if sampleEra == "2017":
            norm_factor = Norm_Map_2017[name_sample]*luminosity_norm


    mytree = root_file[name_sample+"_"+sampleEra].Get("WPiGammaAnalysis/mytree")
    #print "Rootfile: ", root_file[name_sample+"_"+sampleEra] #Verify to get the right rootfile for the processed sample
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

        # if not "Signal" in name_sample: continue

        # if not (name_sample == "ttbar" or name_sample == "Signal"):
        #     continue

        # if name_sample == "ZGTo2LG" or name_sample == "TTGJets" or name_sample == "ttbarlnu" or name_sample == "ttbar":# or name_sample == "WGToLNuG" or:
        #     continue


        ############################################################################
        #                                                                          #
        #------------------------ Access the tree variables -----------------------#
        #                                                                          #
        ############################################################################

        #sampleEra = mytree.runningEra # Define on which year we work
        isMuon = mytree.is_muon
        LepPiOppositeCharge = mytree.LepPiOppositeCharge

        nPV = mytree.nPV

        isSingleMuTrigger_24 = mytree.isSingleMuTrigger_24
        isSingleMuTrigger_50 = mytree.isSingleMuTrigger_50
        if sampleEra == "2017":
            isSingleMuTrigger_27 = mytree.isSingleMuTrigger_27           

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
        
        deltaphi_lep_pi = math.fabs(lep_phi - pi_phi)
        if deltaphi_lep_pi > 3.14:
            deltaphi_lep_pi = 6.28 - deltaphi_lep_pi

        deltaphi_lep_W = math.fabs(lep_phi - W_phi)
        if deltaphi_lep_W > 3.14:
            deltaphi_lep_W = 6.28 - deltaphi_lep_W

        deltaphi_lep_gamma = math.fabs(lep_phi - gamma_phi)
        if deltaphi_lep_gamma > 3.14:
            deltaphi_lep_gamma = 6.28 - deltaphi_lep_gamma

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

        if myWF.post_preselection_cuts(lep_eta,lep_pT,isMuon,LepPiOppositeCharge,sampleEra):
           continue
            
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


        ################ MUON SFs ################

        if isMuon: # Get muon scale factors, which are different for two groups of datasets, and weight them for the respective integrated lumi 
            if sampleEra == "2016":
                isSingleMuTrigger_LOW = isSingleMuTrigger_24
            if sampleEra == "2017":
                isSingleMuTrigger_LOW = isSingleMuTrigger_27

            mu_weight, mu_weight_err = myWF.get_muon_scale(lep_pT,lep_eta,isSingleMuTrigger_LOW,sampleEra)


            if random_mu_SF:
                mu_weight = _Nrandom_for_Gaus_SF.Gaus(mu_weight,mu_weight_err)


        ############## ELECTRON SFs ##############

        else:
            ele_weight, ele_weight_err = myWF.get_ele_scale(lep_pT,ele_etaSC,sampleEra)

            if random_ele_SF:
                ele_weight = _Nrandom_for_Gaus_SF.Gaus(ele_weight,ele_weight_err) 


        ############### PHOTON SFs ###############

        ph_weight, ph_weight_err = myWF.get_photon_scale(gamma_eT,gamma_etaSC,sampleEra)

        
        if random_ph_SF:
            ph_weight = _Nrandom_for_Gaus_SF.Gaus(ph_weight,ph_weight_err)
        

        ############### Multiply weights and SFs for MC. Set weight to 1 for data ###############

        if not "Data" in name_sample:
            MC_Weight = mytree.MC_Weight # Add MC weight        
            PU_Weight = mytree.PU_Weight # Add Pile Up weight        
            Event_Weight = norm_factor*ph_weight*MC_Weight*PU_Weight/math.fabs(MC_Weight) # Just take the sign of the gen weight

            if isMuon:
                Event_Weight = Event_Weight*mu_weight
            else:
                Event_Weight = Event_Weight*ele_weight

            # Correct for the difference in pT of the generated W in Pythia and Madgraph samples
            if sampleEra == "2016" and name_sample == "Signal" and is_signal_Wplus:
                local_Wplus_pT = Wplus_pT
                if Wplus_pT > 300.:
                    local_Wplus_pT = 300.

                Event_Weight = Event_Weight*ttbar_sig_calib.GetBinContent(ttbar_sig_calib.GetXaxis().FindBin(local_Wplus_pT))


            if sampleEra == "2016" and name_sample == "Signal" and is_signal_Wminus:
                local_Wminus_pT = Wminus_pT
                if Wminus_pT > 300.:
                    local_Wminus_pT = 300.

                Event_Weight = Event_Weight*ttbar_sig_calib.GetBinContent(ttbar_sig_calib.GetXaxis().FindBin(local_Wminus_pT))

        else:
            Event_Weight = 1.


        Nevts_per_sample += Event_Weight # Increment the number of events survived in the analyzed sample


        ############################################################################
        #                                                                          #
        #-------------------------- Retrieve BDT output  --------------------------#
        #                                                                          #
        ############################################################################

        if isBDT_with_Wmass:
            BDT_out = myWF.get_BDT_output(pi_pT/Wmass,gamma_eT/Wmass,nBjets_25,lep_pT,piRelIso_05_ch,met,deltaphi_lep_pi,isMuon)  
        else:
            BDT_out = myWF.get_BDT_output(pi_pT,gamma_eT,nBjets_25,lep_pT,piRelIso_05_ch,met,deltaphi_lep_pi,isMuon)

   
        ############################################################################
        #                                                                          #
        #------------------------------- Fill histos ------------------------------#
        #                                                                          #
        ############################################################################

        if (isMuon and BDT_out > -0.1 and BDT_out < BDT_OUT_MU) or (not isMuon and BDT_out > -0.1 and BDT_out < BDT_OUT_ELE): # Alternative cut on BDT score. Two histos will be filled. Eventually, they will be the ratio of Wmass distributions: one with pi_pT and gamma_ET normalized to Wmass, one not
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

        #if not "Signal" in name_sample and not "Data" in name_sample:
        if name_sample == "Data":
            h_PUdistrib.Fill(nPV,Event_Weight)
        

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
            h_base[theSampleName+"h_deltaphi_mu_pi"].Fill(deltaphi_lep_pi,Event_Weight)
            h_base[theSampleName+"h_deltaphi_mu_gamma"].Fill(deltaphi_lep_gamma,Event_Weight)

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
            h_base[theSampleName+"h_deltaphi_ele_pi"].Fill(deltaphi_lep_pi,Event_Weight)
            h_base[theSampleName+"h_deltaphi_ele_gamma"].Fill(deltaphi_lep_gamma,Event_Weight)
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
        #--------------------- Define histo and legend features -------------------#
        #                                                                          #
        ############################################################################


    for idx_histo,hname in enumerate(list_histos):

        # Set to 0 the bins containing negative values, due to negative weights
        hsize = h_base[theSampleName+hname].GetSize() - 2 # GetSize() returns the number of bins +2 (that is + overflow + underflow) 
        for bin in range(1,hsize+1): # The +1 is in order to get the last bin
            bincontent = h_base[theSampleName+hname].GetBinContent(bin)
            if bincontent < 0.:
                h_base[theSampleName+hname].SetBinContent(bin,0.)

        if QCDflag:
            h_base[theSampleName+hname].SetFillColor(12)
        elif name_sample == myWF.sig_samplename:
            h_base[theSampleName+hname].SetLineStyle(2)   #dashed
            h_base[theSampleName+hname].SetLineColor(2)   #red
            h_base[theSampleName+hname].SetLineWidth(4)   #kind of thick
        elif name_sample == myWF.data_samplename:
            h_base[theSampleName+hname].SetMarkerStyle(20)   #point
        else:
            h_base[theSampleName+hname].SetFillColor(colors_mask[theSampleName])


        if idx_histo == 0 and (Nevts_per_sample > 800 or name_sample == myWF.sig_samplename): #Only plot in the legends those samples which give significant contribution
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

        #Here is where histos are added in the THStack. Signal and Data are added later
        if not QCDflag and not name_sample == myWF.sig_samplename and not name_sample == myWF.data_samplename: 
            hs[hname].Add(h_base[theSampleName+hname])

    # if not QCDflag and not name_sample == myWF.sig_samplename and not name_sample == myWF.data_samplename:
    #     idx_sample += 1

print "Finished runnning over samples!"



##################### Plotting all the backgrounds in separate histos for a given distribution ######################

eleeta_canvas = dict()

# for sample_name, sampleEra in zip(samplename_list,sampleEra_list):
for name_sample_withYear in samplename_list:

    ###################################################
    #                                                 #
    #------- Select which era will be plotted --------#
    #                                                 #
    ###################################################

    # if runningEra == 0 and not sampleEra == "2016":
    #     continue
    # if runningEra == 1 and not sampleEra == "2017":
    #     continue
    # if runningEra == 2 and not sampleEra == "2018":
    #     continue
    # if runningEra == 3 and not (sampleEra == "2016" or sampleEra == "2017"):
    #     continue


    if runningEra == 0 and not "2016" in name_sample_withYear:
        continue
    if runningEra == 1 and not "2017" in name_sample_withYear:
        continue
    if runningEra == 2 and not "2018" in name_sample_withYear:
        continue
    if runningEra == 3 and not ("2016" or "2017") in name_sample_withYear:
        continue

    name_sample = name_sample_withYear.split("_")[0]
    sampleEra   = name_sample_withYear.split("_")[1]

    #######################################################################
    
    if "QCD" in name_sample or "Data" in name_sample or "Signal" in name_sample:
        continue
    else:
        eleeta_canvas[name_sample] = ROOT.TCanvas(name_sample,name_sample,200,106,600,600)
        h_base[name_sample+"h_eleeta"].SetTitle(name_sample)
        h_base[name_sample+"h_eleeta"].SetFillColor(colors_mask[name_sample])
        h_base[name_sample+"h_eleeta"].Draw("hist")
        print name_sample, " integral: ", h_base[name_sample+"h_eleeta"].Integral()

        eleeta_canvas[name_sample].SaveAs("plots/" + sampleEra + "/eleeta/h_eleeta_" + name_sample + ".pdf")

# PU_canvas = ROOT.TCanvas(name_sample,name_sample,200,106,600,600)
# print "INTEGRALE!!!!!! ", h_PUdistrib.Integral()
# h_PUdistrib.Scale(1/h_PUdistrib.Integral())
# h_PUdistrib.Draw("hist")
# PU_canvas.SaveAs("plots/h_PUdistrib_Data.pdf")
# PU_file = ROOT.TFile("PUdistrib_Data.root","RECREATE")
# h_PUdistrib.Write("PUdistrib_Data")
# PU_file.Close()

####################################################################################################################
        
    

for idx_histo,hname in enumerate(list_histos):
    hs[hname].Add(h_base["QCD_"+hname])

for hname in list_histos:

    canvas[hname].cd()

    hs[hname].Draw("histo")

    if "h_Wmass_" in hname:
        hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),100.))
    if hname == "h_Wmass":
        hs[hname].SetMaximum(max(hs[hname].GetHistogram().GetMaximum(),140.))
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
    if hname == "h_pieta_sig" or hname == "h_pipt_sig" or hname == "h_gammaeta_sig" or hname == "h_gammaet_sig":
        hs[hname].SetMaximum(0.3)
    if hname == "h_mueta_sig" or hname == "h_mupt_sig" or hname == "h_eleeta_sig" or hname == "h_elept_sig":
        hs[hname].SetMaximum(0.2)


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
        hs[hname].GetXaxis().SetTitle("#Delta#varphi_{#mu-#pi}")

    if hname == "h_deltaphi_ele_pi":
        hs[hname].GetXaxis().SetTitle("#Delta#varphi_{e-#pi}")

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
 
    if runningEra == 0:
        canvas[hname].SaveAs("plots/2016/" + hname + ".pdf")
    if runningEra == 1:
        canvas[hname].SaveAs("plots/2017/" + hname + ".pdf")
    if runningEra == 2:
        canvas[hname].SaveAs("plots/2018/" + hname + ".pdf")
    if runningEra == 3:
        canvas[hname].SaveAs("plots/2016_2017/" + hname + ".pdf")

print "Wmass_mu: ", Wmass_mu.Integral()
print "Wmass_ele: ", Wmass_ele.Integral()
print "Wmass_ratio_mu: ", Wmass_ratio_mu.Integral()
print "Wmass_ratio_ele: ", Wmass_ratio_ele.Integral()


Wmass_mu.Scale(1/Wmass_mu.Integral())
Wmass_ratio_mu.Scale(1/Wmass_ratio_mu.Integral())
Wmass_ele.Scale(1/Wmass_ele.Integral())
Wmass_ratio_ele.Scale(1/Wmass_ratio_ele.Integral())

###################################
canvas_Wmass_mu = ROOT.TCanvas()
ROOT.gStyle.SetOptStat(0)
Wmass_mu.SetMarkerStyle(21)
Wmass_mu.SetTitle(" ")
Wmass_mu.GetXaxis().SetTitle("m_{#pi#gamma} (GeV)")
Wmass_mu.Draw("Pe")
canvas_Wmass_mu.SaveAs("plots/Wmass_mu_standalone.pdf")
canvas_Wmass_ratio_mu = ROOT.TCanvas()
ROOT.gStyle.SetOptStat(0)
Wmass_ratio_mu.SetMarkerStyle(21)
Wmass_ratio_mu.SetTitle(" ")
Wmass_ratio_mu.GetXaxis().SetTitle("m_{#pi#gamma} (GeV)")
Wmass_ratio_mu.Draw("Pe")
canvas_Wmass_ratio_mu.SaveAs("plots/Wmass_ratio_mu_standalone.pdf")
###################################


Wmass_ratio_mu.Divide(Wmass_ratio_mu,Wmass_mu,1.0,1.0,"B")
Wmass_ratio_ele.Divide(Wmass_ratio_ele,Wmass_ele,1.0,1.0,"B")


canvas5 = ROOT.TCanvas()
Wmass_ratio_mu.SetMarkerStyle(21)
Wmass_ratio_mu.SetTitle(" ")
Wmass_ratio_mu.GetXaxis().SetTitle("m_{#pi#gamma} (GeV)")
Wmass_ratio_mu.Draw("Pe")
if isBDT_with_Wmass:
    canvas5.SaveAs("plots/Wmass_ratio_mu_Wmass.pdf")
else:
    canvas5.SaveAs("plots/Wmass_ratio_mu.pdf")


canvas6 = ROOT.TCanvas()
Wmass_ratio_ele.SetMarkerStyle(21)
Wmass_ratio_ele.SetTitle(" ")
Wmass_ratio_ele.GetXaxis().SetTitle("m_{#pi#gamma} (GeV)")
Wmass_ratio_ele.Draw("Pe")
if isBDT_with_Wmass:
    canvas6.SaveAs("plots/Wmass_ratio_ele_Wmass.pdf")
else:
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
ROOT.gStyle.SetOptStat(0)
W_ttbar_hist.SetMarkerStyle(21)
W_ttbar_hist.SetTitle(" ")
W_ttbar_hist.GetXaxis().SetTitle("p_{T}^{W} (GeV)")
W_ttbar_hist.Draw("Pe")
canvas8.SaveAs("plots/ttbar_signal_ratio.pdf")

#-------- Create and write on rootfile with Signal/ttbar ratio --------#

#ttbar_sig_file = ROOT.TFile("ttbar_signal_ratio.root","RECREATE")
#W_ttbar_hist.Write("ttbar_signal_ratio")
#ttbar_sig_file.Close()






print "Number of expected events for ", luminosity_norm, " in fb-1"
print "Number of signal events passed = ", Nsig_passed
print "Number of background events passed = ", Nbkg_passed
print "number of Signal events in muon channel: ", Sevts_mu_SFvariation
print "number of Signal events in electron channel: ", Sevts_ele_SFvariation
print "total number of S evts weighted -mu: ", Sevts_weighted_mu
print "total number of B evts weighted -mu: ", Bevts_weighted_mu
print "total number of S evts weighted -ele: ", Sevts_weighted_ele
print "total number of B evts weighted -ele: ", Bevts_weighted_ele

print "N_WGToLNuG_mu", N_WGToLNuG_mu
print "N_DoubleEMEnriched", N_DoubleEMEnriched
