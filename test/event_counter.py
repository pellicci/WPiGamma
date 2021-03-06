import ROOT
import os
import math
import numpy as np
import copy
from array import array
from Workflow_Handler import Workflow_Handler

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

isBDT_with_Wmass = False # If true, pT(pi) and ET(gamma) in the BDT are normalized to Wmass 
myWF = Workflow_Handler("Signal","Data",isBDT_with_Wmass)

ttbar_sig_calib_file = ROOT.TFile.Open("ttbar_signal_ratio.root")
ttbar_sig_calib = ttbar_sig_calib_file.Get("ttbar_signal_ratio")

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

##Here starts the program
Norm_Map = myWF.get_normalizations_map()

##Get the files and the names of the samples
samplename_list = myWF.get_samples_names()
root_file = myWF.get_root_files()

# Sevts_mu_SFvariation  = 0 # Counters for the number of signal events (weighted) when variating scale factors
# Sevts_ele_SFvariation = 0
# Sevts_tot = 0
# Bevts_tot = 0
# Sevts_weighted_mu = 0
# Bevts_weighted_mu = 0
# Sevts_weighted_ele = 0
# Bevts_weighted_ele = 0
_Nrandom_for_SF = ROOT.TRandom3(44317)
_Nrandom_for_Gaus_SF = ROOT.TRandom3(44329)


event_counter_mu = dict()
event_counter_el = dict()
QCD_sum_mu = 0.
QCD_sum_el = 0.
event_counter_mu_noBDT = dict()
event_counter_el_noBDT = dict()
QCD_sum_mu_noBDT = 0.
QCD_sum_el_noBDT = 0.

# Initializing all the event counters
for name_sample in samplename_list:
    event_counter_mu[name_sample] = 0.
    event_counter_el[name_sample] = 0.
    event_counter_mu_noBDT[name_sample] = 0.
    event_counter_el_noBDT[name_sample] = 0.


for name_sample in samplename_list:

    theSampleName = name_sample

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

        if "Data" in name_sample: continue  #-------------Excluding data-------------#

        #---------- Access the tree variables ----------#

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

        #---------Determine the total event weight---------#
        
        if isMuon: # Get muon scale factors, which are different for two groups of datasets, and weight them for the respective integrated lumi 
            mu_weight_BtoF, mu_weight_BtoF_err = myWF.get_muon_scale_BtoF(lep_pT,lep_eta,isSingleMuTrigger_24)
            mu_weight_GH, mu_weight_GH_err     = myWF.get_muon_scale_GH(lep_pT,lep_eta,isSingleMuTrigger_24)

            # Use a random number to select which muon scale factor to use, depending on the respective lumi fraction
            Nrandom_for_SF = _Nrandom_for_SF.Rndm()

            if Nrandom_for_SF <= (luminosity_BtoF/luminosity_norm): # Access muon SF: B to F

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
            MC_Weight = mytree.MC_Weight # Add MC weight        
            PU_Weight = mytree.PU_Weight # Add Pile Up weight        
            Event_Weight = norm_factor*ph_weight*MC_Weight*PU_Weight/math.fabs(MC_Weight) # Just take the sign of the gen weight

            if isMuon:
                Event_Weight = Event_Weight*mu_weight
            else:
                Event_Weight = Event_Weight*ele_weight

            # Correct for the difference in pT of the generated W in Pythia and Madgraph samples
            if name_sample == "Signal" and is_signal_Wplus:
                local_Wplus_pT = Wplus_pT
                if Wplus_pT > 300.:
                    local_Wplus_pT = 300.

                Event_Weight = Event_Weight*ttbar_sig_calib.GetBinContent(ttbar_sig_calib.GetXaxis().FindBin(local_Wplus_pT))


            if name_sample == "Signal" and is_signal_Wminus:
                local_Wminus_pT = Wminus_pT
                if Wminus_pT > 300.:
                    local_Wminus_pT = 300.

                Event_Weight = Event_Weight*ttbar_sig_calib.GetBinContent(ttbar_sig_calib.GetXaxis().FindBin(local_Wminus_pT))

            # Obtaining the number of sig and bkg events (weighted)
            
            # if "Signal" in name_sample and isMuon:
            #     Sevts_weighted_mu += Event_Weight
            # if not "Signal" in name_sample and isMuon:
            #     Bevts_weighted_mu += Event_Weight
            # if "Signal" in name_sample and not isMuon:
            #     Sevts_weighted_ele += Event_Weight
            # if not "Signal" in name_sample and not isMuon:
            #     Bevts_weighted_ele += Event_Weight
            

        else:
            Event_Weight = 1.
 
        if isMuon and not "Data" in name_sample:
            event_counter_mu_noBDT[theSampleName] += Event_Weight

        if not isMuon and not "Data" in name_sample:
            event_counter_el_noBDT[theSampleName] += Event_Weight


        #---------Retrieve the BDT output----------#

        if isBDT_with_Wmass:
            BDT_out = myWF.get_BDT_output(pi_pT/Wmass,gamma_eT/Wmass,nBjets_25,lep_pT,piRelIso_05_ch,met,deltaphi_lep_pi,isMuon)  
        else:
            BDT_out = myWF.get_BDT_output(pi_pT,gamma_eT,nBjets_25,lep_pT,piRelIso_05_ch,met,deltaphi_lep_pi,isMuon)       

        #---------------------Here's where the BDT selection starts---------------------#
      
        if (isMuon and BDT_out >= 0.220) or (not isMuon and BDT_out >= 0.170):
            if (Wmass >= 65. and Wmass <= 90.): # Count the number of events in the signal region only
            #if (Wmass >= 50. and Wmass <= 100.):

                if isMuon:
                    if not "Data" in name_sample:
                        event_counter_mu[theSampleName] += Event_Weight
                        
                        # if "Signal" in name_sample:
                        #     Sevts_mu_SFvariation += Event_Weight

 
                if not isMuon:# and lep_iso <= 0.35:
                    if not "Data" in name_sample:
                        event_counter_el[theSampleName] += Event_Weight

                        # if "Signal" in name_sample:
                        #     Sevts_ele_SFvariation += Event_Weight


        #-------BDT cut variation------#
        # if (isMuon and BDT_out >= 0.087) or (not isMuon and BDT_out >= 0.072): #Extra block
        #     if (Wmass >= 50. and Wmass <= 100.) and isMuon:
        #         h_base[theSampleName+"h_Wmass_mu_minus"].Fill(Wmass,Event_Weight)
        #     if (Wmass >= 50. and Wmass <= 100.) and not isMuon and lep_iso <= 0.35:
        #         h_base[theSampleName+"h_Wmass_ele_minus"].Fill(Wmass,Event_Weight)

        # if (isMuon and BDT_out >= 0.101) or (not isMuon and BDT_out >= 0.080): #Extra block
        #     if (Wmass >= 50. and Wmass <= 100.) and isMuon:
        #         h_base[theSampleName+"h_Wmass_mu_plus"].Fill(Wmass,Event_Weight)
        #     if (Wmass >= 50. and Wmass <= 100.) and not isMuon and lep_iso <= 0.35:
        #         h_base[theSampleName+"h_Wmass_ele_plus"].Fill(Wmass,Event_Weight)


print "Starting to sum QCD events"

for name_sample in samplename_list:

    if "QCD" in name_sample:
        QCD_sum_mu += event_counter_mu[name_sample]
        QCD_sum_el += event_counter_el[name_sample]
        QCD_sum_mu_noBDT += event_counter_mu_noBDT[name_sample]
        QCD_sum_el_noBDT += event_counter_el_noBDT[name_sample]

print " ##################################### "

for name_sample in samplename_list:
    
    print "N of events in " , name_sample , " in muon channel: ", event_counter_mu[name_sample], " and in electron channel: ", event_counter_el[name_sample]
    print "N of events before BDT in " , name_sample , " in muon channel: ", event_counter_mu_noBDT[name_sample], " and in electron channel: ", event_counter_el_noBDT[name_sample]
    print
    print "######################"
    print

print "N of events in QCD in muon channel: ", QCD_sum_mu, "N of events in QCD in electron channel: ", QCD_sum_el
print "N of events in QCD in muon channel before BDT: ", QCD_sum_mu_noBDT, "N of events in QCD in electron channel before BDT: ", QCD_sum_el_noBDT

