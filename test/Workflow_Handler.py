import ROOT
import math
import os
from array import array


###################################################################################
#                                                                                 #
#----------- Acess scale factor histos for EGamma and Muon corrections -----------#
#                                                                                 #
#        https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2         #
#                                                                                 #
#   https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2016LegacyRereco    #
#                                                                                 #
###################################################################################

#------------------------------- Scale factors 2016 ------------------------------#

el_ID_scale_name  = "scale_factors/2016LegacyReReco_ElectronMVA90_Fall17V2.root"
el_ID_scale_file  = ROOT.TFile(el_ID_scale_name)
el_ID_scale_histo = ROOT.TH2F()
el_ID_scale_histo = el_ID_scale_file.Get("EGamma_SF2D")

ph_ID_scale_name  = "scale_factors/Fall17V2_2016_MVAwp90_photons.root"
ph_ID_scale_file  = ROOT.TFile(ph_ID_scale_name)
ph_ID_scale_histo = ROOT.TH2F()
ph_ID_scale_histo = ph_ID_scale_file.Get("EGamma_SF2D")

ph_pixVeto_scale_name  = "scale_factors/Photon_pixVeto_2D.root"
ph_pixVeto_scale_file  = ROOT.TFile(ph_pixVeto_scale_name)
ph_pixVeto_scale_histo = ROOT.TH2F()
ph_pixVeto_scale_histo = ph_pixVeto_scale_file.Get("Scaling_Factors_HasPix_R9 Inclusive")

mu_ID_scale_name_BCDEF  = "scale_factors/RunBCDEF_SF_ID.root"
mu_ID_scale_file_BCDEF  = ROOT.TFile(mu_ID_scale_name_BCDEF)
mu_ID_scale_histo_BCDEF = ROOT.TH2F()
mu_ID_scale_histo_BCDEF = mu_ID_scale_file_BCDEF.Get("NUM_MediumID_DEN_genTracks_eta_pt")

mu_ID_scale_name_GH  = "scale_factors/RunGH_SF_ID.root"
mu_ID_scale_file_GH  = ROOT.TFile(mu_ID_scale_name_GH)
mu_ID_scale_histo_GH = ROOT.TH2F()
mu_ID_scale_histo_GH = mu_ID_scale_file_GH.Get("NUM_MediumID_DEN_genTracks_eta_pt")

mu_Iso_scale_name_BCDEF  = "scale_factors/RunBCDEF_SF_ISO.root"
mu_Iso_scale_file_BCDEF  = ROOT.TFile(mu_Iso_scale_name_BCDEF)
mu_Iso_scale_histo_BCDEF = ROOT.TH2F()
mu_Iso_scale_histo_BCDEF = mu_Iso_scale_file_BCDEF.Get("NUM_LooseRelIso_DEN_MediumID_eta_pt")

mu_Iso_scale_name_GH  = "scale_factors/RunGH_SF_ISO.root"
mu_Iso_scale_file_GH  = ROOT.TFile(mu_Iso_scale_name_GH)
mu_Iso_scale_histo_GH = ROOT.TH2F()
mu_Iso_scale_histo_GH = mu_Iso_scale_file_GH.Get("NUM_LooseRelIso_DEN_MediumID_eta_pt")

mu_Trigger_scale_name_BCDEF       = "scale_factors/EfficienciesAndSF_RunBtoF.root"
mu_Trigger_scale_file_BCDEF       = ROOT.TFile(mu_Trigger_scale_name_BCDEF)
mu_Trigger_scale_histo_BCDEF_Mu24 = ROOT.TH2F()
mu_Trigger_scale_histo_BCDEF_Mu24 = mu_Trigger_scale_file_BCDEF.Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio")
mu_Trigger_scale_histo_BCDEF_Mu50 = ROOT.TH2F()
mu_Trigger_scale_histo_BCDEF_Mu50 = mu_Trigger_scale_file_BCDEF.Get("Mu50_OR_TkMu50_PtEtaBins/abseta_pt_ratio")

mu_Trigger_scale_name_GH       = "scale_factors/EfficienciesAndSF_Period4.root"
mu_Trigger_scale_file_GH       = ROOT.TFile(mu_Trigger_scale_name_GH)
mu_Trigger_scale_histo_GH_Mu24 = ROOT.TH2F()
mu_Trigger_scale_histo_GH_Mu24 = mu_Trigger_scale_file_GH.Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio")
mu_Trigger_scale_histo_GH_Mu50 = ROOT.TH2F()
mu_Trigger_scale_histo_GH_Mu50 = mu_Trigger_scale_file_GH.Get("Mu50_OR_TkMu50_PtEtaBins/abseta_pt_ratio")


#------------------------------- Scale factors 2017 ------------------------------#

el_ID_scale_name_2017  = "scale_factors/2017_ElectronMVA90.root"
el_ID_scale_file_2017  = ROOT.TFile(el_ID_scale_name_2017)
el_ID_scale_histo_2017 = ROOT.TH2F()
el_ID_scale_histo_2017 = el_ID_scale_file_2017.Get("EGamma_SF2D")

el_reco_scale_name_2017  = "scale_factors/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root"
el_reco_scale_file_2017  = ROOT.TFile(el_reco_scale_name_2017)
el_reco_scale_histo_2017 = ROOT.TH2F()
el_reco_scale_histo_2017 = el_reco_scale_file_2017.Get("EGamma_SF2D")

ph_ID_scale_name_2017  = "scale_factors/2017_PhotonsMVAwp90.root"
ph_ID_scale_file_2017  = ROOT.TFile(ph_ID_scale_name_2017)
ph_ID_scale_histo_2017 = ROOT.TH2F()
ph_ID_scale_histo_2017 = ph_ID_scale_file_2017.Get("EGamma_SF2D")

ph_pixVeto_scale_name_2017  = "scale_factors/PixelSeed_ScaleFactors_2017.root"
ph_pixVeto_scale_file_2017  = ROOT.TFile(ph_pixVeto_scale_name_2017)
ph_pixVeto_scale_histo_2017 = ROOT.TH1F()
ph_pixVeto_scale_histo_2017 = ph_pixVeto_scale_file_2017.Get("MVA_ID")

mu_ID_scale_name_2017  = "scale_factors/RunBCDEF_SF_ID_2017.root"
mu_ID_scale_file_2017  = ROOT.TFile(mu_ID_scale_name_2017)
mu_ID_scale_histo_2017 = ROOT.TH2F()
mu_ID_scale_histo_2017 = mu_ID_scale_file_2017.Get("NUM_MediumID_DEN_genTracks_pt_abseta")

mu_Iso_scale_name_2017  = "scale_factors/RunBCDEF_SF_ISO_2017.root"
mu_Iso_scale_file_2017  = ROOT.TFile(mu_Iso_scale_name_2017)
mu_Iso_scale_histo_2017 = ROOT.TH2F()
mu_Iso_scale_histo_2017 = mu_Iso_scale_file_2017.Get("NUM_LooseRelIso_DEN_MediumID_pt_abseta")

mu_Trigger_scale_name_2017       = "scale_factors/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root"
mu_Trigger_scale_file_2017       = ROOT.TFile(mu_Trigger_scale_name_2017)
mu_Trigger_scale_histo_2017_Mu27 = ROOT.TH2F()
mu_Trigger_scale_histo_2017_Mu27 = mu_Trigger_scale_file_2017.Get("IsoMu27_PtEtaBins/abseta_pt_ratio")
mu_Trigger_scale_histo_2017_Mu50 = ROOT.TH2F()
mu_Trigger_scale_histo_2017_Mu50 = mu_Trigger_scale_file_2017.Get("Mu50_PtEtaBins/abseta_pt_ratio")


###############################################
#                                             #
#------- Arrays and reader for the BDT -------#
#                                             #
###############################################

pi_pT_array           = array('f', [0.])
gamma_eT_array        = array('f', [0.])
nBjets_25_array       = array('f', [0.])
lep_pT_array          = array('f', [0.])
piRelIso_05_array     = array('f', [0.])
met_array             = array('f', [0.])
deltaphi_lep_pi_array = array('f', [0.])
isMuon_array          = array('i', [0])

reader = ROOT.TMVA.Reader("!Color")


class Workflow_Handler:

    def __init__(self,signalname,dataname,isBDT_with_Wmass,runningEra,subprocess="//"):

        # Where the root files are
        self.subprocess = subprocess

        self.data_samplename = dataname  
        self.sig_samplename = signalname
        
        # self.norm_filename = "rootfiles/" + self.subprocess + "Medium/Normalizations_table.txt"
        # self.dir_back_input = "rootfiles/" + self.subprocess + "Medium/backgrounds/"
        # self.sig_filename = "rootfiles/" + self.subprocess + "Medium/signals/" + "WPiGammaAnalysis_" + self.sig_samplename + ".root"

        if runningEra == 0:
            self.norm_filename = "rootfiles/" + self.subprocess + "MC/2016/Normalizations_table.txt"
            self.dir_back_input = "rootfiles/" + self.subprocess + "MC/2016/backgrounds/"
            self.sig_filename = "rootfiles/" + self.subprocess + "MC/2016/signals/" + "WPiGammaAnalysis_" + self.sig_samplename + ".root"

            #---------------- Data ----------------#
            self.dir_data_input = "rootfiles/" + self.subprocess + "data/2016/"


        if runningEra == 1:
            #----------------- MC -----------------#
            self.norm_filename = "rootfiles/" + self.subprocess + "MC/2017/Normalizations_table.txt"
            self.dir_back_input = "rootfiles/" + self.subprocess + "MC/2017/backgrounds/"
            self.sig_filename = "rootfiles/" + self.subprocess + "MC/2017/signals/" + "WPiGammaAnalysis_" + self.sig_samplename + ".root"

            #---------------- Data ----------------#
            self.dir_data_input = "rootfiles/" + self.subprocess + "data/2017/"


        
        ###################################################################################
        #                                                                                 #
        #------------------------ Add BDT variables to the reader ------------------------#
        #                                                                                 #
        ###################################################################################
        
        
        if isBDT_with_Wmass:
            reader.AddVariable("pi_pT/Wmass",pi_pT_array)
            reader.AddVariable("gamma_eT/Wmass",gamma_eT_array)
        else:
            reader.AddVariable("pi_pT",pi_pT_array)
            reader.AddVariable("gamma_eT",gamma_eT_array)

        reader.AddVariable("nBjets_25",nBjets_25_array)
        reader.AddVariable("lep_pT",lep_pT_array)
        reader.AddVariable("piRelIso_05_ch",piRelIso_05_array)
        reader.AddVariable("MET",met_array)
        reader.AddVariable("deltaphi_lep_pi",deltaphi_lep_pi_array)

        if isBDT_with_Wmass:
            reader.BookMVA("BDT_mu","MVA/default/weights/TMVAClassification_BDT.weights_mu_Wmass.xml")# First argument is arbitrary. To be chosen in order to distinguish among methods
            reader.BookMVA("BDT_ele","MVA/default/weights/TMVAClassification_BDT.weights_ele_Wmass.xml")

        else:
            reader.BookMVA("BDT_mu","MVA/default/weights/TMVAClassification_BDT.weights_mu.xml")# First argument is arbitrary. To be chosen in order to distinguish among methods
            reader.BookMVA("BDT_ele","MVA/default/weights/TMVAClassification_BDT.weights_ele.xml")


    ###############################################################################################################################################

    def get_BDT_output(self,pi_pT,gamma_eT,nBjets_25,lep_pT,piRelIso_05_ch,met,deltaphi_lep_pi,isMuon):

        pi_pT_array[0]           = pi_pT
        gamma_eT_array[0]        = gamma_eT
        nBjets_25_array[0]       = nBjets_25
        lep_pT_array[0]          = lep_pT
        piRelIso_05_array[0]     = piRelIso_05_ch
        met_array[0]             = met
        deltaphi_lep_pi_array[0] = deltaphi_lep_pi
        isMuon_array[0]          = int(isMuon)
        
        if isMuon:
            return reader.EvaluateMVA("BDT_mu")
        else:
            return reader.EvaluateMVA("BDT_ele")

    ###############################################################################################################################################

    # get the sample names (MC and data)
    def get_samples_names(self,Add_Signal=True,Add_Data=True):
        list_dirs_bkg = os.listdir(self.dir_back_input)
        samplename_list = []

        for dirname in list_dirs_bkg:
            tmp_samplename = dirname.split("WPiGammaAnalysis_")[1]
            tmp_samplename2 = tmp_samplename.replace(".root","")
            samplename_list.append(tmp_samplename2)

        if Add_Signal:
            samplename_list.append(self.sig_samplename)
        if Add_Data:
            samplename_list.append(self.data_samplename)
        return samplename_list 

    ###############################################################################################################################################

    # get the corresponding root files for background, signal and data sample names
    def get_root_files(self,Add_Signal=True):
        list_dirs = os.listdir(self.dir_back_input)
        list_dirs_data = os.listdir(self.dir_data_input)
        root_file = dict()

        for dirname in list_dirs:
            tmp_samplename = dirname.split("WPiGammaAnalysis_")[1]
            tmp_samplename2 = tmp_samplename.replace(".root","")
            root_file[tmp_samplename2] = ROOT.TFile(self.dir_back_input + dirname)

        for dirname in list_dirs_data:
            tmp_samplename3 = dirname.split("WPiGammaAnalysis_")[1]
            tmp_samplename4 = tmp_samplename3.replace(".root","")
            root_file[tmp_samplename4] = ROOT.TFile(self.dir_data_input + dirname)

        if Add_Signal:
            root_file[self.sig_samplename] = ROOT.TFile(self.sig_filename)

        return root_file

    ###############################################################################################################################################

    def get_normalizations_map(self):
        in_file = open(self.norm_filename,"r")
        norm_map = dict()

        for line in in_file:
            data_norm = line.split()
            norm_map[data_norm[0]] = float(data_norm[1])

        return norm_map

    ###############################################################################################################################################

    def post_preselection_cuts(self, lep_eta, lep_pT, isMuon, LepPiOppositeCharge, runningEra):

        cut = False

        if runningEra == 0:
            if (not isMuon and math.fabs(lep_eta) > 2.4) or (not isMuon and lep_pT < 28.) or (isMuon and lep_pT < 25.) or not LepPiOppositeCharge:
                cut = True

        if runningEra == 1:
            if (not isMuon and math.fabs(lep_eta) > 2.4) or (not isMuon and lep_pT < 28.) or (isMuon and lep_pT < 25.) or not LepPiOppositeCharge:
                cut = True

        return cut

    ###############################################################################################################################################
        
    def get_ele_scale(self, lep_pt, lep_eta, runningEra):

        if runningEra == 0: # Scale factors for 2016

            local_lep_pt = lep_pt
            if local_lep_pt > 150.: # This is because corrections are up to 150 GeV
                local_lep_pt = 150.

            local_lep_eta = lep_eta
            if local_lep_eta >= 2.5:
                local_lep_eta = 2.49
            if local_lep_eta <= -2.5:
                local_lep_eta = -2.49

            scale_factor      = 1.
            scale_factor_reco = 1.
            scale_factor_ID   = 1.
            el_reco_err       = 0.
            el_ID_err         = 0.
            tot_err           = 0.

            scale_factor_ID   = el_ID_scale_histo.GetBinContent( el_ID_scale_histo.GetXaxis().FindBin(local_lep_eta), el_ID_scale_histo.GetYaxis().FindBin(local_lep_pt) )
            el_ID_err         = el_ID_scale_histo.GetBinError( el_ID_scale_histo.GetXaxis().FindBin(local_lep_eta), el_ID_scale_histo.GetYaxis().FindBin(local_lep_pt) )
            
            scale_factor = scale_factor_ID
            tot_err      = el_ID_err


        if runningEra == 1: # Scale factors for 2017

            local_lep_pt = lep_pt
            if local_lep_pt > 499.: # This is because corrections are up to 499 GeV
                local_lep_pt = 499.

            local_lep_eta = lep_eta
            if local_lep_eta >= 2.5:
                local_lep_eta = 2.49
            if local_lep_eta < -2.5: # Corrections reach down to eta = -2.5, but only up to eta = 2.49
                local_lep_eta = -2.5

            scale_factor      = 1.
            scale_factor_reco = 1.
            scale_factor_ID   = 1.
            el_reco_err       = 0.
            el_ID_err         = 0.
            tot_err           = 0.
            
            scale_factor_ID   = el_ID_scale_histo_2017.GetBinContent( el_ID_scale_histo_2017.GetXaxis().FindBin(local_lep_eta), el_ID_scale_histo_2017.GetYaxis().FindBin(local_lep_pt) )
            el_ID_err         = el_ID_scale_histo_2017.GetBinError( el_ID_scale_histo_2017.GetXaxis().FindBin(local_lep_eta), el_ID_scale_histo_2017.GetYaxis().FindBin(local_lep_pt) )

            
            scale_factor_reco   = el_reco_scale_histo_2017.GetBinContent( el_reco_scale_histo_2017.GetXaxis().FindBin(local_lep_eta), el_reco_scale_histo_2017.GetYaxis().FindBin(local_lep_pt) )
            el_reco_err         = el_reco_scale_histo_2017.GetBinError( el_reco_scale_histo_2017.GetXaxis().FindBin(local_lep_eta), el_reco_scale_histo_2017.GetYaxis().FindBin(local_lep_pt) )
            
            scale_factor = scale_factor_ID * scale_factor_reco
            tot_err      = math.sqrt( scale_factor_ID * scale_factor_ID * el_reco_err * el_reco_err + scale_factor_reco * scale_factor_reco * el_ID_err * el_ID_err)


        return scale_factor, tot_err

    ###############################################################################################################################################

    def get_photon_scale(self, ph_pt, ph_eta, runningEra):

        if runningEra == 0: # Scale factors for 2016

            local_ph_pt = ph_pt
            if local_ph_pt > 150.: # This is because corrections are up to 150 GeV
                local_ph_pt = 150.

            local_ph_eta = ph_eta
            if local_ph_eta >= 2.5:
                local_ph_eta = 2.49
            if local_ph_eta <= -2.5:
                local_ph_eta = -2.49
                
            scale_factor         = 1.
            scale_factor_ID      = 1.
            scale_factor_pixVeto = 1.
            ph_ID_err            = 0.
            ph_pixVeto_err       = 0.
            tot_err              = 0.
            
            scale_factor_ID      = ph_ID_scale_histo.GetBinContent( ph_ID_scale_histo.GetXaxis().FindBin(local_ph_eta), ph_ID_scale_histo.GetYaxis().FindBin(local_ph_pt) )
            ph_ID_err            = ph_ID_scale_histo.GetBinError( ph_ID_scale_histo.GetXaxis().FindBin(local_ph_eta), ph_ID_scale_histo.GetYaxis().FindBin(local_ph_pt) )
            
            scale_factor_pixVeto = ph_pixVeto_scale_histo.GetBinContent( ph_pixVeto_scale_histo.GetXaxis().FindBin(math.fabs(local_ph_eta)), ph_pixVeto_scale_histo.GetYaxis().FindBin(local_ph_pt) )
            ph_pixVeto_err       = ph_pixVeto_scale_histo.GetBinError( ph_pixVeto_scale_histo.GetXaxis().FindBin(math.fabs(local_ph_eta)), ph_pixVeto_scale_histo.GetYaxis().FindBin(local_ph_pt) )


        if runningEra == 1:  # Scale factors for 2017

            local_ph_pt = ph_pt
            if local_ph_pt > 499.: # This is because corrections are up to 499 GeV
                local_ph_pt = 499.

            local_ph_eta = ph_eta
            if local_ph_eta >= 2.5:
                local_ph_eta = 2.49
            if local_ph_eta < -2.5: # Corrections reach down to eta = -2.5, but only up to eta = 2.49
                local_ph_eta = -2.5

            scale_factor         = 1.
            scale_factor_ID      = 1.
            scale_factor_pixVeto = 1.
            ph_ID_err            = 0.
            ph_pixVeto_err       = 0.
            tot_err              = 0.
            
            scale_factor_ID      = ph_ID_scale_histo_2017.GetBinContent( ph_ID_scale_histo_2017.GetXaxis().FindBin(local_ph_eta), ph_ID_scale_histo_2017.GetYaxis().FindBin(local_ph_pt) )
            ph_ID_err            = ph_ID_scale_histo_2017.GetBinError( ph_ID_scale_histo_2017.GetXaxis().FindBin(local_ph_eta), ph_ID_scale_histo_2017.GetYaxis().FindBin(local_ph_pt) )
            
            if local_ph_eta <= 1.5: # Get the inclusive R9 SF for EB
                scale_factor_pixVeto = ph_pixVeto_scale_histo_2017.GetBinContent(1)
                ph_pixVeto_err       = ph_pixVeto_scale_histo_2017.GetBinError(1)

            if local_ph_eta > 1.5: # Get the inclusive R9 SF for EE
                scale_factor_pixVeto = ph_pixVeto_scale_histo_2017.GetBinContent(4)
                ph_pixVeto_err       = ph_pixVeto_scale_histo_2017.GetBinError(4)


        scale_factor = scale_factor_ID * scale_factor_pixVeto
        tot_err      = math.sqrt( scale_factor_pixVeto * scale_factor_pixVeto * ph_ID_err * ph_ID_err + scale_factor_ID * scale_factor_ID * ph_pixVeto_err * ph_pixVeto_err )


        return scale_factor, tot_err

    ###############################################################################################################################################

    def get_muon_scale(self, lep_pT, lep_eta, runningEra, Nrandom_for_SF, luminosity_BtoF, luminosity_norm, isSingleMuTrigger_LOW):

        local_lep_pt = lep_pt
        if local_lep_pt >= 120.:  # This is because corrections go up to 120 GeV (excluded)
            local_lep_pt = 119.
        if local_lep_pt < 20.:
            local_lep_pt = 20.

        local_lep_eta = lep_eta
        if local_lep_eta >= 2.4:
            local_lep_eta = 2.39
        if local_lep_eta <= -2.4:
            local_lep_eta = -2.39


        scale_factor            = 1.
        scale_factor_ID         = 1.
        scale_factor_Iso        = 1.
        scale_factor_Tracking   = 1.
        scale_factor_Trigger    = 1.
        mu_ID_err               = 0.
        mu_Iso_err              = 0.
        mu_Trigger_err          = 0.
        tot_err                 = 0.


        ###############################################
        #                                             #
        #------------ Get muon SF for 2016 -----------#
        #                                             #
        ###############################################

        if runningEra == 0:
            
            if Nrandom_for_SF <= (luminosity_BtoF/luminosity_norm): # Access muon SF: B to F

                scale_factor_ID       = mu_ID_scale_histo_BCDEF.GetBinContent( mu_ID_scale_histo_BCDEF.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_ID_scale_histo_BCDEF.GetYaxis().FindBin(local_lep_pt) )
                mu_ID_err             = mu_ID_scale_histo_BCDEF.GetBinError( mu_ID_scale_histo_BCDEF.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_ID_scale_histo_BCDEF.GetYaxis().FindBin(local_lep_pt) )
                scale_factor_Iso      = mu_Iso_scale_histo_BCDEF.GetBinContent( mu_Iso_scale_histo_BCDEF.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Iso_scale_histo_BCDEF.GetYaxis().FindBin(local_lep_pt) )
                mu_Iso_err            = mu_Iso_scale_histo_BCDEF.GetBinError( mu_Iso_scale_histo_BCDEF.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Iso_scale_histo_BCDEF.GetYaxis().FindBin(local_lep_pt) )


                if isSingleMuTrigger_LOW: # Trigger SF go up to higher energies than 120 GeV, so no local_lep_pt is used for those

                    local_lep_pt_forTrigger = lep_pt # Corrections related to this trigger start from 26 GeV
                    if local_lep_pt_forTrigger < 26.:
                        local_lep_pt_forTrigger = 26.

                    scale_factor_Trigger = mu_Trigger_scale_histo_BCDEF_Mu24.GetBinContent( mu_Trigger_scale_histo_BCDEF_Mu24.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_BCDEF_Mu24.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                    mu_Trigger_err       = mu_Trigger_scale_histo_BCDEF_Mu24.GetBinError( mu_Trigger_scale_histo_BCDEF_Mu24.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_BCDEF_Mu24.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                    
                    scale_factor         = scale_factor_ID * scale_factor_Iso * scale_factor_Trigger
                    tot_err              = math.sqrt( scale_factor_Iso * scale_factor_Iso * scale_factor_Trigger * scale_factor_Trigger * mu_ID_err * mu_ID_err + scale_factor_ID * scale_factor_ID * scale_factor_Trigger * scale_factor_Trigger * mu_Iso_err * mu_Iso_err + scale_factor_ID * scale_factor_ID * scale_factor_Iso * scale_factor_Iso * mu_Trigger_err * mu_Trigger_err )

                    return scale_factor, tot_err



                else: # An event can have more than one trigger

                    local_lep_pt_forTrigger = lep_pt # Corrections related to this trigger start from 52 GeV
                    if local_lep_pt_forTrigger < 52.:
                        local_lep_pt_forTrigger = 52.

                    scale_factor_Trigger = mu_Trigger_scale_histo_BCDEF_Mu50.GetBinContent( mu_Trigger_scale_histo_BCDEF_Mu50.GetXaxis().FindBin(math.fabs(lep_eta)), mu_Trigger_scale_histo_BCDEF_Mu50.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                    mu_Trigger_err       = mu_Trigger_scale_histo_BCDEF_Mu50.GetBinError( mu_Trigger_scale_histo_BCDEF_Mu50.GetXaxis().FindBin(math.fabs(lep_eta)), mu_Trigger_scale_histo_BCDEF_Mu50.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                    
                    scale_factor         = scale_factor_ID * scale_factor_Iso * scale_factor_Trigger
                    tot_err              = math.sqrt( scale_factor_Iso * scale_factor_Iso * scale_factor_Trigger * scale_factor_Trigger * mu_ID_err * mu_ID_err + scale_factor_ID * scale_factor_ID * scale_factor_Trigger * scale_factor_Trigger * mu_Iso_err * mu_Iso_err + scale_factor_ID * scale_factor_ID * scale_factor_Iso * scale_factor_Iso * mu_Trigger_err * mu_Trigger_err )

                    return scale_factor, tot_err

            else:

                scale_factor_ID       = mu_ID_scale_histo_GH.GetBinContent( mu_ID_scale_histo_GH.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_ID_scale_histo_GH.GetYaxis().FindBin(local_lep_pt) )
                mu_ID_err             = mu_ID_scale_histo_GH.GetBinError( mu_ID_scale_histo_GH.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_ID_scale_histo_GH.GetYaxis().FindBin(local_lep_pt) )
                scale_factor_Iso      = mu_Iso_scale_histo_GH.GetBinContent( mu_Iso_scale_histo_GH.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Iso_scale_histo_GH.GetYaxis().FindBin(local_lep_pt) )
                mu_Iso_err            = mu_Iso_scale_histo_GH.GetBinError( mu_Iso_scale_histo_GH.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Iso_scale_histo_GH.GetYaxis().FindBin(local_lep_pt) )

                if isSingleMuTrigger_24: # Trigger SF go up to higher energies than 120 GeV, so no local_lep_pt is used for those

                    local_lep_pt_forTrigger = lep_pt 
                    if local_lep_pt_forTrigger < 26.:
                        local_lep_pt_forTrigger = 26.

                    scale_factor_Trigger = mu_Trigger_scale_histo_GH_Mu24.GetBinContent( mu_Trigger_scale_histo_GH_Mu24.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_GH_Mu24.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                    mu_Trigger_err       = mu_Trigger_scale_histo_GH_Mu24.GetBinError( mu_Trigger_scale_histo_GH_Mu24.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_GH_Mu24.GetYaxis().FindBin(local_lep_pt_forTrigger) )

                    scale_factor         = scale_factor_ID * scale_factor_Iso * scale_factor_Trigger
                    tot_err              = math.sqrt( scale_factor_Iso * scale_factor_Iso * scale_factor_Trigger * scale_factor_Trigger * mu_ID_err * mu_ID_err + scale_factor_ID * scale_factor_ID * scale_factor_Trigger * scale_factor_Trigger * mu_Iso_err * mu_Iso_err + scale_factor_ID * scale_factor_ID * scale_factor_Iso * scale_factor_Iso * mu_Trigger_err * mu_Trigger_err )

                    return scale_factor, tot_err

                else: # An event can have more than one trigger

                    local_lep_pt_forTrigger = lep_pt 
                    if local_lep_pt_forTrigger < 52.:
                        local_lep_pt_forTrigger = 52.

                    scale_factor_Trigger = mu_Trigger_scale_histo_GH_Mu50.GetBinContent( mu_Trigger_scale_histo_GH_Mu50.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_GH_Mu50.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                    mu_Trigger_err       = mu_Trigger_scale_histo_GH_Mu50.GetBinError( mu_Trigger_scale_histo_GH_Mu50.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_GH_Mu50.GetYaxis().FindBin(local_lep_pt_forTrigger) )

                    scale_factor         = scale_factor_ID * scale_factor_Iso * scale_factor_Trigger
                    tot_err              = math.sqrt( scale_factor_Iso * scale_factor_Iso * scale_factor_Trigger * scale_factor_Trigger * mu_ID_err * mu_ID_err + scale_factor_ID * scale_factor_ID * scale_factor_Trigger * scale_factor_Trigger * mu_Iso_err * mu_Iso_err + scale_factor_ID * scale_factor_ID * scale_factor_Iso * scale_factor_Iso * mu_Trigger_err * mu_Trigger_err )
                    
                    return scale_factor, tot_err


        ###############################################
        #                                             #
        #------------ Get muon SF for 2017 -----------#
        #                                             #
        ###############################################


        if runningEra == 1: # Get muon SFs for 2017. ID and ISO histos have pT and eta positions inverted wrt 2016 histos
                              
            scale_factor_ID       = mu_ID_scale_histo_2017.GetBinContent( mu_ID_scale_histo_2017.GetYaxis().FindBin(local_lep_pt), mu_ID_scale_histo_2017.GetXaxis().FindBin(math.fabs(local_lep_eta)) )
            mu_ID_err             = mu_ID_scale_histo_2017.GetBinError( mu_ID_scale_histo_2017.GetYaxis().FindBin(local_lep_pt), mu_ID_scale_histo_2017.GetXaxis().FindBin(math.fabs(local_lep_eta)) )
            scale_factor_Iso      = mu_Iso_scale_histo_2017.GetBinContent( mu_Iso_scale_histo_2017.GetYaxis().FindBin(local_lep_pt), mu_Iso_scale_histo_2017.GetXaxis().FindBin(math.fabs(local_lep_eta)) )
            mu_Iso_err            = mu_Iso_scale_histo_2017.GetBinError( mu_Iso_scale_histo_2017.GetYaxis().FindBin(local_lep_pt),mu_Iso_scale_histo_2017.GetXaxis().FindBin(math.fabs(local_lep_eta)) )


            if isSingleMuTrigger_LOW: # Trigger SF go up to higher energies than 120 GeV, so no local_lep_pt is used for those

                local_lep_pt_forTrigger = lep_pt 
                if local_lep_pt_forTrigger < 29.:
                    local_lep_pt_forTrigger = 29.

                scale_factor_Trigger = mu_Trigger_scale_histo_2017_Mu27.GetBinContent( mu_Trigger_scale_histo_2017_Mu27.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_2017_Mu27.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                mu_Trigger_err       = mu_Trigger_scale_histo_2017_Mu27.GetBinError( mu_Trigger_scale_histo_2017_Mu27.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_2017_Mu27.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                
                scale_factor         = scale_factor_ID * scale_factor_Iso * scale_factor_Trigger
                tot_err              = math.sqrt( scale_factor_Iso * scale_factor_Iso * scale_factor_Trigger * scale_factor_Trigger * mu_ID_err * mu_ID_err + scale_factor_ID * scale_factor_ID * scale_factor_Trigger * scale_factor_Trigger * mu_Iso_err * mu_Iso_err + scale_factor_ID * scale_factor_ID * scale_factor_Iso * scale_factor_Iso * mu_Trigger_err * mu_Trigger_err )
                
                return scale_factor, tot_err



            else: # An event can have more than one trigger

                local_lep_pt_forTrigger = lep_pt 
                if local_lep_pt_forTrigger < 52.:
                    local_lep_pt_forTrigger = 52.

                scale_factor_Trigger = mu_Trigger_scale_histo_2017_Mu50.GetBinContent( mu_Trigger_scale_histo_2017_Mu50.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_2017_Mu50.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                mu_Trigger_err       = mu_Trigger_scale_histo_2017_Mu50.GetBinError( mu_Trigger_scale_histo_2017_Mu50.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_2017_Mu50.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                
                scale_factor         = scale_factor_ID * scale_factor_Iso * scale_factor_Trigger
                tot_err              = math.sqrt( scale_factor_Iso * scale_factor_Iso * scale_factor_Trigger * scale_factor_Trigger * mu_ID_err * mu_ID_err + scale_factor_ID * scale_factor_ID * scale_factor_Trigger * scale_factor_Trigger * mu_Iso_err * mu_Iso_err + scale_factor_ID * scale_factor_ID * scale_factor_Iso * scale_factor_Iso * mu_Trigger_err * mu_Trigger_err )
                
                return scale_factor, tot_err



