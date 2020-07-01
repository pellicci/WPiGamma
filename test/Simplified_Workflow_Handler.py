import ROOT
import math
import os
from array import array


###################################################################################
#                                                                                 #
#----------- Acess scale factor histos for EGamma and Muon corrections -----------#
#                                                                                 #
#        https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2         #
#   https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2016LegacyRereco    #
#                                                                                 #
###################################################################################

#------------------------------- Scale factors 2016 ------------------------------#

el_ID_scale_name_2016  = "scale_factors/2016LegacyReReco_ElectronMVA90_Fall17V2.root"
el_ID_scale_file_2016  = ROOT.TFile(el_ID_scale_name_2016)
el_ID_scale_histo_2016 = ROOT.TH2F()
el_ID_scale_histo_2016 = el_ID_scale_file_2016.Get("EGamma_SF2D")

el_reco_scale_name_2016  = "scale_factors/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root"
el_reco_scale_file_2016  = ROOT.TFile(el_reco_scale_name_2016)
el_reco_scale_histo_2016 = ROOT.TH2F()
el_reco_scale_histo_2016 = el_reco_scale_file_2016.Get("EGamma_SF2D")

el_trigger_scale_name_2016  = "scale_factors/egammaTriggerEfficiency_2016.root"
el_trigger_scale_file_2016  = ROOT.TFile(el_trigger_scale_name_2016)
el_trigger_scale_histo_2016 = ROOT.TH2F()
el_trigger_scale_histo_2016 = el_trigger_scale_file_2016.Get("EGamma_SF2D")

ph_ID_scale_name_2016  = "scale_factors/Fall17V2_2016_MVAwp90_photons.root"
ph_ID_scale_file_2016  = ROOT.TFile(ph_ID_scale_name_2016)
ph_ID_scale_histo_2016 = ROOT.TH2F()
ph_ID_scale_histo_2016 = ph_ID_scale_file_2016.Get("EGamma_SF2D")

ph_pixVeto_scale_name_2016  = "scale_factors/Photon_pixVeto_2D_2016.root"
ph_pixVeto_scale_file_2016  = ROOT.TFile(ph_pixVeto_scale_name_2016)
ph_pixVeto_scale_histo_2016 = ROOT.TH2F()
ph_pixVeto_scale_histo_2016 = ph_pixVeto_scale_file_2016.Get("Scaling_Factors_HasPix_R9 Inclusive")

mu_ID_scale_name_BCDEF_2016  = "scale_factors/RunBCDEF_SF_ID_muon_2016.root"
mu_ID_scale_file_BCDEF_2016  = ROOT.TFile(mu_ID_scale_name_BCDEF_2016)
mu_ID_scale_histo_BCDEF_2016 = ROOT.TH2F()
mu_ID_scale_histo_BCDEF_2016 = mu_ID_scale_file_BCDEF_2016.Get("NUM_MediumID_DEN_genTracks_eta_pt")

mu_ID_scale_name_GH_2016  = "scale_factors/RunGH_SF_ID_muon_2016.root"
mu_ID_scale_file_GH_2016  = ROOT.TFile(mu_ID_scale_name_GH_2016)
mu_ID_scale_histo_GH_2016 = ROOT.TH2F()
mu_ID_scale_histo_GH_2016 = mu_ID_scale_file_GH_2016.Get("NUM_MediumID_DEN_genTracks_eta_pt")

mu_Iso_scale_name_BCDEF_2016  = "scale_factors/RunBCDEF_SF_ISO_muon_2016.root"
mu_Iso_scale_file_BCDEF_2016  = ROOT.TFile(mu_Iso_scale_name_BCDEF_2016)
mu_Iso_scale_histo_BCDEF_2016 = ROOT.TH2F()
mu_Iso_scale_histo_BCDEF_2016 = mu_Iso_scale_file_BCDEF_2016.Get("NUM_LooseRelIso_DEN_MediumID_eta_pt")

mu_Iso_scale_name_GH_2016  = "scale_factors/RunGH_SF_ISO_muon_2016.root"
mu_Iso_scale_file_GH_2016  = ROOT.TFile(mu_Iso_scale_name_GH_2016)
mu_Iso_scale_histo_GH_2016 = ROOT.TH2F()
mu_Iso_scale_histo_GH_2016 = mu_Iso_scale_file_GH_2016.Get("NUM_LooseRelIso_DEN_MediumID_eta_pt")

mu_Trigger_scale_name_BCDEF_2016       = "scale_factors/EfficienciesAndSF_RunBtoF_muon_2016.root"
mu_Trigger_scale_file_BCDEF_2016       = ROOT.TFile(mu_Trigger_scale_name_BCDEF_2016)
mu_Trigger_scale_histo_BCDEF_Mu24_2016 = ROOT.TH2F()
mu_Trigger_scale_histo_BCDEF_Mu24_2016 = mu_Trigger_scale_file_BCDEF_2016.Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio")
mu_Trigger_scale_histo_BCDEF_Mu50_2016 = ROOT.TH2F()
mu_Trigger_scale_histo_BCDEF_Mu50_2016 = mu_Trigger_scale_file_BCDEF_2016.Get("Mu50_OR_TkMu50_PtEtaBins/abseta_pt_ratio")

mu_Trigger_scale_name_GH_2016       = "scale_factors/EfficienciesAndSF_Period4_muonTrigger_2016.root"
mu_Trigger_scale_file_GH_2016       = ROOT.TFile(mu_Trigger_scale_name_GH_2016)
mu_Trigger_scale_histo_GH_Mu24_2016 = ROOT.TH2F()
mu_Trigger_scale_histo_GH_Mu24_2016 = mu_Trigger_scale_file_GH_2016.Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio")
mu_Trigger_scale_histo_GH_Mu50_2016 = ROOT.TH2F()
mu_Trigger_scale_histo_GH_Mu50_2016 = mu_Trigger_scale_file_GH_2016.Get("Mu50_OR_TkMu50_PtEtaBins/abseta_pt_ratio")


#------------------------------- Scale factors 2017 ------------------------------#

el_ID_scale_name_2017  = "scale_factors/2017_ElectronMVA90.root"
el_ID_scale_file_2017  = ROOT.TFile(el_ID_scale_name_2017)
el_ID_scale_histo_2017 = ROOT.TH2F()
el_ID_scale_histo_2017 = el_ID_scale_file_2017.Get("EGamma_SF2D")

el_reco_scale_name_2017  = "scale_factors/egammaEffi.txt_EGM2D_runBCDEF_passingRECO_electron_2017.root"
el_reco_scale_file_2017  = ROOT.TFile(el_reco_scale_name_2017)
el_reco_scale_histo_2017 = ROOT.TH2F()
el_reco_scale_histo_2017 = el_reco_scale_file_2017.Get("EGamma_SF2D")

el_trigger_scale_name_2017  = "scale_factors/egammaTriggerEfficiency_2017.root"
el_trigger_scale_file_2017  = ROOT.TFile(el_trigger_scale_name_2017)
el_trigger_scale_histo_2017 = ROOT.TH2F()
el_trigger_scale_histo_2017 = el_trigger_scale_file_2017.Get("EGamma_SF2D")

ph_ID_scale_name_2017  = "scale_factors/2017_PhotonsMVAwp90.root"
ph_ID_scale_file_2017  = ROOT.TFile(ph_ID_scale_name_2017)
ph_ID_scale_histo_2017 = ROOT.TH2F()
ph_ID_scale_histo_2017 = ph_ID_scale_file_2017.Get("EGamma_SF2D")

ph_pixVeto_scale_name_2017  = "scale_factors/PixelSeed_ScaleFactors_2017.root"
ph_pixVeto_scale_file_2017  = ROOT.TFile(ph_pixVeto_scale_name_2017)
ph_pixVeto_scale_histo_2017 = ROOT.TH1F()
ph_pixVeto_scale_histo_2017 = ph_pixVeto_scale_file_2017.Get("MVA_ID")

mu_ID_scale_name_2017  = "scale_factors/RunBCDEF_SF_ID_muon_2017.root"
mu_ID_scale_file_2017  = ROOT.TFile(mu_ID_scale_name_2017)
mu_ID_scale_histo_2017 = ROOT.TH2F()
mu_ID_scale_histo_2017 = mu_ID_scale_file_2017.Get("NUM_MediumID_DEN_genTracks_pt_abseta")

mu_Iso_scale_name_2017  = "scale_factors/RunBCDEF_SF_ISO_muon_2017.root"
mu_Iso_scale_file_2017  = ROOT.TFile(mu_Iso_scale_name_2017)
mu_Iso_scale_histo_2017 = ROOT.TH2F()
mu_Iso_scale_histo_2017 = mu_Iso_scale_file_2017.Get("NUM_LooseRelIso_DEN_MediumID_pt_abseta")

mu_Trigger_scale_name_2017       = "scale_factors/EfficienciesAndSF_RunBtoF_Nov17Nov2017_muonTrigger.root"
mu_Trigger_scale_file_2017       = ROOT.TFile(mu_Trigger_scale_name_2017)
mu_Trigger_scale_histo_2017_Mu27 = ROOT.TH2F()
mu_Trigger_scale_histo_2017_Mu27 = mu_Trigger_scale_file_2017.Get("IsoMu27_PtEtaBins/abseta_pt_ratio")
mu_Trigger_scale_histo_2017_Mu50 = ROOT.TH2F()
mu_Trigger_scale_histo_2017_Mu50 = mu_Trigger_scale_file_2017.Get("Mu50_PtEtaBins/abseta_pt_ratio")


#------------------------------- Scale factors 2018 ------------------------------#

el_ID_scale_name_2018  = "scale_factors/2018_ElectronMVA90.root"
el_ID_scale_file_2018  = ROOT.TFile(el_ID_scale_name_2018)
el_ID_scale_histo_2018 = ROOT.TH2F()
el_ID_scale_histo_2018 = el_ID_scale_file_2018.Get("EGamma_SF2D")

el_reco_scale_name_2018  = "scale_factors/egammaEffi.txt_EGM2D_updatedAll_2018.root"
el_reco_scale_file_2018  = ROOT.TFile(el_reco_scale_name_2018)
el_reco_scale_histo_2018 = ROOT.TH2F()
el_reco_scale_histo_2018 = el_reco_scale_file_2018.Get("EGamma_SF2D")

el_trigger_scale_name_2018  = "scale_factors/egammaTriggerEfficiency_2018.root"
el_trigger_scale_file_2018  = ROOT.TFile(el_trigger_scale_name_2018)
el_trigger_scale_histo_2018 = ROOT.TH2F()
el_trigger_scale_histo_2018 = el_trigger_scale_file_2018.Get("EGamma_SF2D")

ph_ID_scale_name_2018  = "scale_factors/2018_PhotonsMVAwp90.root"
ph_ID_scale_file_2018  = ROOT.TFile(ph_ID_scale_name_2018)
ph_ID_scale_histo_2018 = ROOT.TH2F()
ph_ID_scale_histo_2018 = ph_ID_scale_file_2018.Get("EGamma_SF2D")

ph_pixVeto_scale_name_2018  = "scale_factors/HasPix_2018.root"
ph_pixVeto_scale_file_2018  = ROOT.TFile(ph_pixVeto_scale_name_2018)
ph_pixVeto_scale_histo_2018 = ROOT.TH1F()
ph_pixVeto_scale_histo_2018 = ph_pixVeto_scale_file_2018.Get("eleVeto_SF")

mu_ID_scale_name_2018  = "scale_factors/RunABCD_SF_ID_muon_2018.root"
mu_ID_scale_file_2018  = ROOT.TFile(mu_ID_scale_name_2018)
mu_ID_scale_histo_2018 = ROOT.TH2F()
mu_ID_scale_histo_2018 = mu_ID_scale_file_2018.Get("NUM_MediumID_DEN_TrackerMuons_pt_abseta")

mu_Iso_scale_name_2018  = "scale_factors/RunABCD_SF_ISO_muon_2018.root"
mu_Iso_scale_file_2018  = ROOT.TFile(mu_Iso_scale_name_2018)
mu_Iso_scale_histo_2018 = ROOT.TH2F()
mu_Iso_scale_histo_2018 = mu_Iso_scale_file_2018.Get("NUM_LooseRelIso_DEN_MediumID_pt_abseta")

mu_Trigger_scale_name_BeforeHLTUpdate_2018       = "scale_factors/EfficienciesAndSF_2018Data_BeforeMuonHLTUpdate.root"
mu_Trigger_scale_file_BeforeHLTUpdate_2018       = ROOT.TFile(mu_Trigger_scale_name_BeforeHLTUpdate_2018)
mu_Trigger_scale_histo_BeforeHLTUpdate_2018_Mu24 = ROOT.TH2F()
mu_Trigger_scale_histo_BeforeHLTUpdate_2018_Mu24 = mu_Trigger_scale_file_BeforeHLTUpdate_2018.Get("IsoMu24_PtEtaBins/abseta_pt_ratio")
mu_Trigger_scale_histo_BeforeHLTUpdate_2018_Mu50 = ROOT.TH2F()
mu_Trigger_scale_histo_BeforeHLTUpdate_2018_Mu50 = mu_Trigger_scale_file_BeforeHLTUpdate_2018.Get("Mu50_OR_OldMu100_OR_TkMu100_PtEtaBins/abseta_pt_ratio")

mu_Trigger_scale_name_AfterHLTUpdate_2018       = "scale_factors/EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root"
mu_Trigger_scale_file_AfterHLTUpdate_2018       = ROOT.TFile(mu_Trigger_scale_name_AfterHLTUpdate_2018)
mu_Trigger_scale_histo_AfterHLTUpdate_2018_Mu24 = ROOT.TH2F()
mu_Trigger_scale_histo_AfterHLTUpdate_2018_Mu24 = mu_Trigger_scale_file_AfterHLTUpdate_2018.Get("IsoMu24_PtEtaBins/abseta_pt_ratio")
mu_Trigger_scale_histo_AfterHLTUpdate_2018_Mu50 = ROOT.TH2F()
mu_Trigger_scale_histo_AfterHLTUpdate_2018_Mu50 = mu_Trigger_scale_file_AfterHLTUpdate_2018.Get("Mu50_OR_OldMu100_OR_TkMu100_PtEtaBins/abseta_pt_ratio")


###############################################
#                                             #
#------- Arrays and reader for the BDT -------#
#                                             #
###############################################

pi_pT_array              = array('f', [0.])
gamma_eT_array           = array('f', [0.])
nBjets_25_array          = array('f', [0.])
lep_pT_array             = array('f', [0.])
piRelIso_05_array        = array('f', [0.])
met_array                = array('f', [0.])
#deltaphi_lep_pi_array = array('f', [0.])
#deltaphi_lep_gamma_array = array('f', [0.])
isMuon_array             = array('i', [0])
#mT_lep_met_array         = array('f', [0])

reader = ROOT.TMVA.Reader("!Color")

class Simplified_Workflow_Handler:

    def __init__(self,signalname,dataname,isBDT,isBDT_with_Wmass,runningEra):

        if runningEra == 0:
            year = "2016"
        if runningEra == 1:
            year = "2017"
        if runningEra == 2:
            year = "2018"

        # Where the files are
        self.dir_bkg_input  = "rootfiles/latest_production/MC/backgrounds/"
        self.dir_sig_input  = "rootfiles/latest_production/MC/signals/"
        self.dir_data_input = "rootfiles/latest_production/dataprocess/"

        self.norm_filename_2016 = "rootfiles/latest_production/MC/normalizations/Normalizations_table_2016.txt"
        self.norm_filename_2017 = "rootfiles/latest_production/MC/normalizations/Normalizations_table_2017.txt"
        self.norm_filename_2018 = "rootfiles/latest_production/MC/normalizations/Normalizations_table_2018.txt"

        self.ttbar_sig_calib_file = ROOT.TFile.Open("ttbar_signal_ratio_" + year + ".root")
        self.Pythia_pT_modeling_file = ROOT.TFile.Open("Pythia_pT_modeling_" + year + ".root")
        self.Pythia_polarization_modeling_file = ROOT.TFile.Open("Pythia_polarization_modeling_Pi_W_" + year + ".root")
        
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
        #reader.AddVariable("mT_lep_met",mT_lep_met_array)
        #reader.AddVariable("deltaphi_lep_pi",deltaphi_lep_pi_array)
        #reader.AddVariable("deltaphi_lep_gamma",deltaphi_lep_gamma_array)

        if isBDT_with_Wmass:
            reader.BookMVA("BDT_mu","MVA/default/weights/TMVAClassification_BDT.weights_mu_Wmass.xml")# First argument is arbitrary. To be chosen in order to distinguish among methods
            reader.BookMVA("BDT_ele","MVA/default/weights/TMVAClassification_BDT.weights_ele_Wmass.xml")

        if isBDT and not isBDT_with_Wmass:
            reader.BookMVA("BDT_mu","MVA/default/weights/TMVAClassification_BDT.weights_mu.xml")
            reader.BookMVA("BDT_ele","MVA/default/weights/TMVAClassification_BDT.weights_ele.xml")


    ###############################################################################################################################################

    def get_BDT_output(self,pi_pT,gamma_eT,nBjets_25,lep_pT,piRelIso_05_ch,met,isMuon):

        pi_pT_array[0]           = pi_pT
        gamma_eT_array[0]        = gamma_eT
        nBjets_25_array[0]       = nBjets_25
        lep_pT_array[0]          = lep_pT
        piRelIso_05_array[0]     = piRelIso_05_ch
        met_array[0]             = met
        #mT_lep_met_array[0]      = mT_lep_met
        #deltaphi_lep_pi_array[0] = deltaphi_lep_pi
        #deltaphi_lep_gamma_array[0] = deltaphi_lep_gamma
        isMuon_array[0]          = int(isMuon)
        
        if isMuon:
            return reader.EvaluateMVA("BDT_mu")
        else:
            return reader.EvaluateMVA("BDT_ele")

    ###############################################################################################################################################

    def get_normalizations_map(self, runningEra):
        #list_dirs_norm = os.listdir(self.dir_norm_input)
        if runningEra == 0:
            in_file = open(self.norm_filename_2016,"r")
        if runningEra == 1:
            in_file = open(self.norm_filename_2017,"r")
        if runningEra == 2:
            in_file = open(self.norm_filename_2018,"r")
        norm_map = dict()

        for line in in_file:
            data_norm = line.split()
            norm_map[data_norm[0]] = float(data_norm[1])

        return norm_map

    ###############################################################################################################################################

    def get_ttbar_signal_reweight(self):

        ttbar_sig_calib = self.ttbar_sig_calib_file.Get("ttbar_signal_ratio")

        return ttbar_sig_calib

    ###############################################################################################################################################

    def get_Pythia_pT_modeling_reweight(self,up,W_pT):

        h_Pythia_sig_rescale_up   = self.Pythia_pT_modeling_file.Get("h_genW_pT_up")
        h_Pythia_sig_rescale_down = self.Pythia_pT_modeling_file.Get("h_genW_pT_down")
        Pythia_sig_rescale_up     = h_Pythia_sig_rescale_up.GetBinContent(h_Pythia_sig_rescale_up.FindBin(W_pT))
        Pythia_sig_rescale_down   = h_Pythia_sig_rescale_down.GetBinContent(h_Pythia_sig_rescale_down.FindBin(W_pT))

        if up:
            return Pythia_sig_rescale_up
        else:
            return Pythia_sig_rescale_down
 
    ###############################################################################################################################################

    def get_Pythia_polarization_modeling_reweight(self,sin2,angle_Pi_W):

        h_sin2 = self.Pythia_polarization_modeling_file.Get("h_sin2")
        h_cos  = self.Pythia_polarization_modeling_file.Get("h_cos")
        Pythia_polarization_sin2 = h_sin2.GetBinContent(h_sin2.FindBin(math.cos(angle_Pi_W)))
        Pythia_polarization_cos  = h_cos.GetBinContent(h_cos.FindBin(math.cos(angle_Pi_W)))

        if sin2:
            return Pythia_polarization_sin2
        else:
            return Pythia_polarization_cos
 
    ###############################################################################################################################################

    # get the sample names (MC and data)
    def get_samples_names(self,Add_Signal=True,Add_Data=True):
        list_dirs_bkg  = os.listdir(self.dir_bkg_input)
        list_dirs_sig  = os.listdir(self.dir_sig_input)
        list_dirs_data = os.listdir(self.dir_data_input)
        samplename_list = []

        for dirname in list_dirs_bkg:
            tmp_samplename = dirname.split("WPiGammaAnalysis_")[1].replace(".root","")
            samplename_list.append(tmp_samplename)

        if Add_Signal:
            for dirname in list_dirs_sig:
                tmp_samplename = dirname.split("WPiGammaAnalysis_")[1].replace(".root","")
                samplename_list.append(tmp_samplename)

        if Add_Data:
            for dirname in list_dirs_data:
                tmp_samplename = dirname.split("WPiGammaAnalysis_")[1].replace(".root","")
                samplename_list.append(tmp_samplename)

        return samplename_list

    ###############################################################################################################################################

    # get the corresponding root files for background, signal and data sample names
    def get_root_files(self,Add_Signal=True,Add_Data=True):
        list_dirs_bkg  = os.listdir(self.dir_bkg_input)
        list_dirs_sig  = os.listdir(self.dir_sig_input)
        list_dirs_data = os.listdir(self.dir_data_input)
        root_file = dict()

        for dirname in list_dirs_bkg:
            tmp_samplename = dirname.split("_")[1]
            tmp_sampleEra  = dirname.split("_")[2].replace(".root","") 
            root_file[tmp_samplename+"_"+tmp_sampleEra] = ROOT.TFile(self.dir_bkg_input + dirname)

        if Add_Signal:
            for dirname in list_dirs_sig:
                tmp_samplename = dirname.split("_")[1]
                tmp_sampleEra  = dirname.split("_")[2].replace(".root","")
                root_file[tmp_samplename+"_"+tmp_sampleEra] = ROOT.TFile(self.dir_sig_input + dirname)

        if Add_Data:
            for dirname in list_dirs_data:
                tmp_samplename = dirname.split("_")[1]
                tmp_sampleEra  = dirname.split("_")[2].replace(".root","")
                root_file[tmp_samplename+"_"+tmp_sampleEra] = ROOT.TFile(self.dir_data_input + dirname)

        return root_file

    ###############################################################################################################################################

    def post_preselection_cuts(self, lep_eta, lep_pT, ele_etaSC, gamma_etaSC, isMuon, deltaphi_lep_pi, deltaphi_lep_gamma, isTriggerMatched, runningEra):

        if runningEra == 0:
            if (not isMuon and math.fabs(lep_eta) > 2.4) or (not isMuon and lep_pT < 30.) or (isMuon and deltaphi_lep_pi <= 0.09) or (not isMuon and deltaphi_lep_pi <= 0.05) or (deltaphi_lep_gamma < 0.04) or (isMuon and lep_pT < 25.) or (not isMuon and (math.fabs(ele_etaSC) >= 1.4442 and math.fabs(ele_etaSC) <= 1.566)) or (math.fabs(gamma_etaSC) >= 1.4442 and math.fabs(gamma_etaSC) <= 1.566) or not isTriggerMatched:
                return True

        if runningEra == 1:
            if (not isMuon and math.fabs(lep_eta) > 2.4) or (not isMuon and lep_pT < 33.) or (isMuon and deltaphi_lep_pi <= 0.09) or (not isMuon and deltaphi_lep_pi <= 0.05) or (deltaphi_lep_gamma < 0.04) or (isMuon and lep_pT < 28.) or (not isMuon and (math.fabs(ele_etaSC) >= 1.4442 and math.fabs(ele_etaSC) <= 1.566)) or (math.fabs(gamma_etaSC) >= 1.4442 and math.fabs(gamma_etaSC) <= 1.566) or not isTriggerMatched:
                return True

        if runningEra == 2:
            if (not isMuon and math.fabs(lep_eta) > 2.4) or (not isMuon and lep_pT < 33.) or (isMuon and deltaphi_lep_pi <= 0.09) or (not isMuon and deltaphi_lep_pi <= 0.05) or (deltaphi_lep_gamma < 0.04) or (isMuon and lep_pT < 25.) or (not isMuon and (math.fabs(ele_etaSC) >= 1.4442 and math.fabs(ele_etaSC) <= 1.566)) or (math.fabs(gamma_etaSC) >= 1.4442 and math.fabs(gamma_etaSC) <= 1.566) or not isTriggerMatched:
                return True

        return False

    ###############################################################################################################################################
        
    def get_ele_scale(self, lep_pt, lep_eta, runningEra):

        if runningEra == 0: # Scale factors for 2016

            local_lep_pt = lep_pt
            if local_lep_pt > 150.: # This is because corrections are up to 150 GeV
                local_lep_pt = 150.

            local_lep_pt_4reco = lep_pt
            if local_lep_pt_4reco > 499.:
                local_lep_pt_4reco = 499.

            local_lep_eta = lep_eta
            if local_lep_eta >= 2.5:
                local_lep_eta = 2.49
            if local_lep_eta <= -2.5:
                local_lep_eta = -2.49


            scale_factor_ID   = el_ID_scale_histo_2016.GetBinContent( el_ID_scale_histo_2016.GetXaxis().FindBin(local_lep_eta), el_ID_scale_histo_2016.GetYaxis().FindBin(local_lep_pt) )
            el_ID_err         = el_ID_scale_histo_2016.GetBinError( el_ID_scale_histo_2016.GetXaxis().FindBin(local_lep_eta), el_ID_scale_histo_2016.GetYaxis().FindBin(local_lep_pt) )

            scale_factor_reco   = el_reco_scale_histo_2016.GetBinContent( el_reco_scale_histo_2016.GetXaxis().FindBin(local_lep_eta), el_reco_scale_histo_2016.GetYaxis().FindBin(local_lep_pt_4reco) )
            el_reco_err         = el_reco_scale_histo_2016.GetBinError( el_reco_scale_histo_2016.GetXaxis().FindBin(local_lep_eta), el_reco_scale_histo_2016.GetYaxis().FindBin(local_lep_pt_4reco) )

            scale_factor_trigger   = el_trigger_scale_histo_2016.GetBinContent( el_trigger_scale_histo_2016.GetXaxis().FindBin(local_lep_eta), el_trigger_scale_histo_2016.GetYaxis().FindBin(local_lep_pt_4reco) ) #trigger SFs have the same pT maximum of reco SFs
            el_trigger_err         = el_trigger_scale_histo_2016.GetBinError( el_trigger_scale_histo_2016.GetXaxis().FindBin(local_lep_eta), el_trigger_scale_histo_2016.GetYaxis().FindBin(local_lep_pt_4reco) ) #trigger SFs have the same pT maximum of reco SFs
            
            scale_factor = scale_factor_ID * scale_factor_reco * scale_factor_trigger
            tot_err      = math.sqrt( scale_factor_ID * scale_factor_ID * scale_factor_trigger * scale_factor_trigger * el_reco_err * el_reco_err + scale_factor_reco * scale_factor_reco * scale_factor_trigger * scale_factor_trigger * el_ID_err * el_ID_err + scale_factor_ID * scale_factor_ID * scale_factor_reco * scale_factor_reco * el_trigger_err * el_trigger_err )


        if runningEra == 1: # Scale factors for 2017

            local_lep_pt = lep_pt
            if local_lep_pt > 499.: # This is because corrections are up to 499 GeV
                local_lep_pt = 499.

            local_lep_eta = lep_eta
            if local_lep_eta >= 2.5:
                local_lep_eta = 2.49
            if local_lep_eta < -2.5: # Corrections reach down to eta = -2.5, but only up to eta = 2.49
                local_lep_eta = -2.5

            
            scale_factor_ID   = el_ID_scale_histo_2017.GetBinContent( el_ID_scale_histo_2017.GetXaxis().FindBin(local_lep_eta), el_ID_scale_histo_2017.GetYaxis().FindBin(local_lep_pt) )
            el_ID_err         = el_ID_scale_histo_2017.GetBinError( el_ID_scale_histo_2017.GetXaxis().FindBin(local_lep_eta), el_ID_scale_histo_2017.GetYaxis().FindBin(local_lep_pt) )
            
            scale_factor_reco   = el_reco_scale_histo_2017.GetBinContent( el_reco_scale_histo_2017.GetXaxis().FindBin(local_lep_eta), el_reco_scale_histo_2017.GetYaxis().FindBin(local_lep_pt) )
            el_reco_err         = el_reco_scale_histo_2017.GetBinError( el_reco_scale_histo_2017.GetXaxis().FindBin(local_lep_eta), el_reco_scale_histo_2017.GetYaxis().FindBin(local_lep_pt) )

            scale_factor_trigger   = el_trigger_scale_histo_2017.GetBinContent( el_trigger_scale_histo_2017.GetXaxis().FindBin(local_lep_eta), el_trigger_scale_histo_2017.GetYaxis().FindBin(local_lep_pt) ) #trigger SFs have the same pT maximum of reco SFs
            el_trigger_err         = el_trigger_scale_histo_2017.GetBinError( el_trigger_scale_histo_2017.GetXaxis().FindBin(local_lep_eta), el_trigger_scale_histo_2017.GetYaxis().FindBin(local_lep_pt) ) #trigger SFs have the same pT maximum of reco SFs
            
            scale_factor = scale_factor_ID * scale_factor_reco * scale_factor_trigger
            tot_err      = math.sqrt( scale_factor_ID * scale_factor_ID * scale_factor_trigger * scale_factor_trigger * el_reco_err * el_reco_err + scale_factor_reco * scale_factor_reco * scale_factor_trigger * scale_factor_trigger * el_ID_err * el_ID_err + scale_factor_ID * scale_factor_ID * scale_factor_reco * scale_factor_reco * el_trigger_err * el_trigger_err )



        if runningEra == 2: # Scale factors for 2018

            local_lep_pt = lep_pt
            if local_lep_pt > 499.: # This is because corrections are up to 499 GeV
                local_lep_pt = 499.

            local_lep_eta = lep_eta
            if local_lep_eta >= 2.5:
                local_lep_eta = 2.49
            if local_lep_eta < -2.5: # Corrections reach down to eta = -2.5, but only up to eta = 2.49
                local_lep_eta = -2.5

            
            scale_factor_ID   = el_ID_scale_histo_2018.GetBinContent( el_ID_scale_histo_2018.GetXaxis().FindBin(local_lep_eta), el_ID_scale_histo_2018.GetYaxis().FindBin(local_lep_pt) )
            el_ID_err         = el_ID_scale_histo_2018.GetBinError( el_ID_scale_histo_2018.GetXaxis().FindBin(local_lep_eta), el_ID_scale_histo_2018.GetYaxis().FindBin(local_lep_pt) )

            scale_factor_reco   = el_reco_scale_histo_2018.GetBinContent( el_reco_scale_histo_2018.GetXaxis().FindBin(local_lep_eta), el_reco_scale_histo_2018.GetYaxis().FindBin(local_lep_pt) )
            el_reco_err         = el_reco_scale_histo_2018.GetBinError( el_reco_scale_histo_2018.GetXaxis().FindBin(local_lep_eta), el_reco_scale_histo_2018.GetYaxis().FindBin(local_lep_pt) )

            scale_factor_trigger   = el_trigger_scale_histo_2018.GetBinContent( el_trigger_scale_histo_2018.GetXaxis().FindBin(local_lep_eta), el_trigger_scale_histo_2018.GetYaxis().FindBin(local_lep_pt) ) #trigger SFs have the same pT maximum of reco SFs
            el_trigger_err         = el_trigger_scale_histo_2018.GetBinError( el_trigger_scale_histo_2018.GetXaxis().FindBin(local_lep_eta), el_trigger_scale_histo_2018.GetYaxis().FindBin(local_lep_pt) ) #trigger SFs have the same pT maximum of reco SFs
            
            scale_factor = scale_factor_ID * scale_factor_reco * scale_factor_trigger
            tot_err      = math.sqrt( scale_factor_ID * scale_factor_ID * scale_factor_trigger * scale_factor_trigger * el_reco_err * el_reco_err + scale_factor_reco * scale_factor_reco * scale_factor_trigger * scale_factor_trigger * el_ID_err * el_ID_err + scale_factor_ID * scale_factor_ID * scale_factor_reco * scale_factor_reco * el_trigger_err * el_trigger_err )


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
                
            
            scale_factor_ID      = ph_ID_scale_histo_2016.GetBinContent( ph_ID_scale_histo_2016.GetXaxis().FindBin(local_ph_eta), ph_ID_scale_histo_2016.GetYaxis().FindBin(local_ph_pt) )
            ph_ID_err            = ph_ID_scale_histo_2016.GetBinError( ph_ID_scale_histo_2016.GetXaxis().FindBin(local_ph_eta), ph_ID_scale_histo_2016.GetYaxis().FindBin(local_ph_pt) )
            
            scale_factor_pixVeto = ph_pixVeto_scale_histo_2016.GetBinContent( ph_pixVeto_scale_histo_2016.GetXaxis().FindBin(math.fabs(local_ph_eta)), ph_pixVeto_scale_histo_2016.GetYaxis().FindBin(local_ph_pt) )
            ph_pixVeto_err       = ph_pixVeto_scale_histo_2016.GetBinError( ph_pixVeto_scale_histo_2016.GetXaxis().FindBin(math.fabs(local_ph_eta)), ph_pixVeto_scale_histo_2016.GetYaxis().FindBin(local_ph_pt) )



        if runningEra == 1:  # Scale factors for 2017

            local_ph_pt = ph_pt
            if local_ph_pt > 499.: # This is because corrections are up to 499 GeV
                local_ph_pt = 499.

            local_ph_eta = ph_eta
            if local_ph_eta >= 2.5:
                local_ph_eta = 2.49
            if local_ph_eta < -2.5: # Corrections reach down to eta = -2.5, but only up to eta = 2.49
                local_ph_eta = -2.5

            
            scale_factor_ID      = ph_ID_scale_histo_2017.GetBinContent( ph_ID_scale_histo_2017.GetXaxis().FindBin(local_ph_eta), ph_ID_scale_histo_2017.GetYaxis().FindBin(local_ph_pt) )
            ph_ID_err            = ph_ID_scale_histo_2017.GetBinError( ph_ID_scale_histo_2017.GetXaxis().FindBin(local_ph_eta), ph_ID_scale_histo_2017.GetYaxis().FindBin(local_ph_pt) )
            
            if local_ph_eta <= 1.5: # Get the inclusive R9 SF for EB
                scale_factor_pixVeto = ph_pixVeto_scale_histo_2017.GetBinContent(1)
                ph_pixVeto_err       = ph_pixVeto_scale_histo_2017.GetBinError(1)

            if local_ph_eta > 1.5: # Get the inclusive R9 SF for EE
                scale_factor_pixVeto = ph_pixVeto_scale_histo_2017.GetBinContent(4)
                ph_pixVeto_err       = ph_pixVeto_scale_histo_2017.GetBinError(4)



        if runningEra == 2:  # Scale factors for 2018

            local_ph_pt = ph_pt
            if local_ph_pt > 499.: # This is because corrections are up to 499 GeV
                local_ph_pt = 499.

            local_ph_eta = ph_eta
            if local_ph_eta >= 2.5:
                local_ph_eta = 2.49
            if local_ph_eta < -2.5: # Corrections reach down to eta = -2.5, but only up to eta = 2.49
                local_ph_eta = -2.5

            
            scale_factor_ID      = ph_ID_scale_histo_2018.GetBinContent( ph_ID_scale_histo_2018.GetXaxis().FindBin(local_ph_eta), ph_ID_scale_histo_2018.GetYaxis().FindBin(local_ph_pt) )
            ph_ID_err            = ph_ID_scale_histo_2018.GetBinError( ph_ID_scale_histo_2018.GetXaxis().FindBin(local_ph_eta), ph_ID_scale_histo_2018.GetYaxis().FindBin(local_ph_pt) )

            scale_factor_pixVeto = ph_pixVeto_scale_histo_2018.GetBinContent( ph_pixVeto_scale_histo_2018.GetXaxis().FindBin(local_ph_pt), ph_pixVeto_scale_histo_2018.GetYaxis().FindBin(math.fabs(local_ph_eta)) )
            ph_pixVeto_err       = ph_pixVeto_scale_histo_2018.GetBinError( ph_pixVeto_scale_histo_2018.GetXaxis().FindBin(local_ph_pt), ph_pixVeto_scale_histo_2018.GetYaxis().FindBin(math.fabs(local_ph_eta)) )


        scale_factor = scale_factor_ID * scale_factor_pixVeto
        tot_err      = math.sqrt( scale_factor_pixVeto * scale_factor_pixVeto * ph_ID_err * ph_ID_err + scale_factor_ID * scale_factor_ID * ph_pixVeto_err * ph_pixVeto_err )


        return scale_factor, tot_err

    ###############################################################################################################################################

    def get_muon_scale(self, lep_pt, lep_eta, isSingleMuTrigger_LOW, runningEra):

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


        ###############################################
        #                                             #
        #------------ Get muon SF for 2016 -----------#
        #                                             #
        ###############################################

        if runningEra == 0:

            luminosity_norm = 35.86
            luminosity_BtoF = 19.72 #For 2016
            luminosity_GH   = 16.14 #For 2016

            # Use a random number to select which muon scale factor to use, depending on the respective lumi fraction (only for 2016)
            Nrandom_for_SF = ROOT.TRandom3(44317).Rndm()
            
            if Nrandom_for_SF <= (luminosity_BtoF/luminosity_norm): # Access muon SF: B to F

                scale_factor_ID       = mu_ID_scale_histo_BCDEF_2016.GetBinContent( mu_ID_scale_histo_BCDEF_2016.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_ID_scale_histo_BCDEF_2016.GetYaxis().FindBin(local_lep_pt) )
                mu_ID_err             = mu_ID_scale_histo_BCDEF_2016.GetBinError( mu_ID_scale_histo_BCDEF_2016.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_ID_scale_histo_BCDEF_2016.GetYaxis().FindBin(local_lep_pt) )
                scale_factor_Iso      = mu_Iso_scale_histo_BCDEF_2016.GetBinContent( mu_Iso_scale_histo_BCDEF_2016.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Iso_scale_histo_BCDEF_2016.GetYaxis().FindBin(local_lep_pt) )
                mu_Iso_err            = mu_Iso_scale_histo_BCDEF_2016.GetBinError( mu_Iso_scale_histo_BCDEF_2016.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Iso_scale_histo_BCDEF_2016.GetYaxis().FindBin(local_lep_pt) )


                if isSingleMuTrigger_LOW: # Trigger SF go up to higher energies than 120 GeV, so no local_lep_pt is used for those

                    local_lep_pt_forTrigger = lep_pt # Corrections related to this trigger start from 26 GeV
                    if local_lep_pt_forTrigger < 26.:
                        local_lep_pt_forTrigger = 26.

                    scale_factor_Trigger = mu_Trigger_scale_histo_BCDEF_Mu24_2016.GetBinContent( mu_Trigger_scale_histo_BCDEF_Mu24_2016.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_BCDEF_Mu24_2016.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                    mu_Trigger_err       = mu_Trigger_scale_histo_BCDEF_Mu24_2016.GetBinError( mu_Trigger_scale_histo_BCDEF_Mu24_2016.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_BCDEF_Mu24_2016.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                    
                    scale_factor         = scale_factor_ID * scale_factor_Iso * scale_factor_Trigger
                    tot_err              = math.sqrt( scale_factor_Iso * scale_factor_Iso * scale_factor_Trigger * scale_factor_Trigger * mu_ID_err * mu_ID_err + scale_factor_ID * scale_factor_ID * scale_factor_Trigger * scale_factor_Trigger * mu_Iso_err * mu_Iso_err + scale_factor_ID * scale_factor_ID * scale_factor_Iso * scale_factor_Iso * mu_Trigger_err * mu_Trigger_err )


                else: # An event can have more than one trigger

                    local_lep_pt_forTrigger = lep_pt # Corrections related to this trigger start from 52 GeV
                    if local_lep_pt_forTrigger < 52.:
                        local_lep_pt_forTrigger = 52.

                    scale_factor_Trigger = mu_Trigger_scale_histo_BCDEF_Mu50_2016.GetBinContent( mu_Trigger_scale_histo_BCDEF_Mu50_2016.GetXaxis().FindBin(math.fabs(lep_eta)), mu_Trigger_scale_histo_BCDEF_Mu50_2016.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                    mu_Trigger_err       = mu_Trigger_scale_histo_BCDEF_Mu50_2016.GetBinError( mu_Trigger_scale_histo_BCDEF_Mu50_2016.GetXaxis().FindBin(math.fabs(lep_eta)), mu_Trigger_scale_histo_BCDEF_Mu50_2016.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                    
                    scale_factor         = scale_factor_ID * scale_factor_Iso * scale_factor_Trigger
                    tot_err              = math.sqrt( scale_factor_Iso * scale_factor_Iso * scale_factor_Trigger * scale_factor_Trigger * mu_ID_err * mu_ID_err + scale_factor_ID * scale_factor_ID * scale_factor_Trigger * scale_factor_Trigger * mu_Iso_err * mu_Iso_err + scale_factor_ID * scale_factor_ID * scale_factor_Iso * scale_factor_Iso * mu_Trigger_err * mu_Trigger_err )



            else: #Use GH SFs

                scale_factor_ID       = mu_ID_scale_histo_GH_2016.GetBinContent( mu_ID_scale_histo_GH_2016.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_ID_scale_histo_GH_2016.GetYaxis().FindBin(local_lep_pt) )
                mu_ID_err             = mu_ID_scale_histo_GH_2016.GetBinError( mu_ID_scale_histo_GH_2016.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_ID_scale_histo_GH_2016.GetYaxis().FindBin(local_lep_pt) )
                scale_factor_Iso      = mu_Iso_scale_histo_GH_2016.GetBinContent( mu_Iso_scale_histo_GH_2016.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Iso_scale_histo_GH_2016.GetYaxis().FindBin(local_lep_pt) )
                mu_Iso_err            = mu_Iso_scale_histo_GH_2016.GetBinError( mu_Iso_scale_histo_GH_2016.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Iso_scale_histo_GH_2016.GetYaxis().FindBin(local_lep_pt) )

                if isSingleMuTrigger_LOW: # Trigger SF go up to higher energies than 120 GeV, so no local_lep_pt is used for those

                    local_lep_pt_forTrigger = lep_pt 
                    if local_lep_pt_forTrigger < 26.:
                        local_lep_pt_forTrigger = 26.

                    scale_factor_Trigger = mu_Trigger_scale_histo_GH_Mu24_2016.GetBinContent( mu_Trigger_scale_histo_GH_Mu24_2016.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_GH_Mu24_2016.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                    mu_Trigger_err       = mu_Trigger_scale_histo_GH_Mu24_2016.GetBinError( mu_Trigger_scale_histo_GH_Mu24_2016.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_GH_Mu24_2016.GetYaxis().FindBin(local_lep_pt_forTrigger) )

                    scale_factor         = scale_factor_ID * scale_factor_Iso * scale_factor_Trigger
                    tot_err              = math.sqrt( scale_factor_Iso * scale_factor_Iso * scale_factor_Trigger * scale_factor_Trigger * mu_ID_err * mu_ID_err + scale_factor_ID * scale_factor_ID * scale_factor_Trigger * scale_factor_Trigger * mu_Iso_err * mu_Iso_err + scale_factor_ID * scale_factor_ID * scale_factor_Iso * scale_factor_Iso * mu_Trigger_err * mu_Trigger_err )


                else: # An event can have more than one trigger

                    local_lep_pt_forTrigger = lep_pt 
                    if local_lep_pt_forTrigger < 52.:
                        local_lep_pt_forTrigger = 52.

                    scale_factor_Trigger = mu_Trigger_scale_histo_GH_Mu50_2016.GetBinContent( mu_Trigger_scale_histo_GH_Mu50_2016.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_GH_Mu50_2016.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                    mu_Trigger_err       = mu_Trigger_scale_histo_GH_Mu50_2016.GetBinError( mu_Trigger_scale_histo_GH_Mu50_2016.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_GH_Mu50_2016.GetYaxis().FindBin(local_lep_pt_forTrigger) )

                    scale_factor         = scale_factor_ID * scale_factor_Iso * scale_factor_Trigger
                    tot_err              = math.sqrt( scale_factor_Iso * scale_factor_Iso * scale_factor_Trigger * scale_factor_Trigger * mu_ID_err * mu_ID_err + scale_factor_ID * scale_factor_ID * scale_factor_Trigger * scale_factor_Trigger * mu_Iso_err * mu_Iso_err + scale_factor_ID * scale_factor_ID * scale_factor_Iso * scale_factor_Iso * mu_Trigger_err * mu_Trigger_err )
                    


        ###############################################
        #                                             #
        #------------ Get muon SF for 2017 -----------#
        #                                             #
        ###############################################


        if runningEra == 1: # Get muon SFs for 2017. ID and ISO histos have pT and eta positions inverted wrt 2016 histos
                              
            scale_factor_ID       = mu_ID_scale_histo_2017.GetBinContent( mu_ID_scale_histo_2017.GetXaxis().FindBin(local_lep_pt), mu_ID_scale_histo_2017.GetYaxis().FindBin(math.fabs(local_lep_eta)) )
            mu_ID_err             = mu_ID_scale_histo_2017.GetBinError( mu_ID_scale_histo_2017.GetXaxis().FindBin(local_lep_pt), mu_ID_scale_histo_2017.GetYaxis().FindBin(math.fabs(local_lep_eta)) )
            scale_factor_Iso      = mu_Iso_scale_histo_2017.GetBinContent( mu_Iso_scale_histo_2017.GetXaxis().FindBin(local_lep_pt), mu_Iso_scale_histo_2017.GetYaxis().FindBin(math.fabs(local_lep_eta)) )
            mu_Iso_err            = mu_Iso_scale_histo_2017.GetBinError( mu_Iso_scale_histo_2017.GetXaxis().FindBin(local_lep_pt),mu_Iso_scale_histo_2017.GetYaxis().FindBin(math.fabs(local_lep_eta)) )


            if isSingleMuTrigger_LOW: # Trigger SF go up to higher energies than 120 GeV, so no local_lep_pt is used for those

                local_lep_pt_forTrigger = lep_pt 
                if local_lep_pt_forTrigger < 29.:
                    local_lep_pt_forTrigger = 29.

                scale_factor_Trigger = mu_Trigger_scale_histo_2017_Mu27.GetBinContent( mu_Trigger_scale_histo_2017_Mu27.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_2017_Mu27.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                mu_Trigger_err       = mu_Trigger_scale_histo_2017_Mu27.GetBinError( mu_Trigger_scale_histo_2017_Mu27.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_2017_Mu27.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                

            else: # An event can have more than one trigger

                local_lep_pt_forTrigger = lep_pt 
                if local_lep_pt_forTrigger < 52.:
                    local_lep_pt_forTrigger = 52.

                scale_factor_Trigger = mu_Trigger_scale_histo_2017_Mu50.GetBinContent( mu_Trigger_scale_histo_2017_Mu50.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_2017_Mu50.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                mu_Trigger_err       = mu_Trigger_scale_histo_2017_Mu50.GetBinError( mu_Trigger_scale_histo_2017_Mu50.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_2017_Mu50.GetYaxis().FindBin(local_lep_pt_forTrigger) )

                
            scale_factor         = scale_factor_ID * scale_factor_Iso * scale_factor_Trigger
            tot_err              = math.sqrt( scale_factor_Iso * scale_factor_Iso * scale_factor_Trigger * scale_factor_Trigger * mu_ID_err * mu_ID_err + scale_factor_ID * scale_factor_ID * scale_factor_Trigger * scale_factor_Trigger * mu_Iso_err * mu_Iso_err + scale_factor_ID * scale_factor_ID * scale_factor_Iso * scale_factor_Iso * mu_Trigger_err * mu_Trigger_err )



        ###############################################
        #                                             #
        #------------ Get muon SF for 2018 -----------#
        #                                             #
        ###############################################

        if runningEra == 2:

            luminosity_norm = 59.69
            luminosity_BeforeHLTUpdate = 8.98
            luminosity_AfterHLTUpdate = 50.71

            # Use a random number to select which muon scale factor to use for trigger, depending on the respective lumi fraction (only for 2018)
            Nrandom_for_SF = ROOT.TRandom3(92562).Rndm()
            
            scale_factor_ID       = mu_ID_scale_histo_2018.GetBinContent( mu_ID_scale_histo_2018.GetXaxis().FindBin(local_lep_pt), mu_ID_scale_histo_2018.GetYaxis().FindBin(math.fabs(local_lep_eta)) )
            mu_ID_err             = mu_ID_scale_histo_2018.GetBinError( mu_ID_scale_histo_2018.GetXaxis().FindBin(local_lep_pt), mu_ID_scale_histo_2018.GetYaxis().FindBin(math.fabs(local_lep_eta)) )
            scale_factor_Iso      = mu_Iso_scale_histo_2018.GetBinContent( mu_Iso_scale_histo_2018.GetXaxis().FindBin(local_lep_pt), mu_Iso_scale_histo_2018.GetYaxis().FindBin(math.fabs(local_lep_eta)) )
            mu_Iso_err            = mu_Iso_scale_histo_2018.GetBinError( mu_Iso_scale_histo_2018.GetXaxis().FindBin(local_lep_pt), mu_Iso_scale_histo_2018.GetYaxis().FindBin(math.fabs(local_lep_eta)) )

            if Nrandom_for_SF <= (luminosity_BeforeHLTUpdate/luminosity_norm): # Access muon Trigger SF: Before HLT Update

                if isSingleMuTrigger_LOW: # Trigger SF go up to higher energies than 120 GeV, so no local_lep_pt is used for those

                    local_lep_pt_forTrigger = lep_pt # Corrections related to this trigger start from 26 GeV
                    if local_lep_pt_forTrigger < 26.:
                        local_lep_pt_forTrigger = 26.

                    scale_factor_Trigger = mu_Trigger_scale_histo_BeforeHLTUpdate_2018_Mu24.GetBinContent( mu_Trigger_scale_histo_BeforeHLTUpdate_2018_Mu24.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_BeforeHLTUpdate_2018_Mu24.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                    mu_Trigger_err       = mu_Trigger_scale_histo_BeforeHLTUpdate_2018_Mu24.GetBinError( mu_Trigger_scale_histo_BeforeHLTUpdate_2018_Mu24.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_BeforeHLTUpdate_2018_Mu24.GetYaxis().FindBin(local_lep_pt_forTrigger) )


                else: # Trigger SF go up to higher energies than 120 GeV, so no local_lep_pt is used for those

                    local_lep_pt_forTrigger = lep_pt # Corrections related to this trigger start from 52 GeV
                    if local_lep_pt_forTrigger < 52.:
                        local_lep_pt_forTrigger = 52.
                
                    scale_factor_Trigger = mu_Trigger_scale_histo_BeforeHLTUpdate_2018_Mu50.GetBinContent( mu_Trigger_scale_histo_BeforeHLTUpdate_2018_Mu50.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_BeforeHLTUpdate_2018_Mu50.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                    mu_Trigger_err       = mu_Trigger_scale_histo_BeforeHLTUpdate_2018_Mu50.GetBinError( mu_Trigger_scale_histo_BeforeHLTUpdate_2018_Mu50.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_BeforeHLTUpdate_2018_Mu50.GetYaxis().FindBin(local_lep_pt_forTrigger) )


            else: # Access muon Trigger SF: After HLT Update

                if isSingleMuTrigger_LOW: # Trigger SF go up to higher energies than 120 GeV, so no local_lep_pt is used for those

                    local_lep_pt_forTrigger = lep_pt # Corrections related to this trigger start from 26 GeV
                    if local_lep_pt_forTrigger < 26.:
                        local_lep_pt_forTrigger = 26.

                    scale_factor_Trigger = mu_Trigger_scale_histo_AfterHLTUpdate_2018_Mu24.GetBinContent( mu_Trigger_scale_histo_AfterHLTUpdate_2018_Mu24.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_AfterHLTUpdate_2018_Mu24.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                    mu_Trigger_err       = mu_Trigger_scale_histo_AfterHLTUpdate_2018_Mu24.GetBinError( mu_Trigger_scale_histo_AfterHLTUpdate_2018_Mu24.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_AfterHLTUpdate_2018_Mu24.GetYaxis().FindBin(local_lep_pt_forTrigger) )


                else: # Trigger SF go up to higher energies than 120 GeV, so no local_lep_pt is used for those

                    local_lep_pt_forTrigger = lep_pt # Corrections related to this trigger start from 52 GeV
                    if local_lep_pt_forTrigger < 52.:
                        local_lep_pt_forTrigger = 52.
                
                    scale_factor_Trigger = mu_Trigger_scale_histo_AfterHLTUpdate_2018_Mu50.GetBinContent( mu_Trigger_scale_histo_AfterHLTUpdate_2018_Mu50.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_AfterHLTUpdate_2018_Mu50.GetYaxis().FindBin(local_lep_pt_forTrigger) )
                    mu_Trigger_err       = mu_Trigger_scale_histo_AfterHLTUpdate_2018_Mu50.GetBinError( mu_Trigger_scale_histo_AfterHLTUpdate_2018_Mu50.GetXaxis().FindBin(math.fabs(local_lep_eta)), mu_Trigger_scale_histo_AfterHLTUpdate_2018_Mu50.GetYaxis().FindBin(local_lep_pt_forTrigger) )


            scale_factor         = scale_factor_ID * scale_factor_Iso * scale_factor_Trigger
            tot_err              = math.sqrt( scale_factor_Iso * scale_factor_Iso * scale_factor_Trigger * scale_factor_Trigger * mu_ID_err * mu_ID_err + scale_factor_ID * scale_factor_ID * scale_factor_Trigger * scale_factor_Trigger * mu_Iso_err * mu_Iso_err + scale_factor_ID * scale_factor_ID * scale_factor_Iso * scale_factor_Iso * mu_Trigger_err * mu_Trigger_err )


        return scale_factor, tot_err
