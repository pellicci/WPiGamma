import ROOT
import math
import os
from array import array

#------Acessing scale factor histos for EGamma and Muon corrections --------#

eg_reco_scale_name  = "scale_factors/Electron_reco_2D.root"
eg_reco_scale_file  = ROOT.TFile(eg_reco_scale_name)
eg_reco_scale_histo = ROOT.TH2F()
eg_reco_scale_histo = eg_reco_scale_file.Get("EGamma_SF2D")

eg_ID_scale_name  = "scale_factors/Electron_ID_2D.root"
eg_ID_scale_file  = ROOT.TFile(eg_ID_scale_name)
eg_ID_scale_histo = ROOT.TH2F()
eg_ID_scale_histo = eg_ID_scale_file.Get("EGamma_SF2D")

ph_ID_scale_name  = "scale_factors/Photon_ID_2D.root"
ph_ID_scale_file  = ROOT.TFile(ph_ID_scale_name)
ph_ID_scale_histo = ROOT.TH2F()
ph_ID_scale_histo = ph_ID_scale_file.Get("EGamma_SF2D")

ph_pixVeto_scale_name  = "scale_factors/Photon_pixVeto_2D.root"
ph_pixVeto_scale_file  = ROOT.TFile(ph_pixVeto_scale_name)
ph_pixVeto_scale_histo = ROOT.TH2F()
ph_pixVeto_scale_histo = ph_pixVeto_scale_file.Get("Scaling_Factors_HasPix_R9 Inclusive")

mu_ID_scale_name_BCDEF  = "scale_factors/EfficienciesAndSF_BCDEF.root"
mu_ID_scale_file_BCDEF  = ROOT.TFile(mu_ID_scale_name_BCDEF)
mu_ID_scale_histo_BCDEF = ROOT.TH2F()
mu_ID_scale_histo_BCDEF = mu_ID_scale_file_BCDEF.Get("MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio")

mu_ID_scale_name_GH  = "scale_factors/EfficienciesAndSF_GH.root"
mu_ID_scale_file_GH  = ROOT.TFile(mu_ID_scale_name_GH)
mu_ID_scale_histo_GH = ROOT.TH2F()
mu_ID_scale_histo_GH = mu_ID_scale_file_GH.Get("MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio")

mu_Iso_scale_name_BCDEF  = "scale_factors/EfficienciesAndSF_Iso_BCDEF.root"
mu_Iso_scale_file_BCDEF  = ROOT.TFile(mu_Iso_scale_name_BCDEF)
mu_Iso_scale_histo_BCDEF = ROOT.TH2F()
mu_Iso_scale_histo_BCDEF = mu_Iso_scale_file_BCDEF.Get("LooseISO_MediumID_pt_eta/abseta_pt_ratio")

mu_Iso_scale_name_GH  = "scale_factors/EfficienciesAndSF_Iso_GH.root"
mu_Iso_scale_file_GH  = ROOT.TFile(mu_Iso_scale_name_GH)
mu_Iso_scale_histo_GH = ROOT.TH2F()
mu_Iso_scale_histo_GH = mu_Iso_scale_file_GH.Get("LooseISO_MediumID_pt_eta/abseta_pt_ratio")

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

mu_Tracking_scale_name_BCDEFGH  = "scale_factors/Tracking_EfficienciesAndSF_BCDEFGH.root"
mu_Tracking_scale_file_BCDEFGH  = ROOT.TFile(mu_Tracking_scale_name_BCDEFGH)
mu_Tracking_scale_graph_BCDEFGH = ROOT.TGraphAsymmErrors()
mu_Tracking_scale_graph_BCDEFGH = mu_Tracking_scale_file_BCDEFGH.Get("ratio_eff_aeta_dr030e030_corr")

#----------Some arrays for the BDT----------#
pi_pT_array           = array('f', [0.])
gamma_eT_array        = array('f', [0.])
nBjets_25_array       = array('f', [0.])
deltaphi_lep_pi_array = array('f', [0.])
lep_pT_array          = array('f', [0.])
piRelIso_05_array     = array('f', [0.])
met_array             = array('f', [0.])
isMuon_array          = array('i', [0])

reader = ROOT.TMVA.Reader("!Color")


class Workflow_Handler:

    def __init__(self,signalname,dataname,isMedium,subprocess="//"):

        ##Where the root files are
        self.subprocess = subprocess

        self.data_samplename = dataname  
        self.sig_samplename = signalname
        #self.isMedium = isMedium
        
        #----Medium selection----#
        if isMedium:
            #self.norm_filename = "rootfiles/" + self.subprocess + "Medium/Normalizations_table.txt"
            self.norm_filename = "rootfiles/" + self.subprocess + "Medium_AfterFix/Normalizations_table.txt"
            #self.dir_back_input = "rootfiles/" + self.subprocess + "Medium/backgrounds/"
            self.dir_back_input = "rootfiles/" + self.subprocess + "Medium_AfterFix/backgrounds/"
            #self.sig_filename = "rootfiles/" + self.subprocess + "Medium/signals/" + "WPiGammaAnalysis_" + self.sig_samplename + ".root"
            self.sig_filename = "rootfiles/" + self.subprocess + "Medium_AfterFix/signals/" + "WPiGammaAnalysis_" + self.sig_samplename + ".root"

        #----Tight selection----#
        if not isMedium:
            self.norm_filename = "rootfiles/" + self.subprocess + "Tight/Normalizations_table.txt" 
            self.dir_back_input = "rootfiles/" + self.subprocess + "Tight/backgrounds/" 
            self.sig_filename = "rootfiles/" + self.subprocess + "Tight/signals/" + "WPiGammaAnalysis_" + self.sig_samplename + ".root"

        #----Data----#
        #self.dir_data_input = "rootfiles/" + self.subprocess + "data/"
        self.dir_data_input = "rootfiles/" + self.subprocess + "data_AfterFix/"

                
        reader.AddVariable("pi_pT",pi_pT_array)
        reader.AddVariable("gamma_eT",gamma_eT_array)
        reader.AddVariable("nBjets_25",nBjets_25_array)
        #reader.AddVariable("deltaphi_lep_pi",deltaphi_lep_pi_array)
        reader.AddVariable("lep_pT",lep_pT_array)
        reader.AddVariable("piRelIso_05",piRelIso_05_array)
        #reader.AddSpectator("isMuon",isMuon_array)
        reader.AddVariable("MET",met_array)

        reader.BookMVA("BDT_mu","MVA/weights/TMVAClassification_BDT.weights_mu.xml")#The first argument is arbitrary. To be chosen in order to distinguish among methods
        reader.BookMVA("BDT_ele","MVA/weights/TMVAClassification_BDT.weights_ele.xml")
        
    def get_ele_scale(self, ele_pt, ele_eta):
        #This is because corrections are up to 150 GeV
        local_ele_pt = ele_pt
        if local_ele_pt > 150.:
            local_ele_pt = 150.

        scale_factor      = 1.
        scale_factor_reco = 1.
        scale_factor_ID   = 1.
        eg_reco_err       = 0.
        eg_ID_err         = 0.
        tot_err           = 0.

        scale_factor_reco = eg_reco_scale_histo.GetBinContent( eg_reco_scale_histo.GetXaxis().FindBin(ele_eta), eg_reco_scale_histo.GetYaxis().FindBin(local_ele_pt) )
        eg_reco_err       = eg_reco_scale_histo.GetBinError( eg_reco_scale_histo.GetXaxis().FindBin(ele_eta), eg_reco_scale_histo.GetYaxis().FindBin(local_ele_pt) )
        scale_factor_ID   = eg_ID_scale_histo.GetBinContent( eg_ID_scale_histo.GetXaxis().FindBin(ele_eta), eg_ID_scale_histo.GetYaxis().FindBin(local_ele_pt) )
        eg_ID_err         = eg_ID_scale_histo.GetBinError( eg_ID_scale_histo.GetXaxis().FindBin(ele_eta), eg_ID_scale_histo.GetYaxis().FindBin(local_ele_pt) )

        scale_factor = scale_factor_reco * scale_factor_ID
        tot_err      = math.sqrt( scale_factor_reco * scale_factor_reco * eg_ID_err * eg_ID_err + scale_factor_ID * scale_factor_ID * eg_reco_err * eg_reco_err )

        return scale_factor, tot_err

    def get_photon_scale(self, ph_pt, ph_eta):
        #This is because corrections are up to 150 GeV
        local_ph_pt = ph_pt
        if local_ph_pt > 150.:
            local_ph_pt = 150.

        scale_factor         = 1.
        scale_factor_ID      = 1.
        scale_factor_pixVeto = 1.
        ph_ID_err            = 0.
        ph_pixVeto_err       = 0.
        tot_err              = 0.

        scale_factor_ID      = ph_ID_scale_histo.GetBinContent( ph_ID_scale_histo.GetXaxis().FindBin(ph_eta), ph_ID_scale_histo.GetYaxis().FindBin(local_ph_pt) )
        ph_ID_err            = ph_ID_scale_histo.GetBinError( ph_ID_scale_histo.GetXaxis().FindBin(ph_eta), ph_ID_scale_histo.GetYaxis().FindBin(local_ph_pt) )
        scale_factor_pixVeto = ph_pixVeto_scale_histo.GetBinContent( ph_pixVeto_scale_histo.GetXaxis().FindBin(math.fabs(ph_eta)), ph_pixVeto_scale_histo.GetYaxis().FindBin(local_ph_pt) )
        ph_pixVeto_err       = ph_pixVeto_scale_histo.GetBinError( ph_pixVeto_scale_histo.GetXaxis().FindBin(ph_eta), ph_pixVeto_scale_histo.GetYaxis().FindBin(local_ph_pt) )

        scale_factor = scale_factor_ID * scale_factor_pixVeto
        tot_err      = math.sqrt( scale_factor_pixVeto * scale_factor_pixVeto * ph_ID_err * ph_ID_err + scale_factor_ID * scale_factor_ID * ph_pixVeto_err * ph_pixVeto_err )

        return scale_factor, tot_err

    def get_muon_scale_BtoF(self, lep_pt, lep_eta, isSingleMuTrigger_24):

        if lep_pt > 120. : # This is because corrections go up to 120 GeV (excluded)
            lep_pt = 119.

        scale_factor            = 1.
        scale_factor_ID         = 1.
        scale_factor_Iso        = 1.
        scale_factor_Tracking   = 1.
        scale_factor_Trigger    = 1.
        mu_ID_err               = 0.
        mu_Iso_err              = 0.
        mu_Trigger_err          = 0.
        tot_err                 = 0.
        

        scale_factor_ID       = mu_ID_scale_histo_BCDEF.GetBinContent( mu_ID_scale_histo_BCDEF.GetXaxis().FindBin(math.fabs(lep_eta)), mu_ID_scale_histo_BCDEF.GetYaxis().FindBin(lep_pt) )
        mu_ID_err             = mu_ID_scale_histo_BCDEF.GetBinError( mu_ID_scale_histo_BCDEF.GetXaxis().FindBin(math.fabs(lep_eta)), mu_ID_scale_histo_BCDEF.GetYaxis().FindBin(lep_pt) )
        scale_factor_Iso      = mu_Iso_scale_histo_BCDEF.GetBinContent( mu_Iso_scale_histo_BCDEF.GetXaxis().FindBin(math.fabs(lep_eta)), mu_Iso_scale_histo_BCDEF.GetYaxis().FindBin(lep_pt) )
        mu_Iso_err            = mu_Iso_scale_histo_BCDEF.GetBinError( mu_Iso_scale_histo_BCDEF.GetXaxis().FindBin(math.fabs(lep_eta)), mu_Iso_scale_histo_BCDEF.GetYaxis().FindBin(lep_pt) )
        scale_factor_Tracking = mu_Tracking_scale_graph_BCDEFGH.Eval(math.fabs(lep_eta))

        if isSingleMuTrigger_24:
            scale_factor_Trigger = mu_Trigger_scale_histo_BCDEF_Mu24.GetBinContent( mu_Trigger_scale_histo_BCDEF_Mu24.GetXaxis().FindBin(math.fabs(lep_eta)), mu_Trigger_scale_histo_BCDEF_Mu24.GetYaxis().FindBin(lep_pt) )
            mu_Trigger_err       = mu_Trigger_scale_histo_BCDEF_Mu24.GetBinError( mu_Trigger_scale_histo_BCDEF_Mu24.GetXaxis().FindBin(math.fabs(lep_eta)), mu_Trigger_scale_histo_BCDEF_Mu24.GetYaxis().FindBin(lep_pt) )

            scale_factor         = scale_factor_ID * scale_factor_Iso * scale_factor_Tracking * scale_factor_Trigger
            tot_err              = math.sqrt( scale_factor_Iso * scale_factor_Iso * scale_factor_Trigger * scale_factor_Trigger * scale_factor_Tracking * scale_factor_Tracking * mu_ID_err * mu_ID_err + scale_factor_ID * scale_factor_ID * scale_factor_Trigger * scale_factor_Trigger * scale_factor_Tracking * scale_factor_Tracking * mu_Iso_err * mu_Iso_err + scale_factor_ID * scale_factor_ID * scale_factor_Iso * scale_factor_Iso * scale_factor_Tracking * scale_factor_Tracking * mu_Trigger_err * mu_Trigger_err )

            return scale_factor, tot_err

        else: # An event can have more than one trigger
            scale_factor_Trigger = mu_Trigger_scale_histo_BCDEF_Mu50.GetBinContent( mu_Trigger_scale_histo_BCDEF_Mu50.GetXaxis().FindBin(math.fabs(lep_eta)), mu_Trigger_scale_histo_BCDEF_Mu50.GetYaxis().FindBin(lep_pt) )
            mu_Trigger_err       = mu_Trigger_scale_histo_BCDEF_Mu50.GetBinError( mu_Trigger_scale_histo_BCDEF_Mu50.GetXaxis().FindBin(math.fabs(lep_eta)), mu_Trigger_scale_histo_BCDEF_Mu50.GetYaxis().FindBin(lep_pt) )

            scale_factor         = scale_factor_ID * scale_factor_Iso * scale_factor_Tracking * scale_factor_Trigger
            tot_err              = math.sqrt( scale_factor_Iso * scale_factor_Iso * scale_factor_Trigger * scale_factor_Trigger * scale_factor_Tracking * scale_factor_Tracking * mu_ID_err * mu_ID_err + scale_factor_ID * scale_factor_ID * scale_factor_Trigger * scale_factor_Trigger * scale_factor_Tracking * scale_factor_Tracking * mu_Iso_err * mu_Iso_err + scale_factor_ID * scale_factor_ID * scale_factor_Iso * scale_factor_Iso * scale_factor_Tracking * scale_factor_Tracking * mu_Trigger_err * mu_Trigger_err )

            return scale_factor, tot_err

    def get_muon_scale_GH(self, lep_pt, lep_eta, isSingleMuTrigger_24):

        if lep_pt > 120.: # This is because corrections go up to 120 GeV (excluded)
            lep_pt = 119.

        scale_factor            = 1.
        scale_factor_ID         = 1.
        scale_factor_Iso        = 1.
        scale_factor_Tracking   = 1.
        scale_factor_Trigger    = 1.
        mu_ID_err               = 0.
        mu_Iso_err              = 0.
        mu_Trigger_err          = 0.
        tot_err                 = 0.
                              
        scale_factor_ID       = mu_ID_scale_histo_GH.GetBinContent( mu_ID_scale_histo_GH.GetXaxis().FindBin(math.fabs(lep_eta)), mu_ID_scale_histo_GH.GetYaxis().FindBin(lep_pt) )
        mu_ID_err             = mu_ID_scale_histo_GH.GetBinError( mu_ID_scale_histo_GH.GetXaxis().FindBin(math.fabs(lep_eta)), mu_ID_scale_histo_GH.GetYaxis().FindBin(lep_pt) )
        scale_factor_Iso      = mu_Iso_scale_histo_GH.GetBinContent( mu_Iso_scale_histo_GH.GetXaxis().FindBin(math.fabs(lep_eta)), mu_Iso_scale_histo_GH.GetYaxis().FindBin(lep_pt) )
        mu_Iso_err            = mu_Iso_scale_histo_GH.GetBinError( mu_Iso_scale_histo_GH.GetXaxis().FindBin(math.fabs(lep_eta)), mu_Iso_scale_histo_GH.GetYaxis().FindBin(lep_pt) )
        scale_factor_Tracking = mu_Tracking_scale_graph_BCDEFGH.Eval(math.fabs(lep_eta))

        if isSingleMuTrigger_24:
            scale_factor_Trigger = mu_Trigger_scale_histo_GH_Mu24.GetBinContent( mu_Trigger_scale_histo_GH_Mu24.GetXaxis().FindBin(math.fabs(lep_eta)), mu_Trigger_scale_histo_GH_Mu24.GetYaxis().FindBin(lep_pt) )
            mu_Trigger_err       = mu_Trigger_scale_histo_GH_Mu24.GetBinError( mu_Trigger_scale_histo_GH_Mu24.GetXaxis().FindBin(math.fabs(lep_eta)), mu_Trigger_scale_histo_GH_Mu24.GetYaxis().FindBin(lep_pt) )

            scale_factor         = scale_factor_ID * scale_factor_Iso * scale_factor_Tracking * scale_factor_Trigger
            tot_err              = math.sqrt( scale_factor_Iso * scale_factor_Iso * scale_factor_Trigger * scale_factor_Trigger * scale_factor_Tracking * scale_factor_Tracking * mu_ID_err * mu_ID_err + scale_factor_ID * scale_factor_ID * scale_factor_Trigger * scale_factor_Trigger * scale_factor_Tracking * scale_factor_Tracking * mu_Iso_err * mu_Iso_err + scale_factor_ID * scale_factor_ID * scale_factor_Iso * scale_factor_Iso * scale_factor_Tracking * scale_factor_Tracking * mu_Trigger_err * mu_Trigger_err )

            return scale_factor, tot_err

        else: # An event can have more than one trigger
            scale_factor_Trigger = mu_Trigger_scale_histo_GH_Mu50.GetBinContent( mu_Trigger_scale_histo_GH_Mu50.GetXaxis().FindBin(math.fabs(lep_eta)), mu_Trigger_scale_histo_GH_Mu50.GetYaxis().FindBin(lep_pt) )
            mu_Trigger_err       = mu_Trigger_scale_histo_GH_Mu50.GetBinError( mu_Trigger_scale_histo_GH_Mu50.GetXaxis().FindBin(math.fabs(lep_eta)), mu_Trigger_scale_histo_GH_Mu50.GetYaxis().FindBin(lep_pt) )

            scale_factor         = scale_factor_ID * scale_factor_Iso * scale_factor_Tracking * scale_factor_Trigger
            tot_err              = math.sqrt( scale_factor_Iso * scale_factor_Iso * scale_factor_Trigger * scale_factor_Trigger * scale_factor_Tracking * scale_factor_Tracking * mu_ID_err * mu_ID_err + scale_factor_ID * scale_factor_ID * scale_factor_Trigger * scale_factor_Trigger * scale_factor_Tracking * scale_factor_Tracking * mu_Iso_err * mu_Iso_err + scale_factor_ID * scale_factor_ID * scale_factor_Iso * scale_factor_Iso * scale_factor_Tracking * scale_factor_Tracking * mu_Trigger_err * mu_Trigger_err )

            return scale_factor, tot_err
    """
    def get_muon_scale_tracking_BtoH(self, lep_eta):

        scale_factor = 1.


        tot_err      = 
        
        return scale_factor
    """
    def get_normalizations_map(self):
        in_file = open(self.norm_filename,"r")
        norm_map = dict()

        for line in in_file:
            data_norm = line.split()
            norm_map[data_norm[0]] = float(data_norm[1])

        return norm_map

    def get_BDT_output(self,pi_pT,gamma_eT,nBjets_25,deltaphi_lep_pi,lep_pT,piRelIso_05,isMuon):

        pi_pT_array[0] = pi_pT
        gamma_eT_array[0] = gamma_eT
        nBjets_25_array[0] = nBjets_25
        deltaphi_lep_pi_array[0] = deltaphi_lep_pi
        lep_pT_array[0] = lep_pT
        piRelIso_05_array[0] = piRelIso_05
        isMuon_array[0] = int(isMuon)
        
        if isMuon:
            return reader.EvaluateMVA("BDT_mu")
        else:
            return reader.EvaluateMVA("BDT_ele")

   # get the sample names and signal names
    def get_samples_names(self, Add_Signal=True,Add_Data=True):
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

   # get data sample names
    def get_dataSample_names(self):
        list_dirs_data = os.listdir(self.dir_data_input)
        dataName_list = []

        for dirname in list_dirs_data:
            tmp_dataName = dirname.split("WPiGammaAnalysis_")[1]
            tmp_dataName2 = tmp_dataName.replace(".root","")
            dataName_list.append(tmp_dataName2)

        return dataName_list

  # get the corresponding root files for the background, signal and data sample names
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
