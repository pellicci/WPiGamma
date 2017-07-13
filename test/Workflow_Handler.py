import ROOT
import math
import os

class Workflow_Handler:

    def __init__(self,signalname, subprocess="//"):

        ##Where the root files are
        self.subprocess = subprocess

        #self.norm_filename = "rootfiles/" + self.subprocess + "Tight/Normalizations_table.txt"
        self.norm_filename = "rootfiles/" + self.subprocess + "Medium/Normalizations_table.txt"
        #self.dir_back_input = "rootfiles/" + self.subprocess + "Tight/backgrounds/"
        self.dir_back_input = "rootfiles/" + self.subprocess + "Medium/backgrounds/"
        self.dir_data_input = "rootfiles/" + self.subprocess + "data/"      

        self.sig_samplename = signalname
        #self.sig_filename = "rootfiles/" + self.subprocess + "Tight/signals/" + "WPiGammaAnalysis_" + self.sig_samplename + ".root"
        self.sig_filename = "rootfiles/" + self.subprocess + "Medium/signals/" + "WPiGammaAnalysis_" + self.sig_samplename + ".root"

        eg_reco_scale_name = "scale_factors/Electron_reco_2D.root"
        eg_reco_scale_file = ROOT.TFile(eg_reco_scale_name)
        self.eg_reco_scale_histo = eg_reco_scale_file.Get("EGamma_SF2D")

        eg_ID_scale_name = "scale_factors/Electron_ID_2D.root"
        eg_ID_scale_file = ROOT.TFile(eg_ID_scale_name)
        self.eg_ID_scale_histo = eg_ID_scale_file.Get("EGamma_SF2D")

        ph_ID_scale_name = "scale_factors/Photon_ID_2D.root"
        ph_ID_scale_file = ROOT.TFile(ph_ID_scale_name)
        self.ph_ID_scale_histo = ph_ID_scale_file.Get("EGamma_SF2D")

        ph_pixVeto_scale_name = "scale_factors/Photon_pixVeto_2D.root"
        ph_pixVeto_scale_file = ROOT.TFile(ph_pixVeto_scale_name)
        self.ph_pixVeto_scale_histo = ph_pixVeto_scale_file.Get("Scaling_Factors_HasPix_R9 Inclusive")
        
    def get_ele_scale(self, ele_pt, ele_eta):
        #This is because corrections are up to 200 GeV
        local_ele_pt = ele_pt
        if local_ele_pt >= 200.:
            local_ele_pt = 199.

        scale_factor = 1.
        scale_factor = scale_factor * eg_reco_scale_histo.GetBinContent( eg_reco_scale_histo.GetXaxis().FindBin(ele_eta), eg_reco_scale_histo.GetYaxis().FindBin(local_ele_pt) )
        scale_factor = scale_factor * eg_ID_scale_histo.GetBinContent( eg_ID_scale_histo.GetXaxis().FindBin(ele_eta), eg_ID_scale_histo.GetYaxis().FindBin(local_ele_pt) )

        return scale_factor

    def get_photon_scale(self, ph_pt, ph_eta):
        #This is because corrections are up to 200 GeV
        local_ph_pt = ph_pt
        if local_ph_pt >= 200.:
            local_ph_pt = 199.

        scale_factor = 1.
        scale_factor = scale_factor * ph_ID_scale_histo.GetBinContent( ph_ID_scale_histo.GetXaxis().FindBin(ph_eta), ph_ID_scale_histo.GetYaxis().FindBin(local_ph_pt) )
        scale_factor = scale_factor * ph_pixVeto_scale_histo.GetBinContent( ph_pixVeto_scale_histo.GetXaxis().FindBin(math.abs(ph_eta)), ph_pixVeto_scale_histo.GetYaxis().FindBin(local_ph_pt) )

        return scale_factor
    
    def get_normalizations_map(self):
        in_file = open(self.norm_filename,"r")
        norm_map = dict()

        for line in in_file:
            data_norm = line.split()
            norm_map[data_norm[0]] = float(data_norm[1])

        return norm_map

   # get the sample names and signal names
    def get_samples_names(self, Add_Signal=True):
        list_dirs_bkg = os.listdir(self.dir_back_input)
        samplename_list = []

        for dirname in list_dirs_bkg:
            tmp_samplename = dirname.split("WPiGammaAnalysis_")[1]
            tmp_samplename2 = tmp_samplename.replace(".root","")
            samplename_list.append(tmp_samplename2)

        if Add_Signal:
            samplename_list.append(self.sig_samplename)
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

   # get the corresponding root files for the background and signal sample names
    def get_root_files(self,Add_Signal=True):
        list_dirs = os.listdir(self.dir_back_input)
        root_file = dict()

        for dirname in list_dirs:
            tmp_samplename = dirname.split("WPiGammaAnalysis_")[1]
            tmp_samplename2 = tmp_samplename.replace(".root","")
            root_file[tmp_samplename2] = ROOT.TFile(self.dir_back_input + dirname)

        if Add_Signal:
            root_file[self.sig_samplename] = ROOT.TFile(self.sig_filename)
        return root_file


