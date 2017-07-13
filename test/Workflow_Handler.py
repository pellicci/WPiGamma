import ROOT
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


