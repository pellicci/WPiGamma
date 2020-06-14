import ROOT
import os
import math
import numpy as np
import shutil

source = "plots/latest_production/2016_2017_2018/"

source_list = os.listdir(source)

destination = "/afs/cern.ch/user/r/rselvati/www/WPiGamma/InterestingVariables/BDT/09_06_2020/2016_2017_2018/"

for files in source_list:
    fullpath = os.path.join(source, files)
    if files.startswith("h_") and files.endswith("eta.pdf"):
        shutil.copy(fullpath,destination)
    if files.startswith("h_") and files.endswith("pt.pdf"):
        shutil.copy(fullpath,destination)
    if files == "h_gammaet.pdf": 
        shutil.copy(fullpath,destination)
    # if files == "h_ele_gamma_InvMass.pdf": 
    #     shutil.copy(fullpath,destination)
    # if files == "h_deltaphi_ele_gamma.pdf": 
    #     shutil.copy(fullpath,destination)
    # if files == "h_mu_gamma_InvMass.pdf": 
    #     shutil.copy(fullpath,destination)
    # if files == "h_deltaphi_mu_gamma.pdf": 
    #     shutil.copy(fullpath,destination)
    if files.startswith("h_met") and (files.endswith("_mu.pdf") or files.endswith("_ele.pdf")):
        shutil.copy(fullpath,destination)
    if files.startswith("h_nPV") and (files.endswith("_mu.pdf") or files.endswith("_ele.pdf")):
        shutil.copy(fullpath,destination)
    if files == "h_Wmass.pdf": 
        shutil.copy(fullpath,destination)
    if files == "h_Wmass_flag_mu.pdf": 
        shutil.copy(fullpath,destination)
    if files == "h_Wmass_flag_ele.pdf": 
        shutil.copy(fullpath,destination)
    if files.startswith("Wmass_ratio_") and files.endswith(".pdf"):
        shutil.copy(fullpath,destination)
    if files == "h_piRelIso_05_mu_ch.pdf": 
        shutil.copy(fullpath,destination)
    if files == "h_piRelIso_05_ele_ch.pdf": 
        shutil.copy(fullpath,destination)
    if files == "h_piRelIso_05_ch.pdf": 
        shutil.copy(fullpath,destination)
    if files == "h_deltaR_mu_gamma.pdf": 
        shutil.copy(fullpath,destination)
    if files == "h_deltaR_ele_gamma.pdf": 
        shutil.copy(fullpath,destination)
    if files == "h_nBjets_25.pdf": 
        shutil.copy(fullpath,destination)
    if files == "h_mu_met_mT.pdf": 
        shutil.copy(fullpath,destination)
    if files == "h_ele_met_mT.pdf": 
        shutil.copy(fullpath,destination)

print "All histos copied!"

