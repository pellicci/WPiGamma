import ROOT
import argparse

p = argparse.ArgumentParser(description='Select whether to fill the histograms after pre-selection or after BDT')
p.add_argument('runningEra_option', help='Type <<0>> for 2016, <<1>> for 2017, <<2>> for 2018, <<3>> for combination 2016+2017')
p.add_argument('inputfile_option', help='Provide input file name')
p.add_argument('outputfile_option', help='Provide output file name')
args = p.parse_args()

runningEra = int(args.runningEra_option)

input_filename = args.inputfile_option
output_filename = args.outputfile_option

h_base  = dict()

list_histos = ["h_mupt", "h_elept", "h_pipt", "h_gammaet", "h_mueta", "h_eleeta","h_pieta","h_gammaeta", "h_nBjets_25","h_deltaphi_mu_pi","h_deltaphi_ele_pi","h_deltaphi_mu_W","h_deltaphi_ele_W","h_deltaeta_mu_pi","h_deltaeta_ele_pi","h_Wmass","h_Wmass_flag_mu","h_Wmass_flag_ele","h_mu_gamma_InvMass","h_ele_gamma_InvMass","h_piIso_05_mu","h_piRelIso_05_mu_ch","h_piRelIso_05_mu","h_piIso_05_ele","h_piRelIso_05_ele_ch","h_piRelIso_05_ele","h_met_mu","h_met_ele","h_met_puppi","h_Wmass_ratio_mu","h_Wmass_ratio_ele","h_nPV_mu","h_nPV_ele","h_deltaphi_mu_gamma","h_deltaphi_ele_gamma"]

