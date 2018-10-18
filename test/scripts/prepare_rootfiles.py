import ROOT
import os
import subprocess
import argparse


#---------------------------------#
p = argparse.ArgumentParser(description='Select whether to download MC or data')
p.add_argument('isData_option', help='Type <<MC>> or <<data>>')
args = p.parse_args()

# Switch from muon to electron channel
if args.isData_option == "MC":
    isData = False
if args.isData_option == "data":
    isData = True
#---------------------------------#

if not isData:
    dir_input = "crab_projects/samples_Medium/"
    dir_output_bkg = "rootfiles/Medium/backgrounds/"
    dir_output_sig = "rootfiles/Medium/signals/"

if isData:
    dir_input = "crab_projects/dataprocess/"
    dir_output_data = "rootfiles/data/"


list_dirs = os.listdir(dir_input)

complementary_samples_list = ["ttbarWlnu","ttbarZlnu","DY_10_50","DY_50","QCD_HT200to300","QCD_HT300to500","QCD_HT500to700","QCD_HT700to1000","QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf","WZ","WGToLNuG","TTGJets","ZGTo2LG"]


if not os.path.exists("rootfiles"):
    os.makedirs("rootfiles")

if not isData and not os.path.exists(dir_output_bkg):
    os.makedirs(dir_output_bkg)

if not isData and not os.path.exists(dir_output_sig):
    os.makedirs(dir_output_sig)

if isData and not os.path.exists(dir_output_data):
    os.makedirs(dir_output_data)

for dirname in list_dirs:

    print "Processing sample dir " + dirname
    
    n_jobs_command = "crab status -d " + dir_input + dirname + " | grep status: " + "| awk " + """'{split($0,array,"/") ; print array[2]}'""" + "| sed 's/.$//'"
    n_jobs = subprocess.check_output(n_jobs_command, shell=True)

    print "Number of jobs to be retrieved: ", n_jobs

    if n_jobs <= 500:
        crab_command = "crab getoutput -d " + dir_input + dirname
        os.system(crab_command)
    elif (n_jobs > 500 and n_jobs <= 1000):
        crab_command = "crab getoutput -d " + dir_input + dirname + " --jobids 1-500"
        os.system(crab_command)
        crab_command_1 = "crab getoutput -d " + dir_input + dirname + " --jobids 501-" + n_jobs
        os.system(crab_command_1)
    else:
        crab_command = "crab getoutput -d " + dir_input + dirname + " --jobids 1-500"
        os.system(crab_command)
        crab_command_1 = "crab getoutput -d " + dir_input + dirname + " --jobids 501-1000"
        os.system(crab_command_1)
        crab_command_2 = "crab getoutput -d " + dir_input + dirname + " --jobids 1001-" + n_jobs
        os.system(crab_command_2)

    samplename = dirname.split("crab_WPiGammaAnalysis_") #--which means "dirname"-"crab_WPiGammaAnalysys_"

    if "Signal" in dirname:
        hadd_command = "hadd -f " + dir_output_sig + "/WPiGammaAnalysis_" + samplename[1] + ".root " + dir_input + dirname + "/results/*.root"
    elif isData:
        hadd_command = "hadd -f " + dir_output_data + "/WPiGammaAnalysis_" + samplename[1] + ".root " + dir_input + dirname + "/results/*.root"
    else:
        hadd_command = "hadd -f " + dir_output_bkg + "/WPiGammaAnalysis_" + samplename[1] + ".root " + dir_input + dirname + "/results/*.root"

    os.system(hadd_command)


# Now add samples with different names but same xsec
if not isData:
    list_signals = os.listdir(dir_output_sig)
    if len(list_signals) > 1:
        hadd_command = "hadd -f " + dir_output_sig + "/WPiGammaAnalysis_Signal.root " + dir_output_sig + "/WPiGammaAnalysis_Signal_*.root "
        rm_command = "rm -rf " + dir_output_sig + "/WPiGammaAnalysis_Signal_*.root "

        os.system(hadd_command)
        os.system(rm_command)

    for sample in complementary_samples_list:
        hadd_command = "hadd -f " + dir_output_bkg + "/WPiGammaAnalysis_" + sample + ".root " + dir_output_bkg + "/WPiGammaAnalysis_" + sample + "_*.root "
        rm_command = "rm -rf " + dir_output_bkg + "/WPiGammaAnalysis_" + sample + "_*.root "

        os.system(hadd_command)
        os.system(rm_command)        


if isData:
    list_data = os.listdir(dir_output_data)
    if len(list_data) > 1:
        hadd_command = "hadd -f " + dir_output_data + "/WPiGammaAnalysis_Data.root " + dir_output_data + "/WPiGammaAnalysis_*.root "
        rm_command = "rm -rf " + dir_output_data + "/WPiGammaAnalysis_Single*.root "

        os.system(hadd_command)
        os.system(rm_command)

print "All done!"
