import ROOT
import os
import subprocess
import argparse


#---------------------------------#
p = argparse.ArgumentParser(description='Select whether to download MC or data')
p.add_argument('isData_option', help='Type <<MC>> or <<data>>')
p.add_argument('year_option', help='Type <<2016>>, <<2017>> or <<2018>>')
p.add_argument('do_option', help='Type <<doDownload>>, <<doHadd>> or <<doBoth>>')
args = p.parse_args()

# Switch from muon to electron channel
if args.isData_option == "MC":
    isData = False
if args.isData_option == "data":
    isData = True

doDownload = False
doHadd = False
doBoth = False
if args.do_option == "doDownload":
    doDownload = True
if args.do_option == "doHadd":
    doHadd = True
if args.do_option == "doBoth":
    doBoth = True

year = args.year_option
#---------------------------------#

if not isData:
    dir_input = "crab_projects/samples_MC_" + year + "/"
    dir_output_bkg = "rootfiles/latest_production/MC/backgrounds/"
    dir_output_sig = "rootfiles/latest_production/MC/signals/" 
    #dir_output_bkg = "rootfiles/medium_production/MC/backgrounds/"
    #dir_output_sig = "rootfiles/medium_production/MC/signals/" 
    #dir_output_bkg = "rootfiles/bTagEfficiency_medium/MC/backgrounds/"
    #dir_output_sig = "rootfiles/bTagEfficiency_medium/MC/signals/" 

if isData:
    dir_input = "crab_projects/samples_data_" + year + "/"
    dir_output_data = "rootfiles/latest_production/dataprocess/"
    #dir_output_data = "rootfiles/medium_production/dataprocess/"


list_dirs = os.listdir(dir_input)

complementary_samples_list_2016 = ["ttbarWlnu","ttbarZlnu","WJetsToLNu","DY10to50","QCDHT200to300","QCDHT300to500","QCDHT500to700","QCDHT700to1000","QCDHT1000to1500","QCDHT1500to2000","QCDHT2000toInf","WZ","TTGJets"]

complementary_samples_list_2017 = ["WJetsToLNu1J","WJetsToLNu2J","DY50","TTGJets"]

complementary_samples_list_2018 = ["ttbarZQQ","DY10to50","DY50"]

if year == "2016":
    complementary_samples_list = complementary_samples_list_2016
if year == "2017":
    complementary_samples_list = complementary_samples_list_2017
if year == "2018":
    complementary_samples_list = complementary_samples_list_2018


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

    if doDownload or doBoth:
    
        n_jobs_command = "crab status -d " + dir_input + dirname + " | grep status: " + "| awk " + """'{split($0,array,"/") ; print array[2]}'""" + "| sed 's/.$//'"
        n_jobs = int(subprocess.check_output(n_jobs_command, shell=True))

        print "Number of jobs to be retrieved: ", n_jobs

        if n_jobs <= 500:
            crab_command = "crab getoutput -d " + dir_input + dirname
            os.system(crab_command)
        elif (n_jobs > 500 and n_jobs <= 1000):
            crab_command = "crab getoutput -d " + dir_input + dirname + " --jobids 1-500"
            os.system(crab_command)
            crab_command_1 = "crab getoutput -d " + dir_input + dirname + " --jobids 501-" + str(n_jobs) # Because it is impossible to concatenate str and int objects
            os.system(crab_command_1)
        elif (n_jobs > 1000 and n_jobs <= 1500):
            crab_command = "crab getoutput -d " + dir_input + dirname + " --jobids 1-500"
            os.system(crab_command)
            crab_command_1 = "crab getoutput -d " + dir_input + dirname + " --jobids 501-1000"
            os.system(crab_command_1)
            crab_command_2 = "crab getoutput -d " + dir_input + dirname + " --jobids 1001" + str(n_jobs) 
            os.system(crab_command_2)
        elif (n_jobs > 1500):
            crab_command = "crab getoutput -d " + dir_input + dirname + " --jobids 1-500"
            os.system(crab_command)
            crab_command_1 = "crab getoutput -d " + dir_input + dirname + " --jobids 501-1000"
            os.system(crab_command_1)
            crab_command_2 = "crab getoutput -d " + dir_input + dirname + " --jobids 1001-1500" 
            os.system(crab_command_2)
            crab_command_3 = "crab getoutput -d " + dir_input + dirname + " --jobids 1501-" + str(n_jobs) # Because it is impossible to concatenate str and int objects
            os.system(crab_command_3)

    if doHadd or doBoth:

        samplename = dirname.split("crab_" + year + "_WPiGammaAnalysis_") #--which means "dirname"-"crab_2017_WPiGammaAnalysys_"

        if "Signal" in dirname:
            hadd_command = "hadd -f " + dir_output_sig + "WPiGammaAnalysis_" + samplename[1] + "_" + year + ".root " + dir_input + dirname + "/results/*.root"
        elif isData:
            hadd_command = "hadd -f " + dir_output_data + "WPiGammaAnalysis_" + samplename[1] + "_" + year + ".root " + dir_input + dirname + "/results/*.root"
        else:
            hadd_command = "hadd -f " + dir_output_bkg + "WPiGammaAnalysis_" + samplename[1] + "_" + year + ".root " + dir_input + dirname + "/results/*.root"
            
        os.system(hadd_command)

if doHadd or doBoth:
    # Now add samples with different names but same xsec
    if not isData:
        list_signals = os.listdir(dir_output_sig)
        if len(list_signals) > 1:
            hadd_command = "hadd -f " + dir_output_sig + "WPiGammaAnalysis_Signal_" + year + ".root " + dir_output_sig + "WPiGammaAnalysis_Signal_*_" + year + ".root "
            rm_command = "rm -rf " + dir_output_sig + "WPiGammaAnalysis_Signal_*" + "_" + year + ".root "
            
            os.system(hadd_command)
            os.system(rm_command)

        for sample in complementary_samples_list:
            hadd_command = "hadd -f " + dir_output_bkg + "WPiGammaAnalysis_" + sample + "_" + year + ".root " + dir_output_bkg + "WPiGammaAnalysis_" + sample + "_?_" + year + ".root "
            rm_command = "rm -rf " + dir_output_bkg + "WPiGammaAnalysis_" + sample + "_?_" + year + ".root "
                
            os.system(hadd_command)
            os.system(rm_command)        
            
            
    if isData:
        list_data = os.listdir(dir_output_data)
        if len(list_data) > 1:
            hadd_command = "hadd -f " + dir_output_data + "WPiGammaAnalysis_Data_" + year + ".root " + dir_output_data + "WPiGammaAnalysis_*_" + year + ".root "
            rm_command = "rm -rf " + dir_output_data + "WPiGammaAnalysis_Single*.root "
            
            os.system(hadd_command)
            os.system(rm_command)

if doDownload:
    print "Download done!"
if doHadd:
    print "Hadd done!"
if doBoth:
    print "All done!"
