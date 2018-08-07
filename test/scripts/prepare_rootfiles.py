import ROOT
import os
import subprocess

isData = False ##---------switch from DATA to MC and vice versa---------##

if not isData:
    dir_input = "crab_projects/samples_Medium/"
    dir_output_bkg = "rootfiles/Medium/backgrounds/"
    #dir_output_bkg = "rootfiles/Medium_AfterFix/backgrounds/"
    dir_output_sig = "rootfiles/Medium/signals/"
    #dir_output_sig = "rootfiles/Medium_AfterFix/signals/"

if isData:
    dir_input = "crab_projects/dataprocess/"
    dir_output_data = "rootfiles/data/"
    #dir_output_data = "rootfiles/data_AfterFix/"

list_dirs = os.listdir(dir_input)

WGToLNuG_samples = 0
TTGJets_samples = 0

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

    if "WGToLNuG" in dirname:
        WGToLNuG_samples += 1

    if "TTGJets" in dirname:
        TTGJets_samples += 1

if not isData:
    list_signals = os.listdir(dir_output_sig)
    if len(list_signals) > 1:
        hadd_command = "hadd -f " + dir_output_sig + "/WPiGammaAnalysis_Signal.root " + dir_output_sig + "/WPiGammaAnalysis_Signal_*.root "
        rm_command = "rm -rf " + dir_output_sig + "/WPiGammaAnalysis_Signal_*.root "

        os.system(hadd_command)
        os.system(rm_command)

    if WGToLNuG_samples > 1:
        hadd_command = "hadd -f " + dir_output_bkg + "/WPiGammaAnalysis_WGToLNuG.root " + dir_output_bkg + "/WPiGammaAnalysis_WGToLNuG_ext*.root "
        rm_command = "rm -rf " + dir_output_bkg + "/WPiGammaAnalysis_WGToLNuG_ext*.root "

        os.system(hadd_command)
        os.system(rm_command)

    if TTGJets_samples > 1:
        hadd_command = "hadd -f " + dir_output_bkg + "/WPiGammaAnalysis_TTGJets.root " + dir_output_bkg + "/WPiGammaAnalysis_TTGJets_*.root "
        rm_command = "rm -rf " + dir_output_bkg + "/WPiGammaAnalysis_WGToLNuG_*.root "

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
