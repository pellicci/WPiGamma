import ROOT
import os

dir_input = "crab_projects/samples/"
dir_output_bkg = "rootfiles/backgrounds/"
dir_output_sig = "rootfiles/signals/"

list_dirs = os.listdir(dir_input)

if not os.path.exists("rootfiles"):
    os.makedirs("rootfiles")

if not os.path.exists(dir_output_bkg):
    os.makedirs(dir_output_bkg)

if not os.path.exists(dir_output_sig):
    os.makedirs(dir_output_sig)

for dirname in list_dirs:

    print "Processing sample dir " + dirname
    crab_command = "crab getoutput -d " + dir_input + dirname
    os.system(crab_command)

    samplename = dirname.split("crab_WPiGammaAnalysis_")

    if "Signal" in dirname:
        hadd_command = "hadd -f " + dir_output_sig + "/WPiGammaAnalysis_" + samplename[1] + ".root " + dir_input + dirname + "/results/*.root"
    else:
        hadd_command = "hadd -f " + dir_output_bkg + "/WPiGammaAnalysis_" + samplename[1] + ".root " + dir_input + dirname + "/results/*.root"

    os.system(hadd_command)

list_signals = os.listdir(dir_output_sig)
if len(list_signals) > 1:
    hadd_command = "hadd -f " + dir_output_sig + "/WPiGammaAnalysis_Signal.root " + dir_output_sig + "/WPiGammaAnalysis_Signal_*.root "
    rm_command = "rm -rf " + dir_output_sig + "/WPiGammaAnalysis_Signal_*.root "

    os.system(hadd_command)
    os.system(rm_command)

print "All done!"
