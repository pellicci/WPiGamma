import os

###All normalizations are provided to 1fb-1 of lumi in these tables

dir_input = "crab_projects/samples/"
list_dirs = os.listdir(dir_input)

if not os.path.exists("rootfiles"):
    os.makedirs("rootfiles")

output_filename = "rootfiles/Normalizations_table.txt"

##These are in pb
def get_xsec_fromsample(samplename):
    
    if samplename == "ttbar":
        return 831.76

    if samplename == "ttbarlnu":
        return 87.31

    if samplename == "ttbarWQQ":
        return 0.4062

    if samplename == "ttbarWlnu":
        return 0.2043
                  
    if samplename == "ttbarZQQ":
        return 0.5297 

    if samplename == "SingleTop_tW":
        return 35.85

    if samplename == "SingleAntiTop_tW":
        return 35.85

    if samplename == "WJetsToLNu":
        return 20508.9*3.

    if samplename == "DY_10_50":
        return 18610.0

    if samplename == "DY_50":
        return 5765.0

    if samplename == "QCD_HT100to200":
        return 27540000.0 

    if samplename == "QCD_HT200to300_1":
        return 1717000.0

    if samplename == "QCD_HT200to300_2":
        return 1717000.0

    if samplename == "QCD_HT300to500_1":
        return 351300.0

    if samplename == "QCD_HT300to500_2":
        return 351300.0

    if samplename == "QCD_HT500to700_1":
        return 31630.0

    if samplename == "QCD_HT500to700_2":
        return 31630.0

    if samplename == "QCD_HT700to1000_1":
        return 6802.0

    if samplename == "QCD_HT700to1000_2":
        return 6802.0

    if samplename == "QCD_HT1000to1500_1":
        return 1206.0

    if samplename == "QCD_HT1000to1500_2":
        return 1206.0

    if samplename == "QCD_HT1500to2000_1":
        return 120.4 

    if samplename == "QCD_HT1500to2000_2":
        return 120.4

    if samplename == "QCD_HT2000toInf_1":
        return 25.25

    if samplename == "QCD_HT2000toInf_2":
        return 25.25

    if samplename == "ZZ":
        return 8.16

    if samplename == "WW":
        return 51.723

    if samplename == "WZ":
        return 47.13

    if "Signal" in samplename:
        return 186444.*0.0000001       #cross section taken from 1603.09222, BR assumed 1*10-7

##Now starts the program

out_file = open(output_filename,"w")

for dirname in list_dirs:
    samplename = dirname.split("crab_WPiGammaAnalysis_")[1]
    print "Processing sample dir " + dirname
    crab_command = "crab report -d " + dir_input + dirname + " | grep read"
    print crab_command

    event_string = os.popen(crab_command).read()
    print "event string: ", event_string
    number_events = float((event_string.split())[4])
    print "No. of events processed = " + str (number_events) + "\n"
    xsection = float(get_xsec_fromsample(samplename))
    print "crossection = ", xsection
    if number_events == 0:
        scale_factor = 0.
        print "NUMBER OF EVENTS RETRIEVED = 0. SCALE FACTOR SET TO 0"
    else:
        scale_factor = float(xsection*1000./number_events)
        print "scale_factor = ", scale_factor
        write_string = samplename + " " + str(scale_factor) + "\n"
        print "Output Norm = ", write_string
        out_file.write(write_string)

print "All done!"
