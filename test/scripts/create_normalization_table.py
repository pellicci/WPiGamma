import os

###All normalizations are provided to 1fb-1 of lumi in these tables

dir_input = "crab_projects/samples_Medium/"
list_dirs = os.listdir(dir_input)

if not os.path.exists("rootfiles"):
    os.makedirs("rootfiles")

#output_filename = "rootfiles/Tight/Normalizations_table.txt"
output_filename = "rootfiles/Medium/Normalizations_table.txt"

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

    if samplename == "ttbarZlnu":
        return 0.2529 

    if samplename == "SingleTop_tW":
        return 35.85

    if samplename == "SingleAntiTop_tW":
        return 35.85

    if samplename == "WJetsToLNu":
        return 20508.9*3.

    if samplename == "DY_10_50":
        return 18610.0

    if samplename == "DY_50":
        return 1921.8*3.

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

    if samplename == "GammaJets_20_40":
        return 137751.

    if samplename == "GammaJets_40_Inf":
        return 16792.

    if samplename == "GammaJets_20_Inf":
        return 154500.

    #if samplename == "WGToLNuG":
    #    return 489.
    
    if "Signal" in samplename:
        #cross section taken from https://arxiv.org/pdf/1611.04040.pdf, BR assumed 10-6, last factor 2 is because we have two possible final states (one for W+ and one for W-)
        return 831.76*0.1086*2.*0.000001*2.

##Now starts the program

out_file = open(output_filename,"w")

signal_events_cumul = 0.

for dirname in list_dirs:

    samplename = dirname.split("crab_WPiGammaAnalysis_")[1]
    print "Processing sample dir " + dirname
    crab_command = "crab report -d " + dir_input + dirname + " | grep read"
    print crab_command

    event_string = os.popen(crab_command).read()
    print "event string: ", event_string
    number_events = float((event_string.split())[4])
    print "No. of events processed = " + str (number_events) + "\n"
    
    #Treat signal differently to account for different samples with same xsec
    if "Signal" in samplename:
        signal_events_cumul = signal_events_cumul + number_events
        continue

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

if signal_events_cumul > 0.:
    xsection = float(get_xsec_fromsample("Signal"))
    scale_factor = float(xsection*1000./signal_events_cumul)
    print "Signal scale_factor = ", scale_factor
    write_string = "Signal" + " " + str(scale_factor) + "\n"
    print "Output Norm = ", write_string
    out_file.write(write_string)

print "All done!"
