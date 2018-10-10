import os

###All normalizations are provided to 1fb-1 of lumi in these tables


##These are in pb
def get_xsec_fromsample(samplename):
    
    if samplename == "ttbar":
        return 831.76 - 87.31

    if samplename == "ttbarlnu":
        return 87.31

    if samplename == "ttbarWQQ":
        return 0.4062

    if "ttbarWlnu" in samplename:
        return 0.2043
                  
    if samplename == "ttbarZQQ":
        return 0.5297 

    if "ttbarZlnu" in samplename:
        return 0.2529 

    if samplename == "SingleTop_tW":
        return 35.85

    if samplename == "SingleAntiTop_tW":
        return 35.85

    if samplename == "WJetsToLNu":
        return 20508.9*3.

    if "DY_10_50" in samplename:
        return 18610.0

    if "DY_50" in samplename:
        return 1921.8*3.

    if samplename == "QCD_HT100to200":
        return 27540000.0 

    if "QCD_HT200to300" in samplename:
        return 1717000.0

    if "QCD_HT300to500" in samplename:
        return 351300.0

    if "QCD_HT500to700" in samplename:
        return 31630.0

    if "QCD_HT700to1000" in samplename:
        return 6802.0

    if "QCD_HT1000to1500" in samplename:
        return 1206.0

    if "QCD_HT1500to2000" in samplename:
        return 120.4 

    if "QCD_HT2000toInf" in samplename:
        return 25.25

    if samplename == "ZZ":
        return 8.16

    if samplename == "WW":
        return 51.723

    if "WZ" in samplename:
        return 47.13

    if samplename == "GammaJets_20_40":
        return 219.2

    if samplename == "GammaJets_40_Inf":
        return 862.4

    if samplename == "GammaJets_20_Inf":
        return 3255.0

    if samplename == "QCD_DoubleEMEnriched_30to40":
        return 22180.0

    if samplename == "QCD_DoubleEMEnriched_30toInf":
        return 247000.0

    if samplename == "QCD_DoubleEMEnriched_40toInf":
        return 113100.0

    if "TTGJets" in samplename:
        return 3.795
    
    # if "WGToLNuG" in samplename:
    #     return 489.

    # if "ZGTo2LG" in samplename:
    #     return 117.864
    
    if "Signal" in samplename:
        #cross section taken from https://arxiv.org/pdf/1611.04040.pdf, BR assumed 10-6, last factor 2 is because we have two possible final states (one for W+ and one for W-)
        return 831.76*0.1086*2.*0.000001*2.

##Now starts the program

def main():

    dir_input = "crab_projects/samples_Medium/"
    list_dirs = os.listdir(dir_input)

    if not os.path.exists("rootfiles"):
        os.makedirs("rootfiles")
        

    output_filename = "rootfiles/Medium/Normalizations_table.txt"

    out_file = open(output_filename,"w")


    complementary_samples_list = ["ttbarWlnu","ttbarZlnu","DY_10_50","DY_50","QCD_HT200to300","QCD_HT300to500","QCD_HT500to700","QCD_HT700to1000","QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf","WZ","TTGJets"]#"WGToLNuG","ZGTo2LG"

    events_cumul = dict()

    for sample in complementary_samples_list:
        events_cumul[sample] = 0.

    #Counter to treat signal differently: 2 samples, same xsec
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
        
        #Treat differently different samples with same xsec

        if "ttbarWlnu" in samplename:
            events_cumul["ttbarWlnu"] = events_cumul["ttbarWlnu"] + number_events
            print "events_cumul: ", events_cumul["ttbarWlnu"]
            continue

        if "ttbarZlnu" in samplename:
            events_cumul["ttbarZlnu"] = events_cumul["ttbarZlnu"] + number_events
            print "events_cumul: ", events_cumul["ttbarZlnu"]
            continue

        if "DY_10_50" in samplename:
            events_cumul["DY_10_50"] = events_cumul["DY_10_50"] + number_events
            print "events_cumul: ", events_cumul["DY_10_50"]
            continue

        if "DY_50" in samplename:
            events_cumul["DY_50"] = events_cumul["DY_50"] + number_events
            print "events_cumul: ", events_cumul["DY_50"]
            continue

        if "QCD_HT200to300" in samplename:
            events_cumul["QCD_HT200to300"] = events_cumul["QCD_HT200to300"] + number_events
            continue

        if "QCD_HT300to500" in samplename:
            events_cumul["QCD_HT300to500"] = events_cumul["QCD_HT300to500"] + number_events
            continue

        if "QCD_HT500to700" in samplename:
            events_cumul["QCD_HT500to700"] = events_cumul["QCD_HT500to700"] + number_events
            continue

        if "QCD_HT700to1000" in samplename:
            events_cumul["QCD_HT700to1000"] = events_cumul["QCD_HT700to1000"] + number_events
            continue

        if "QCD_HT1000to1500" in samplename:
            events_cumul["QCD_HT1000to1500"] = events_cumul["QCD_HT1000to1500"] + number_events
            continue

        if "QCD_HT1500to2000" in samplename:
            events_cumul["QCD_HT1500to2000"] = events_cumul["QCD_HT1500to2000"] + number_events
            continue

        if "QCD_HT2000toInf" in samplename:
            events_cumul["QCD_HT2000toInf"] = events_cumul["QCD_HT2000toInf"] + number_events
            continue

        if "WZ" in samplename:
            events_cumul["WZ"] = events_cumul["WZ"] + number_events
            continue

        if "WGToLNuG" in samplename:
            events_cumul["WGToLNuG"] = events_cumul["WGToLNuG"] + number_events
            continue

        if "TTGJets" in samplename:
            events_cumul["TTGJets"] = events_cumul["TTGJets"] + number_events
            continue

        if "ZGTo2LG" in samplename:
            events_cumul["ZGTo2LG"] = events_cumul["ZGTo2LG"] + number_events
            continue

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
        print "Signal scale factor = ", scale_factor
        write_string = "Signal" + " " + str(scale_factor) + "\n"
        print "Output Norm = ", write_string
        out_file.write(write_string)

    for sample in complementary_samples_list:
        xsection = float(get_xsec_fromsample(sample))
        scale_factor = float(xsection*1000./events_cumul[sample])
        print sample + " scale factor = ", scale_factor
        write_string = sample + " " + str(scale_factor) + "\n"
        print "Output Norm = ", write_string
        out_file.write(write_string)        

            
    print "All done!"

if __name__ == "__main__":
    main()
