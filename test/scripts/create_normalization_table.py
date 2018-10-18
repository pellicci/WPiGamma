###All normalizations are provided to 1fb-1 of lumi in these tables
import os
import sys

#xsec lookup in pb
secs_table = dict()
secs_table["ttbar"] = 831.76 - 87.31
secs_table["ttbarlnu"] = 87.31
secs_table["ttbarWQQ"] = 0.4062
secs_table["ttbarWlnu"] = 0.2043
secs_table["ttbarZQQ"] = 0.5297
secs_table["ttbarZlnu"] = 0.2529
secs_table["SingleTop_tW"] = 35.85
secs_table["SingleAntiTop_tW"] = 35.85
secs_table["WJetsToLNu"] = 20508.9*3.
secs_table["DY_10_50"] = 18610.0
secs_table["DY_50"] = 1921.8*3.
secs_table["QCD_HT100to200"] = 27540000.0
secs_table["QCD_HT200to300"] = 1717000.0
secs_table["QCD_HT300to500"] = 351300.0
secs_table["QCD_HT500to700"] = 31630.0
secs_table["QCD_HT700to1000"] = 6802.0
secs_table["QCD_HT1000to1500"] = 1206.0
secs_table["QCD_HT1500to2000"] = 120.4
secs_table["QCD_HT2000toInf"] = 25.25
secs_table["ZZ"] = 8.16
secs_table["WW"] = 51.723
secs_table["WZ"] = 47.13
secs_table["GammaJets_20_40"] = 219.2
secs_table["GammaJets_40_Inf"] = 862.4
secs_table["GammaJets_20_Inf"] = 3255.0
secs_table["QCD_DoubleEMEnriched_30to40"] = 22180.0
secs_table["QCD_DoubleEMEnriched_30toInf"] = 247000.0
secs_table["QCD_DoubleEMEnriched_40toInf"] = 113100.0
secs_table["TTGJets"] = 3.795
secs_table["WGToLNuG"] = 489.0
secs_table["ZGTo2LG"] = 117.864
#cross section taken from https://arxiv.org/pdf/1611.04040.pdf, BR assumed 10-6, last factor 2 is because we have two possible final states (one for W+ and one for W-)
secs_table["Signal"] = 831.76*0.1086*2.*0.000001*2.

complementary_samples_list = ["ttbarWlnu","ttbarZlnu","DY_10_50","DY_50","QCD_HT200to300","QCD_HT300to500","QCD_HT500to700","QCD_HT700to1000","QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf","WZ","WGToLNuG","TTGJets","ZGTo2LG","Signal"]

##Now starts the program
def main():

    dir_input = "crab_projects/samples_Medium/"
    list_dirs = os.listdir(dir_input)

    if not os.path.exists("rootfiles"):
        os.makedirs("rootfiles")

    output_filename = "rootfiles/Medium/Normalizations_table.txt"
    out_file = open(output_filename,"w")

    events_cumul = dict()
    for sample in complementary_samples_list:
        events_cumul[sample] = 0.

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
        for name_complement in complementary_samples_list:
            if name_complement in samplename :
                events_cumul[name_complement] = events_cumul[name_complement] + number_events
                print "events_cumul: ", events_cumul[name_complement]
                continue

        events_cumul[samplename] = number_events

    for sample,event_count in events_cumul.iteritems():
        if event_count == 0:
            scale_factor = 0.
            print "NUMBER OF EVENTS RETRIEVED = 0. SCALE FACTOR SET TO 0"
        else:
            xsection = secs_table[sample]
            print "Cross section = ", xsection
            scale_factor = float(xsection*1000./events_cumul[sample])
            print sample + " scale factor = ", scale_factor

        write_string = sample + " " + str(scale_factor) + "\n"
        print "Output Norm = ", write_string
        out_file.write(write_string)        
            
print "All done!"

if __name__ == "__main__":
    main()
