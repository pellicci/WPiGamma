###All normalizations are provided to 1fb-1 of lumi in these tables
import os
import sys
import argparse

#---------------------------------#
p = argparse.ArgumentParser(description='Select whether to download MC or data')
p.add_argument('year_option', help='Type <<2016>> or <<2017>>')
args = p.parse_args()

year = args.year_option
#---------------------------------#

#######################################
#                                     #
#--------------- 2016 ----------------#
#                                     #
#######################################

#xsec lookup in pb
secs_table_2016 = dict()
secs_table_2016["ttbar"] = 831.76 - 88.29
secs_table_2016["ttbarlnu"] = 88.29 #Controllare questo conto
secs_table_2016["ttbarWQQ"] = 0.405
secs_table_2016["ttbarWlnu"] = 0.2001
secs_table_2016["ttbarZQQ"] = 0.5297
secs_table_2016["ttbarZlnu"] = 0.2529
secs_table_2016["SingleTop_tW"] = 35.09
secs_table_2016["SingleAntiTop_tW"] = 35.09
secs_table_2016["WJetsToLNu"] = 60430.0
secs_table_2016["DY_10_50"] = 18810.0
secs_table_2016["DY_50"] = 4963.0
secs_table_2016["QCD_HT100to200"] = 27540000.0
secs_table_2016["QCD_HT200to300"] = 1717000.0
secs_table_2016["QCD_HT300to500"] = 351300.0
secs_table_2016["QCD_HT500to700"] = 31630.0
secs_table_2016["QCD_HT700to1000"] = 6802.0
secs_table_2016["QCD_HT1000to1500"] = 1206.0
secs_table_2016["QCD_HT1500to2000"] = 120.4
secs_table_2016["QCD_HT2000toInf"] = 25.25
#secs_table_2016["ZZ"] = 6.912
secs_table_2016["WW"] = 45.2
secs_table_2016["WZ"] = 23.43
secs_table_2016["GammaJets_20_40"] = 219.2
secs_table_2016["GammaJets_40_Inf"] = 862.4
secs_table_2016["GammaJets_20_Inf"] = 3255.0
secs_table_2016["QCD_DoubleEMEnriched_30to40"] = 22180.0
secs_table_2016["QCD_DoubleEMEnriched_30toInf"] = 247000.0
secs_table_2016["QCD_DoubleEMEnriched_40toInf"] = 113100.0
secs_table_2016["WGToLNuG"] = 510.6
secs_table_2016["TTGJets"] = 3.795
secs_table_2016["ZGTo2LG"] = 123.8
secs_table_2016["Signal"] = 831.76*0.1086*2.*0.000001*2. #cross section taken from https://arxiv.org/pdf/1611.04040.pdf, BR assumed 10-6, last factor 2 is because we have two possible final states (one for W+ and one for W-)

#fraction of negative-weighted events in NLO samples (2016)
frac_table_2016 = dict()
frac_table_2016["ttbar"] = 0.
frac_table_2016["ttbarlnu"] = 0.
frac_table_2016["ttbarWQQ"] = 0.2426
frac_table_2016["ttbarWlnu"] = 0.2419
frac_table_2016["ttbarZQQ"] = 0.2656
frac_table_2016["ttbarZlnu"] = 0.2687
frac_table_2016["SingleTop_tW"] = 0.
frac_table_2016["SingleAntiTop_tW"] = 0.
frac_table_2016["WJetsToLNu"] = 0.1581
frac_table_2016["DY_10_50"] = 0.1367
frac_table_2016["DY_50"] = 0.
frac_table_2016["QCD_HT100to200"] = 0.
frac_table_2016["QCD_HT200to300"] = 0.
frac_table_2016["QCD_HT300to500"] = 0.
frac_table_2016["QCD_HT500to700"] = 0.
frac_table_2016["QCD_HT700to1000"] = 0.
frac_table_2016["QCD_HT1000to1500"] = 0.
frac_table_2016["QCD_HT1500to2000"] = 0.
frac_table_2016["QCD_HT2000toInf"] = 0.
#frac_table_2016["ZZ"] = 0.1894
frac_table_2016["WW"] = 0.
frac_table_2016["WZ"] = 0.
frac_table_2016["GammaJets_20_40"] = 0.
frac_table_2016["GammaJets_40_Inf"] = 0.
frac_table_2016["GammaJets_20_Inf"] = 0.
frac_table_2016["QCD_DoubleEMEnriched_30to40"] = 0.
frac_table_2016["QCD_DoubleEMEnriched_30toInf"] = 0.
frac_table_2016["QCD_DoubleEMEnriched_40toInf"] = 0.
frac_table_2016["WGToLNuG"] = 0.1799
frac_table_2016["TTGJets"] = 0.3381
frac_table_2016["ZGTo2LG"] = 0.1584
frac_table_2016["Signal"] = 0.


#######################################
#                                     #
#--------------- 2017 ----------------#
#                                     #
#######################################


secs_table_2017 = dict()
secs_table_2017["ttbar_ToHadronic"] = 377.96
secs_table_2017["ttbar_ToSemiLeptonic"] = 365.34 # accounting for the 2 possible charge signs of the W
secs_table_2017["ttbarlnu"] = 88.29 #NNLO-2017
secs_table_2017["ttbarWQQ"] = 0.4316
secs_table_2017["ttbarWlnu"] = 0.2149
secs_table_2017["ttbarZQQ"] = 0.5104
secs_table_2017["ttbarZlnu"] = 0.2432
secs_table_2017["SingleTop_tW"] = 34.91
secs_table_2017["SingleAntiTop_tW"] = 34.91
secs_table_2017["WJetsToLNu"] = 60430.0 #oppure usare il valore LO che e' intorno a 52k?
secs_table_2017["DY_10_50"] = 18810.0
#secs_table_2017["DY_50"] = 4963.0 #LO 2017
secs_table_2017["DY_50"] = 6529.0 #amcatnlo 2017
secs_table_2017["QCD_HT100to200"] = 23700000.0
secs_table_2017["QCD_HT200to300"] = 1547000.0
secs_table_2017["QCD_HT300to500"] = 322600.0
secs_table_2017["QCD_HT500to700"] = 29980.0
secs_table_2017["QCD_HT700to1000"] = 6334.0
secs_table_2017["QCD_HT1000to1500"] = 1088.0
secs_table_2017["QCD_HT1500to2000"] = 99.11
secs_table_2017["QCD_HT2000toInf"] = 20.23
#secs_table_2017["ZZ"] = 6.912
secs_table_2017["WW"] = 47.73
secs_table_2017["WZ"] = 27.6
secs_table_2017["GammaJets_20_40"] = 232.8
secs_table_2017["GammaJets_40_Inf"] = 872.8
secs_table_2017["GammaJets_20_Inf"] = 3164.0
secs_table_2017["QCD_DoubleEMEnriched_30to40"] = 24750.0
secs_table_2017["QCD_DoubleEMEnriched_30toInf"] = 242700.0
secs_table_2017["QCD_DoubleEMEnriched_40toInf"] = 117400.0
secs_table_2017["WGToLNuG"] = 510.6
secs_table_2017["TTGJets"] = 4.078
secs_table_2017["ZGTo2LG"] = 123.8
secs_table_2017["Signal"] = 831.76*0.1086*2.*0.000001*2. #cross section taken from https://arxiv.org/pdf/1611.04040.pdf, BR assumed 10-6, last factor 2 is because we have two possible final states (one for W+ and one for W-)

#fraction of negative-weighted events in NLO samples (2017)
frac_table_2017 = dict()
frac_table_2017["ttbar_ToHadronic"] = 0.
frac_table_2017["ttbar_ToSemiLeptonic"] = 0.
frac_table_2017["ttbarlnu"] = 0.
frac_table_2017["ttbarWQQ"] = 0.228
frac_table_2017["ttbarWlnu"] = 0.2268
frac_table_2017["ttbarZQQ"] = 0.2645
frac_table_2017["ttbarZlnu"] = 0.2652
frac_table_2017["SingleTop_tW"] = 0.003758
frac_table_2017["SingleAntiTop_tW"] = 0.0034
frac_table_2017["WJetsToLNu"] = 0.0004079
frac_table_2017["DY_10_50"] = 0.
frac_table_2017["DY_50"] = 0.1624
frac_table_2017["QCD_HT100to200"] = 0.
frac_table_2017["QCD_HT200to300"] = 0.000598
frac_table_2017["QCD_HT300to500"] = 0.0009162
frac_table_2017["QCD_HT500to700"] = 0.001485
frac_table_2017["QCD_HT700to1000"] = 0.002061
frac_table_2017["QCD_HT1000to1500"] = 0.003427
frac_table_2017["QCD_HT1500to2000"] = 0.005569
frac_table_2017["QCD_HT2000toInf"] = 0.009878
#frac_table_2017["ZZ"] = 0.1894
frac_table_2017["WW"] = 0.001755
frac_table_2017["WZ"] = 0.
frac_table_2017["GammaJets_20_40"] = 0.
frac_table_2017["GammaJets_40_Inf"] = 0.
frac_table_2017["GammaJets_20_Inf"] = 0.
frac_table_2017["QCD_DoubleEMEnriched_30to40"] = 0.
frac_table_2017["QCD_DoubleEMEnriched_30toInf"] = 0.
frac_table_2017["QCD_DoubleEMEnriched_40toInf"] = 0.
frac_table_2017["WGToLNuG"] = 0. #The sample is not in XSDB, for some reason, but it is madgraph and not amcatnlo like it was in 2016
frac_table_2017["TTGJets"] = 0.3049
frac_table_2017["ZGTo2LG"] = 0.1584
frac_table_2017["Signal"] = 0.



secs_table = dict()
frac_table = dict()

complementary_samples_list_2016 = ["ttbarWlnu","ttbarZlnu","DY_10_50","DY_50","QCD_HT200to300","QCD_HT300to500","QCD_HT500to700","QCD_HT700to1000","QCD_HT1000to1500","QCD_HT1500to2000","QCD_HT2000toInf","WZ","WGToLNuG","TTGJets","ZGTo2LG","Signal"]

complementary_samples_list_2017 = ["WJetsToLNu","DY_50","TTGJets","Signal"]

if year == "2016":
    complementary_samples_list = complementary_samples_list_2016
    secs_table = secs_table_2016
    frac_table = frac_table_2016

if year == "2017":
    complementary_samples_list = complementary_samples_list_2017
    secs_table = secs_table_2017
    frac_table = frac_table_2017

##Now starts the program
def main():

    dir_input = "crab_projects/samples_MC_" + year + "/"
    list_dirs = os.listdir(dir_input)

    if not os.path.exists("rootfiles"):
        os.makedirs("rootfiles")

    output_filename = "rootfiles/latest_production/MC/" + year + "/"
    out_file = open(output_filename,"w")

    events_cumul = dict()
    for sample in complementary_samples_list:
        events_cumul[sample] = 0.

    for dirname in list_dirs:

        is_in_complementary_sample_list = False
        
        samplename = dirname.split("crab_WPiGammaAnalysis_")[1]
        print "Processing sample dir " + dirname
        crab_command = "crab report -d " + dir_input + dirname + " | grep read"
        print crab_command

        event_string = os.popen(crab_command).read()
        print "event string: ", event_string
        number_events = float((event_string.split())[4])
        print "No. of events processed = " + str(number_events) + "\n"
        # number_events = number_events*(1-frac_table[samplename])
        # print "No. of events with positive weights = " + str(number_events) + "\n"
        
        #Treat differently different samples with same xsec
        for name_complement in complementary_samples_list:
            if name_complement in samplename :
                events_cumul[name_complement] = events_cumul[name_complement] + number_events*(1-frac_table[name_complement])
                print "events_cumul: ", events_cumul[name_complement]
                is_in_complementary_sample_list = True
                continue
        if not is_in_complementary_sample_list:
            events_cumul[samplename] = number_events*(1-frac_table[samplename])

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
