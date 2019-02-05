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

if isData:
    dir_input = "crab_projects/dataprocess/"
else:
    dir_input = "crab_projects/samples_Medium/"

list_dirs = os.listdir(dir_input)

for dirname in list_dirs:

    print "Verifying sample dir " + dirname
    
    #n_jobs_total_command  = "crab status -d " + dir_input + dirname + " | grep status: " + "| awk " + """'{split($0,array,"/") ; print array[2]}'""" + "| sed 's/.$//'"
    n_jobs_finished_command = "crab status -d " + dir_input + dirname + " | grep finished " + "| awk " + """'{split($0,array,"/") ; print array[1]}'""" + "| awk " + """'{split($0,array,"(") ; print array[2]}'"""
    n_rootfiles_in_dir_command = "ls " + dir_input + dirname + "/results/*.root | wc -l"

    n_jobs_finished = int(subprocess.check_output(n_jobs_finished_command, shell=True))

    n_rootfiles_in_dir = int(subprocess.check_output(n_rootfiles_in_dir_command, shell=True))


    if not n_jobs_finished == n_rootfiles_in_dir:
        print "!!!! N JOBS IN STATUS FINISHED != N ROOTFILES IN DIRECTORY !!!!"
        print "n jobs finished: ",  n_jobs_finished
        print "n rootfiles: ",  n_rootfiles_in_dir

    print "/////////////////"
