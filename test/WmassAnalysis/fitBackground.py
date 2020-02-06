#This program fits both lepton samples at the same time

import ROOT
import math
import argparse

ROOT.gROOT.ProcessLineSync(".L dCB/RooDoubleCBFast.cc+")

#---------------------------------#
p = argparse.ArgumentParser(description='Select whether to fit muon or electron channels, and for which year')
p.add_argument('isMuon_option', help='Type <<muon>> or <<electron>>')
p.add_argument('runningEra_option', help='Type <<0>> for 2016, <<1>> for 2017, <<2>> for 2018, <<3>> for summing 2016+2017+2018')
args = p.parse_args()

# Switch from muon to electron channel
isMuon = False

if args.isMuon_option == "muon":
    isMuon = True

runningEra = int(args.runningEra_option)

useChebychev = True

#---------------------------------#

################################################################
#                                                              #
#------------------------ Instructions ------------------------#
#                                                              #
################################################################


################################################################
#                                                              #
#------------------- Define the observable --------------------#
#                                                              #
################################################################

Wmass = ROOT.RooRealVar("Wmass","m_{#pi#gamma}",50.,100.,"GeV/c^{2}")
Wmass.setRange("LowSideband",55.,65.)
Wmass.setRange("HighSideband",90.,95.)

################################################################
#                                                              #
#--------------------- Retrieve the sample --------------------#
#                                                              #
################################################################

#fInput = ROOT.TFile("Tree_input_massfit_Data_" + str(runningEra) + ".root")
fInput = ROOT.TFile("Tree_input_massfit_Data_" + str(runningEra) + ".root")

fInput.cd()

mytree = fInput.Get("minitree")

################################################################
#                                                              #
#-------------------- Define the categories -------------------#
#                                                              #
################################################################

#Define the mu/ele category for 2016 and 2017
Categorization = ROOT.RooCategory("Categorization","Categorization") #Signal categories are required to compute the number of events in the signal regions, so that the restricted CR can contain similar statistic
Categorization.defineType("MuonCR",0)
Categorization.defineType("MuonSignal",1)
Categorization.defineType("ElectronCR",2)
Categorization.defineType("ElectronSignal",3)

################################################################
#                                                              #
#----------------------- Get BDT output -----------------------#
#                                                              #
################################################################

#Support variables
BDT_out = ROOT.RooRealVar("BDT_out","Output of BDT",-1.,1.)

#Create the RooDataSet. No need to import weight for signal only analysis
data_initial = ROOT.RooDataSet("data","data", ROOT.RooArgSet(Wmass,Categorization,BDT_out), ROOT.RooFit.Import(mytree))

#Check the number of events contained in the signal regions before reducing the dataset
if runningEra == 0 and isMuon:
    print "number of events mu 2016  - SR: ", data_initial.sumEntries("Categorization==1")
    data = data_initial.reduce("(BDT_out > 0.135 && Categorization==Categorization::MuonCR)")
    print "number of events mu 2016  - CR: ", data.sumEntries("Categorization==0")
if runningEra == 0 and not isMuon:
    print "number of events ele 2016 - SR: ", data_initial.sumEntries("Categorization==3")
    data = data_initial.reduce("(BDT_out > 0.052 && Categorization==Categorization::ElectronCR)")
    print "number of events ele 2016 - CR: ", data.sumEntries("Categorization==2")
if runningEra == 1 and isMuon:
    print "number of events mu 2017  - SR: ", data_initial.sumEntries("Categorization==5")
    data = data_initial.reduce("(BDT_out > 0.160 && Categorization==Categorization::MuonCR)")
    print "number of events mu 2017  - CR: ", data.sumEntries("Categorization==4")
if runningEra == 1 and not isMuon:
    print "number of events ele 2017 - SR: ", data_initial.sumEntries("Categorization==7")
    data = data_initial.reduce("(BDT_out > 0.110 && Categorization==Categorization::ElectronCR)")
    print "number of events ele 2017 - CR: ", data.sumEntries("Categorization==6")

if runningEra == 3 and isMuon:
    print "number of events mu - SR: ", data_initial.sumEntries("Categorization==1")
    data = data_initial.reduce("BDT_out > 0.150 && Categorization==Categorization::MuonCR")
    print "number of events mu - CR: ", data.sumEntries("Categorization==0")
if runningEra == 3 and not isMuon:
    print "number of events ele - SR: ", data_initial.sumEntries("Categorization==3")
    data = data_initial.reduce("BDT_out > 0.090 && Categorization==Categorization::ElectronCR")
    print "number of events ele - CR: ", data.sumEntries("Categorization==2")

print "Using ", data.numEntries(), " events to fit"


################################################################
#                                                              #
#---------------- Variables for bkg description ---------------#
#                                                              #
################################################################

#Now describe the background

#Parameters for muon Chebychev
a0_mu = ROOT.RooRealVar("a0_mu","a0_mu",0.15,-5.,5.)
a1_mu = ROOT.RooRealVar("a1_mu","a1_mu",-0.7,-5.,5.)
a2_mu = ROOT.RooRealVar("a2_mu","a2_mu",-0.2,-5.,5.)
a3_mu = ROOT.RooRealVar("a3_mu","a3_mu",-0.15,-5.,5.)
a4_mu = ROOT.RooRealVar("a4_mu","a4_mu",0.15,-5.,5.)
a5_mu = ROOT.RooRealVar("a5_mu","a5_mu",-0.7,-5.,5.)
a6_mu = ROOT.RooRealVar("a6_mu","a6_mu",-0.2,-5.,5.)
a7_mu = ROOT.RooRealVar("a7_mu","a7_mu",-0.15,-5.,5.)
# a8_mu = ROOT.RooRealVar("a8_mu","a8_mu",0.,-5.,5.)
# a5_mu = ROOT.RooRealVar("a5_mu","a5_mu",0.1,-5.,5.)
# a6_mu = ROOT.RooRealVar("a6_mu","a6_mu",0.1,-5.,5.)

#Parameters for electron Chebychev
a0_el = ROOT.RooRealVar("a0_el","a0_el",0.1,-5.,5.)
a1_el = ROOT.RooRealVar("a1_el","a1_el",-0.6,-5.,5.)
a2_el = ROOT.RooRealVar("a2_el","a2_el",-0.05,-5.,5.)
a3_el = ROOT.RooRealVar("a3_el","a3_el",-0.2,-5.,5.)
a4_el = ROOT.RooRealVar("a4_el","a4_el",0.1,-5.,5.)
a5_el = ROOT.RooRealVar("a5_el","a5_el",-0.6,-5.,5.)
a6_el = ROOT.RooRealVar("a6_el","a6_el",-0.05,-5.,5.)
a7_el = ROOT.RooRealVar("a7_el","a7_el",-0.2,-5.,5.)
a8_el = ROOT.RooRealVar("a8_el","a8_el",0.,-5.,5.)

#Parameters for muon Bernstein
b0_mu = ROOT.RooRealVar("b0_mu","b0_mu",0.06,0.,2.0)
b1_mu = ROOT.RooRealVar("b1_mu","b1_mu",0.09,0.,3.8)
b2_mu = ROOT.RooRealVar("b2_mu","b2_mu",7.,0.,10.)
b3_mu = ROOT.RooRealVar("b3_mu","b3_mu",3.8,0.,7.5)
b4_mu = ROOT.RooRealVar("b4_mu","b4_mu",1.7,0.,5)
b5_mu = ROOT.RooRealVar("b5_mu","b5_mu",3.7,0.,5.)
b6_mu = ROOT.RooRealVar("b6_mu","b6_mu",6.9,0.,11.)
b7_mu = ROOT.RooRealVar("b7_mu","b7_mu",3.8,0.,7.5)
#b8_mu = ROOT.RooRealVar("b8_mu","b8_mu",5.,0.,10.)

#Parameters for electron Bernstein
b0_el = ROOT.RooRealVar("b0_el","b0_el",2.5,0.,4.)
b1_el = ROOT.RooRealVar("b1_el","b1_el",1.5,0.,4.)
b2_el = ROOT.RooRealVar("b2_el","b2_el",7.,0.,10.)
b3_el = ROOT.RooRealVar("b3_el","b3_el",4.5,0.,8.)
b4_el = ROOT.RooRealVar("b4_el","b4_el",0.5,0.,3.)
b5_el = ROOT.RooRealVar("b5_el","b5_el",0.5,0.,3.)
b6_el = ROOT.RooRealVar("b6_el","b6_el",0.6,0.,2.)
b7_el = ROOT.RooRealVar("b7_el","b7_el",3.,0.,6.)
b8_el = ROOT.RooRealVar("b8_el","b8_el",2.,0.,5.)
# b9_el = ROOT.RooRealVar("b9_el","b9_el",5.,0.,10.)


backPDF_cheb_mu = ROOT.RooChebychev("backPDF_cheb_mu","backPDF_cheb_mu",Wmass,ROOT.RooArgList(a0_mu))#,a1_mu))#,a2_mu))#,a3_mu))#,a4_mu))#,a5_mu,a6_mu))
backPDF_cheb_el = ROOT.RooChebychev("backPDF_cheb_el","backPDF_cheb_el",Wmass,ROOT.RooArgList(a0_el))#,a1_el))#,a2_el,a3_el))#,a4_el))#,a6_el))


backPDF_bern_mu = ROOT.RooBernstein("backPDF_bern_mu","backPDF_bern_mu",Wmass,ROOT.RooArgList(b0_mu,b1_mu))#,b2_mu))#,b3_mu))#,b4_mu)) #,b5_mu,b6_mu))
backPDF_bern_el = ROOT.RooBernstein("backPDF_bern_el","backPDF_bern_el",Wmass,ROOT.RooArgList(b0_el,b1_el))#,b2_el))#,b3_el,b4_el))#,b5_el))#,b6_el))



if runningEra == 0 and useChebychev and isMuon:
    backPDF = backPDF_cheb_mu
if runningEra == 0 and useChebychev and not isMuon:
    backPDF = backPDF_cheb_el
if runningEra == 0 and not useChebychev and isMuon:
    backPDF = backPDF_bern_mu
if runningEra == 0 and not useChebychev and not isMuon:
    backPDF = backPDF_bern_el

if runningEra == 1 and useChebychev and isMuon:
    backPDF = backPDF_cheb_mu
if runningEra == 1 and useChebychev and not isMuon:
    backPDF = backPDF_cheb_el
if runningEra == 1 and not useChebychev and isMuon:
    backPDF = backPDF_bern_mu
if runningEra == 1 and not useChebychev and not isMuon:
    backPDF = backPDF_bern_el

if runningEra == 3 and useChebychev and isMuon:
    backPDF = backPDF_cheb_mu
if runningEra == 3 and useChebychev and not isMuon:
    backPDF = backPDF_cheb_el
if runningEra == 3 and not useChebychev and isMuon:
    backPDF = backPDF_bern_mu
if runningEra == 3 and not useChebychev and not isMuon:
    backPDF = backPDF_bern_el



################################################################
#                                                              #
#---------------------- Fit (and F-Test) ----------------------#
#                                                              #
################################################################

result_dataFit = backPDF.fitTo(data,ROOT.RooFit.Extended(0), ROOT.RooFit.NumCPU(4), ROOT.RooFit.Save() ) 
#For the signal region, I want the fit to be extended (Poisson fluctuation of unmber of events) to take into account that the total number of events is the sum of signal and background events.
#Either I do this, or I use a fraction frac*Nbkg+(1-frac)*Nsig, which will become a parameter of the fit and will have a Gaussian behavior (whilst the extended fit preserves the natural Poisson behavior)

print "minNll = ", result_dataFit.minNll()
print "2Delta_minNll = ", 2*(5417.81160974-result_dataFit.minNll()) # If 2*(NLL(N)-NLL(N+1)) > 3.85 -> N+1 is significant improvement


################################################################
#                                                              #
#-------------------------- Plotting --------------------------#
#                                                              #
################################################################

#First plot the signal region

xframe = Wmass.frame(55.,95.,15)
xframe.SetTitle(" ")
xframe.SetTitleOffset(1.4,"y")
xframe.SetMaximum(60)

data.plotOn(xframe, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
backPDF.plotOn(xframe)
    
canvas = ROOT.TCanvas()
xframe.Draw()

# Save the plot
if isMuon:
    canvas.SaveAs("plots/fitBackground_muon.pdf")
else:
    canvas.SaveAs("plots/fitBackground_electron.pdf")

################################################################
#                                                              #
#------------------------ Write to file -----------------------#
#                                                              #
################################################################
        
#Save the fit into a root file
fOutput = ROOT.TFile("fitBackground.root","RECREATE")

fOutput.cd()

workspace_bkg_out = ROOT.RooWorkspace("workspace_bkg")
getattr(workspace_bkg_out,'import')(backPDF_cheb_mu)
getattr(workspace_bkg_out,'import')(backPDF_cheb_el)
getattr(workspace_bkg_out,'import')(backPDF_bern_el)
getattr(workspace_bkg_out,'import')(backPDF_bern_mu)


workspace_bkg_out.Write()

fOutput.Close()

raw_input()
