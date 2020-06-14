#This program fits both lepton samples at the same time

import ROOT
import math

ROOT.gROOT.ProcessLineSync(".L dCB/RooDoubleCBFast.cc+")

useChebychev = True #Choose to use Chebychev or Exponential to fit the background

################################################################
#                                                              #
#------------------- Define the observable --------------------#
#                                                              #
################################################################

Wmass = ROOT.RooRealVar("Wmass","m_{#pi#gamma}",50.,100.,"GeV/c^{2}")

################################################################
#                                                              #
#--------------------- Retrieve the sample --------------------#
#                                                              #
################################################################

fInput = ROOT.TFile("Tree_input_massfit_Data_3.root")

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
print "number of events mu - SR: ", data_initial.sumEntries("Categorization==1")
print "number of events ele - SR: ", data_initial.sumEntries("Categorization==3")
data = data_initial.reduce("(BDT_out > 0.206 && Categorization==Categorization::MuonCR) || (BDT_out > 0.209 && Categorization==Categorization::ElectronCR)")
print "number of events mu - CR: ", data.sumEntries("Categorization==0")
print "number of events ele - CR: ", data.sumEntries("Categorization==2")

print "Using ", data.numEntries(), " events to fit"


################################################################
#                                                              #
#---------------- Variables for bkg description ---------------#
#                                                              #
################################################################

#Parameters for Chebychev
a0_bkg = ROOT.RooRealVar("a0_bkg","a0_bkg",0.29,-1.,1.)
a1_bkg = ROOT.RooRealVar("a1_bkg","a1_bkg",-0.7,-5.,5.)
a2_bkg = ROOT.RooRealVar("a2_bkg","a2_bkg",-0.2,-5.,5.)
a3_bkg = ROOT.RooRealVar("a3_bkg","a3_bkg",-0.15,-5.,5.)

#Parameters for exponential
b0_bkg = ROOT.RooRealVar("b0_bkg","b0_bkg",0.001,0.,0.01)

backPDF_cheb = ROOT.RooChebychev("backPDF_cheb","backPDF_cheb",Wmass,ROOT.RooArgList(a0_bkg))#,a1_bkg))#,a2_bkg))#,a3_bkg))

backPDF_exp = ROOT.RooExponential("backPDF_exp","backPDF_exp",Wmass,b0_bkg)

if useChebychev:
    backPDF = backPDF_cheb
else:
    backPDF = backPDF_exp


################################################################
#                                                              #
#---------------------- Fit (and F-Test) ----------------------#
#                                                              #
################################################################

result_dataFit = backPDF.fitTo(data,ROOT.RooFit.Extended(0), ROOT.RooFit.NumCPU(4), ROOT.RooFit.Save() ) 
#For the signal region, I want the fit to be extended (Poisson fluctuation of unmber of events) to take into account that the total number of events is the sum of signal and background events.
#Either I do this, or I use a fraction frac*Nbkg+(1-frac)*Nsig, which will become a parameter of the fit and will have a Gaussian behavior (whilst the extended fit preserves the natural Poisson behavior)

print "minNll = ", result_dataFit.minNll()
print "2Delta_minNll = ", 2*(4106.73001635-result_dataFit.minNll()) # If 2*(NLL(N)-NLL(N+1)) > 3.85 -> N+1 is significant improvement

################################################################
#                                                              #
#-------------------------- Plotting --------------------------#
#                                                              #
################################################################

xframe = Wmass.frame(50.,100.,10)
xframe.SetTitle(" ")
xframe.SetTitleOffset(1.4,"y")
xframe.SetMaximum(60)

data.plotOn(xframe, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
backPDF.plotOn(xframe)
    
canvas = ROOT.TCanvas()
xframe.Draw()

# Save the plot
canvas.SaveAs("plots/fitBackground.pdf")

################################################################
#                                                              #
#------------------------ Write to file -----------------------#
#                                                              #
################################################################
        
#Save the fit into a root file
fOutput = ROOT.TFile("fitBackground.root","RECREATE")

fOutput.cd()

workspace_bkg_out = ROOT.RooWorkspace("workspace_bkg")
getattr(workspace_bkg_out,'import')(backPDF_cheb)
getattr(workspace_bkg_out,'import')(backPDF_exp)

workspace_bkg_out.Write()

fOutput.Close()

raw_input()
