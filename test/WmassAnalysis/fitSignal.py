#This code only fits the signal lineshape from MC and saves it into a RooWorkspace in a ROOT file

import ROOT
import argparse

ROOT.gROOT.ProcessLineSync(".L dCB/RooDoubleCBFast.cc+")

#---------------------------------#
p = argparse.ArgumentParser(description='Select whether to fit data or MC')
p.add_argument('runningEra_option', help='Type <<0>> for 2016, <<1>> for 2017, <<2>> for 2016+2017, <<3>> for 2016+2017+2018')
args = p.parse_args()

runningEra = args.runningEra_option
#---------------------------------#

#Define the observable
Wmass = ROOT.RooRealVar("Wmass","m_{#pi#gamma}",50.,100.,"GeV")

#Retrive the sample
fInput = ROOT.TFile("Tree_input_massfit_MC_" + runningEra + ".root")
fInput.cd()

mytree = fInput.Get("minitree")

#Define the event weight
weight = ROOT.RooRealVar("weight","The event weight",0.,10.)

#Define the signal category
isSignal = ROOT.RooCategory("isSignal","isSignal")
isSignal.defineType("Signal",1)
isSignal.defineType("Background",0)

#Define the mu/ele category
Categorization = ROOT.RooCategory("Categorization","Categorization")
Categorization.defineType("MuonCR",0)
Categorization.defineType("MuonSignal",1)
Categorization.defineType("ElectronCR",2)
Categorization.defineType("ElectronSignal",3)

Categorization.setRange("SignalRegion","MuonSignal,ElectronSignal") #The signal region used for the actual fit, i.e. the sum of the signal regions (defined by BDT cut) of muon and electron channels
Categorization.setRange("SignalRegion_mu","MuonSignal") #Only select events in the muon channel
Categorization.setRange("SignalRegion_ele","ElectronSignal") #Only select events in the lepton channel

#Create the RooDataSet, selecting only signal events (isSignal==1)
sample = ROOT.RooDataSet("sample","sample", ROOT.RooArgSet(Wmass,isSignal,weight,Categorization), ROOT.RooFit.Import(mytree),ROOT.RooFit.Cut("isSignal==1") ,ROOT.RooFit.WeightVar("weight"))

#Skim the signal according to the chosen range
data_Signal = sample.reduce(ROOT.RooFit.CutRange("SignalRegion"))
#data_Signal = sample.reduce(ROOT.RooFit.CutRange("SignalRegion_mu"))
#data_Signal = sample.reduce(ROOT.RooFit.CutRange("SignalRegion_ele"))

print "Using ", data_Signal.numEntries(), " events to fit the signal shape"

#Define the signal lineshape
Gauss_pole = ROOT.RooRealVar("Gauss_pole","The gaussian pole", 81.,70.,90.)
Gauss_sigma = ROOT.RooRealVar("Gauss_sigma","The gaussian sigma",5.,0.1,10.)
Gauss_W = ROOT.RooGaussian("Gauss_W","The Gaussian",Wmass,Gauss_pole,Gauss_sigma)

fracSig_prime = ROOT.RooRealVar("fracSig_prime","Partial fraction",0.5,0.,1.)
fracSig = ROOT.RooRealVar("fracSig","Fraction",0.5,0.,1.)

#Second, the resolution part
dCB_pole  = ROOT.RooRealVar("dCB_pole", "Double CB pole", 80.,75.,90.)
dCB_width = ROOT.RooRealVar("dCB_width", "Double CB width",1.,0.01,10.)
dCB_aL    = ROOT.RooRealVar("dCB_aL", "Double CB alpha left", 3., 0.1, 50.)
dCB_aR    = ROOT.RooRealVar("dCB_aR", "Double CB alpha right", 1., 0.1, 50.)
dCB_nL    = ROOT.RooRealVar("dCB_nL", "Double CB n left", 3., 0.1, 50.)
dCB_nR    = ROOT.RooRealVar("dCB_nR", "Double CB n right", 1., 0.1, 50.)
dCB       = ROOT.RooDoubleCBFast("dCB", "Double Crystal Ball", Wmass, dCB_pole, dCB_width, dCB_aL, dCB_nL, dCB_aR, dCB_nR)

totSignal = ROOT.RooAddPdf("totSignal","Total signal PDF",ROOT.RooArgList(dCB,Gauss_W),ROOT.RooArgList(fracSig))

totSignal.fitTo(data_Signal, ROOT.RooFit.SumW2Error(0))

#Make the plots
massplot = Wmass.frame()
massplot.SetTitle(" ")
massplot.SetTitleOffset(1.5,"y")
data_Signal.plotOn(massplot)
totSignal.plotOn(massplot)

canvas = ROOT.TCanvas()
canvas.cd()
massplot.Draw()
canvas.SaveAs("plots/SignalFit.pdf")

#Create the workspace and import the signal PDF
workspace = ROOT.RooWorkspace("myworkspace")
getattr(workspace,'import')(totSignal)

fOutput = ROOT.TFile("Signal_model_" + runningEra + ".root","RECREATE")
fOutput.cd()
workspace.Write()
massplot.Write("massplot")
h_signal_mu.Write()
fOutput.Close()

raw_input()
