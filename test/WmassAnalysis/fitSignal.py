#This code only fits the signal lineshape from MC and saves it into a RooWorkspace in a ROOT file

import ROOT

ROOT.gROOT.ProcessLineSync(".L dCB/RooDoubleCBFast.cc+")

#Define the observable
Wmass = ROOT.RooRealVar("Wmass","m_{#pi#gamma}",50.,100.,"GeV")

#Retrive the sample
fInput = ROOT.TFile("Tree_input_massfit_MC.root")
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

Categorization.setRange("SignalRegion","MuonSignal,ElectronSignal")


#Create the RooDataSet. No need to import weight for signal only analysis
sample = ROOT.RooDataSet("sample","sample", ROOT.RooArgSet(Wmass,isSignal,weight,Categorization), ROOT.RooFit.Import(mytree),ROOT.RooFit.Cut("isSignal==1"))

#Skim the signal only
data_Signal = sample.reduce(ROOT.RooFit.CutRange("SignalRegion"))


print "Using ", data_Signal.numEntries(), " events to fit the signal shape"

#Define the signal lineshape
Gauss_pole = ROOT.RooRealVar("Gauss_pole","The gaussian pole", 73.,70.,80.)
Gauss_sigma = ROOT.RooRealVar("Gauss_sigma","The gaussian sigma",4,0.1,10.)
Gauss_W = ROOT.RooGaussian("Gauss_W","The Gaussian",Wmass,Gauss_pole,Gauss_sigma)


Gauss_pole_2 = ROOT.RooRealVar("Gauss_pole_2","The second gaussian pole", 80.,70.,90.)
Gauss_sigma_2 = ROOT.RooRealVar("Gauss_sigma_2","The second gaussian sigma",5.,0.1,10.)
Gauss_W_2 = ROOT.RooGaussian("Gauss_W_2","The second Gaussian",Wmass,Gauss_pole_2,Gauss_sigma_2)


fracSig_prime = ROOT.RooRealVar("fracSig_prime","Partial fraction",0.5,0.,1.)
fracSig = ROOT.RooRealVar("fracSig","Fraction",0.5,0.,1.)


#Second the resolution part
dCB_pole  = ROOT.RooRealVar("dCB_pole", "Double CB pole", 80.,75.,90.)
dCB_width = ROOT.RooRealVar("dCB_width", "Double CB width",1.,0.01,10.)
dCB_aL    = ROOT.RooRealVar("dCB_aL", "Double CB alpha left", 3., 0.1, 50.)
dCB_aR    = ROOT.RooRealVar("dCB_aR", "Double CB alpha right", 1., 0.1, 50.)
dCB_nL    = ROOT.RooRealVar("dCB_nL", "Double CB n left", 3., 0.1, 50.)
dCB_nR    = ROOT.RooRealVar("dCB_nR", "Double CB n right", 1., 0.1, 50.)
dCB       = ROOT.RooDoubleCBFast("dCB", "Double Crystal Ball", Wmass, dCB_pole, dCB_width, dCB_aL, dCB_nL, dCB_aR, dCB_nR)

#totSignal_prime = ROOT.RooAddPdf("totSignal_prime","Partial signal PDF",ROOT.RooArgList(dCB,Gauss_W),ROOT.RooArgList(fracSig_prime))
#totSignal = ROOT.RooAddPdf("totSignal","Total signal PDF",ROOT.RooArgList(totSignal_prime,Gauss_W_2),ROOT.RooArgList(fracSig))
totSignal = ROOT.RooAddPdf("totSignal","Total signal PDF",ROOT.RooArgList(dCB,Gauss_W_2),ROOT.RooArgList(fracSig))


totSignal.fitTo(data_Signal)

#Make the plots
massplot = Wmass.frame()
massplot.SetTitle(" ")
massplot.SetTitleOffset(1.5,"y")
data_Signal.plotOn(massplot)
totSignal.plotOn(massplot)

#totSignal.paramOn(massplot,ROOT.RooFit.Layout(0.55))
#data_Signal.statOn(massplot,ROOT.RooFit.Layout(0.55,0.99,0.8))

canvas = ROOT.TCanvas()
canvas.cd()
massplot.Draw()
canvas.SaveAs("plots/SignalFit.pdf")


#Create the workspace and import the signal PDF
workspace = ROOT.RooWorkspace("myworkspace")
getattr(workspace,'import')(totSignal)


fOutput = ROOT.TFile("Signal_model.root","RECREATE")
fOutput.cd()
workspace.Write()
fOutput.Close()

raw_input()
