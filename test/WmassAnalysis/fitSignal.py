
#This code only fits the signal lineshape from MC and saves it into a RooWorkspace in a ROOT file

import ROOT

#Define the observable
Wmass = ROOT.RooRealVar("Wmass","m_{#pi#gamma}",50.,100.,"GeV/c^{2}")

#Retrive the sample
fInput = ROOT.TFile("Tree_MC.root")
fInput.cd()

mytree = fInput.Get("minitree")

#Define the event weight
weight = ROOT.RooRealVar("weight","The event weight",0.,10.)

#Define the signal category
isSignal = ROOT.RooCategory("isSignal","isSignal")
isSignal.defineType("Signal",1)
isSignal.defineType("Background",0)

#Define the mu/ele category
isMuon = ROOT.RooCategory("isMuon","isMuon")
isMuon.defineType("Muon",1)
isMuon.defineType("Electron",0)

#Create the RooDataSet. No need to import weight for signal only analysis
sample = ROOT.RooDataSet("sample","sample", ROOT.RooArgSet(Wmass,isSignal,weight), ROOT.RooFit.Import(mytree))

#Skim the signal only
data_Signal = sample.reduce("isSignal==1")

print "Using ", data_Signal.numEntries(), " events to fit the signal shape"

#Define the signal lineshape
Gauss_pole = ROOT.RooRealVar("Gauss_pole","The gaussian pole", 75.,70.,80.)
Gauss_sigma = ROOT.RooRealVar("Gauss_sigma","The gaussian sigma",8.,0.1,10.)
Gauss_W = ROOT.RooGaussian("Gauss_W","The Gaussian",Wmass,Gauss_pole,Gauss_sigma)
fracSig = ROOT.RooRealVar("fracSig","Fraction",0.5,0.,1.)

#Second the resolution part
W_resol_pole  = ROOT.RooRealVar("W_resol_pole","The center of the resolution",80.,60.,90.)
W_resol_width = ROOT.RooRealVar("W_resol_width","The width of resolution",1.,0.01,10.)
W_resol_alpha = ROOT.RooRealVar("W_resol_alpha","The alpha of resolution",1.,-10.,10.)
W_resol_n     = ROOT.RooRealVar("W_resol_n","The n of resolution",2.,-10.,10.)
CBshape = ROOT.RooCBShape("CBshape","The resolution function",Wmass,W_resol_pole,W_resol_width,W_resol_alpha,W_resol_n)

totSignal = ROOT.RooAddPdf("totSignal","Total signal PDF",ROOT.RooArgList(CBshape,Gauss_W),ROOT.RooArgList(fracSig))

totSignal.fitTo(data_Signal)

massplot = Wmass.frame()
data_Signal.plotOn(massplot)
totSignal.plotOn(massplot)

canvas = ROOT.TCanvas()
canvas.cd()
massplot.Draw()
canvas.SaveAs("SignalFit.pdf")


workspace = ROOT.RooWorkspace("myworkspace")
getattr(workspace,'import')(totSignal)

fOutput = ROOT.TFile("Signal_model.root","RECREATE")
fOutput.cd()
workspace.Write()
fOutput.Close()
