
#This program only fits a single lepton sample at a time

import ROOT

#Define the observable
Wmass = ROOT.RooRealVar("Wmass","#pi-#gamma invariant mass",40.,120.)

#Retrive the sample
fInput = ROOT.TFile("Tree_MC.root")
fInput.cd()

mytree = fInput.Get("mytree")

#Define the signal category
isSignal = ROOT.RooCategory("isSignal","isSignal")
isSignal.defineType("Signal",1)
isSignal.defineType("Background",0)

#Define the mu/ele category
isMuon = ROOT.RooCategory("isMuon","isMuon")
isMuon.defineType("Muon",1)
isMuon.defineType("Electron",0)

#Define the event weight
weight = ROOT.RooRealVar("weight","The event weight",0.,10.)

#Create the RooDataSet. No need to import weight for signal only analysis
data = ROOT.RooDataSet("data","data", ROOT.RooArgSet(Wmass,isSignal,isMuon,weight), ROOT.RooFit.Import(mytree), ROOT.RooFit.WeightVar("weight"))

#Skim one lepton sample only
data_lep = data.reduce("isMuon==isMuon::Muon")
#data_lep = data.reduce("isMuon==isMuon::Electron")

print "Using ", data_lep.numEntries(), " events to fit the signal shape"

#Import the signal
fInput_sigmodel = ROOT.TFile("Signal_model.root")
fInput_sigmodel.cd()

workspace = fInput_sigmodel.Get("myworkspace")

totSignal = workspace.pdf("totSignal")

#Now describe the background
