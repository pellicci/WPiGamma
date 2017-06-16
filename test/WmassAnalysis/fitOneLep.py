
#This program only fits a single lepton sample at a time

import ROOT

#Define the observable
Wmass = ROOT.RooRealVar("Wmass","#pi-#gamma invariant mass",40.,120.)

#Retrive the sample
fInput = ROOT.TFile("Tree_MC.root")
fInput.cd()

mytree = fInput.Get("minitree")

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
data_lep = data.reduce("isMuon==0")
#data_lep = data.reduce("isMuon==isMuon::Electron")

print "Using ", data_lep.numEntries(), " events to fit the lepton shape"

#Import the signal
fInput_sigmodel = ROOT.TFile("Signal_model.root")
fInput_sigmodel.cd()

workspace = fInput_sigmodel.Get("myworkspace")

#Fix the signal parametrization
workspace.var("W_resol_pole").setConstant(1)
workspace.var("W_resol_width").setConstant(1)
workspace.var("W_resol_alpha").setConstant(1)
workspace.var("W_resol_n").setConstant(1)

totSignal = workspace.pdf("totSignal")

#Now describe the background
a0 = ROOT.RooRealVar("a0","a0",1.,-20.,20.)
a1 = ROOT.RooRealVar("a1","a1",0.,-20.,20.)
a2 = ROOT.RooRealVar("a2","a2",0.,-20.,20.)

backPDF = ROOT.RooChebychev("backPDF","backPDF",Wmass,ROOT.RooArgList(a0,a1,a2))

#Now fit signal + background

#Define signal as sigma*BR*eff*lumi
#Nsig = ROOT.RooRealVar("Nsig","Nsig",10.,0.,1000.)
W_cross_sec  = ROOT.RooRealVar("W_cross_sec","W_cross_sec",20.64/(2.*10.86) * 1000. ) #Taken from the ATLAS paper, in pb
W_pigamma_BR = ROOT.RooRealVar("W_pigamma_BR","W_pigamma_BR",0.1,0.,0.5) # The parameter of interest
W_eff        = ROOT.RooRealVar("W_eff","W_eff",708.*2./36900.) #For now, just the raw MC passed/generated number
lumi         = ROOT.RooRealVar("lumi","lumi",36.46 * 1000.) #In pb

Nsig = ROOT.RooFormulaVar("Nsig","@0*@1*@2*@3", ROOT.RooArgList(W_cross_sec,W_pigamma_BR,W_eff,lumi))
Nbkg = ROOT.RooRealVar("Nbkg","Nbkg",1000.,0.,3000.)

totPDF = ROOT.RooAddPdf("totPDF","Total PDF",ROOT.RooArgList(totSignal,backPDF),ROOT.RooArgList(Nsig,Nbkg))

totPDF.fitTo(data_lep,ROOT.RooFit.Extended(1), ROOT.RooFit.SumW2Error(0) )

xframe = Wmass.frame()
data_lep.plotOn(xframe)
totPDF.plotOn(xframe)

canvas = ROOT.TCanvas()
canvas.cd()
xframe.Draw()
canvas.SaveAs("fitOneLep.png")
