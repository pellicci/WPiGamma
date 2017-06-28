
#This program only fits a single lepton sample at a time

import ROOT

#Define the observable
Wmass = ROOT.RooRealVar("Wmass","#pi-#gamma invariant mass",40.,120.)

#Retrive the sample
fInput = ROOT.TFile("Tree_MC.root")
fInput.cd()

mytree = fInput.Get("minitree")

#Define the signal category
#isSignal = ROOT.RooCategory("isSignal","isSignal")
#isSignal.defineType("Signal",1)
#isSignal.defineType("Background",0)

#Define the mu/ele category
isMuon = ROOT.RooCategory("isMuon","isMuon")
isMuon.defineType("Muon",1)
isMuon.defineType("Electron",0)

#Define the event weight
weight = ROOT.RooRealVar("weight","The event weight",0.,10.)

#Create the RooDataSet. No need to import weight for signal only analysis
data = ROOT.RooDataSet("data","data", ROOT.RooArgSet(Wmass,isMuon,weight), ROOT.RooFit.Import(mytree), ROOT.RooFit.WeightVar("weight"))

print "Using ", data.numEntries(), " events to fit"

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

#First the muon
a0_mu = ROOT.RooRealVar("a0_mu","a0_mu",0.1,-5.,5.)
a1_mu = ROOT.RooRealVar("a1_mu","a1_mu",-0.1,-5.,5.)
backPDF_mu = ROOT.RooChebychev("backPDF_mu","backPDF_mu",Wmass,ROOT.RooArgList(a0_mu,a1_mu))

#Then the electron
a0_el = ROOT.RooRealVar("a0_el","a0_el",0.1,-5.,5.)
a1_el = ROOT.RooRealVar("a1_el","a1_el",-0.1,-5.,5.)
backPDF_el = ROOT.RooChebychev("backPDF_el","backPDF_el",Wmass,ROOT.RooArgList(a0_el,a1_el))


#Now fit signal + background

#Define signal as sigma*BR*eff*lumi
W_cross_sec  = ROOT.RooRealVar("W_cross_sec","W_cross_sec",20.64/(2.*10.86) * 1000. ) #Taken from the ATLAS paper, in pb
lumi         = ROOT.RooRealVar("lumi","lumi",36.46 * 1000.) #In pb
W_pigamma_BR = ROOT.RooRealVar("W_pigamma_BR","W_pigamma_BR",0.00001,0.,0.1) # The parameter of interest

W_eff_mu     = ROOT.RooRealVar("W_eff_mu","W_eff_mu",708.*2./36900.) #For now, just the raw MC passed/generated number
W_eff_el     = ROOT.RooRealVar("W_eff_el","W_eff_el",1326.*2./36900.) #For now, just the raw MC passed/generated number

Nsig_mu = ROOT.RooFormulaVar("Nsig_mu","@0*@1*@2*@3", ROOT.RooArgList(W_cross_sec,W_pigamma_BR,W_eff_mu,lumi))
Nsig_el = ROOT.RooFormulaVar("Nsig_el","@0*@1*@2*@3", ROOT.RooArgList(W_cross_sec,W_pigamma_BR,W_eff_el,lumi))

Nbkg_mu = ROOT.RooRealVar("Nbkg_mu","Nbkg_mu",100.,0.,500.)
Nbkg_el = ROOT.RooRealVar("Nbkg_el","Nbkg_el",160.,0.,1000.)

totPDF_mu = ROOT.RooAddPdf("totPDF_mu","Total PDF for the mu channel",ROOT.RooArgList(totSignal,backPDF_mu),ROOT.RooArgList(Nsig_mu,Nbkg_mu))
totPDF_el = ROOT.RooAddPdf("totPDF_el","Total PDF for the el channel",ROOT.RooArgList(totSignal,backPDF_el),ROOT.RooArgList(Nsig_el,Nbkg_el))

totPDF = ROOT.RooSimultaneous("totPDF","The total PDF",isMuon)
totPDF.addPdf(totPDF_mu,"Muon")
totPDF.addPdf(totPDF_el,"Electron")

totPDF.fitTo(data,ROOT.RooFit.Extended(1), ROOT.RooFit.SumW2Error(1) )

xframe_mu = Wmass.frame()
data.plotOn(xframe_mu, ROOT.RooFit.Cut("isMuon==1"))
totPDF.plotOn(xframe_mu, ROOT.RooFit.Slice(isMuon,"Muon"), ROOT.RooFit.ProjWData(data))

xframe_el = Wmass.frame()
data.plotOn(xframe_el, ROOT.RooFit.Cut("isMuon==0"))
totPDF.plotOn(xframe_el, ROOT.RooFit.Slice(isMuon,"Electron"), ROOT.RooFit.ProjWData(data))


canvas = ROOT.TCanvas()
canvas.Divide(2,1)
canvas.cd(1)
xframe_mu.Draw()
canvas.cd(2)
xframe_el.Draw()
canvas.SaveAs("fitAllLep.png")

#Save the fit into a root file
fOutput = ROOT.TFile("fitAllLep.root","RECREATE")
fOutput.cd()

workspace = ROOT.RooWorkspace("workspace")
getattr(workspace,'import')(data)
getattr(workspace,'import')(totPDF)

workspace.Write()

fOutput.Close()
