
#This program only fits a single lepton sample at a time

import ROOT
import math

#Define the observable
Wmass = ROOT.RooRealVar("Wmass","#pi-#gamma invariant mass",50.,110.)

#Retrive the sample
fInput = ROOT.TFile("Tree_MC.root")
fInput.cd()

mytree = fInput.Get("minitree")

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

#First the cross section, with a modifier for systematics
#CMS ttbar measurement/W->lnu BR (it is measured with both W in lnu), in pb
#http://cms-results.web.cern.ch/cms-results/public-results/publications/TOP-16-005/index.html
W_cross_sec       = ROOT.RooRealVar("W_cross_sec","W_cross_sec", 815./0.1086 )
W_cross_sec_ratio = ROOT.RooRealVar("W_cross_sec_ratio","W_cross_sec_ratio",1.,0.,2.)
glb_W_xsec        = ROOT.RooRealVar("glb_W_xsec","glb_W_xsec",1.,0.,3.)
W_xsec_syst       = ROOT.RooRealVar("W_xsec_syst","W_xsec_syst",43./815.)
gauss_W_xsec      = ROOT.RooGaussian("gauss_W_xsec","gauss_W_xsec",glb_W_xsec,W_cross_sec_ratio,W_xsec_syst)

#Represent the luminosity with a modifier for systematics. For 2016, 2.5% systematics
lumi       = ROOT.RooRealVar("lumi","lumi",36.46 * 1000.) #In pb
lumi_ratio = ROOT.RooRealVar("lumi_ratio","lumi_ratio",1.,0.,2.)
glb_lumi   = ROOT.RooRealVar("glb_lumi","glb_lumi",1.,0.,3.)
lumi_syst  = ROOT.RooRealVar("lumi_syst","lumi_syst",0.025)
gauss_lumi = ROOT.RooGaussian("gauss_lumi","gauss_lumi",glb_lumi,lumi_ratio,lumi_syst) 

W_eff_mu     = ROOT.RooRealVar("W_eff_mu","W_eff_mu",708.*2./36900.) #For now, just the raw MC passed/generated number
eff_mu_ratio = ROOT.RooRealVar("eff_mu_ratio","eff_mu_ratio",1.,0.,2.)
glb_eff_mu   = ROOT.RooRealVar("glb_eff_mu","glb_eff_mu",1.,0.,3.)
eff_mu_syst  = ROOT.RooRealVar("eff_mu_syst","eff_mu_syst", math.sqrt(1./708. + 2./36900.)*708.*2./36900.)
gauss_eff_mu = ROOT.RooGaussian("gauss_eff_mu","gauss_eff_mu",glb_eff_mu,eff_mu_ratio,eff_mu_syst) 

W_eff_el     = ROOT.RooRealVar("W_eff_el","W_eff_el",1326.*2./36900.) #For now, just the raw MC passed/generated number
eff_el_ratio = ROOT.RooRealVar("eff_el_ratio","eff_el_ratio",1.,0.,2.)
glb_eff_el   = ROOT.RooRealVar("glb_eff_el","glb_eff_el",1.,0.,3.)
eff_el_syst  = ROOT.RooRealVar("eff_el_syst","eff_el_syst", math.sqrt(1./708. + 2./36900.)*708.*2./36900.)
gauss_eff_el = ROOT.RooGaussian("gauss_eff_el","gauss_eff_el",glb_eff_el,eff_el_ratio,eff_el_syst) 

W_pigamma_BR = ROOT.RooRealVar("W_pigamma_BR","W_pigamma_BR",0.00001,0.,0.1) # The parameter of interest

glb_W_xsec.setConstant(1)
glb_lumi.setConstant(1)
glb_eff_mu.setConstant(1)
glb_eff_el.setConstant(1)

Nsig_mu = ROOT.RooFormulaVar("Nsig_mu","@0*@1*@2*@3*@4*@5", ROOT.RooArgList(W_cross_sec,W_pigamma_BR,W_eff_mu,lumi,W_cross_sec_ratio,lumi_ratio,eff_mu_ratio))
Nsig_el = ROOT.RooFormulaVar("Nsig_el","@0*@1*@2*@3*@4*@5", ROOT.RooArgList(W_cross_sec,W_pigamma_BR,W_eff_el,lumi,W_cross_sec_ratio,lumi_ratio,eff_el_ratio))

Nbkg_mu = ROOT.RooRealVar("Nbkg_mu","Nbkg_mu",100.,0.,500.)
Nbkg_el = ROOT.RooRealVar("Nbkg_el","Nbkg_el",160.,0.,1000.)

totPDF_mu_unconstr = ROOT.RooAddPdf("totPDF_mu_unconstr","Total PDF for the mu channel",ROOT.RooArgList(totSignal,backPDF_mu),ROOT.RooArgList(Nsig_mu,Nbkg_mu))
totPDF_el_unconstr = ROOT.RooAddPdf("totPDF_el_unconstr","Total PDF for the el channel",ROOT.RooArgList(totSignal,backPDF_el),ROOT.RooArgList(Nsig_el,Nbkg_el))

totPDF_mu = ROOT.RooProdPdf("totPDF_mu","totPDF_mu", ROOT.RooArgList(totPDF_mu_unconstr,gauss_lumi,gauss_W_xsec,gauss_eff_mu))
totPDF_el = ROOT.RooProdPdf("totPDF_el","totPDF_el", ROOT.RooArgList(totPDF_el_unconstr,gauss_lumi,gauss_W_xsec,gauss_eff_el))

totPDF = ROOT.RooSimultaneous("totPDF","The total PDF",isMuon)
totPDF.addPdf(totPDF_mu,"Muon")
totPDF.addPdf(totPDF_el,"Electron")

constrained_params = ROOT.RooArgSet()
constrained_params.add(lumi_ratio)
constrained_params.add(W_cross_sec_ratio)
constrained_params.add(eff_mu_ratio)
constrained_params.add(eff_el_ratio)

totPDF.fitTo(data,ROOT.RooFit.Extended(1), ROOT.RooFit.SumW2Error(1), ROOT.RooFit.NumCPU(4), ROOT.RooFit.Constrain(constrained_params) )

xframe_mu = Wmass.frame(50)
data.plotOn(xframe_mu, ROOT.RooFit.Cut("isMuon==1"))
totPDF.plotOn(xframe_mu, ROOT.RooFit.Slice(isMuon,"Muon"), ROOT.RooFit.ProjWData(data))

xframe_el = Wmass.frame(50)
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
