#This program fits both lepton samples at the same time

import ROOT
import math

ROOT.gROOT.ProcessLineSync(".L dCB/RooDoubleCBFast.cc+")

#Define if working on MC or DATA
isData = False

#--------some bools for scale factor systematics----------#

random_mu_SF  = False #-------if True, muon scale factors are sampled from a Gaussian
random_ele_SF = False #-------if True, electron scale factors are sampled from a Gaussian
random_ph_SF  = False #-------if True, photon scale factors are sampled from a Gaussian

#---------------------------------------------------------#

#Define the observable
Wmass = ROOT.RooRealVar("Wmass","m_{#pi#gamma}",50.,100.,"GeV/c^{2}")

#Retrive the sample
if isData:
    fInput = ROOT.TFile("Tree_Data.root")
else:
    if random_mu_SF:
        fInput = ROOT.TFile("Tree_MC_muSF.root")
    elif random_ele_SF:
        fInput = ROOT.TFile("Tree_MC_eleSF.root")
    elif random_ph_SF:
        fInput = ROOT.TFile("Tree_MC_phSF.root")
    else:
        fInput = ROOT.TFile("Tree_input_massfit_MC.root")
fInput.cd()

mytree = fInput.Get("minitree")

#Define the mu/ele category
isMuon = ROOT.RooCategory("isMuon","isMuon")
isMuon.defineType("Muon",1)
isMuon.defineType("Electron",0)

isSignalRegion = ROOT.RooCategory("isSignalRegion","isSignalRegion")
isSignalRegion.defineType("Signal",1)
isSignalRegion.defineType("CR",0)

#Define the event weight
weight = ROOT.RooRealVar("weight","The event weight",0.,10.)

#Create the RooDataSet. No need to import weight for signal only analysis
if isData:
    data_initial = ROOT.RooDataSet("data","data", ROOT.RooArgSet(Wmass,isMuon), ROOT.RooFit.Import(mytree))
else:
    data_initial = ROOT.RooDataSet("data","data", ROOT.RooArgSet(Wmass,isMuon,weight,isSignalRegion), ROOT.RooFit.Import(mytree), ROOT.RooFit.WeightVar("weight"))

data = data_initial.reduce("isSignalRegion == 1")

    
print "Using ", data.numEntries(), " events to fit"

#Import the signal
if random_mu_SF:
    fInput_sigmodel = ROOT.TFile("Signal_model_muSF.root")
elif random_ele_SF:
    fInput_sigmodel = ROOT.TFile("Signal_model_eleSF.root")
elif random_ph_SF:
    fInput_sigmodel = ROOT.TFile("Signal_model_phSF.root")
else:
    fInput_sigmodel = ROOT.TFile("Signal_model.root")
fInput_sigmodel.cd()

workspace = fInput_sigmodel.Get("myworkspace")

#Fix the signal parametrization
# workspace.var("W_resol_pole").setConstant(1)
# #workspace.var("W_resol_width").setConstant(1)
# workspace.var("W_resol_alpha").setConstant(1)
# workspace.var("W_resol_n").setConstant(1)
# workspace.var("Gauss_pole").setConstant(1)
# workspace.var("Gauss_sigma").setConstant(1)
# workspace.var("fracSig").setConstant(1)
workspace.var("dCB_pole").setConstant(1)
#workspace.var("dCB_width").setConstant(1)
workspace.var("dCB_aL").setConstant(1)
workspace.var("dCB_nL").setConstant(1)
workspace.var("dCB_aR").setConstant(1)
workspace.var("dCB_nR").setConstant(1)
workspace.var("Gauss_pole").setConstant(1)
workspace.var("Gauss_sigma").setConstant(1)
workspace.var("Gauss_pole_2").setConstant(1)
workspace.var("Gauss_sigma_2").setConstant(1)
workspace.var("fracSig_prime").setConstant(1)
workspace.var("fracSig").setConstant(1)

totSignal = workspace.pdf("totSignal")

#Now describe the background

#First the muon
a0_mu = ROOT.RooRealVar("a0_mu","a0_mu",-0.1,-2.,2.)
a1_mu = ROOT.RooRealVar("a1_mu","a1_mu",-0.1,-2.,2.)
# a2_mu = ROOT.RooRealVar("a2_mu","a2_mu",-0.1,-2.,2.)

# a0_mu = ROOT.RooRealVar("a0_mu","a0_mu",2,0.,10.)
# a1_mu = ROOT.RooRealVar("a1_mu","a1_mu",2,0.,10.)
# a2_mu = ROOT.RooRealVar("a2_mu","a2_mu",2,0.,10.)
backPDF_mu = ROOT.RooChebychev("backPDF_mu","backPDF_mu",Wmass,ROOT.RooArgList(a0_mu,a1_mu))

#Then the electron
a0_el = ROOT.RooRealVar("a0_el","a0_el",0.3,-1.,1.)
a1_el = ROOT.RooRealVar("a1_el","a1_el",-0.3,-1.,1.)
# a2_el = ROOT.RooRealVar("a2_el","a2_el",-0.1,-2.,2.)

# a0_el = ROOT.RooRealVar("a0_el","a0_el",2,0.,10.)
# a1_el = ROOT.RooRealVar("a1_el","a1_el",3,0.,10.)
# a2_el = ROOT.RooRealVar("a2_el","a2_el",2,0.,10.)
backPDF_el = ROOT.RooChebychev("backPDF_el","backPDF_el",Wmass,ROOT.RooArgList(a0_el,a1_el))

#exp_factor_mu = ROOT.RooRealVar("exp_factor_mu","exp_factor_mu", 0.5,-1.,1.)
#backPDF_mu = ROOT.RooExponential("backPDF_mu","backPDF_mu",Wmass,exp_factor_mu)

#exp_factor_el = ROOT.RooRealVar("exp_factor_el","exp_factor_el", 0.5,-1.,1.)
#backPDF_el = ROOT.RooExponential("backPDF_el","backPDF_el",Wmass,exp_factor_el)

#Gaussian distribution of W resolution width for systematics
# W_resol_width = workspace.var("W_resol_width")
dCB_width = workspace.var("dCB_width")
dCB_width_constr = ROOT.RooRealVar("dCB_width_constr","dCB_width_constr",dCB_width.getVal())
dCB_width_err = ROOT.RooRealVar("dCB_width_err","dCB_width_err",dCB_width.getError())
gauss_W_resol = ROOT.RooGaussian("gauss_W_resol","gauss_W_resol",dCB_width,dCB_width_constr,dCB_width_err)

#Now fit signal + background

#Define signal as sigma*BR*eff*lumi

#First the cross section, with a modifier for systematics
#CMS ttbar measurement/W->lnu BR (it is measured with both W in lnu), in pb
#http://cms-results.web.cern.ch/cms-results/public-results/publications/TOP-16-005/index.html
glb_W_xsec    = ROOT.RooRealVar("glb_W_xsec","glb_W_xsec", 2.*815.*0.1086, 0., 1000.)
W_xsec_constr = ROOT.RooRealVar("W_xsec_constr","W_x_sec_constr", 2.*815.*0.1086, 0., 1000.)
W_xsec_syst   = ROOT.RooRealVar("W_xsec_syst","W_xsec_syst",43.*2.*0.1086)
gauss_W_xsec  = ROOT.RooGaussian("gauss_W_xsec","gauss_W_xsec",glb_W_xsec,W_xsec_constr,W_xsec_syst)

#Represent the luminosity with a modifier for systematics. For 2016, 2.5% systematics
glb_lumi    = ROOT.RooRealVar("glb_lumi","glb_lumi",35.86 * 1000., 0., 50000.) #In pb
lumi_constr = ROOT.RooRealVar("lumi_constr","lumi_constr", 35.86 * 1000., 0., 50000.)
lumi_syst   = ROOT.RooRealVar("lumi_syst","lumi_syst", 35.86*0.025*1000.)
gauss_lumi  = ROOT.RooGaussian("gauss_lumi","gauss_lumi",glb_lumi,lumi_constr,lumi_syst) 


#Now the efficiency
totsig = 107810.  #total number of signal events
totmu = 7516.     #total number of signal muon events
totel = 6141.     #total number of signal electron events

glb_eff_mu    = ROOT.RooRealVar("glb_eff_mu","glb_eff_mu",totmu*2./totsig, 0., 1.) #For now, just the raw MC passed/generated number
eff_mu_constr = ROOT.RooRealVar("eff_mu_constr","eff_mu_constr", totmu*2./totsig, 0., 1.)
eff_mu_syst   = ROOT.RooRealVar("eff_mu_syst","eff_mu_syst", 4*totmu*(totsig-2*totmu)/(totsig*totsig*totsig))
gauss_eff_mu  = ROOT.RooGaussian("gauss_eff_mu","gauss_eff_mu",glb_eff_mu,eff_mu_constr,eff_mu_syst) 

glb_eff_el    = ROOT.RooRealVar("glb_eff_el","glb_eff_el", totel*2./totsig, 0., 1.) #For now, just the raw MC passed/generated number
eff_el_constr = ROOT.RooRealVar("eff_el_constr","eff_el_constr",totel*2./totsig, 0., 1.)
eff_el_syst   = ROOT.RooRealVar("eff_el_syst","eff_el_syst",  4*totel*(totsig-2*totel)/(totsig*totsig*totsig))
gauss_eff_el  = ROOT.RooGaussian("gauss_eff_el","gauss_eff_el",glb_eff_el,eff_el_constr,eff_el_syst) 

W_pigamma_BR = ROOT.RooRealVar("W_pigamma_BR","W_pigamma_BR",0.000001,0.,0.01) # The parameter of interest

#Introducing systematic connected to the background parametrization
eta = ROOT.RooRealVar("eta","eta", 1.,0.0001,3.)

glb_bkg_param    = ROOT.RooRealVar("glb_bkg_param","glb_bkg_param", 1., 0.0001, 3.)
bkg_param_syst   = ROOT.RooRealVar("bkg_param_syst","bkg_param_syst",0.42)
gauss_bkg_param  = ROOT.RooGaussian("gauss_bkg_param","gauss_bkg_param",glb_bkg_param,eta,bkg_param_syst)

glb_W_xsec.setConstant(1)
glb_bkg_param.setConstant(1)
glb_lumi.setConstant(1)
glb_eff_mu.setConstant(1)
glb_eff_el.setConstant(1)
#dCB_width_constr.setConstant(1)

Nsig_mu = ROOT.RooFormulaVar("Nsig_mu","@0*@1*@2*@3*@4", ROOT.RooArgList(W_pigamma_BR, W_xsec_constr,lumi_constr,eff_mu_constr,eta))
Nsig_el = ROOT.RooFormulaVar("Nsig_el","@0*@1*@2*@3*@4", ROOT.RooArgList(W_pigamma_BR, W_xsec_constr,lumi_constr,eff_el_constr,eta))

Nbkg_mu = ROOT.RooRealVar("Nbkg_mu","Nbkg_mu",125.,10.,300.)
Nbkg_el = ROOT.RooRealVar("Nbkg_el","Nbkg_el",150.,10.,300.)

totPDF_mu_unconstr = ROOT.RooAddPdf("totPDF_mu_unconstr","Total PDF for the mu channel",ROOT.RooArgList(totSignal,backPDF_mu),ROOT.RooArgList(Nsig_mu,Nbkg_mu))
totPDF_el_unconstr = ROOT.RooAddPdf("totPDF_el_unconstr","Total PDF for the el channel",ROOT.RooArgList(totSignal,backPDF_el),ROOT.RooArgList(Nsig_el,Nbkg_el))

totPDF_mu = ROOT.RooProdPdf("totPDF_mu","totPDF_mu", ROOT.RooArgList(totPDF_mu_unconstr,gauss_lumi,gauss_W_xsec,gauss_eff_mu,gauss_W_resol,gauss_bkg_param))
totPDF_el = ROOT.RooProdPdf("totPDF_el","totPDF_el", ROOT.RooArgList(totPDF_el_unconstr,gauss_lumi,gauss_W_xsec,gauss_eff_el,gauss_W_resol,gauss_bkg_param))

totPDF = ROOT.RooSimultaneous("totPDF","The total PDF",isMuon)
totPDF.addPdf(totPDF_mu,"Muon")
totPDF.addPdf(totPDF_el,"Electron")

constrained_params = ROOT.RooArgSet()
constrained_params.add(dCB_width)
constrained_params.add(eta)
constrained_params.add(W_xsec_constr)
constrained_params.add(lumi_constr)
constrained_params.add(eff_mu_constr)
constrained_params.add(eff_el_constr)

if isData:
    totPDF.fitTo(data,ROOT.RooFit.Extended(1), ROOT.RooFit.NumCPU(4), ROOT.RooFit.Constrain(constrained_params) )
else:
    totPDF.fitTo(data,ROOT.RooFit.Extended(1), ROOT.RooFit.SumW2Error(0), ROOT.RooFit.NumCPU(4), ROOT.RooFit.Constrain(constrained_params) )

xframe_mu = Wmass.frame(10)
xframe_mu.SetTitle("Fit to m_{#pi-#gamma} for the #mu channel")
xframe_mu.SetMaximum(30)
xframe_mu.SetTitleOffset(1.4,"y")
data.plotOn(xframe_mu, ROOT.RooFit.Cut("isMuon==1"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
totPDF.plotOn(xframe_mu, ROOT.RooFit.Slice(isMuon,"Muon"), ROOT.RooFit.ProjWData(data))

xframe_el = Wmass.frame(10)
xframe_el.SetTitle("Fit to m_{#pi-#gamma} for the e channel")
xframe_el.SetMaximum(30)
xframe_el.SetTitleOffset(1.4,"y")
data.plotOn(xframe_el, ROOT.RooFit.Cut("isMuon==0"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
totPDF.plotOn(xframe_el, ROOT.RooFit.Slice(isMuon,"Electron"), ROOT.RooFit.ProjWData(data))

canvas = ROOT.TCanvas()
canvas.Divide(2,1)
canvas.cd(1)
xframe_mu.Draw()
canvas.cd(2)
xframe_el.Draw()

if isData:
    canvas.SaveAs("plots/fitData.pdf")
else:
    if random_mu_SF:
        canvas.SaveAs("plots/fitMC_muSF.pdf")
    elif random_ele_SF:
        canvas.SaveAs("plots/fitMC_eleSF.pdf")
    elif random_ph_SF:
        canvas.SaveAs("plots/fitMC_phSF.pdf")
    else:
        canvas.SaveAs("plots/fitMC.pdf")

#Save the fit into a root file
if isData:
    fOutput = ROOT.TFile("fitData.root","RECREATE")
else:
    if random_mu_SF:
        fOutput = ROOT.TFile("fitMC_muSF.root","RECREATE")
    elif random_ele_SF:
        fOutput = ROOT.TFile("fitMC_eleSF.root","RECREATE")
    elif random_ph_SF:
        fOutput = ROOT.TFile("fitMC_phSF.root","RECREATE")
    else:
        fOutput = ROOT.TFile("fitMC.root","RECREATE")
fOutput.cd()

workspace = ROOT.RooWorkspace("workspace")
getattr(workspace,'import')(data)
getattr(workspace,'import')(totPDF)

workspace.Write()

fOutput.Close()

raw_input()
