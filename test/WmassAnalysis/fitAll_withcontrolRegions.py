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
    fInput = ROOT.TFile("Tree_input_massfit_Data.root")
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
Categorization = ROOT.RooCategory("Categorization","Categorization")
Categorization.defineType("MuonCR",0)
Categorization.defineType("MuonSignal",1)
Categorization.defineType("ElectronCR",2)
Categorization.defineType("ElectronSignal",3)

isSignal = ROOT.RooCategory("isSignal","isSignal")
isSignal.defineType("Signal",1)
isSignal.defineType("Background",0)

#Define the event weight
weight = ROOT.RooRealVar("weight","The event weight",0.,10.)

#Support variables
BDT_out = ROOT.RooRealVar("BDT_out","Output of BDT",-1.,1.)

#Create the RooDataSet. No need to import weight for signal only analysis

if isData:
    data_initial = ROOT.RooDataSet("data","data", ROOT.RooArgSet(Wmass,Categorization,BDT_out), ROOT.RooFit.Import(mytree))
else:
    data_initial = ROOT.RooDataSet("data","data", ROOT.RooArgSet(Wmass,Categorization,weight,BDT_out,isSignal), ROOT.RooFit.Import(mytree), ROOT.RooFit.WeightVar("weight"))

data = data_initial.reduce("BDT_out > -0.1")

# prova = ROOT.TFile("prova.root")
# prova.cd()
# data = prova.Get("totPDFData")

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
a0_mu = ROOT.RooRealVar("a0_mu","a0_mu",0.15,-5.,5.)
a1_mu = ROOT.RooRealVar("a1_mu","a1_mu",-0.7,-5.,5.)
a2_mu = ROOT.RooRealVar("a2_mu","a2_mu",-0.2,-5.,5.)
a3_mu = ROOT.RooRealVar("a3_mu","a3_mu",-0.15,-5.,5.)
a4_mu = ROOT.RooRealVar("a4_mu","a4_mu",-0.15,-5.,5.)
a5_mu = ROOT.RooRealVar("a5_mu","a5_mu",0.1,-5.,5.)
a6_mu = ROOT.RooRealVar("a6_mu","a6_mu",0.1,-5.,5.)

#a0_mu = ROOT.RooRealVar("a0_mu","a0_mu",0.5,0.,1.)
#a1_mu = ROOT.RooRealVar("a1_mu","a1_mu",0.7,0.,1.)
#a2_mu = ROOT.RooRealVar("a2_mu","a2_mu",0.7,0.,1.)
backPDF_mu = ROOT.RooChebychev("backPDF_mu","backPDF_mu",Wmass,ROOT.RooArgList(a0_mu,a1_mu,a2_mu,a3_mu)) #,a4_mu)) #,a5_mu,a6_mu))
#backPDF_mu = ROOT.RooBernstein("backPDF_mu","backPDF_mu",Wmass,ROOT.RooArgList(a0_mu,a1_mu,a2_mu)) #,a3_mu)) #,a4_mu)) #,a5_mu,a6_mu))

#Then the electron
a0_el = ROOT.RooRealVar("a0_el","a0_el",0.1,-5.,5.)
a1_el = ROOT.RooRealVar("a1_el","a1_el",-0.6,-5.,5.)
a2_el = ROOT.RooRealVar("a2_el","a2_el",-0.05,-5.,5.)
a3_el = ROOT.RooRealVar("a3_el","a3_el",-0.2,-5.,5.)
a4_el = ROOT.RooRealVar("a4_el","a4_el",0.1,-5.,5.)
a5_el = ROOT.RooRealVar("a5_el","a5_el",0.1,-5.,5.)
a6_el = ROOT.RooRealVar("a6_el","a6_el",0.1,-5.,5.)

#a0_el = ROOT.RooRealVar("a0_el","a0_el",0.5,0.,1.)
#a1_el = ROOT.RooRealVar("a1_el","a1_el",0.9,0.,1.)
#a2_el = ROOT.RooRealVar("a2_el","a2_el",0.7,0.,1.)
backPDF_el = ROOT.RooChebychev("backPDF_el","backPDF_el",Wmass,ROOT.RooArgList(a0_el,a1_el,a2_el))#,a3_el))#a4_el)) #,a5_el,a6_el))
#backPDF_el = ROOT.RooBernstein("backPDF_el","backPDF_el",Wmass,ROOT.RooArgList(a0_el,a1_el,a2_el)) #,a3_el)) #,a4_el)) #,a5_el,a6_el))

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
#totmu = 8913.     #total number of signal muon events
#totel = 6293.     #total number of signal electron events

totsig = 107810.  #total number of signal events
totmu = 9139.     #total number of signal muon events
totel = 6212.     #total number of signal electron events

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
bkg_param_syst   = ROOT.RooRealVar("bkg_param_syst","bkg_param_syst",0.27)
gauss_bkg_param  = ROOT.RooGaussian("gauss_bkg_param","gauss_bkg_param",glb_bkg_param,eta,bkg_param_syst)

glb_W_xsec.setConstant(1)
glb_bkg_param.setConstant(1)
glb_lumi.setConstant(1)
glb_eff_mu.setConstant(1)
glb_eff_el.setConstant(1)
#dCB_width_constr.setConstant(1)

#This is for the signal region
Nsig_mu = ROOT.RooFormulaVar("Nsig_mu","@0*@1*@2*@3*@4", ROOT.RooArgList(W_pigamma_BR, W_xsec_constr,lumi_constr,eff_mu_constr,eta))
Nsig_el = ROOT.RooFormulaVar("Nsig_el","@0*@1*@2*@3*@4", ROOT.RooArgList(W_pigamma_BR, W_xsec_constr,lumi_constr,eff_el_constr,eta))

Nbkg_mu = ROOT.RooRealVar("Nbkg_mu","Nbkg_mu",125.,1.,500.)
Nbkg_el = ROOT.RooRealVar("Nbkg_el","Nbkg_el",150.,1.,500.)

totPDF_mu_unconstr = ROOT.RooAddPdf("totPDF_mu_unconstr","Total PDF for the mu channel",ROOT.RooArgList(totSignal,backPDF_mu),ROOT.RooArgList(Nsig_mu,Nbkg_mu))
totPDF_el_unconstr = ROOT.RooAddPdf("totPDF_el_unconstr","Total PDF for the el channel",ROOT.RooArgList(totSignal,backPDF_el),ROOT.RooArgList(Nsig_el,Nbkg_el))

totPDF_mu = ROOT.RooProdPdf("totPDF_mu","totPDF_mu", ROOT.RooArgList(totPDF_mu_unconstr,gauss_lumi,gauss_W_xsec,gauss_eff_mu,gauss_W_resol,gauss_bkg_param))
totPDF_el = ROOT.RooProdPdf("totPDF_el","totPDF_el", ROOT.RooArgList(totPDF_el_unconstr,gauss_lumi,gauss_W_xsec,gauss_eff_el,gauss_W_resol,gauss_bkg_param))

#Now for the CR
Nbkg_mu_CR = ROOT.RooRealVar("Nbkg_mu_CR","Nbkg_mu_CR",100.,1.,40000.)
Nbkg_el_CR = ROOT.RooRealVar("Nbkg_el_CR","Nbkg_el_CR",100.,1.,40000.)

totPDF_mu_CR = ROOT.RooExtendPdf("totPDF_mu_CR","Background PDF for the CR for mu",backPDF_mu,Nbkg_mu_CR)
totPDF_el_CR = ROOT.RooExtendPdf("totPDF_el_CR","Background PDF for the CR for el",backPDF_el,Nbkg_el_CR)

#totPDF_mu_unconstr_CR = ROOT.RooExtendPdf("totPDF_mu_unconstr_CR","Background PDF for the CR for mu",backPDF_mu,Nbkg_mu_CR)
#totPDF_el_unconstr_CR = ROOT.RooExtendPdf("totPDF_el_unconstr_CR","Background PDF for the CR for el",backPDF_el,Nbkg_el_CR)

#totPDF_mu_CR = ROOT.RooProdPdf("totPDF_mu_CR","totPDF_mu_CR", ROOT.RooArgList(totPDF_mu_unconstr_CR,gauss_bkg_param))
#totPDF_el_CR = ROOT.RooProdPdf("totPDF_el_CR","totPDF_el_CR", ROOT.RooArgList(totPDF_el_unconstr_CR,gauss_bkg_param))

#Create the global simultaneous PDF
totPDF = ROOT.RooSimultaneous("totPDF","The total PDF",Categorization)
if not isData:
    totPDF.addPdf(totPDF_mu_CR,"MuonCR")
    totPDF.addPdf(totPDF_mu,"MuonSignal")
    totPDF.addPdf(totPDF_el_CR,"ElectronCR")
    totPDF.addPdf(totPDF_el,"ElectronSignal")
else: #Only fit the CR in case of data (before unblinding)
    totPDF.addPdf(totPDF_mu_CR,"MuonCR")
    totPDF.addPdf(totPDF_el_CR,"ElectronCR")

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
    totPDF.fitTo(data,ROOT.RooFit.Extended(1), ROOT.RooFit.SumW2Error(0), ROOT.RooFit.NumCPU(2), ROOT.RooFit.Constrain(constrained_params) )

#First plot the signal region
if not isData:
    xframe_mu = Wmass.frame(55.,95.,15)
    #xframe_mu.SetTitle("Fit to m_{#pi-#gamma} for the #mu channel")
    xframe_mu.SetTitle(" ")
    xframe_mu.SetMaximum(30)
    xframe_mu.SetTitleOffset(1.4,"y")
    data.plotOn(xframe_mu, ROOT.RooFit.Cut("Categorization==1"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    totPDF.plotOn(xframe_mu, ROOT.RooFit.Slice(Categorization,"MuonSignal"), ROOT.RooFit.ProjWData(data))
    
    xframe_el = Wmass.frame(55.,95.,15)
    #xframe_el.SetTitle("Fit to m_{#pi-#gamma} for the e channel")
    xframe_el.SetTitle(" ")
    xframe_el.SetMaximum(30)
    xframe_el.SetTitleOffset(1.4,"y")
    data.plotOn(xframe_el, ROOT.RooFit.Cut("Categorization==3"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    totPDF.plotOn(xframe_el, ROOT.RooFit.Slice(Categorization,"ElectronSignal"), ROOT.RooFit.ProjWData(data))
    
    canvas = ROOT.TCanvas()
    canvas.Divide(2,1)
    canvas.cd(1)
    xframe_mu.Draw()
    canvas.cd(2)
    xframe_el.Draw()

# if isData:
#     canvas.SaveAs("plots/fitData_signalR.pdf")
if not isData:
    if random_mu_SF:
        canvas.SaveAs("plots/fitMC_muSF.pdf")
    elif random_ele_SF:
        canvas.SaveAs("plots/fitMC_eleSF.pdf")
    elif random_ph_SF:
        canvas.SaveAs("plots/fitMC_phSF.pdf")
    else:
        canvas.SaveAs("plots/fitMC_signalR.pdf")

#Now plot the CR
xframe_mu_CR = Wmass.frame(55.,95.,15)
#xframe_mu_CR.SetTitle("Fit to m_{#pi-#gamma} for the #mu channel - Control Region")
xframe_mu_CR.SetTitle(" ")
xframe_mu_CR.SetMaximum(30)
xframe_mu_CR.SetTitleOffset(1.4,"y")
data.plotOn(xframe_mu_CR, ROOT.RooFit.Cut("Categorization==0"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
totPDF.plotOn(xframe_mu_CR, ROOT.RooFit.Slice(Categorization,"MuonCR"), ROOT.RooFit.ProjWData(data))

xframe_el_CR = Wmass.frame(55.,95.,15)
#xframe_el_CR.SetTitle("Fit to m_{#pi-#gamma} for the e channel - Control Region")
xframe_el_CR.SetTitle(" ")
xframe_el_CR.SetMaximum(30)
xframe_el_CR.SetTitleOffset(1.4,"y")
data.plotOn(xframe_el_CR, ROOT.RooFit.Cut("Categorization==2"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
totPDF.plotOn(xframe_el_CR, ROOT.RooFit.Slice(Categorization,"ElectronCR"), ROOT.RooFit.ProjWData(data))

canvas_CR = ROOT.TCanvas()
canvas_CR.Divide(2,1)
canvas_CR.cd(1)
xframe_mu_CR.Draw()
canvas_CR.cd(2)
xframe_el_CR.Draw()

if isData:
    canvas_CR.SaveAs("plots/fitData_CR.pdf")
else:
    canvas_CR.SaveAs("plots/fitMC_CR.pdf")
        
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

# W_pigamma_BR.setVal(0.)
# newdata = totPDF.generate(ROOT.RooArgSet(Categorization),7921,ROOT.RooFit.ProtoData(data))
# prova = ROOT.TFile("prova.root","RECREATE")
# prova.cd()
# newdata.Write()
# prova.Close()


raw_input()
