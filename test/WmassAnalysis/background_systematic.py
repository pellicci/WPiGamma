#This program fits both lepton samples at the same time

import ROOT
import math

#Define the observable
Wmass = ROOT.RooRealVar("Wmass","m_{#pi#gamma}",50.,100.,"GeV/c^{2}")


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

backPDF_mu = ROOT.RooChebychev("backPDF_mu","backPDF_mu",Wmass,ROOT.RooArgList(a0_mu,a1_mu))

#Alternative PDF for muon channel
a2_mu = ROOT.RooRealVar("a2_mu","a2_mu",2.,0.,3.)
a3_mu = ROOT.RooRealVar("a3_mu","a3_mu",0.7,0.,1.)
a4_mu = ROOT.RooRealVar("a4_mu","a4_mu",2.,0.,4.)
a5_mu = ROOT.RooRealVar("a5_mu","a5_mu",1.,0.,1.)

backPDF_mu_alt = ROOT.RooBernstein("backPDF_mu_alt","backPDF_mu_alt",Wmass,ROOT.RooArgList(a2_mu,a3_mu,a4_mu,a5_mu))

#Then the electron
a0_el = ROOT.RooRealVar("a0_el","a0_el",-0.1,-2.,2.)
a1_el = ROOT.RooRealVar("a1_el","a1_el",-0.1,-2.,2.)

backPDF_el = ROOT.RooChebychev("backPDF_el","backPDF_el",Wmass,ROOT.RooArgList(a0_el,a1_el))

#Aternative PDF for electron channel
a2_el = ROOT.RooRealVar("a2_el","a2_el",0.6,0.,1.)
a3_el = ROOT.RooRealVar("a3_el","a3_el",1.5,0.,5.)
a4_el = ROOT.RooRealVar("a4_el","a4_el",1.,0.,4.)
a5_el = ROOT.RooRealVar("a5_el","a5_el",0.3,0.,2.)

backPDF_el_alt = ROOT.RooBernstein("backPDF_el_alt","backPDF_el_alt",Wmass,ROOT.RooArgList(a2_el,a3_el,a4_el,a5_el))


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

#Systematic connected to scale factors

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

glb_W_xsec.setConstant(1)
glb_lumi.setConstant(1)
glb_eff_mu.setConstant(1)
glb_eff_el.setConstant(1)
#dCB_width_constr.setConstant(1)

Nsig_mu = ROOT.RooFormulaVar("Nsig_mu","@0*@1*@2*@3", ROOT.RooArgList(W_pigamma_BR, W_xsec_constr,lumi_constr,eff_mu_constr))
Nsig_el = ROOT.RooFormulaVar("Nsig_el","@0*@1*@2*@3", ROOT.RooArgList(W_pigamma_BR, W_xsec_constr,lumi_constr,eff_el_constr))

Nbkg_mu = ROOT.RooRealVar("Nbkg_mu","Nbkg_mu",125.,10.,300.)
Nbkg_el = ROOT.RooRealVar("Nbkg_el","Nbkg_el",150.,10.,310.)

totPDF_mu_unconstr = ROOT.RooAddPdf("totPDF_mu_unconstr","Total PDF for the mu channel",ROOT.RooArgList(totSignal,backPDF_mu),ROOT.RooArgList(Nsig_mu,Nbkg_mu))
totPDF_el_unconstr = ROOT.RooAddPdf("totPDF_el_unconstr","Total PDF for the el channel",ROOT.RooArgList(totSignal,backPDF_el),ROOT.RooArgList(Nsig_el,Nbkg_el))

totPDF_mu = ROOT.RooProdPdf("totPDF_mu","totPDF_mu", ROOT.RooArgList(totPDF_mu_unconstr,gauss_lumi,gauss_W_xsec,gauss_eff_mu,gauss_W_resol))
totPDF_el = ROOT.RooProdPdf("totPDF_el","totPDF_el", ROOT.RooArgList(totPDF_el_unconstr,gauss_lumi,gauss_W_xsec,gauss_eff_el,gauss_W_resol))

totPDF = ROOT.RooSimultaneous("totPDF","The total PDF",isMuon)
totPDF.addPdf(totPDF_mu,"Muon")
totPDF.addPdf(totPDF_el,"Electron")

#-----------Alternative total PDFs-------------#
totPDF_mu_unconstr_alt = ROOT.RooAddPdf("totPDF_mu_unconstr_alt","Total PDF for the mu channel - alternative",ROOT.RooArgList(totSignal,backPDF_mu_alt),ROOT.RooArgList(Nsig_mu,Nbkg_mu))
totPDF_el_unconstr_alt = ROOT.RooAddPdf("totPDF_el_unconstr_alt","Total PDF for the el channel - alternative",ROOT.RooArgList(totSignal,backPDF_el_alt),ROOT.RooArgList(Nsig_el,Nbkg_el))

totPDF_mu_alt = ROOT.RooProdPdf("totPDF_mu_alt","totPDF_mu_alt", ROOT.RooArgList(totPDF_mu_unconstr_alt,gauss_lumi,gauss_W_xsec,gauss_eff_mu,gauss_W_resol))
totPDF_el_alt = ROOT.RooProdPdf("totPDF_el_alt","totPDF_el_alt", ROOT.RooArgList(totPDF_el_unconstr_alt,gauss_lumi,gauss_W_xsec,gauss_eff_el,gauss_W_resol))

totPDF_alt = ROOT.RooSimultaneous("totPDF_alt","The total PDF - alternative",isMuon)
totPDF_alt.addPdf(totPDF_mu_alt,"Muon")
totPDF_alt.addPdf(totPDF_el_alt,"Electron")
#----------------------------------------------#


constrained_params = ROOT.RooArgSet()
constrained_params.add(dCB_width)
constrained_params.add(W_xsec_constr)
constrained_params.add(lumi_constr)
constrained_params.add(eff_mu_constr)
constrained_params.add(eff_el_constr)


totPDF.fitTo(data,ROOT.RooFit.Extended(1), ROOT.RooFit.SumW2Error(0), ROOT.RooFit.NumCPU(4), ROOT.RooFit.Constrain(constrained_params) )


# xframe_mu = Wmass.frame(10)
# xframe_mu.SetTitle("Fit to m_{#pi-#gamma} for the #mu channel")
# xframe_mu.SetMaximum(30)
# xframe_mu.SetTitleOffset(1.4,"y")
# data.plotOn(xframe_mu, ROOT.RooFit.Cut("isMuon==1"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
# totPDF.plotOn(xframe_mu, ROOT.RooFit.Slice(isMuon,"Muon"), ROOT.RooFit.ProjWData(data))

# xframe_el = Wmass.frame(10)
# xframe_el.SetTitle("Fit to m_{#pi-#gamma} for the e channel")
# xframe_el.SetMaximum(30)
# xframe_el.SetTitleOffset(1.4,"y")
# data.plotOn(xframe_el, ROOT.RooFit.Cut("isMuon==0"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
# totPDF.plotOn(xframe_el, ROOT.RooFit.Slice(isMuon,"Electron"), ROOT.RooFit.ProjWData(data))

# canvas = ROOT.TCanvas()
# canvas.Divide(2,1)
# canvas.cd(1)
# xframe_mu.Draw()
# canvas.cd(2)
# xframe_el.Draw()

# canvas.SaveAs("plots/fitMC.pdf")


#---------------------------Enable event generation according to fit PDFs----------------------------#
dataset_mu = totPDF_mu.generate(ROOT.RooArgSet(Wmass), 125, ROOT.RooFit.Extended(0))
dataset_el = totPDF_el.generate(ROOT.RooArgSet(Wmass), 290, ROOT.RooFit.Extended(0))
gen_data = ROOT.RooDataSet("gen_data", "comb dataset", ROOT.RooArgSet(Wmass), ROOT.RooFit.Index(isMuon), ROOT.RooFit.Import("Muon", dataset_mu), ROOT.RooFit.Import("Electron", dataset_el))
#----------------------------------------------------------------------------------------------------#

W_pigamma_BR.setRange(-0.01,0.01)

totPDF.fitTo(gen_data,ROOT.RooFit.Extended(1), ROOT.RooFit.SumW2Error(0), ROOT.RooFit.NumCPU(4), ROOT.RooFit.Constrain(constrained_params) )

BR_1 = W_pigamma_BR.getVal()
Sigma_BR_1 = W_pigamma_BR.getError()


xframe_mu = Wmass.frame(10)
xframe_mu.SetTitle("Fit to m_{#pi-#gamma} for the #mu channel")
xframe_mu.SetMaximum(30)
xframe_mu.SetTitleOffset(1.4,"y")
gen_data.plotOn(xframe_mu, ROOT.RooFit.Cut("isMuon==1"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
totPDF.plotOn(xframe_mu, ROOT.RooFit.Slice(isMuon,"Muon"), ROOT.RooFit.ProjWData(data))

xframe_el = Wmass.frame(10)
xframe_el.SetTitle("Fit to m_{#pi-#gamma} for the e channel")
xframe_el.SetMaximum(30)
xframe_el.SetTitleOffset(1.4,"y")
gen_data.plotOn(xframe_el, ROOT.RooFit.Cut("isMuon==0"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
totPDF.plotOn(xframe_el, ROOT.RooFit.Slice(isMuon,"Electron"), ROOT.RooFit.ProjWData(data))

canvas = ROOT.TCanvas()
canvas.Divide(2,1)
canvas.cd(1)
xframe_mu.Draw()
canvas.cd(2)
xframe_el.Draw()

canvas.SaveAs("plots/RefitMC.pdf")

totPDF_alt.fitTo(gen_data,ROOT.RooFit.Extended(1), ROOT.RooFit.SumW2Error(0), ROOT.RooFit.NumCPU(4), ROOT.RooFit.Constrain(constrained_params) )

BR_2 = W_pigamma_BR.getVal()

xframe_mu_alt = Wmass.frame(10)
xframe_mu_alt.SetTitle("Alternative fit to m_{#pi-#gamma} for the #mu channel")
xframe_mu_alt.SetMaximum(30)
xframe_mu_alt.SetTitleOffset(1.4,"y")
gen_data.plotOn(xframe_mu_alt, ROOT.RooFit.Cut("isMuon==1"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
totPDF_alt.plotOn(xframe_mu_alt, ROOT.RooFit.Slice(isMuon,"Muon"), ROOT.RooFit.ProjWData(data))

xframe_el_alt = Wmass.frame(10)
xframe_el_alt.SetTitle("Alternative fit to m_{#pi-#gamma} for the e channel")
xframe_el_alt.SetMaximum(30)
xframe_el_alt.SetTitleOffset(1.4,"y")
gen_data.plotOn(xframe_el_alt, ROOT.RooFit.Cut("isMuon==0"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
totPDF_alt.plotOn(xframe_el_alt, ROOT.RooFit.Slice(isMuon,"Electron"), ROOT.RooFit.ProjWData(data))

canvas_alt = ROOT.TCanvas()
canvas_alt.Divide(2,1)
canvas_alt.cd(1)
xframe_mu_alt.Draw()
canvas_alt.cd(2)
xframe_el_alt.Draw()

canvas_alt.SaveAs("plots/ReFitMC_alt.pdf")

#Save the fit into a root file

fOutput = ROOT.TFile("ReFitMC.root","RECREATE")
# fOutput = ROOT.TFile("fitMC_gen.root","RECREATE")
fOutput.cd()

# workspace = ROOT.RooWorkspace("workspace")
# getattr(workspace,'import')(data)
# # getattr(workspace,'import')(gen_data)
# getattr(workspace,'import')(totPDF)

# workspace.Write()

fOutput.Close()

print "BR_1: ", BR_1, "  BR_2: ", BR_2, "  Sigma_BR_1: ", Sigma_BR_1, "  pull: ", math.fabs(BR_2-BR_1)/Sigma_BR_1

raw_input()
