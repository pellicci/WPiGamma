#This program fits both lepton samples at the same time

import ROOT
import math
import argparse

ROOT.gROOT.ProcessLineSync(".L dCB/RooDoubleCBFast.cc+")

#---------------------------------#
p = argparse.ArgumentParser(description='Select whether to fit data or MC')
p.add_argument('isData_option', help='Type <<data>> or <<MC>>')
p.add_argument('runningEra_option', help='Type <<0>> for 2016, <<1>> for 2017, <<2>> for 2016+2017')
args = p.parse_args()

isData = False
# Switch from MC to data channel
if args.isData_option == "data":
    isData = True

runningEra = int(args.runningEra_option)
#---------------------------------#

################################################################
#                                                              #
#------------------------ Instructions ------------------------#
#                                                              #
################################################################

########------- Fit to restricted CR to determine the best background parametrization -------########
# useChebychev = 0 will use Chebychev PDFs to describe the backround (they can be replaced by other PDFs in the " if useChebychev == 0:" statement);
# restricted_CR = True will perform the fit on restricted CRs, cuttin on a different BDT output according to channel and year
# isAlternativeBkgDescription = False will make sure the Categorization has also MuonCR and ElectronCR for the two years

########------- Determine the systematic on bkg parametrization -------########
# the useChebychev variable has to be put equal to 0, then to 1, 2, 3, and 4. This way, the PDFs used to describe the background will be changed one by one (per channel and per year), and it will be possible to obtain the value of W_pigamma_BR for each case
# restricted_CR = False will make sure the fit is performed on the signal regions
# isAlternativeBkgDescription = True will allow W_pigamma_BR to float negative, and will set the value of the systematic on the bkg parametrization to be negligible during the fit

########------- Fit to data signal regions -------########
# useChebychev = 0 will use the nominal background parametrization
# restricted_CR = False will perform the fit (and the plotting) on the signal regions
# isAlternativeBkgDescription = False will make sure the W_pigamma_BR parameter is positive definite and the systematic on background parametrization is introduced

useChebychev = 0 # 0: use Chebychev for ALL the bkg PDFs. 1: use Bernstein for mu 2016. 2: use Bernstein for ele 2016. 3: use Bernstein for mu 2017. 4: use Bernstein for ele 2017
restricted_CR = False # If True, data will be reduced to have the content of the restricted control regions
isAlternativeBkgDescription = False # To be used when trying to fit with alternative bkg description, in order to estimate a systematic. If True, it will allow W_pigamma_BR to float negative. Moreover, it will use Signal+Background in the totPDF, so that the fit to the restricted CRs will contain also the POI BR, which will be used to calculate the pull and hence to estimate the systematic


################################################################
#                                                              #
#------------------- Define the observable --------------------#
#                                                              #
################################################################

Wmass = ROOT.RooRealVar("Wmass","m_{#pi#gamma}",50.,100.,"GeV/c^{2}")
Wmass.setRange("LowSideband",55.,65.)
Wmass.setRange("HighSideband",90.,95.)

################################################################
#                                                              #
#--------------------- Retrieve the sample --------------------#
#                                                              #
################################################################

if isData:
    fInput = ROOT.TFile("Tree_input_massfit_Data_" + str(runningEra) + ".root")
else:
    fInput = ROOT.TFile("Tree_input_massfit_MC_ " + str(runningEra) + ".root")

fInput.cd()

mytree = fInput.Get("minitree")

################################################################
#                                                              #
#-------------------- Define the categories -------------------#
#                                                              #
################################################################

#Define the mu/ele category for 2016 and 2017
Categorization = ROOT.RooCategory("Categorization","Categorization")
if runningEra == 0:
    Categorization.defineType("MuonSignal_2016",1)
    Categorization.defineType("ElectronSignal_2016",3)
if runningEra == 1:
    Categorization.defineType("MuonSignal_2017",5)
    Categorization.defineType("ElectronSignal_2017",7)
if runningEra == 2:
    Categorization.defineType("MuonSignal_2016",1)
    Categorization.defineType("ElectronSignal_2016",3)
    Categorization.defineType("MuonSignal_2017",5)
    Categorization.defineType("ElectronSignal_2017",7)
if restricted_CR and not isAlternativeBkgDescription:
    Categorization.defineType("MuonCR_2016",0)
    Categorization.defineType("ElectronCR_2016",2)
    Categorization.defineType("MuonCR_2017",4)
    Categorization.defineType("ElectronCR_2017",6)

isSignal = ROOT.RooCategory("isSignal","isSignal")
isSignal.defineType("Signal",1)
isSignal.defineType("Background",0)

isMuon = ROOT.RooCategory("isMuon","isMuon")
isMuon.defineType("Muon",1)
isMuon.defineType("Electron",0)

################################################################
#                                                              #
#--------------- Get event weight and BDT output --------------#
#                                                              #
################################################################

minWeight_inMC = -50. # Set the minimum weight an event can have when processing MC. This is useful to remove spikes which cannot be fitted
maxWeight_inMC = 50.  # Set the maximum weight an event can have when processing MC. This is useful to remove spikes which cannot be fitted

#Define the event weight
weight = ROOT.RooRealVar("weight","The event weight",minWeight_inMC,maxWeight_inMC)

#Support variables
BDT_out = ROOT.RooRealVar("BDT_out","Output of BDT",-1.,1.)

#Create the RooDataSet. No need to import weight for signal only analysis

if isData:
    data_initial = ROOT.RooDataSet("data","data", ROOT.RooArgSet(Wmass,Categorization,BDT_out,isMuon), ROOT.RooFit.Import(mytree))
else:
    data_initial = ROOT.RooDataSet("data","data", ROOT.RooArgSet(Wmass,Categorization,weight,BDT_out,isSignal,isMuon), ROOT.RooFit.Import(mytree), ROOT.RooFit.WeightVar("weight"))

if restricted_CR:
    if runningEra == 0:
        data = data_initial.reduce("(BDT_out > 0.155 && isMuon==isMuon::Muon) || (BDT_out > 0.107 && isMuon==isMuon::Electron)")
    if runningEra == 1:
        data = data_initial.reduce("(BDT_out > 0.184 && isMuon==isMuon::Muon) || (BDT_out > 0.122 && isMuon==isMuon::Electron)")
else:
    data = data_initial

print "number of events mu 2016  - SR: ", data.sumEntries("Categorization==1")
print "number of events ele 2016 - SR: ", data.sumEntries("Categorization==3")
print "number of events mu 2017  - SR: ", data.sumEntries("Categorization==5")
print "number of events ele 2017 - SR: ", data.sumEntries("Categorization==7")
if restricted_CR:
    print "number of events mu 2016  - CR: ", data.sumEntries("Categorization==0")
    print "number of events ele 2016 - CR: ", data.sumEntries("Categorization==2")
    print "number of events mu 2017  - CR: ", data.sumEntries("Categorization==4")
    print "number of events ele 2017 - CR: ", data.sumEntries("Categorization==6")

print "Using ", data.numEntries(), " events to fit"

################################################################
#                                                              #
#--------------- Import variables from workspace --------------#
#                                                              #
################################################################

#Import the signal
fInput_sigmodel = ROOT.TFile("Signal_model.root")
fInput_sigmodel.cd()

workspace = fInput_sigmodel.Get("myworkspace")

#Fix the signal parametrization
workspace.var("dCB_pole").setConstant(1)
workspace.var("dCB_aL").setConstant(1)
workspace.var("dCB_nL").setConstant(1)
workspace.var("dCB_aR").setConstant(1)
workspace.var("dCB_nR").setConstant(1)
workspace.var("Gauss_pole_2").setConstant(1)
workspace.var("Gauss_sigma_2").setConstant(1)
workspace.var("fracSig").setConstant(1)

totSignal = workspace.pdf("totSignal")


################################################################
#                                                              #
#---------------- Variables for bkg description ---------------#
#                                                              #
################################################################

#Now describe the background

fInput_bkg = ROOT.TFile("fitBackground.root")
fInput_bkg.cd()

workspace_bkg = fInput_bkg.Get("workspace_bkg")

if useChebychev == 0:
    backPDF_mu_2016 = workspace_bkg.pdf("backPDF_cheb_mu_2016")
    backPDF_el_2016 = workspace_bkg.pdf("backPDF_cheb_el_2016")
    backPDF_mu_2017 = workspace_bkg.pdf("backPDF_cheb_mu_2017")
    backPDF_el_2017 = workspace_bkg.pdf("backPDF_cheb_el_2017")

elif useChebychev == 1: #Use Bernstein for muon channel 2016 (calculation of systematic on background parametrization)
    # workspace_bkg.var("b0_mu_2016").setRange(0.,0.005)
    # workspace_bkg.var("b0_mu_2016").setVal(0.0005)
    backPDF_mu_2016 = workspace_bkg.pdf("backPDF_bern_mu_2016")
    backPDF_el_2016 = workspace_bkg.pdf("backPDF_cheb_el_2016")
    backPDF_mu_2017 = workspace_bkg.pdf("backPDF_cheb_mu_2017")
    backPDF_el_2017 = workspace_bkg.pdf("backPDF_cheb_el_2017")

elif useChebychev == 2: #Use Bernstein for electron channel 2016 (calculation of systematic on background parametrization)
    backPDF_mu_2016 = workspace_bkg.pdf("backPDF_cheb_mu_2016")
    backPDF_el_2016 = workspace_bkg.pdf("backPDF_bern_el_2016")
    backPDF_mu_2017 = workspace_bkg.pdf("backPDF_cheb_mu_2017")
    backPDF_el_2017 = workspace_bkg.pdf("backPDF_cheb_el_2017")

elif useChebychev == 3: #Use Bernstein for muon channel 2017 (calculation of systematic on background parametrization)
    backPDF_mu_2016 = workspace_bkg.pdf("backPDF_cheb_mu_2016")
    backPDF_el_2016 = workspace_bkg.pdf("backPDF_cheb_el_2016")
    backPDF_mu_2017 = workspace_bkg.pdf("backPDF_bern_mu_2017")
    backPDF_el_2017 = workspace_bkg.pdf("backPDF_cheb_el_2017")

elif useChebychev == 4: #Use Bernstein for electron channel 2017 (calculation of systematic on background parametrization)
    backPDF_mu_2016 = workspace_bkg.pdf("backPDF_cheb_mu_2016")
    backPDF_el_2016 = workspace_bkg.pdf("backPDF_cheb_el_2016")
    backPDF_mu_2017 = workspace_bkg.pdf("backPDF_cheb_mu_2017")
    backPDF_el_2017 = workspace_bkg.pdf("backPDF_bern_el_2017")

################################################################
#                                                              #
#-------------- Systematic on W width resolution --------------#
#                                                              #
################################################################

#Gaussian distribution of W resolution width for systematics
# W_resol_width = workspace.var("W_resol_width")
dCB_width = workspace.var("dCB_width")
dCB_width_constr = ROOT.RooRealVar("dCB_width_constr","dCB_width_constr",dCB_width.getVal())
dCB_width_err = ROOT.RooRealVar("dCB_width_err","dCB_width_err",dCB_width.getError())
gauss_W_resol = ROOT.RooGaussian("gauss_W_resol","gauss_W_resol",dCB_width,dCB_width_constr,dCB_width_err)


################################################################
#                                                              #
#------------------ Systematic on ttbar xsec ------------------#
#                                                              #
################################################################

#First the cross section, with a modifier for systematics
#CMS ttbar measurement/W->lnu BR (it is measured with both W in lnu), in pb
#http://cms-results.web.cern.ch/cms-results/public-results/publications/TOP-16-005/index.html
glb_W_xsec    = ROOT.RooRealVar("glb_W_xsec","glb_W_xsec", 2.*831.76*0.1086, 0., 1000.)
W_xsec_constr = ROOT.RooRealVar("W_xsec_constr","W_x_sec_constr", 2.*831.76*0.1086, 0., 1000.)
W_xsec_syst   = ROOT.RooRealVar("W_xsec_syst","W_xsec_syst",43.*2.*0.1086)
gauss_W_xsec  = ROOT.RooGaussian("gauss_W_xsec","gauss_W_xsec",glb_W_xsec,W_xsec_constr,W_xsec_syst)


################################################################
#                                                              #
#------------------ Systematic on luminosity ------------------#
#                                                              #
################################################################

#Represent the luminosity with a modifier for systematics. For 2016 (35.86 fb-1): 2.5% systematic. For 2017 (41.529 fb-1): 2.3% systematic

glb_lumi_2016    = ROOT.RooRealVar("glb_lumi_2016","glb_lumi_2016",35.86*1000., 0., 50000.) #In pb
lumi_constr_2016 = ROOT.RooRealVar("lumi_constr_2016","lumi_constr_2016", 35.86*1000., 0., 50000.)
lumi_syst_2016   = ROOT.RooRealVar("lumi_syst_2016","lumi_syst_2016", 35.86*0.025*1000.)
gauss_lumi_2016  = ROOT.RooGaussian("gauss_lumi_2016","gauss_lumi_2016",glb_lumi_2016,lumi_constr_2016,lumi_syst_2016) 

glb_lumi_2017    = ROOT.RooRealVar("glb_lumi_2017","glb_lumi_2017",41.529*1000., 0., 50000.) #In pb
lumi_constr_2017 = ROOT.RooRealVar("lumi_constr_2017","lumi_constr_2017", 41.529*1000., 0., 50000.)
lumi_syst_2017   = ROOT.RooRealVar("lumi_syst_2017","lumi_syst_2017", 41.529*0.023*1000.)
gauss_lumi_2017  = ROOT.RooGaussian("gauss_lumi_2017","gauss_lumi_2017",glb_lumi_2017,lumi_constr_2017,lumi_syst_2017) 


################################################################
#                                                              #
#------------------ Systematic on efficiency ------------------#
#                                                              #
################################################################

totsig_2016 = 100000. #total number of signal events in 2016
totmu_2016  = 7114.   #total number of signal muon events in 2016
totel_2016  = 5036.   #total number of signal electron events in 2016

totsig_2017 = 80000.  #total number of signal events in 2017
totmu_2017  = 4787.   #total number of signal muon events in 2017
totel_2017  = 3854.   #total number of signal electron events in 2017

glb_eff_mu_2016    = ROOT.RooRealVar("glb_eff_mu_2016","glb_eff_mu_2016",totmu_2016*2./totsig_2016, 0., 1.)
eff_mu_constr_2016 = ROOT.RooRealVar("eff_mu_constr_2016","eff_mu_constr_2016", totmu_2016*2./totsig_2016, 0., 1.)
eff_mu_syst_2016   = ROOT.RooRealVar("eff_mu_syst_2016","eff_mu_syst_2016", 4*totmu_2016*(totsig_2016-2*totmu_2016)/(totsig_2016*totsig_2016*totsig_2016))
gauss_eff_mu_2016  = ROOT.RooGaussian("gauss_eff_mu_2016","gauss_eff_mu_2016",glb_eff_mu_2016,eff_mu_constr_2016,eff_mu_syst_2016) 

glb_eff_el_2016    = ROOT.RooRealVar("glb_eff_el_2016","glb_eff_el_2016", totel_2016*2./totsig_2016, 0., 1.)
eff_el_constr_2016 = ROOT.RooRealVar("eff_el_constr_2016","eff_el_constr_2016",totel_2016*2./totsig_2016, 0., 1.)
eff_el_syst_2016   = ROOT.RooRealVar("eff_el_syst_2016","eff_el_syst_2016", 4*totel_2016*(totsig_2016-2*totel_2016)/(totsig_2016*totsig_2016*totsig_2016))
gauss_eff_el_2016  = ROOT.RooGaussian("gauss_eff_el_2016","gauss_eff_el_2016",glb_eff_el_2016,eff_el_constr_2016,eff_el_syst_2016) 

glb_eff_mu_2017    = ROOT.RooRealVar("glb_eff_mu_2017","glb_eff_mu_2017",totmu_2017*2./totsig_2017, 0., 1.)
eff_mu_constr_2017 = ROOT.RooRealVar("eff_mu_constr_2017","eff_mu_constr_2017", totmu_2017*2./totsig_2017, 0., 1.)
eff_mu_syst_2017   = ROOT.RooRealVar("eff_mu_syst_2017","eff_mu_syst_2017", 4*totmu_2017*(totsig_2017-2*totmu_2017)/(totsig_2017*totsig_2017*totsig_2017))
gauss_eff_mu_2017  = ROOT.RooGaussian("gauss_eff_mu_2017","gauss_eff_mu_2017",glb_eff_mu_2017,eff_mu_constr_2017,eff_mu_syst_2017) 

glb_eff_el_2017    = ROOT.RooRealVar("glb_eff_el_2017","glb_eff_el_2017", totel_2017*2./totsig_2017, 0., 1.)
eff_el_constr_2017 = ROOT.RooRealVar("eff_el_constr_2017","eff_el_constr_2017",totel_2017*2./totsig_2017, 0., 1.)
eff_el_syst_2017   = ROOT.RooRealVar("eff_el_syst_2017","eff_el_syst_2017", 4*totel_2017*(totsig_2017-2*totel_2017)/(totsig_2017*totsig_2017*totsig_2017))
gauss_eff_el_2017  = ROOT.RooGaussian("gauss_eff_el_2017","gauss_eff_el_2017",glb_eff_el_2017,eff_el_constr_2017,eff_el_syst_2017) 


################################################################
#                                                              #
#------------- Systematic on bkg parametrization --------------#
#                                                              #
################################################################

eta_mu_2016 = ROOT.RooRealVar("eta_mu_2016","eta_mu_2016", 1.,0.0001,3.)
eta_el_2016 = ROOT.RooRealVar("eta_el_2016","eta_el_2016", 1.,0.0001,3.)
eta_mu_2017 = ROOT.RooRealVar("eta_mu_2017","eta_mu_2017", 1.,0.0001,3.)
eta_el_2017 = ROOT.RooRealVar("eta_el_2017","eta_el_2017", 1.,0.0001,3.)

if not isAlternativeBkgDescription:
    bkg_syst_mu_2016 = 0.004 #Value of the systematic to use in the fit of the signal regions
    bkg_syst_el_2016 = 0.002
    bkg_syst_mu_2017 = 0.001
    bkg_syst_el_2017 = 0.001
else:
    bkg_syst_mu_2016 = 0.0001 #In the AlternativeBkgDescription mode, one is supposed to not know yet this systematic. 0.0001 is a custom small value
    bkg_syst_el_2016 = 0.0001
    bkg_syst_mu_2017 = 0.0001
    bkg_syst_el_2017 = 0.0001

glb_bkg_param_mu_2016    = ROOT.RooRealVar("glb_bkg_param_mu_2016","glb_bkg_param_mu_2016", 1., 0.0001, 3.)
bkg_param_syst_mu_2016   = ROOT.RooRealVar("bkg_param_syst_mu_2016","bkg_param_syst_mu_2016",bkg_syst_mu_2016)
gauss_bkg_param_mu_2016  = ROOT.RooGaussian("gauss_bkg_param_mu_2016","gauss_bkg_param_2016_mu",glb_bkg_param_mu_2016,eta_mu_2016,bkg_param_syst_mu_2016)

glb_bkg_param_el_2016    = ROOT.RooRealVar("glb_bkg_param_el_2016","glb_bkg_param_el_2016", 1., 0.0001, 3.)
bkg_param_syst_el_2016   = ROOT.RooRealVar("bkg_param_syst_el_2016","bkg_param_syst_el_2016",bkg_syst_el_2016)
gauss_bkg_param_el_2016  = ROOT.RooGaussian("gauss_bkg_param_el_2016","gauss_bkg_param_2016_el",glb_bkg_param_el_2016,eta_el_2016,bkg_param_syst_el_2016)

glb_bkg_param_mu_2017    = ROOT.RooRealVar("glb_bkg_param_mu_2017","glb_bkg_param_mu_2017", 1., 0.0001, 3.)
bkg_param_syst_mu_2017   = ROOT.RooRealVar("bkg_param_syst_mu_2017","bkg_param_syst_mu_2017",bkg_syst_mu_2017)
gauss_bkg_param_mu_2017  = ROOT.RooGaussian("gauss_bkg_param_mu_2017","gauss_bkg_param_2017_mu",glb_bkg_param_mu_2017,eta_mu_2017,bkg_param_syst_mu_2017)

glb_bkg_param_el_2017    = ROOT.RooRealVar("glb_bkg_param_el_2017","glb_bkg_param_el_2017", 1., 0.0001, 3.)
bkg_param_syst_el_2017   = ROOT.RooRealVar("bkg_param_syst_el_2017","bkg_param_syst_el_2017",bkg_syst_el_2017)
gauss_bkg_param_el_2017  = ROOT.RooGaussian("gauss_bkg_param_el_2017","gauss_bkg_param_2017_el",glb_bkg_param_el_2017,eta_el_2017,bkg_param_syst_el_2017)

################################################################
#                                                              #
#----------------------- Define the POI -----------------------#
#                                                              #
################################################################

if isAlternativeBkgDescription:
    W_pigamma_BR = ROOT.RooRealVar("W_pigamma_BR","W_pigamma_BR",0.000001,-0.1,0.01) # The parameter of interest can go negative when we try the alternative bkg description
else:
    W_pigamma_BR = ROOT.RooRealVar("W_pigamma_BR","W_pigamma_BR",0.000001,0.,0.01) # The parameter of interest

W_pigamma_BR_blind = ROOT.RooUnblindOffset("W_pigamma_BR_blind","W_pigamma_BR_blind","aSeedString",0.000001,W_pigamma_BR)


################################################################
#                                                              #
#------------------ Set a few things constant -----------------#
#                                                              #
################################################################

glb_W_xsec.setConstant(1)
glb_bkg_param_mu_2016.setConstant(1)
glb_bkg_param_el_2016.setConstant(1)
glb_bkg_param_mu_2017.setConstant(1)
glb_bkg_param_el_2017.setConstant(1)
glb_lumi_2016.setConstant(1)
glb_lumi_2017.setConstant(1)
glb_eff_mu_2016.setConstant(1)
glb_eff_el_2016.setConstant(1)
glb_eff_mu_2017.setConstant(1)
glb_eff_el_2017.setConstant(1)


################################################################
#                                                              #
#------------------- Define PDFs for the fit ------------------#
#                                                              #
################################################################

#This is for the signal region
# Nsig_mu = ROOT.RooFormulaVar("Nsig_mu","@0*@1*@2*@3*@4", ROOT.RooArgList(W_pigamma_BR, W_xsec_constr,lumi_constr,eff_mu_constr,eta))
# Nsig_el = ROOT.RooFormulaVar("Nsig_el","@0*@1*@2*@3*@4", ROOT.RooArgList(W_pigamma_BR, W_xsec_constr,lumi_constr,eff_el_constr,eta))
Nsig_mu_2016 = ROOT.RooFormulaVar("Nsig_mu_2016","@0*@1*@2*@3*@4", ROOT.RooArgList(W_pigamma_BR_blind, W_xsec_constr, lumi_constr_2016, eff_mu_constr_2016, eta_mu_2016))
Nsig_el_2016 = ROOT.RooFormulaVar("Nsig_el_2016","@0*@1*@2*@3*@4", ROOT.RooArgList(W_pigamma_BR_blind, W_xsec_constr, lumi_constr_2016, eff_el_constr_2016, eta_el_2016))

Nsig_mu_2017 = ROOT.RooFormulaVar("Nsig_mu_2017","@0*@1*@2*@3*@4", ROOT.RooArgList(W_pigamma_BR_blind, W_xsec_constr, lumi_constr_2017, eff_mu_constr_2017, eta_mu_2017))
Nsig_el_2017 = ROOT.RooFormulaVar("Nsig_el_2017","@0*@1*@2*@3*@4", ROOT.RooArgList(W_pigamma_BR_blind, W_xsec_constr, lumi_constr_2017, eff_el_constr_2017, eta_el_2017))


Nbkg_mu_2016 = ROOT.RooRealVar("Nbkg_mu_2016","Nbkg_mu_2016",150.,1.,1000.)
Nbkg_el_2016 = ROOT.RooRealVar("Nbkg_el_2016","Nbkg_el_2016",150.,1.,1000.)

Nbkg_mu_2017 = ROOT.RooRealVar("Nbkg_mu_2017","Nbkg_mu_2017",150.,1.,1000.)
Nbkg_el_2017 = ROOT.RooRealVar("Nbkg_el_2017","Nbkg_el_2017",150.,1.,1000.)


totPDF_mu_unconstr_2016 = ROOT.RooAddPdf("totPDF_mu_unconstr_2016","Total PDF for the mu channel (2016)",ROOT.RooArgList(totSignal,backPDF_mu_2016),ROOT.RooArgList(Nsig_mu_2016,Nbkg_mu_2016))
totPDF_el_unconstr_2016 = ROOT.RooAddPdf("totPDF_el_unconstr_2016","Total PDF for the el channel (2016)",ROOT.RooArgList(totSignal,backPDF_el_2016),ROOT.RooArgList(Nsig_el_2016,Nbkg_el_2016))

totPDF_mu_unconstr_2017 = ROOT.RooAddPdf("totPDF_mu_unconstr_2017","Total PDF for the mu channel (2017)",ROOT.RooArgList(totSignal,backPDF_mu_2017),ROOT.RooArgList(Nsig_mu_2017,Nbkg_mu_2017))
totPDF_el_unconstr_2017 = ROOT.RooAddPdf("totPDF_el_unconstr_2017","Total PDF for the el channel (2017)",ROOT.RooArgList(totSignal,backPDF_el_2017),ROOT.RooArgList(Nsig_el_2017,Nbkg_el_2017))


totPDF_mu_2016 = ROOT.RooProdPdf("totPDF_mu_2016","totPDF_mu_2016", ROOT.RooArgList(totPDF_mu_unconstr_2016,gauss_lumi_2016,gauss_W_xsec,gauss_eff_mu_2016,gauss_W_resol,gauss_bkg_param_mu_2016))
totPDF_el_2016 = ROOT.RooProdPdf("totPDF_el_2016","totPDF_el_2016", ROOT.RooArgList(totPDF_el_unconstr_2016,gauss_lumi_2016,gauss_W_xsec,gauss_eff_el_2016,gauss_W_resol,gauss_bkg_param_el_2016))

totPDF_mu_2017 = ROOT.RooProdPdf("totPDF_mu_2017","totPDF_mu_2017", ROOT.RooArgList(totPDF_mu_unconstr_2017,gauss_lumi_2017,gauss_W_xsec,gauss_eff_mu_2017,gauss_W_resol,gauss_bkg_param_mu_2017))
totPDF_el_2017 = ROOT.RooProdPdf("totPDF_el_2017","totPDF_el_2017", ROOT.RooArgList(totPDF_el_unconstr_2017,gauss_lumi_2017,gauss_W_xsec,gauss_eff_el_2017,gauss_W_resol,gauss_bkg_param_el_2017))


################################################################
#                                                              #
#------------- Create the global simultaneous PDF -------------#
#                                                              #
################################################################

totPDF = ROOT.RooSimultaneous("totPDF","The total PDF",Categorization)
constrained_params = ROOT.RooArgSet()

if runningEra == 0 and restricted_CR: #Fit on 2016 restricted CR with background-only PDF (background description)
    totPDF.addPdf(backPDF_mu_2016,"MuonCR_2016")
    totPDF.addPdf(backPDF_el_2016,"ElectronCR_2016")
if runningEra == 0 and not restricted_CR: #Fit on 2016 signal region with background+signal PDF
    totPDF.addPdf(totPDF_mu_2016,"MuonSignal_2016")
    totPDF.addPdf(totPDF_el_2016,"ElectronSignal_2016")
    constrained_params.add(dCB_width)
    constrained_params.add(eta_mu_2016)
    constrained_params.add(eta_el_2016)
    constrained_params.add(W_xsec_constr)
    constrained_params.add(lumi_constr_2016)
    constrained_params.add(eff_mu_constr_2016)
    constrained_params.add(eff_el_constr_2016)

if runningEra == 1 and restricted_CR: #Fit on 2017 restricted CR with background-only PDF (background description)
    totPDF.addPdf(backPDF_mu_2017,"MuonCR_2017")
    totPDF.addPdf(backPDF_el_2017,"ElectronCR_2017")
if runningEra == 1 and not restricted_CR: #Fit on 2017 signal region with background+signal PDF
    totPDF.addPdf(totPDF_mu_2017,"MuonSignal_2017")
    totPDF.addPdf(totPDF_el_2017,"ElectronSignal_2017")
    constrained_params.add(dCB_width)
    constrained_params.add(eta_mu_2017)
    constrained_params.add(eta_el_2017)
    constrained_params.add(W_xsec_constr)
    constrained_params.add(lumi_constr_2017)
    constrained_params.add(eff_mu_constr_2017)
    constrained_params.add(eff_el_constr_2017)

if runningEra == 2: #Fit on 2016+2017 signal regions
    totPDF.addPdf(totPDF_mu_2016,"MuonSignal_2016")
    totPDF.addPdf(totPDF_el_2016,"ElectronSignal_2016")
    totPDF.addPdf(totPDF_mu_2017,"MuonSignal_2017")
    totPDF.addPdf(totPDF_el_2017,"ElectronSignal_2017")
    constrained_params.add(dCB_width)
    constrained_params.add(eta_mu_2016)
    constrained_params.add(eta_el_2016)
    constrained_params.add(eta_mu_2017)
    constrained_params.add(eta_el_2017)
    constrained_params.add(W_xsec_constr)
    constrained_params.add(lumi_constr_2016)
    constrained_params.add(lumi_constr_2017)
    constrained_params.add(eff_mu_constr_2016)
    constrained_params.add(eff_el_constr_2016)
    constrained_params.add(eff_mu_constr_2017)
    constrained_params.add(eff_el_constr_2017)


################################################################
#                                                              #
#---------------------- Fit (and F-Test) ----------------------#
#                                                              #
################################################################

if isData:
    if restricted_CR:
        result_dataFit = totPDF.fitTo(data,ROOT.RooFit.Extended(0), ROOT.RooFit.NumCPU(4), ROOT.RooFit.Constrain(constrained_params), ROOT.RooFit.Save() )
    else:
        result_dataFit = totPDF.fitTo(data,ROOT.RooFit.Extended(1), ROOT.RooFit.NumCPU(4), ROOT.RooFit.Constrain(constrained_params), ROOT.RooFit.Save() )#For the signal region, I want the fit to be extended (Poisson fluctuation of unmber of events) to take into account that the total number of events is the sum of signal and background events. Either I do this, or I use a fraction frac*Nbkg+(1-frac)*Nsig, which will become a parameter of the fit and will have a Gaussian behavior (whilst the extended fit preserves the natural Poisson behavior)

    print "minNll = ", result_dataFit.minNll()
    print "2Delta_minNll = ", 2*(3250.96600396-result_dataFit.minNll()) # If 2*(NLL(N)-NLL(N+1)) > 3.85 -> N+1 is significant improvement
else:
    totPDF.fitTo(data,ROOT.RooFit.Extended(1), ROOT.RooFit.SumW2Error(0), ROOT.RooFit.NumCPU(4), ROOT.RooFit.Constrain(constrained_params) )


################################################################
#                                                              #
#-------------------------- Plotting --------------------------#
#                                                              #
################################################################

#First plot the signal region

xframe_mu_2016 = Wmass.frame(55.,95.,15)
xframe_mu_2016.SetTitle(" ")
xframe_mu_2016.SetTitleOffset(1.4,"y")
xframe_mu_2016.SetMaximum(60)

xframe_el_2016 = Wmass.frame(55.,95.,15)
xframe_el_2016.SetTitle(" ")
xframe_el_2016.SetTitleOffset(1.4,"y")
xframe_el_2016.SetMaximum(60)

xframe_mu_2017 = Wmass.frame(55.,95.,15)
xframe_mu_2017.SetTitle(" ")
xframe_mu_2017.SetTitleOffset(1.4,"y")
xframe_mu_2017.SetMaximum(60)

xframe_el_2017 = Wmass.frame(55.,95.,15)
xframe_el_2017.SetTitle(" ")
xframe_el_2017.SetTitleOffset(1.4,"y")
xframe_el_2017.SetMaximum(60)


if runningEra == 0 and restricted_CR:
    data.plotOn(xframe_mu_2016, ROOT.RooFit.Cut("Categorization==0"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    data.plotOn(xframe_el_2016, ROOT.RooFit.Cut("Categorization==2"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    totPDF.plotOn(xframe_mu_2016, ROOT.RooFit.Slice(Categorization,"MuonCR_2016"), ROOT.RooFit.ProjWData(data))
    totPDF.plotOn(xframe_el_2016, ROOT.RooFit.Slice(Categorization,"ElectronCR_2016"), ROOT.RooFit.ProjWData(data))

if runningEra == 1 and restricted_CR:
    data.plotOn(xframe_mu_2017, ROOT.RooFit.Cut("Categorization==4"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    data.plotOn(xframe_el_2017, ROOT.RooFit.Cut("Categorization==6"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    totPDF.plotOn(xframe_mu_2017, ROOT.RooFit.Slice(Categorization,"MuonCR_2017"), ROOT.RooFit.ProjWData(data))  
    totPDF.plotOn(xframe_el_2017, ROOT.RooFit.Slice(Categorization,"ElectronCR_2017"), ROOT.RooFit.ProjWData(data))  

#################################################

if runningEra == 0 and not restricted_CR:
    data_reduced = data.reduce("Wmass < 65. || Wmass > 90.")
    data_reduced.plotOn(xframe_mu_2016, ROOT.RooFit.Cut("Categorization==1"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    totPDF.plotOn(xframe_mu_2016, ROOT.RooFit.Range("LowSideband,HighSideband"), ROOT.RooFit.Slice(Categorization,"MuonSignal_2016"), ROOT.RooFit.ProjWData(data))
    data_reduced.plotOn(xframe_el_2016, ROOT.RooFit.Cut("Categorization==3"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    totPDF.plotOn(xframe_el_2016, ROOT.RooFit.Range("LowSideband,HighSideband"), ROOT.RooFit.Slice(Categorization,"ElectronSignal_2016"), ROOT.RooFit.ProjWData(data))
if runningEra == 1 and not restricted_CR:
    data_reduced = data.reduce("Wmass < 65. || Wmass > 90.")
    data_reduced.plotOn(xframe_mu_2017, ROOT.RooFit.Cut("Categorization==5"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    totPDF.plotOn(xframe_mu_2017, ROOT.RooFit.Range("LowSideband,HighSideband"), ROOT.RooFit.Slice(Categorization,"MuonSignal_2017"), ROOT.RooFit.ProjWData(data))
    data_reduced.plotOn(xframe_el_2017, ROOT.RooFit.Cut("Categorization==7"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    totPDF.plotOn(xframe_el_2017, ROOT.RooFit.Range("LowSideband,HighSideband"), ROOT.RooFit.Slice(Categorization,"ElectronSignal_2017"), ROOT.RooFit.ProjWData(data))

#################################################

if runningEra == 2:
    #Exclude the control regions
    data_reduced = data.reduce("Wmass < 65. || Wmass > 90.")
    data_reduced.plotOn(xframe_mu_2016, ROOT.RooFit.Cut("Categorization==1"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    totPDF.plotOn(xframe_mu_2016, ROOT.RooFit.Range("LowSideband,HighSideband"), ROOT.RooFit.Slice(Categorization,"MuonSignal_2016"), ROOT.RooFit.ProjWData(data))
    data_reduced.plotOn(xframe_el_2016, ROOT.RooFit.Cut("Categorization==3"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    totPDF.plotOn(xframe_el_2016, ROOT.RooFit.Range("LowSideband,HighSideband"), ROOT.RooFit.Slice(Categorization,"ElectronSignal_2016"), ROOT.RooFit.ProjWData(data))
    data_reduced.plotOn(xframe_mu_2017, ROOT.RooFit.Cut("Categorization==5"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    totPDF.plotOn(xframe_mu_2017, ROOT.RooFit.Range("LowSideband,HighSideband"), ROOT.RooFit.Slice(Categorization,"MuonSignal_2017"), ROOT.RooFit.ProjWData(data))
    data_reduced.plotOn(xframe_el_2017, ROOT.RooFit.Cut("Categorization==7"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    totPDF.plotOn(xframe_el_2017, ROOT.RooFit.Range("LowSideband,HighSideband"), ROOT.RooFit.Slice(Categorization,"ElectronSignal_2017"), ROOT.RooFit.ProjWData(data))

    # #Do not exclude the control regions
    # data.plotOn(xframe_mu_2016, ROOT.RooFit.Cut("Categorization==1"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    # totPDF.plotOn(xframe_mu_2016, ROOT.RooFit.Slice(Categorization,"MuonSignal_2016"), ROOT.RooFit.ProjWData(data))
    # data.plotOn(xframe_el_2016, ROOT.RooFit.Cut("Categorization==3"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    # totPDF.plotOn(xframe_el_2016, ROOT.RooFit.Slice(Categorization,"ElectronSignal_2016"), ROOT.RooFit.ProjWData(data))
    # data.plotOn(xframe_mu_2017, ROOT.RooFit.Cut("Categorization==5"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    # totPDF.plotOn(xframe_mu_2017, ROOT.RooFit.Slice(Categorization,"MuonSignal_2017"), ROOT.RooFit.ProjWData(data))
    # data.plotOn(xframe_el_2017, ROOT.RooFit.Cut("Categorization==7"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
    # totPDF.plotOn(xframe_el_2017, ROOT.RooFit.Slice(Categorization,"ElectronSignal_2017"), ROOT.RooFit.ProjWData(data))

    
canvas = ROOT.TCanvas()
canvas.Divide(2,2)
canvas.cd(1)
xframe_mu_2016.Draw()
canvas.cd(2)
xframe_el_2016.Draw()
canvas.cd(3)
xframe_mu_2017.Draw()
canvas.cd(4)
xframe_el_2017.Draw()

# Save the plot
if isData:
    canvas.SaveAs("plots/fitData_signalR.pdf")
else:
 canvas.SaveAs("plots/fitMC_signalR.pdf")

# #Now plot the CR
# xframe_mu_CR = Wmass.frame(55.,95.,15)
# xframe_mu_CR.SetTitle(" ")
# xframe_mu_CR.SetTitleOffset(1.4,"y")
# data.plotOn(xframe_mu_CR, ROOT.RooFit.Cut("Categorization==0"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
# totPDF.plotOn(xframe_mu_CR, ROOT.RooFit.Slice(Categorization,"MuonCR"), ROOT.RooFit.ProjWData(data))

# xframe_el_CR = Wmass.frame(55.,95.,15)
# xframe_el_CR.SetTitle(" ")
# xframe_el_CR.SetTitleOffset(1.4,"y")
# data.plotOn(xframe_el_CR, ROOT.RooFit.Cut("Categorization==2"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
# totPDF.plotOn(xframe_el_CR, ROOT.RooFit.Slice(Categorization,"ElectronCR"), ROOT.RooFit.ProjWData(data))

# canvas_CR = ROOT.TCanvas()
# canvas_CR.Divide(2,1)
# canvas_CR.cd(1)
# xframe_mu_CR.Draw()
# canvas_CR.cd(2)
# xframe_el_CR.Draw()

# if isData:
#     canvas_CR.SaveAs("plots/fitData_CR_restricted.pdf")
# else:
#     canvas_CR.SaveAs("plots/fitMC_CR.pdf")


################################################################
#                                                              #
#------------------------ Write to file -----------------------#
#                                                              #
################################################################
        
#Save the fit into a root file
if isData:
    fOutput = ROOT.TFile("fitData.root","RECREATE")
else:
    fOutput = ROOT.TFile("fitMC.root","RECREATE")

fOutput.cd()

workspace_out = ROOT.RooWorkspace("workspace")
getattr(workspace_out,'import')(data)
getattr(workspace_out,'import')(totPDF)

workspace_out.Write()

fOutput.Close()

raw_input()
