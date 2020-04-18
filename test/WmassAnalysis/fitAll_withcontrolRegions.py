#This program fits both lepton samples at the same time

import ROOT
import math
import argparse

ROOT.gROOT.ProcessLineSync(".L dCB/RooDoubleCBFast.cc+")

#---------------------------------#
p = argparse.ArgumentParser(description='Select whether to fit data or MC')
p.add_argument('isData_option', help='Type <<data>> or <<MC>>')
#p.add_argument('runningEra_option', help='Type <<0>> for 2016, <<1>> for 2017, <<2>> for 2016+2017, <<3>> for 2016+2017+2018')
args = p.parse_args()

isData = False
# Switch from MC to data channel
if args.isData_option == "data":
    isData = True

#runningEra = int(args.runningEra_option)
#---------------------------------#

################################################################
#                                                              #
#------------------------ Instructions ------------------------#
#                                                              #
################################################################

selectBkgFunction = 0 # 0: use Chebychev for ALL the bkg PDFs. 1: use alternative PDF for muon channel. 2: use alternative PDF for electron channel
suppressBkgSystematic = False # To be used when trying to fit with alternative bkg description, in order to estimate a systematic. If True, it will allow W_pigamma_BR to float negative. Moreover, it will use Signal+Background in the totPDF, so that the fit to the restricted CRs will contain also the POI BR, which will be used to calculate the pull and hence to estimate the systematic
suppressAllSystematics = False


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
    fInput = ROOT.TFile("Tree_input_massfit_Data_3.root")
else:
    fInput = ROOT.TFile("Tree_input_massfit_MC_3.root")

fInput.cd()

mytree = fInput.Get("minitree")

################################################################
#                                                              #
#-------------------- Define the categories -------------------#
#                                                              #
################################################################

#Define the mu/ele category for 2016 and 2017
Categorization = ROOT.RooCategory("Categorization","Categorization")
Categorization.defineType("MuonSignal",1)
Categorization.defineType("ElectronSignal",3)

isSignal = ROOT.RooCategory("isSignal","isSignal")
isSignal.defineType("Signal",1)
isSignal.defineType("Background",0)

################################################################
#                                                              #
#--------------- Get event weight and BDT output --------------#
#                                                              #
################################################################

minWeight_inMC = -50. # Set the minimum weight an event can have when processing MC. This is useful to remove spikes which cannot be fitted
maxWeight_inMC = 50.  # Set the maximum weight an event can have when processing MC. This is useful to remove spikes which cannot be fitted

#Define the event weight
if not isData:
    weight = ROOT.RooRealVar("weight","The event weight",minWeight_inMC,maxWeight_inMC)

#Import dataset
if isData:
    data_initial = ROOT.RooDataSet("data","data", ROOT.RooArgSet(Wmass,Categorization), ROOT.RooFit.Import(mytree))
else:
    data_initial = ROOT.RooDataSet("data","data", ROOT.RooArgSet(Wmass,Categorization,weight,isSignal), ROOT.RooFit.Import(mytree), ROOT.RooFit.WeightVar("weight"))


data = data_initial.reduce("(Categorization==Categorization::MuonSignal) || (Categorization==Categorization::ElectronSignal)")
print "number of events mu - SR: ", data.sumEntries("Categorization==1")
print "number of events ele - SR: ", data.sumEntries("Categorization==3")

print "Using ", data.numEntries(), " events to fit"

################################################################
#                                                              #
#--------------- Import variables from workspace --------------#
#                                                              #
################################################################

#Import the signal
fInput_sigmodel = ROOT.TFile("Signal_model_3.root")
fInput_sigmodel.cd()

workspace = fInput_sigmodel.Get("myworkspace")

#Fix the signal parametrization
#workspace.var("dCB_pole").setConstant(1)
workspace.var("dCB_aL").setConstant(1)
workspace.var("dCB_nL").setConstant(1)
workspace.var("dCB_aR").setConstant(1)
workspace.var("dCB_nR").setConstant(1)
workspace.var("Gauss_pole").setConstant(1)
workspace.var("Gauss_sigma").setConstant(1)
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

if selectBkgFunction == 0:

    backPDF_mu = workspace_bkg.pdf("backPDF_cheb_mu")
    backPDF_el = workspace_bkg.pdf("backPDF_cheb_el")

elif selectBkgFunction == 1: #Use Exponential for muon channel (calculation of systematic on background parametrization)

    backPDF_mu = workspace_bkg.pdf("backPDF_exp_mu")
    backPDF_el = workspace_bkg.pdf("backPDF_cheb_el")

elif selectBkgFunction == 2: #Use Exponential for electron channel (calculation of systematic on background parametrization)

    backPDF_mu = workspace_bkg.pdf("backPDF_cheb_mu")
    backPDF_el = workspace_bkg.pdf("backPDF_exp_el")


################################################################
#                                                              #
#-------------- Systematic on W pole resolution ---------------#
#                                                              #
################################################################

#Gaussian distribution of W pole resolution for systematics
dCB_pole           = workspace.var("dCB_pole")
dCB_pole_constr    = ROOT.RooRealVar("dCB_pole_constr","dCB_pole_constr",dCB_pole.getVal())
dCB_pole_err       = ROOT.RooRealVar("dCB_pole_err","dCB_pole_err",dCB_pole.getError())
gauss_W_pole_resol = ROOT.RooGaussian("gauss_W_pole_resol","gauss_W_pole_resol",dCB_pole,dCB_pole_constr,dCB_pole_err)

################################################################
#                                                              #
#-------------- Systematic on W width resolution --------------#
#                                                              #
################################################################

#Gaussian distribution of W width resolution for systematics
# W_resol_width = workspace.var("W_resol_width")
dCB_width        = workspace.var("dCB_width")
dCB_width_constr = ROOT.RooRealVar("dCB_width_constr","dCB_width_constr",dCB_width.getVal())
dCB_width_err    = ROOT.RooRealVar("dCB_width_err","dCB_width_err",dCB_width.getError())
gauss_W_resol    = ROOT.RooGaussian("gauss_W_resol","gauss_W_resol",dCB_width,dCB_width_constr,dCB_width_err)

################################################################
#                                                              #
#------------------ Systematic on ttbar xsec ------------------#
#                                                              #
################################################################
#First the cross section, with a modifier for systematics
#CMS ttbar measurement/W->lnu BR (it is measured with both W in lnu), in pb
#http://cms-results.web.cern.ch/cms-results/public-results/publications/TOP-16-005/index.html
W_xsec_nominal  = ROOT.RooRealVar("W_xsec_nominal","W_xsec_nominal", 2.*831.76*0.1086)
W_xsec_syst     = ROOT.RooRealVar("W_xsec_syst","W_xsec_syst",43.*2.*0.1086)
W_xsec_kappa    = ROOT.RooRealVar("W_xsec_kappa","W_xsec_kappa",1.+(W_xsec_syst.getVal()/W_xsec_nominal.getVal()))
W_xsec_beta     = ROOT.RooRealVar("W_xsec_beta","W_xsec_beta",0.,-5.,5.)
W_xsec          = ROOT.RooFormulaVar("W_xsec","@0 * pow(@1,@2)",ROOT.RooArgList(W_xsec_nominal,W_xsec_kappa,W_xsec_beta))
glb_W_xsec      = ROOT.RooRealVar("glb_W_xsec","glb_W_xsec", 0., -5., 5.)
one             = ROOT.RooRealVar("one","one",1.)
gauss_W_xsec    = ROOT.RooGaussian("gauss_W_xsec","gauss_W_xsec",glb_W_xsec,W_xsec_beta,one)

################################################################
#                                                              #
#------------------ Systematic on luminosity ------------------#
#                                                              #
################################################################

#Represent the luminosity with a modifier for systematics. For 2016 (35.86 fb-1): 2.5% systematic. For 2017 (41.529 fb-1): 2.3% systematic

lumi_nominal_2016 = ROOT.RooRealVar("lumi_nominal_2016","lumi_nominal_2016",35.86*1000.) #In pb
lumi_syst_2016    = ROOT.RooRealVar("lumi_syst_2016","lumi_syst_2016", 35.86*0.025*1000.)
lumi_kappa_2016   = ROOT.RooRealVar("lumi_kappa_2016","lumi_kappa_2016",1.+(lumi_syst_2016.getVal()/lumi_nominal_2016.getVal()))
lumi_beta_2016    = ROOT.RooRealVar("lumi_beta_2016","lumi_beta_2016",0.,-5.,5.)
lumi_2016         = ROOT.RooFormulaVar("lumi_2016","@0 * pow(@1,@2)",ROOT.RooArgList(lumi_nominal_2016,lumi_kappa_2016,lumi_beta_2016))
glb_lumi_2016     = ROOT.RooRealVar("glb_lumi_2016","glb_lumi_2016",0.,-5.,5.)
gauss_lumi_2016   = ROOT.RooGaussian("gauss_lumi_2016","gauss_lumi_2016",glb_lumi_2016,lumi_beta_2016,one)


lumi_nominal_2017 = ROOT.RooRealVar("lumi_nominal_2017","lumi_nominal_2017",41.53*1000.) #In pb
lumi_syst_2017    = ROOT.RooRealVar("lumi_syst_2017","lumi_syst_2017", 41.53*0.023*1000.)
lumi_kappa_2017   = ROOT.RooRealVar("lumi_kappa_2017","lumi_kappa_2017",1.+(lumi_syst_2017.getVal()/lumi_nominal_2017.getVal()))
lumi_beta_2017    = ROOT.RooRealVar("lumi_beta_2017","lumi_beta_2017",0.,-5.,5.)
lumi_2017         = ROOT.RooFormulaVar("lumi_2017","@0 * pow(@1,@2)",ROOT.RooArgList(lumi_nominal_2017,lumi_kappa_2017,lumi_beta_2017))
glb_lumi_2017     = ROOT.RooRealVar("glb_lumi_2017","glb_lumi_2017",0.,-5.,5.)
gauss_lumi_2017   = ROOT.RooGaussian("gauss_lumi_2017","gauss_lumi_2017",glb_lumi_2017,lumi_beta_2017,one) 


lumi_nominal_2018 = ROOT.RooRealVar("lumi_nominal_2018","lumi_nominal_2018",59.69*1000.) #In pb
lumi_syst_2018    = ROOT.RooRealVar("lumi_syst_2018","lumi_syst_2018", 59.69*0.025*1000.)
lumi_kappa_2018   = ROOT.RooRealVar("lumi_kappa_2018","lumi_kappa_2018",1.+(lumi_syst_2018.getVal()/lumi_nominal_2018.getVal()))
lumi_beta_2018    = ROOT.RooRealVar("lumi_beta_2018","lumi_beta_2018",0.,-5.,5.)
lumi_2018         = ROOT.RooFormulaVar("lumi_2018","@0 * pow(@1,@2)",ROOT.RooArgList(lumi_nominal_2018,lumi_kappa_2018,lumi_beta_2018))
glb_lumi_2018     = ROOT.RooRealVar("glb_lumi_2018","glb_lumi_2018",0.,-5.,5.)
gauss_lumi_2018   = ROOT.RooGaussian("gauss_lumi_2018","gauss_lumi_2018",glb_lumi_2018,lumi_beta_2018,one) 


################################################################
#                                                              #
#------------------ Systematic on efficiency ------------------#
#                                                              #
################################################################

totsig_2016 = 100000. #total number of signal events in 2016
totmu_2016  = 6046.   #total number of signal muon events in 2016
totel_2016  = 3686.   #total number of signal electron events in 2016

totsig_2017 = 80000.  #total number of signal events in 2017
totmu_2017  = 4310.   #total number of signal muon events in 2017
totel_2017  = 2934.   #total number of signal electron events in 2017

totsig_2018 = 79820.  #total number of signal events in 2018
totmu_2018  = 4373.   #total number of signal muon events in 2018
totel_2018  = 2815.   #total number of signal electron events in 2018

eff_mu_nominal_2016 = ROOT.RooRealVar("eff_mu_nominal_2016","eff_mu_nominal_2016",totmu_2016*2./totsig_2016)
binom_eff_mu_2016   = (4*totmu_2016*(totsig_2016-2*totmu_2016)/(totsig_2016*totsig_2016*totsig_2016))*(4*totmu_2016*(totsig_2016-2*totmu_2016)/(totsig_2016*totsig_2016*totsig_2016))/(eff_mu_nominal_2016.getVal()*eff_mu_nominal_2016.getVal()) #It will be summed in quadrature to the BDT systematic
BDT_syst_mu_2016    = 0.01*0.01 #It will be summed in quadrature to the binomial uncertainty
Pythia_syst_mu_2016 = 0.02*0.02 #It will be summed in quadrature to the binomial uncertainty
eff_mu_syst_2016    = ROOT.RooRealVar("eff_mu_syst_2016","eff_mu_syst_2016", math.sqrt(binom_eff_mu_2016+BDT_syst_mu_2016+Pythia_syst_mu_2016))
eff_mu_kappa_2016   = ROOT.RooRealVar("eff_mu_kappa_2016","lumi_kappa_2016",1.+eff_mu_syst_2016.getVal())
eff_mu_beta_2016    = ROOT.RooRealVar("eff_mu_beta_2016","eff_mu_beta_2016",0.,-5.,5.)
eff_mu_2016         = ROOT.RooFormulaVar("eff_mu_2016","@0 * pow(@1,@2)",ROOT.RooArgList(eff_mu_nominal_2016,eff_mu_kappa_2016,eff_mu_beta_2016))
glb_eff_mu_2016     = ROOT.RooRealVar("glb_eff_mu_2016","glb_eff_mu_2016",0.,-5.,5.)
gauss_eff_mu_2016   = ROOT.RooGaussian("gauss_eff_mu_2016","gauss_eff_mu_2016",glb_eff_mu_2016,eff_mu_beta_2016,one) 


eff_el_nominal_2016 = ROOT.RooRealVar("eff_el_nominal_2016","eff_el_nominal_2016",totel_2016*2./totsig_2016)
binom_eff_el_2016   = (4*totel_2016*(totsig_2016-2*totel_2016)/(totsig_2016*totsig_2016*totsig_2016))*(4*totel_2016*(totsig_2016-2*totel_2016)/(totsig_2016*totsig_2016*totsig_2016))/(eff_el_nominal_2016.getVal()*eff_el_nominal_2016.getVal()) #It will be summed in quadrature to the BDT systematic
BDT_syst_el_2016    = 0.01*0.01 #It will be summed in quadrature to the binomial uncertainty
Pythia_syst_el_2016 = 0.02*0.02 #It will be summed in quadrature to the binomial uncertainty
SF_syst_el_2016     = 0.0125*0.0125 #It will be summed in quadrature to the binomial uncertainty
eff_el_syst_2016    = ROOT.RooRealVar("eff_el_syst_2016","eff_el_syst_2016", math.sqrt(binom_eff_el_2016+BDT_syst_el_2016+Pythia_syst_el_2016+SF_syst_el_2016))
eff_el_kappa_2016   = ROOT.RooRealVar("eff_el_kappa_2016","lumi_kappa_2016",1.+eff_el_syst_2016.getVal())
eff_el_beta_2016    = ROOT.RooRealVar("eff_el_beta_2016","eff_el_beta_2016",0.,-5.,5.)
eff_el_2016         = ROOT.RooFormulaVar("eff_el_2016","@0 * pow(@1,@2)",ROOT.RooArgList(eff_el_nominal_2016,eff_el_kappa_2016,eff_el_beta_2016))
glb_eff_el_2016     = ROOT.RooRealVar("glb_eff_el_2016","glb_eff_el_2016",0.,-5.,5.)
gauss_eff_el_2016   = ROOT.RooGaussian("gauss_eff_el_2016","gauss_eff_el_2016",glb_eff_el_2016,eff_el_beta_2016,one)


eff_mu_nominal_2017 = ROOT.RooRealVar("eff_mu_nominal_2017","eff_mu_nominal_2017",totmu_2017*2./totsig_2017)
binom_eff_mu_2017   = (4*totmu_2017*(totsig_2017-2*totmu_2017)/(totsig_2017*totsig_2017*totsig_2017))*(4*totmu_2017*(totsig_2017-2*totmu_2017)/(totsig_2017*totsig_2017*totsig_2017))/(eff_mu_nominal_2017.getVal()*eff_mu_nominal_2017.getVal()) #It will be summed in quadrature to the BDT systematic
BDT_syst_mu_2017    = 0.01*0.01 #It will be summed in quadrature to the binomial uncertainty
Pythia_syst_mu_2017 = 0.02*0.02 #It will be summed in quadrature to the binomial uncertainty
eff_mu_syst_2017    = ROOT.RooRealVar("eff_mu_syst_2017","eff_mu_syst_2017", math.sqrt(binom_eff_mu_2017+BDT_syst_mu_2017+Pythia_syst_mu_2017))
eff_mu_kappa_2017   = ROOT.RooRealVar("eff_mu_kappa_2017","lumi_kappa_2017",1.+eff_mu_syst_2017.getVal())
eff_mu_beta_2017    = ROOT.RooRealVar("eff_mu_beta_2017","eff_mu_beta_2017",0.,-5.,5.)
eff_mu_2017         = ROOT.RooFormulaVar("eff_mu_2017","@0 * pow(@1,@2)",ROOT.RooArgList(eff_mu_nominal_2017,eff_mu_kappa_2017,eff_mu_beta_2017))
glb_eff_mu_2017     = ROOT.RooRealVar("glb_eff_mu_2017","glb_eff_mu_2017",0.,-5.,5.)
gauss_eff_mu_2017   = ROOT.RooGaussian("gauss_eff_mu_2017","gauss_eff_mu_2017",glb_eff_mu_2017,eff_mu_beta_2017,one)


eff_el_nominal_2017 = ROOT.RooRealVar("eff_el_nominal_2017","eff_el_nominal_2017",totel_2017*2./totsig_2017)
binom_eff_el_2017   = (4*totel_2017*(totsig_2017-2*totel_2017)/(totsig_2017*totsig_2017*totsig_2017))*(4*totel_2017*(totsig_2017-2*totel_2017)/(totsig_2017*totsig_2017*totsig_2017))/(eff_el_nominal_2017.getVal()*eff_el_nominal_2017.getVal()) #It will be summed in quadrature to the BDT systematic
BDT_syst_el_2017    = 0.01*0.01 #It will be summed in quadrature to the binomial uncertainty
Pythia_syst_el_2017 = 0.02*0.02 #It will be summed in quadrature to the binomial uncertainty
SF_syst_el_2017     = 0.0125*0.0125 #It will be summed in quadrature to the binomial uncertainty
eff_el_syst_2017    = ROOT.RooRealVar("eff_el_syst_2017","eff_el_syst_2017", math.sqrt(binom_eff_el_2017+BDT_syst_el_2017+Pythia_syst_el_2017+SF_syst_el_2017))
eff_el_kappa_2017   = ROOT.RooRealVar("eff_el_kappa_2017","lumi_kappa_2017",1.+eff_el_syst_2017.getVal())
eff_el_beta_2017    = ROOT.RooRealVar("eff_el_beta_2017","eff_el_beta_2017",0.,-5.,5.)
eff_el_2017         = ROOT.RooFormulaVar("eff_el_2017","@0 * pow(@1,@2)",ROOT.RooArgList(eff_el_nominal_2017,eff_el_kappa_2017,eff_el_beta_2017))
glb_eff_el_2017     = ROOT.RooRealVar("glb_eff_el_2017","glb_eff_el_2017",0.,-5.,5.)
gauss_eff_el_2017   = ROOT.RooGaussian("gauss_eff_el_2017","gauss_eff_el_2017",glb_eff_el_2017,eff_el_beta_2017,one)


eff_mu_nominal_2018 = ROOT.RooRealVar("eff_mu_nominal_2018","eff_mu_nominal_2018",totmu_2018*2./totsig_2018)
binom_eff_mu_2018   = (4*totmu_2018*(totsig_2018-2*totmu_2018)/(totsig_2018*totsig_2018*totsig_2018))*(4*totmu_2018*(totsig_2018-2*totmu_2018)/(totsig_2018*totsig_2018*totsig_2018))/(eff_mu_nominal_2018.getVal()*eff_mu_nominal_2018.getVal()) #It will be summed in quadrature to the BDT systematic
BDT_syst_mu_2018    = 0.01*0.01 #It will be summed in quadrature to the binomial uncertainty
Pythia_syst_mu_2018 = 0.02*0.02 #It will be summed in quadrature to the binomial uncertainty
eff_mu_syst_2018    = ROOT.RooRealVar("eff_mu_syst_2018","eff_mu_syst_2018", math.sqrt(binom_eff_mu_2018+BDT_syst_mu_2018+Pythia_syst_mu_2018))
eff_mu_kappa_2018   = ROOT.RooRealVar("eff_mu_kappa_2018","lumi_kappa_2018",1.+eff_mu_syst_2018.getVal())
eff_mu_beta_2018    = ROOT.RooRealVar("eff_mu_beta_2018","eff_mu_beta_2018",0.,-5.,5.)
eff_mu_2018         = ROOT.RooFormulaVar("eff_mu_2018","@0 * pow(@1,@2)",ROOT.RooArgList(eff_mu_nominal_2018,eff_mu_kappa_2018,eff_mu_beta_2018))
glb_eff_mu_2018     = ROOT.RooRealVar("glb_eff_mu_2018","glb_eff_mu_2018",0.,-5.,5.)
gauss_eff_mu_2018   = ROOT.RooGaussian("gauss_eff_mu_2018","gauss_eff_mu_2018",glb_eff_mu_2018,eff_mu_beta_2018,one)


eff_el_nominal_2018 = ROOT.RooRealVar("eff_el_nominal_2018","eff_el_nominal_2018",totel_2018*2./totsig_2018)
binom_eff_el_2018   = (4*totel_2018*(totsig_2018-2*totel_2018)/(totsig_2018*totsig_2018*totsig_2018))*(4*totel_2018*(totsig_2018-2*totel_2018)/(totsig_2018*totsig_2018*totsig_2018))/(eff_el_nominal_2018.getVal()*eff_el_nominal_2018.getVal()) #It will be summed in quadrature to the BDT systematic
BDT_syst_el_2018    = 0.01*0.01 #It will be summed in quadrature to the binomial uncertainty
Pythia_syst_el_2018 = 0.02*0.02 #It will be summed in quadrature to the binomial uncertainty
SF_syst_el_2018     = 0.0125*0.0125 #It will be summed in quadrature to the binomial uncertainty
eff_el_syst_2018    = ROOT.RooRealVar("eff_el_syst_2018","eff_el_syst_2018", math.sqrt(binom_eff_el_2018+BDT_syst_el_2018+Pythia_syst_el_2018+SF_syst_el_2018))
eff_el_kappa_2018   = ROOT.RooRealVar("eff_el_kappa_2018","lumi_kappa_2018",1.+eff_el_syst_2018.getVal())
eff_el_beta_2018    = ROOT.RooRealVar("eff_el_beta_2018","eff_el_beta_2018",0.,-5.,5.)
eff_el_2018         = ROOT.RooFormulaVar("eff_el_2018","@0 * pow(@1,@2)",ROOT.RooArgList(eff_el_nominal_2018,eff_el_kappa_2018,eff_el_beta_2018))
glb_eff_el_2018     = ROOT.RooRealVar("glb_eff_el_2018","glb_eff_el_2018",0.,-5.,5.)
gauss_eff_el_2018   = ROOT.RooGaussian("gauss_eff_el_2018","gauss_eff_el_2018",glb_eff_el_2018,eff_el_beta_2018,one) 

print eff_mu_syst_2016.getVal()
print eff_el_syst_2016.getVal()
print eff_mu_syst_2017.getVal()
print eff_el_syst_2017.getVal()
print eff_mu_syst_2018.getVal()
print eff_el_syst_2018.getVal()

################################################################
#                                                              #
#------------- Systematic on bkg parametrization --------------#
#                                                              #
################################################################

if not suppressBkgSystematic:
    bkg_syst_mu = 0.139
    bkg_syst_el = 0.049
else:
    bkg_syst_mu = 0.0001
    bkg_syst_el = 0.0001

glb_bkg_param_mu   = ROOT.RooRealVar("glb_bkg_param_mu","glb_bkg_param_mu", 0.,-5.,-5.)
bkg_param_kappa_mu = ROOT.RooRealVar("bkg_param_kappa_mu","bkg_param_kappa_mu",1.+bkg_syst_mu)
bkg_param_beta_mu  = ROOT.RooRealVar("bkg_param_beta_mu","bkg_param_beta_mu",0.,-5.,5.)
eta_mu             = ROOT.RooFormulaVar("eta_mu","@0 * pow(@1,@2)",ROOT.RooArgList(one,bkg_param_kappa_mu,bkg_param_beta_mu))
gauss_bkg_param_mu = ROOT.RooGaussian("gauss_bkg_param_mu","gauss_bkg_param_mu",glb_bkg_param_mu,bkg_param_beta_mu,one)

glb_bkg_param_el   = ROOT.RooRealVar("glb_bkg_param_el","glb_bkg_param_el", 0.,-5.,-5.)
bkg_param_kappa_el = ROOT.RooRealVar("bkg_param_kappa_el","bkg_param_kappa_el",1.+bkg_syst_el)
bkg_param_beta_el  = ROOT.RooRealVar("bkg_param_beta_el","bkg_param_beta_el",0.,-5.,5.)
eta_el             = ROOT.RooFormulaVar("eta_el","@0 * pow(@1,@2)",ROOT.RooArgList(one,bkg_param_kappa_el,bkg_param_beta_el))
gauss_bkg_param_el = ROOT.RooGaussian("gauss_bkg_param_el","gauss_bkg_param_el",glb_bkg_param_el,bkg_param_beta_el,one)


################################################################
#                                                              #
#----------------------- Define the POI -----------------------#
#                                                              #
################################################################

if suppressBkgSystematic:
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
glb_bkg_param_mu.setConstant(1) 
glb_bkg_param_el.setConstant(1)
glb_lumi_2016.setConstant(1)
glb_lumi_2017.setConstant(1)
glb_lumi_2018.setConstant(1)
glb_eff_mu_2016.setConstant(1)
glb_eff_el_2016.setConstant(1)
glb_eff_mu_2017.setConstant(1)
glb_eff_el_2017.setConstant(1)
glb_eff_mu_2018.setConstant(1)
glb_eff_el_2018.setConstant(1)

################################################################
#                                                              #
#------------------- Define PDFs for the fit ------------------#
#                                                              #
################################################################
Nsig_mu = ROOT.RooFormulaVar("Nsig_mu","@0*@1*@2*(@3*@4+@5*@6+@7*@8)", ROOT.RooArgList(W_pigamma_BR, W_xsec, eta_mu, lumi_2016, eff_mu_2016, lumi_2017, eff_mu_2017, lumi_2018, eff_mu_2018))

Nsig_el = ROOT.RooFormulaVar("Nsig_el","@0*@1*@2*(@3*@4+@5*@6+@7*@8)", ROOT.RooArgList(W_pigamma_BR, W_xsec, eta_el, lumi_2016, eff_el_2016, lumi_2017, eff_el_2017, lumi_2018, eff_el_2018))

Nbkg_mu = ROOT.RooRealVar("Nbkg_mu","Nbkg_mu",900.,300.,3000.)
Nbkg_el = ROOT.RooRealVar("Nbkg_el","Nbkg_el",900.,200.,3000.)

#Per qualche motivo, il numero di argomenti della RooArgList sembra essere limitato. Se non uso gauss_lumi_product, ce n'e' uno uno di troppo in totPDF
gauss_lumi_product = ROOT.RooProdPdf("gauss_lumi_product","gauss_lumi_product",ROOT.RooArgList(gauss_lumi_2016,gauss_lumi_2017,gauss_lumi_2018))
gauss_eff_mu_product = ROOT.RooProdPdf("gauss_eff_mu_product","gauss_eff_mu_product",ROOT.RooArgList(gauss_eff_mu_2016,gauss_eff_mu_2017,gauss_eff_mu_2018))
gauss_eff_el_product = ROOT.RooProdPdf("gauss_eff_el_product","gauss_eff_el_product",ROOT.RooArgList(gauss_eff_el_2016,gauss_eff_el_2017,gauss_eff_el_2018))


totPDF_mu_unconstr = ROOT.RooAddPdf("totPDF_mu_unconstr","Total PDF for the mu channel",ROOT.RooArgList(totSignal,backPDF_mu),ROOT.RooArgList(Nsig_mu,Nbkg_mu))
totPDF_el_unconstr = ROOT.RooAddPdf("totPDF_el_unconstr","Total PDF for the el channel",ROOT.RooArgList(totSignal,backPDF_el),ROOT.RooArgList(Nsig_el,Nbkg_el))


totPDF_mu = ROOT.RooProdPdf("totPDF_mu","totPDF_mu", ROOT.RooArgList(totPDF_mu_unconstr,gauss_lumi_product,gauss_W_xsec,gauss_eff_mu_product,gauss_W_resol,gauss_W_pole_resol,gauss_bkg_param_mu))
totPDF_el = ROOT.RooProdPdf("totPDF_el","totPDF_el", ROOT.RooArgList(totPDF_el_unconstr,gauss_lumi_product,gauss_W_xsec,gauss_eff_el_product,gauss_W_resol,gauss_W_pole_resol,gauss_bkg_param_el))


################################################################
#                                                              #
#------------- Create the global simultaneous PDF -------------#
#                                                              #
################################################################

totPDF = ROOT.RooSimultaneous("totPDF","The total PDF",Categorization)

if not suppressAllSystematics:
    constrained_params = ROOT.RooArgSet()
    totPDF.addPdf(totPDF_mu,"MuonSignal")
    totPDF.addPdf(totPDF_el,"ElectronSignal")
    constrained_params.add(dCB_width)
    constrained_params.add(dCB_pole)
    #constrained_params.add(eta_mu)
    #constrained_params.add(eta_el)
else:
    totPDF.addPdf(totPDF_mu_unconstr,"MuonSignal")
    totPDF.addPdf(totPDF_el_unconstr,"ElectronSignal")
    dCB_width.setConstant(1)
    dCB_pole.setConstant(1)
    W_xsec_beta.setConstant(1)
    lumi_beta_2016.setConstant(1)
    lumi_beta_2017.setConstant(1)
    lumi_beta_2018.setConstant(1)
    eff_mu_beta_2016.setConstant(1)
    eff_el_beta_2016.setConstant(1)
    eff_mu_beta_2017.setConstant(1)
    eff_el_beta_2017.setConstant(1)
    eff_mu_beta_2018.setConstant(1)
    eff_el_beta_2018.setConstant(1)
    bkg_param_beta_mu.setConstant(1)
    bkg_param_beta_el.setConstant(1)
    #eta_mu.setConstant(1)
    #eta_el.setConstant(1)

################################################################
#                                                              #
#---------------------- Fit (and F-Test) ----------------------#
#                                                              #
################################################################

if isData:
    if not suppressAllSystematics:
        result_dataFit = totPDF.fitTo(data,ROOT.RooFit.Extended(1), ROOT.RooFit.Constrain(constrained_params), ROOT.RooFit.NumCPU(4), ROOT.RooFit.Save() )#For the signal region, I want the fit to be extended (Poisson fluctuation of unmber of events) to take into account that the total number of events is the sum of signal and background events. Either I do this, or I use a fraction frac*Nbkg+(1-frac)*Nsig, which will become a parameter of the fit and will have a Gaussian behavior (whilst the extended fit preserves the natural Poisson behavior)
    else:
        result_dataFit = totPDF.fitTo(data,ROOT.RooFit.Extended(1), ROOT.RooFit.NumCPU(4), ROOT.RooFit.Save() )
else:
    totPDF.fitTo(data,ROOT.RooFit.Extended(1), ROOT.RooFit.SumW2Error(0), ROOT.RooFit.NumCPU(4), ROOT.RooFit.Constrain(constrained_params) )


################################################################
#                                                              #
#-------------------------- Plotting --------------------------#
#                                                              #
################################################################

xframe_mu = Wmass.frame(55.,95.,15)
xframe_mu.SetTitle(" ")
xframe_mu.SetTitleOffset(1.4,"y")
xframe_mu.SetMaximum(90)

xframe_el = Wmass.frame(55.,95.,15)
xframe_el.SetTitle(" ")
xframe_el.SetTitleOffset(1.4,"y")
xframe_el.SetMaximum(90)

#################################################

#Exclude the control regions
data_reduced = data.reduce("Wmass < 65. || Wmass > 90.")

data_reduced.plotOn(xframe_mu, ROOT.RooFit.Cut("Categorization==1"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
totPDF.plotOn(xframe_mu, ROOT.RooFit.Range("LowSideband,HighSideband"), ROOT.RooFit.Slice(Categorization,"MuonSignal"), ROOT.RooFit.ProjWData(data))

data_reduced.plotOn(xframe_el, ROOT.RooFit.Cut("Categorization==3"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
totPDF.plotOn(xframe_el, ROOT.RooFit.Range("LowSideband,HighSideband"), ROOT.RooFit.Slice(Categorization,"ElectronSignal"), ROOT.RooFit.ProjWData(data))

#DRAW ON CANVAS
canvas = ROOT.TCanvas()
canvas.Divide(2,1)
canvas.cd(1)
xframe_mu.Draw()
canvas.cd(2)
xframe_el.Draw()

# Save the plot
if isData:
    canvas.SaveAs("plots/fitData_signalR.pdf")
else:
    canvas.SaveAs("plots/fitMC_signalR.pdf")

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
del workspace_out

raw_input()
