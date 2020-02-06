#This program fits both lepton samples at the same time

import ROOT
import math
import argparse

ROOT.gROOT.ProcessLineSync(".L dCB/RooDoubleCBFast.cc+")

#---------------------------------#
p = argparse.ArgumentParser(description='Select whether to fit data or MC')
p.add_argument('isData_option', help='Type <<data>> or <<MC>>')
p.add_argument('runningEra_option', help='Type <<0>> for 2016, <<1>> for 2017, <<2>> for 2016+2017, <<3>> for 2016+2017+2018')
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

selectBkgFunction = 0 # 0: use Chebychev for ALL the bkg PDFs. 1: use Bernstein for mu 2016. 2: use Bernstein for ele 2016. 3: use Bernstein for mu 2017. 4: use Bernstein for ele 2017
suppressBkgSystematic = False # To be used when trying to fit with alternative bkg description, in order to estimate a systematic. If True, it will allow W_pigamma_BR to float negative. Moreover, it will use Signal+Background in the totPDF, so that the fit to the restricted CRs will contain also the POI BR, which will be used to calculate the pull and hence to estimate the systematic


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
Categorization.defineType("MuonSignal",1)
Categorization.defineType("ElectronSignal",3)

isSignal = ROOT.RooCategory("isSignal","isSignal")
isSignal.defineType("Signal",1)
isSignal.defineType("Background",0)

# isMuon = ROOT.RooCategory("isMuon","isMuon")
# isMuon.defineType("Muon",1)
# isMuon.defineType("Electron",0)

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
fInput_sigmodel = ROOT.TFile("Signal_model_" + str(runningEra) + ".root")
fInput_sigmodel.cd()

workspace = fInput_sigmodel.Get("myworkspace")

#Fix the signal parametrization
workspace.var("dCB_pole").setConstant(1)
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

elif selectBkgFunction == 1: #Use Bernstein for muon channel (calculation of systematic on background parametrization)

    backPDF_mu = workspace_bkg.pdf("backPDF_bern_mu")
    backPDF_el = workspace_bkg.pdf("backPDF_cheb_el")

elif selectBkgFunction == 2: #Use Bernstein for electron channel (calculation of systematic on background parametrization)

    backPDF_mu = workspace_bkg.pdf("backPDF_cheb_mu")
    backPDF_el = workspace_bkg.pdf("backPDF_bern_el")


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

glb_lumi_2018    = ROOT.RooRealVar("glb_lumi_2018","glb_lumi_2018",59.69*1000., 0., 65000.) #In pb
lumi_constr_2018 = ROOT.RooRealVar("lumi_constr_2018","lumi_constr_2018", 59.69*1000., 0., 65000.)
lumi_syst_2018   = ROOT.RooRealVar("lumi_syst_2018","lumi_syst_2018", 59.69*0.025*1000.)
gauss_lumi_2018  = ROOT.RooGaussian("gauss_lumi_2018","gauss_lumi_2018",glb_lumi_2018,lumi_constr_2018,lumi_syst_2018) 


################################################################
#                                                              #
#------------------ Systematic on efficiency ------------------#
#                                                              #
################################################################

totsig_2016 = 100000. #total number of signal events in 2016
totmu_2016  = 7425.   #total number of signal muon events in 2016
totel_2016  = 5284.   #total number of signal electron events in 2016

totsig_2017 = 80000.  #total number of signal events in 2017
totmu_2017  = 5485.   #total number of signal muon events in 2017
totel_2017  = 4167.   #total number of signal electron events in 2017

totsig_2018 = 79820.  #total number of signal events in 2018
totmu_2018  = 5654.   #total number of signal muon events in 2018
totel_2018  = 4086.   #total number of signal electron events in 2018

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

glb_eff_mu_2018    = ROOT.RooRealVar("glb_eff_mu_2018","glb_eff_mu_2018",totmu_2018*2./totsig_2018, 0., 1.)
eff_mu_constr_2018 = ROOT.RooRealVar("eff_mu_constr_2018","eff_mu_constr_2018", totmu_2018*2./totsig_2018, 0., 1.)
eff_mu_syst_2018   = ROOT.RooRealVar("eff_mu_syst_2018","eff_mu_syst_2018", 4*totmu_2018*(totsig_2018-2*totmu_2018)/(totsig_2018*totsig_2018*totsig_2018))
gauss_eff_mu_2018  = ROOT.RooGaussian("gauss_eff_mu_2018","gauss_eff_mu_2018",glb_eff_mu_2018,eff_mu_constr_2018,eff_mu_syst_2018) 

glb_eff_el_2018    = ROOT.RooRealVar("glb_eff_el_2018","glb_eff_el_2018", totel_2018*2./totsig_2018, 0., 1.)
eff_el_constr_2018 = ROOT.RooRealVar("eff_el_constr_2018","eff_el_constr_2018",totel_2018*2./totsig_2018, 0., 1.)
eff_el_syst_2018   = ROOT.RooRealVar("eff_el_syst_2018","eff_el_syst_2018", 4*totel_2018*(totsig_2018-2*totel_2018)/(totsig_2018*totsig_2018*totsig_2018))
gauss_eff_el_2018  = ROOT.RooGaussian("gauss_eff_el_2018","gauss_eff_el_2018",glb_eff_el_2018,eff_el_constr_2018,eff_el_syst_2018) 


################################################################
#                                                              #
#------------- Systematic on bkg parametrization --------------#
#                                                              #
################################################################

eta_mu_2016 = ROOT.RooRealVar("eta_mu_2016","eta_mu_2016", 1.,0.0001,3.)
eta_el_2016 = ROOT.RooRealVar("eta_el_2016","eta_el_2016", 1.,0.0001,3.)
eta_mu_2017 = ROOT.RooRealVar("eta_mu_2017","eta_mu_2017", 1.,0.0001,3.)
eta_el_2017 = ROOT.RooRealVar("eta_el_2017","eta_el_2017", 1.,0.0001,3.)
eta_mu_2018 = ROOT.RooRealVar("eta_mu_2018","eta_mu_2018", 1.,0.0001,3.)
eta_el_2018 = ROOT.RooRealVar("eta_el_2018","eta_el_2018", 1.,0.0001,3.)

#SUM OF THE YEARS
eta_mu_2016_2017_2018 = ROOT.RooRealVar("eta_mu_2016_2017_2018","eta_mu_2016_2017_2018", 1.,0.0001,3.)
eta_el_2016_2017_2018 = ROOT.RooRealVar("eta_el_2016_2017_2018","eta_el_2016_2017_2018", 1.,0.0001,3.)

if not suppressBkgSystematic:
    bkg_syst_mu_2016 = 0.004 #Value of the systematic to use in the fit of the signal regions
    bkg_syst_el_2016 = 0.002
    bkg_syst_mu_2017 = 0.001
    bkg_syst_el_2017 = 0.001
    bkg_syst_mu_2018 = 0.001
    bkg_syst_el_2018 = 0.001

    #SUM OF THE YEARS
    bkg_syst_mu_2016_2017_2018 = 0.003
    bkg_syst_el_2016_2017_2018 = 0.003

else:
    bkg_syst_mu_2016 = 0.0001 #In the suppressBkgSystematic mode, one is supposed to not know yet this systematic. 0.0001 is a custom small value
    bkg_syst_el_2016 = 0.0001
    bkg_syst_mu_2017 = 0.0001
    bkg_syst_el_2017 = 0.0001
    bkg_syst_mu_2018 = 0.0001
    bkg_syst_el_2018 = 0.0001
    
    #SUM OF THE YEARS
    bkg_syst_mu_2016_2017_2018 = 0.0001
    bkg_syst_el_2016_2017_2018 = 0.0001

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

glb_bkg_param_mu_2018    = ROOT.RooRealVar("glb_bkg_param_mu_2018","glb_bkg_param_mu_2018", 1., 0.0001, 3.)
bkg_param_syst_mu_2018   = ROOT.RooRealVar("bkg_param_syst_mu_2018","bkg_param_syst_mu_2018",bkg_syst_mu_2018)
gauss_bkg_param_mu_2018  = ROOT.RooGaussian("gauss_bkg_param_mu_2018","gauss_bkg_param_2018_mu",glb_bkg_param_mu_2018,eta_mu_2018,bkg_param_syst_mu_2018)

glb_bkg_param_el_2018    = ROOT.RooRealVar("glb_bkg_param_el_2018","glb_bkg_param_el_2018", 1., 0.0001, 3.)
bkg_param_syst_el_2018   = ROOT.RooRealVar("bkg_param_syst_el_2018","bkg_param_syst_el_2018",bkg_syst_el_2018)
gauss_bkg_param_el_2018  = ROOT.RooGaussian("gauss_bkg_param_el_2018","gauss_bkg_param_2018_el",glb_bkg_param_el_2018,eta_el_2018,bkg_param_syst_el_2018)

#SUM OF THE YEARS
glb_bkg_param_mu_2016_2017_2018    = ROOT.RooRealVar("glb_bkg_param_mu_2016_2017_2018","glb_bkg_param_mu_2016_2017_2018", 1., 0.0001, 3.)
bkg_param_syst_mu_2016_2017_2018   = ROOT.RooRealVar("bkg_param_syst_mu_2016_2017_2018","bkg_param_syst_mu_2016_2017_2018",bkg_syst_mu_2016_2017_2018)
gauss_bkg_param_mu_2016_2017_2018  = ROOT.RooGaussian("gauss_bkg_param_mu_2016_2017_2018","gauss_bkg_param_2016_2017_2018_mu",glb_bkg_param_mu_2016_2017_2018,eta_mu_2016_2017_2018,bkg_param_syst_mu_2016_2017_2018)

glb_bkg_param_el_2016_2017_2018    = ROOT.RooRealVar("glb_bkg_param_el_2016_2017_2018","glb_bkg_param_el_2016_2017_2018", 1., 0.0001, 3.)
bkg_param_syst_el_2016_2017_2018   = ROOT.RooRealVar("bkg_param_syst_el_2016_2017_2018","bkg_param_syst_el_2016_2017_2018",bkg_syst_el_2016_2017_2018)
gauss_bkg_param_el_2016_2017_2018  = ROOT.RooGaussian("gauss_bkg_param_el_2016_2017_2018","gauss_bkg_param_2016_2017_2018_el",glb_bkg_param_el_2016_2017_2018,eta_el_2016_2017_2018,bkg_param_syst_el_2016_2017_2018)

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

if runningEra == 0:
    glb_bkg_param_mu_2016.setConstant(1)
    glb_bkg_param_el_2016.setConstant(1)
    glb_lumi_2016.setConstant(1)
    glb_eff_mu_2016.setConstant(1)
    glb_eff_el_2016.setConstant(1)

elif runningEra == 1:
    glb_bkg_param_mu_2017.setConstant(1)
    glb_bkg_param_el_2017.setConstant(1)
    glb_lumi_2017.setConstant(1)
    glb_eff_mu_2017.setConstant(1)
    glb_eff_el_2017.setConstant(1)

elif runningEra == 2:
    glb_bkg_param_mu_2018.setConstant(1)
    glb_bkg_param_el_2018.setConstant(1)
    glb_lumi_2018.setConstant(1)
    glb_eff_mu_2018.setConstant(1)
    glb_eff_el_2018.setConstant(1)

elif runningEra == 3:
    glb_bkg_param_mu_2016_2017_2018.setConstant(1)
    glb_bkg_param_el_2016_2017_2018.setConstant(1)
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

#This is for the signal region
# Nsig_mu = ROOT.RooFormulaVar("Nsig_mu","@0*@1*@2*@3*@4", ROOT.RooArgList(W_pigamma_BR, W_xsec_constr,lumi_constr,eff_mu_constr,eta))
# Nsig_el = ROOT.RooFormulaVar("Nsig_el","@0*@1*@2*@3*@4", ROOT.RooArgList(W_pigamma_BR, W_xsec_constr,lumi_constr,eff_el_constr,eta))
Nsig_mu_2016 = ROOT.RooFormulaVar("Nsig_mu_2016","@0*@1*@2*@3*@4", ROOT.RooArgList(W_pigamma_BR_blind, W_xsec_constr, lumi_constr_2016, eff_mu_constr_2016, eta_mu_2016))
Nsig_el_2016 = ROOT.RooFormulaVar("Nsig_el_2016","@0*@1*@2*@3*@4", ROOT.RooArgList(W_pigamma_BR_blind, W_xsec_constr, lumi_constr_2016, eff_el_constr_2016, eta_el_2016))

Nsig_mu_2017 = ROOT.RooFormulaVar("Nsig_mu_2017","@0*@1*@2*@3*@4", ROOT.RooArgList(W_pigamma_BR_blind, W_xsec_constr, lumi_constr_2017, eff_mu_constr_2017, eta_mu_2017))
Nsig_el_2017 = ROOT.RooFormulaVar("Nsig_el_2017","@0*@1*@2*@3*@4", ROOT.RooArgList(W_pigamma_BR_blind, W_xsec_constr, lumi_constr_2017, eff_el_constr_2017, eta_el_2017))

Nsig_mu_2018 = ROOT.RooFormulaVar("Nsig_mu_2018","@0*@1*@2*@3*@4", ROOT.RooArgList(W_pigamma_BR_blind, W_xsec_constr, lumi_constr_2018, eff_mu_constr_2018, eta_mu_2018))
Nsig_el_2018 = ROOT.RooFormulaVar("Nsig_el_2018","@0*@1*@2*@3*@4", ROOT.RooArgList(W_pigamma_BR_blind, W_xsec_constr, lumi_constr_2018, eff_el_constr_2018, eta_el_2018))

#SUM OF YEARS
Nsig_mu_2016_2017_2018 = ROOT.RooFormulaVar("Nsig_mu_2016_2017_2018","@0*@1*@2*(@3*@4+@5*@6+@7*@8)", ROOT.RooArgList(W_pigamma_BR_blind, W_xsec_constr, eta_mu_2016_2017_2018, lumi_constr_2016, eff_mu_constr_2016, lumi_constr_2017, eff_mu_constr_2017, lumi_constr_2018, eff_mu_constr_2018))

Nsig_el_2016_2017_2018 = ROOT.RooFormulaVar("Nsig_el_2016_2017_2018","@0*@1*@2*(@3*@4+@5*@6+@7*@8)", ROOT.RooArgList(W_pigamma_BR_blind, W_xsec_constr, eta_el_2016_2017_2018, lumi_constr_2016, eff_el_constr_2016, lumi_constr_2017, eff_el_constr_2017, lumi_constr_2018, eff_el_constr_2018))


Nbkg_mu_2016 = ROOT.RooRealVar("Nbkg_mu_2016","Nbkg_mu_2016",300.,100.,1000.)
Nbkg_el_2016 = ROOT.RooRealVar("Nbkg_el_2016","Nbkg_el_2016",300.,100.,1000.)

Nbkg_mu_2017 = ROOT.RooRealVar("Nbkg_mu_2017","Nbkg_mu_2017",300.,100.,1000.)
Nbkg_el_2017 = ROOT.RooRealVar("Nbkg_el_2017","Nbkg_el_2017",300.,100.,1000.)

Nbkg_mu_2018 = ROOT.RooRealVar("Nbkg_mu_2017","Nbkg_mu_2017",300.,100.,1000.)
Nbkg_el_2018 = ROOT.RooRealVar("Nbkg_el_2017","Nbkg_el_2017",300.,100.,1000.)

#SUM OF YEARS
Nbkg_mu_2016_2017_2018 = ROOT.RooRealVar("Nbkg_mu_2016_2017_2018","Nbkg_mu_2016_2017_2018",900.,300.,3000.)
Nbkg_el_2016_2017_2018 = ROOT.RooRealVar("Nbkg_el_2016_2017_2018","Nbkg_el_2016_2017_2018",900.,300.,3000.)


if runningEra == 0:
    Nsig_mu = Nsig_mu_2016
    Nsig_el = Nsig_el_2016
    Nbkg_mu = Nbkg_mu_2016
    Nbkg_el = Nbkg_el_2016
    gauss_lumi_product = gauss_lumi_2016
    gauss_eff_mu_product = gauss_eff_mu_2016
    gauss_eff_el_product = gauss_eff_el_2016
    gauss_bkg_param_mu = gauss_bkg_param_mu_2016
    gauss_bkg_param_el = gauss_bkg_param_el_2016
elif runningEra == 1:
    Nsig_mu = Nsig_mu_2017
    Nsig_el = Nsig_el_2017
    Nbkg_mu = Nbkg_mu_2017
    Nbkg_el = Nbkg_el_2017
    gauss_lumi_product = gauss_lumi_2017
    gauss_eff_mu_product = gauss_eff_mu_2017
    gauss_eff_el_product = gauss_eff_el_2017
    gauss_bkg_param_mu = gauss_bkg_param_mu_2017
    gauss_bkg_param_el = gauss_bkg_param_el_2017
elif runningEra == 2:
    Nsig_mu = Nsig_mu_2018
    Nsig_el = Nsig_el_2018
    Nbkg_mu = Nbkg_mu_2018
    Nbkg_el = Nbkg_el_2018
    gauss_lumi_product = gauss_lumi_2018
    gauss_eff_mu_product = gauss_eff_mu_2018
    gauss_eff_el_product = gauss_eff_el_2018
    gauss_bkg_param_mu = gauss_bkg_param_mu_2018
    gauss_bkg_param_el = gauss_bkg_param_el_2018
elif runningEra == 3:
    Nsig_mu = Nsig_mu_2016_2017_2018
    Nsig_el = Nsig_el_2016_2017_2018
    Nbkg_mu = Nbkg_mu_2016_2017_2018
    Nbkg_el = Nbkg_el_2016_2017_2018
    #Per qualche motivo, il numero di argomenti della RooArgList sembra essere limitato. Se non uso gauss_lumi_product, ce n'e' uno uno di troppo in totPDF
    gauss_lumi_product = ROOT.RooProdPdf("gauss_lumi_product","gauss_lumi_product",ROOT.RooArgList(gauss_lumi_2016,gauss_lumi_2017,gauss_lumi_2018))
    gauss_eff_mu_product = ROOT.RooProdPdf("gauss_eff_mu_product","gauss_eff_mu_product",ROOT.RooArgList(gauss_eff_mu_2016,gauss_eff_mu_2017,gauss_eff_mu_2018))
    gauss_eff_el_product = ROOT.RooProdPdf("gauss_eff_el_product","gauss_eff_el_product",ROOT.RooArgList(gauss_eff_el_2016,gauss_eff_el_2017,gauss_eff_el_2018))
    gauss_bkg_param_mu = gauss_bkg_param_mu_2016_2017_2018
    gauss_bkg_param_el = gauss_bkg_param_el_2016_2017_2018



totPDF_mu_unconstr = ROOT.RooAddPdf("totPDF_mu_unconstr","Total PDF for the mu channel",ROOT.RooArgList(totSignal,backPDF_mu),ROOT.RooArgList(Nsig_mu,Nbkg_mu))
totPDF_el_unconstr = ROOT.RooAddPdf("totPDF_el_unconstr","Total PDF for the el channel",ROOT.RooArgList(totSignal,backPDF_el),ROOT.RooArgList(Nsig_el,Nbkg_el))


totPDF_mu = ROOT.RooProdPdf("totPDF_mu","totPDF_mu", ROOT.RooArgList(totPDF_mu_unconstr,gauss_lumi_product,gauss_W_xsec,gauss_eff_mu_product,gauss_W_resol,gauss_bkg_param_mu))
totPDF_el = ROOT.RooProdPdf("totPDF_el","totPDF_el", ROOT.RooArgList(totPDF_el_unconstr,gauss_lumi_product,gauss_W_xsec,gauss_eff_el_product,gauss_W_resol,gauss_bkg_param_el))

# totPDF_mu = ROOT.RooProdPdf("totPDF_mu","totPDF_mu", ROOT.RooArgList(totPDF_mu_unconstr,gauss_lumi_product,gauss_W_xsec,gauss_eff_mu_2016,gauss_eff_mu_2017,gauss_eff_mu_2018,gauss_W_resol,gauss_bkg_param_mu))
# totPDF_el = ROOT.RooProdPdf("totPDF_el","totPDF_el", ROOT.RooArgList(totPDF_el_unconstr,gauss_lumi_product,gauss_W_xsec,gauss_eff_el_2016,gauss_eff_el_2017,gauss_eff_el_2018,gauss_W_resol,gauss_bkg_param_el))

#FIXME

#Per qualche motivo, il numero di argomenti della RooArgList sembra essere limitato. Qui ce ne sta uno di troppo
#totPDF_mu_2016_2017_2018 = ROOT.RooProdPdf("totPDF_mu_2016_2017_2018","totPDF_mu_2016_2017_2018", ROOT.RooArgList(totPDF_mu_unconstr_2016_2017_2018,gauss_lumi_2016,gauss_lumi_2017,gauss_lumi_2018,gauss_W_xsec,gauss_eff_mu_2016,gauss_eff_mu_2017,gauss_eff_mu_2018,gauss_W_resol,gauss_bkg_param_mu_2016_2017_2018))
# totPDF_mu_2016_2017_2018 = ROOT.RooProdPdf("totPDF_mu_2016_2017_2018","totPDF_mu_2016_2017_2018", ROOT.RooArgList(totPDF_mu_unconstr_2016_,lumi_product,gauss_W_xsec,gauss_eff_mu_2016,gauss_eff_mu_2017,gauss_eff_mu_2018,gauss_W_resol,gauss_bkg_param_mu_2016_2017_2018))
# totPDF_el_2016_2017_2018 = ROOT.RooProdPdf("totPDF_el_2016_2017_2018","totPDF_el_2016_2017_2018", ROOT.RooArgList(totPDF_el_unconstr_2016_2017_2018,lumi_product,gauss_W_xsec,gauss_eff_el_2016,gauss_eff_el_2017,gauss_eff_el_2018,gauss_W_resol,gauss_bkg_param_el_2016_2017_2018))


################################################################
#                                                              #
#------------- Create the global simultaneous PDF -------------#
#                                                              #
################################################################

totPDF = ROOT.RooSimultaneous("totPDF","The total PDF",Categorization)
constrained_params = ROOT.RooArgSet()

if runningEra == 0: #Fit on 2016 signal region with background+signal PDF
    totPDF.addPdf(totPDF_mu,"MuonSignal")
    totPDF.addPdf(totPDF_el,"ElectronSignal")
    constrained_params.add(dCB_width)
    constrained_params.add(eta_mu_2016)
    constrained_params.add(eta_el_2016)
    constrained_params.add(W_xsec_constr)
    constrained_params.add(lumi_constr_2016)
    constrained_params.add(eff_mu_constr_2016)
    constrained_params.add(eff_el_constr_2016)

if runningEra == 1: #Fit on 2017 signal region with background+signal PDF
    totPDF.addPdf(totPDF_mu,"MuonSignal")
    totPDF.addPdf(totPDF_el,"ElectronSignal")
    constrained_params.add(dCB_width)
    constrained_params.add(eta_mu_2017)
    constrained_params.add(eta_el_2017)
    constrained_params.add(W_xsec_constr)
    constrained_params.add(lumi_constr_2017)
    constrained_params.add(eff_mu_constr_2017)
    constrained_params.add(eff_el_constr_2017)

if runningEra == 2: #Fit on 2018 signal regions
    totPDF.addPdf(totPDF_mu,"MuonSignal")
    totPDF.addPdf(totPDF_el,"ElectronSignal")
    constrained_params.add(dCB_width)
    constrained_params.add(eta_mu_2018)
    constrained_params.add(eta_el_2018)
    constrained_params.add(W_xsec_constr)
    constrained_params.add(lumi_constr_2018)
    constrained_params.add(eff_mu_constr_2018)
    constrained_params.add(eff_el_constr_2018)

#SUM OF YEARS
if runningEra == 3: #Fit on 2016+2017+2018 signal regions, summing the years but not the channels
    totPDF.addPdf(totPDF_mu,"MuonSignal")
    totPDF.addPdf(totPDF_el,"ElectronSignal")
    constrained_params.add(dCB_width)
    constrained_params.add(eta_mu_2016_2017_2018)
    constrained_params.add(eta_el_2016_2017_2018)
    constrained_params.add(W_xsec_constr)
    constrained_params.add(lumi_constr_2016)
    constrained_params.add(lumi_constr_2017)
    constrained_params.add(lumi_constr_2018)
    constrained_params.add(eff_mu_constr_2016)
    constrained_params.add(eff_el_constr_2016)
    constrained_params.add(eff_mu_constr_2017)
    constrained_params.add(eff_el_constr_2017)
    constrained_params.add(eff_mu_constr_2018)
    constrained_params.add(eff_el_constr_2018)


################################################################
#                                                              #
#---------------------- Fit (and F-Test) ----------------------#
#                                                              #
################################################################

if isData:
    # W_pigamma_BR.setVal(0.000006)
    # W_pigamma_BR.setConstant(1)
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


canvas = ROOT.TCanvas()
# canvas.Divide(2,2)
# canvas.cd(1)
# xframe_mu_2016.Draw()
# canvas.cd(2)
# xframe_el_2016.Draw()
# canvas.cd(3)
# xframe_mu_2017.Draw()
# canvas.cd(4)
# xframe_el_2017.Draw()

#FIXME
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
