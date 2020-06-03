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
selectSigShift = 0
suppressSigSystematic = False

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
    data_initial = ROOT.RooDataSet("data","data", ROOT.RooArgSet(Wmass,Categorization,weight), ROOT.RooFit.Import(mytree), ROOT.RooFit.WeightVar("weight"))

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

workspace_sig = fInput_sigmodel.Get("myworkspace")

dCB_pole         = workspace_sig.var("dCB_pole")
dCB_pole_nominal = dCB_pole.getVal()
dCB_pole_err     = dCB_pole.getError()

dCB_width         = workspace_sig.var("dCB_width")
dCB_width_nominal = dCB_width.getVal()
dCB_width_err     = dCB_width.getError()

totSignal = workspace_sig.pdf("totSignal")

#Fix the signal parametrization
workspace_sig.var("dCB_aL").setConstant(1)
workspace_sig.var("dCB_nL").setConstant(1)
workspace_sig.var("dCB_aR").setConstant(1)
workspace_sig.var("dCB_nR").setConstant(1)
workspace_sig.var("Gauss_pole").setConstant(1)
workspace_sig.var("Gauss_sigma").setConstant(1)
workspace_sig.var("fracSig").setConstant(1)
if selectSigShift == 0:
    workspace_sig.var("dCB_pole").setConstant(1)
    workspace_sig.var("dCB_width").setConstant(1)
if selectSigShift == 1:
    print "old value of the pole: ", dCB_pole_nominal, "  dCB_pole_err: ", dCB_pole_err
    print "old value of the width: ", dCB_width_nominal, "  dCB_width_err: ", dCB_width_err
    workspace_sig.var("dCB_pole").setVal(dCB_pole_nominal+dCB_pole_err)
    workspace_sig.var("dCB_width").setVal(dCB_width_nominal+dCB_width_err)
    workspace_sig.var("dCB_pole").setConstant(1)
    workspace_sig.var("dCB_width").setConstant(1)
    print "new value of the pole: ", workspace_sig.var("dCB_pole").getVal()
    print "new value of the width: ", workspace_sig.var("dCB_width").getVal()
if selectSigShift == 2:
    workspace_sig.var("dCB_pole").setVal(dCB_pole_nominal-dCB_pole_err)
    workspace_sig.var("dCB_width").setVal(dCB_width_nominal+dCB_width_err)
    workspace_sig.var("dCB_pole").setConstant(1)
    workspace_sig.var("dCB_width").setConstant(1)
if selectSigShift == 3:
    workspace_sig.var("dCB_pole").setVal(dCB_pole_nominal+dCB_pole_err)
    workspace_sig.var("dCB_width").setVal(dCB_width_nominal-dCB_width_err)
    workspace_sig.var("dCB_pole").setConstant(1)
    workspace_sig.var("dCB_width").setConstant(1)
if selectSigShift == 4:
    workspace_sig.var("dCB_pole").setVal(dCB_pole_nominal-dCB_pole_err)
    workspace_sig.var("dCB_width").setVal(dCB_width_nominal-dCB_width_err)
    workspace_sig.var("dCB_pole").setConstant(1)
    workspace_sig.var("dCB_width").setConstant(1)


################################################################
#                                                              #
#---------------- Variables for bkg description ---------------#
#                                                              #
################################################################

#Now describe the background
fInput_bkgmodel = ROOT.TFile("fitBackground.root")
fInput_bkgmodel.cd()

workspace_bkg = fInput_bkgmodel.Get("workspace_bkg")

if selectBkgFunction == 0:

    backPDF = workspace_bkg.pdf("backPDF_cheb")

elif selectBkgFunction == 1: #Use Exponential for muon channel (calculation of systematic on background parametrization)

    backPDF = workspace_bkg.pdf("backPDF_exp")


################################################################
#                                                              #
#------------------ Systematic on ttbar xsec ------------------#
#                                                              #
################################################################
#First the cross section, with a modifier for systematics
#CMS ttbar measurement/W->lnu BR (it is measured with both W in lnu), in pb
#http://cms-results.web.cern.ch/cms-results/public-results/publications/TOP-16-005/index.html
W_xsec_nominal = 2.*2.*815.*0.1086 #The two factors 2 account for the possible charge signs of the Ws and for the two leptonic decay channels of the tag W
W_xsec_syst    = (43.*2.*2.*0.1086)/W_xsec_nominal

################################################################
#                                                              #
#------------------ Systematic on luminosity ------------------#
#                                                              #
################################################################

#Represent the luminosity with a modifier for systematics. For 2016 (35.86 fb-1): 2.5% systematic. For 2017 (41.529 fb-1): 2.3% systematic
lumi_2016 = 35.86*1000.
lumi_2017 = 41.53*1000.
lumi_2018 = 59.69*1000.
lumi_syst_2016 = 0.025
lumi_syst_2017 = 0.023 
lumi_syst_2018 = 0.025
lumi_nominal = lumi_2016 + lumi_2017 + lumi_2018 #In pb
lumi_syst    = (lumi_syst_2016*lumi_2016 + lumi_syst_2017*lumi_2017 + lumi_syst_2018*lumi_2018)/lumi_nominal

################################################################
#                                                              #
#------------------ Systematic on efficiency ------------------#
#                                                              #
################################################################

totsig_2016 = 100000. #total number of signal events in 2016
totmu_2016 = 6046.
totel_2016 = 3686.
tot_2016  = totmu_2016 + totel_2016

totsig_2017 = 80000.  #total number of signal events in 2017
totmu_2017 = 4310.
totel_2017 = 2934.
tot_2017  = totmu_2017 + totel_2017

totsig_2018 = 79820.  #total number of signal events in 2018
totmu_2018 = 4373.
totel_2018 = 2815.
tot_2018  = totmu_2018 + totel_2018

eff_nominal_2016 = tot_2016*2./totsig_2016
binom_eff_2016   = (4*tot_2016*(totsig_2016-2*tot_2016)/(totsig_2016*totsig_2016*totsig_2016))*(4*tot_2016*(totsig_2016-2*tot_2016)/(totsig_2016*totsig_2016*totsig_2016))/(eff_nominal_2016*eff_nominal_2016) #It will be summed in quadrature to the BDT systematic
BDT_syst_2016    = ((0.01*totmu_2016 + 0.01 * totel_2016)/tot_2016)**2 #It will be summed in quadrature to the binomial uncertainty
Pythia_syst_pT_2016 = ((0.02*totmu_2016 + 0.02*totel_2016)/tot_2016)**2 #It will be summed in quadrature to the binomial uncertainty
Pythia_syst_angle_2016 = ((0.05*totmu_2016 + 0.05*totel_2016)/tot_2016)**2 #It will be summed in quadrature to the binomial uncertainty
eff_syst_2016    = math.sqrt(binom_eff_2016+BDT_syst_2016+Pythia_syst_pT_2016+Pythia_syst_angle_2016)

eff_nominal_2017 = tot_2017*2./totsig_2017
binom_eff_2017   = (4*tot_2017*(totsig_2017-2*tot_2017)/(totsig_2017*totsig_2017*totsig_2017))*(4*tot_2017*(totsig_2017-2*tot_2017)/(totsig_2017*totsig_2017*totsig_2017))/(eff_nominal_2017*eff_nominal_2017) #It will be summed in quadrature to the BDT systematic
BDT_syst_2017    = ((0.01*totmu_2017 + 0.01 * totel_2017)/tot_2017)**2 #It will be summed in quadrature to the binomial uncertainty
Pythia_syst_pT_2017 = ((0.02*totmu_2017 + 0.02*totel_2017)/tot_2017)**2 #It will be summed in quadrature to the binomial uncertainty
Pythia_syst_angle_2017 = ((0.05*totmu_2017 + 0.05*totel_2017)/tot_2017)**2 #It will be summed in quadrature to the binomial uncertainty
eff_syst_2017    = math.sqrt(binom_eff_2017+BDT_syst_2017+Pythia_syst_pT_2017+Pythia_syst_angle_2017)

eff_nominal_2018 = tot_2018*2./totsig_2018
binom_eff_2018   = (4*tot_2018*(totsig_2018-2*tot_2018)/(totsig_2018*totsig_2018*totsig_2018))*(4*tot_2018*(totsig_2018-2*tot_2018)/(totsig_2018*totsig_2018*totsig_2018))/(eff_nominal_2018*eff_nominal_2018) #It will be summed in quadrature to the BDT systematic
BDT_syst_2018    = ((0.01*totmu_2018 + 0.01 * totel_2018)/tot_2018)**2 #It will be summed in quadrature to the binomial uncertainty
Pythia_syst_pT_2018 = ((0.02*totmu_2018 + 0.02*totel_2018)/tot_2018)**2 #It will be summed in quadrature to the binomial uncertainty
Pythia_syst_angle_2018 = ((0.05*totmu_2018 + 0.05*totel_2018)/tot_2018)**2 #It will be summed in quadrature to the binomial uncertainty
eff_syst_2018    = math.sqrt(binom_eff_2018+BDT_syst_2018+Pythia_syst_pT_2018+Pythia_syst_angle_2018)

efflumi_2016_nominal = eff_nominal_2016*lumi_2016
efflumi_2017_nominal = eff_nominal_2017*lumi_2017
efflumi_2018_nominal = eff_nominal_2018*lumi_2018

efflumi_2016_syst = (eff_syst_2016+lumi_syst_2016)*efflumi_2016_nominal
efflumi_2017_syst = (eff_syst_2017+lumi_syst_2017)*efflumi_2017_nominal
efflumi_2018_syst = (eff_syst_2018+lumi_syst_2018)*efflumi_2018_nominal

efflumi_nominal = efflumi_2016_nominal + efflumi_2017_nominal + efflumi_2018_nominal
efflumi_syst = math.sqrt( efflumi_2016_syst**2. + efflumi_2017_syst**2. + efflumi_2018_syst**2. )/efflumi_nominal

################################################################
#                                                              #
#------------- Systematic on bkg parametrization --------------#
#                                                              #
################################################################

if not suppressSigSystematic:
    sig_syst = 0.011
else:
    sig_syst = 0.0001

if not suppressBkgSystematic:
    bkg_syst = 0.198 
else:
    bkg_syst = 0.0001

eta_nominal = 1.
eta_sig_syst = sig_syst
eta_bkg_syst = bkg_syst

################################################################
#                                                              #
#----------------------- Define the POI -----------------------#
#                                                              #
################################################################

if suppressBkgSystematic or suppressSigSystematic:
    W_pigamma_BR = ROOT.RooRealVar("W_pigamma_BR","W_pigamma_BR",0.000005,-0.0001,0.01) # The parameter of interest can go negative when we try the alternative bkg description
else:
    W_pigamma_BR = ROOT.RooRealVar("W_pigamma_BR","W_pigamma_BR",0.000001,0.,0.01) # The parameter of interest

W_pigamma_BR_blind = ROOT.RooUnblindOffset("W_pigamma_BR_blind","W_pigamma_BR_blind","aSeedString",0.000001,W_pigamma_BR)


################################################################
#                                                              #
#------------------- Define PDFs for the fit ------------------#
#                                                              #
################################################################

Nsig_multiplier_value = W_xsec_nominal * efflumi_nominal * eta_nominal
Nsig_multiplier_value_err = W_xsec_syst + efflumi_syst + eta_sig_syst + eta_bkg_syst

Multi_param_nominal = ROOT.RooRealVar("Multi_param_nominal","Multi_param_nominal", Nsig_multiplier_value)
glb_Multi_param   = ROOT.RooRealVar("glb_Multi_param","glb_Multi_param", 0.,-5.,5.)
Multi_param_kappa = ROOT.RooRealVar("Multi_param_kappa","Multi_param_kappa",1.+Nsig_multiplier_value_err)
Multi_param_beta  = ROOT.RooRealVar("Multi_param_beta","Multi_param_beta",0.,-5.,5.)
Nsig_multiplier   = ROOT.RooFormulaVar("Nsig_multiplier","@0 * pow(@1,@2)",ROOT.RooArgList(Multi_param_nominal,Multi_param_kappa,Multi_param_beta))
one               = ROOT.RooRealVar("one","one",1.)
gauss_Multi_param = ROOT.RooGaussian("gauss_Multi_param","gauss_Multi_param",glb_Multi_param,Multi_param_beta,one)

Nsig = ROOT.RooFormulaVar("Nsig_mu","@0*@1", ROOT.RooArgList(W_pigamma_BR, Nsig_multiplier))

Nbkg = ROOT.RooRealVar("Nbkg","Nbkg",900.,100.,5000.)

totPDF_unconstr = ROOT.RooAddPdf("totPDF_unconstr","Total PDF",ROOT.RooArgList(totSignal,backPDF),ROOT.RooArgList(Nsig,Nbkg))

totPDF = ROOT.RooProdPdf("totPDF","totPDF", ROOT.RooArgList(totPDF_unconstr,gauss_Multi_param))

################################################################
#                                                              #
#------------------ Set a few things constant -----------------#
#                                                              #
################################################################
glb_Multi_param.setConstant(1) 

if suppressAllSystematics:
    Multi_param_beta.setConstant(1)

################################################################
#                                                              #
#---------------------- Fit (and F-Test) ----------------------#
#                                                              #
################################################################

constrained_params = ROOT.RooArgSet()
constrained_params.add(Multi_param_beta)

if isData:
    if not suppressAllSystematics:
        result_dataFit = totPDF.fitTo(data, ROOT.RooFit.Extended(1), ROOT.RooFit.Constrain(constrained_params), ROOT.RooFit.Save() )#For the signal region, I want the fit to be extended (Poisson fluctuation of unmber of events) to take into account that the total number of events is the sum of signal and background events. Either I do this, or I use a fraction frac*Nbkg+(1-frac)*Nsig, which will become a parameter of the fit and will have a Gaussian behavior (whilst the extended fit preserves the natural Poisson behavior)
    else:
        result_dataFit = totPDF.fitTo(data, ROOT.RooFit.Extended(1), ROOT.RooFit.Save() )
else:
    totPDF.fitTo(data,ROOT.RooFit.Extended(1), ROOT.RooFit.SumW2Error(0), ROOT.RooFit.NumCPU(4) )

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


################################################################
#                                                              #
#-------------------------- Plotting --------------------------#
#                                                              #
################################################################

xframe = Wmass.frame(55.,95.,15)
xframe.SetTitle(" ")
xframe.SetTitleOffset(1.4,"y")

#################################################

#Exclude the control regions
data_reduced = data.reduce("Wmass < 65. || Wmass > 90.")

data_reduced.plotOn(xframe, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
totPDF.plotOn(xframe, ROOT.RooFit.Range("LowSideband,HighSideband"))


#DRAW ON CANVAS
canvas = ROOT.TCanvas()
canvas.cd()
xframe.Draw()


# Save the plot
if isData:
    canvas.SaveAs("plots/fitData_signalR.pdf")
else:
    canvas.SaveAs("plots/fitMC_signalR.pdf")

del workspace_out

raw_input()
