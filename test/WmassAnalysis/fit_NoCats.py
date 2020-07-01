#This program fits both lepton samples at the same time

import ROOT
import math

ROOT.gROOT.ProcessLineSync(".L dCB/RooDoubleCBFast.cc+")

################################################################
#                                                              #
#------------------------ Instructions ------------------------#
#                                                              #
################################################################
selectSigShift = 0 #Zero is the nominal case: signal parameters fixed to their central values from fit on MC signal. Numbers from 1 to 4 are plus/minus 1sigma
suppressSigSystematic = False

selectBkgFunction = 0 # 0: use Chebychev for ALL the bkg PDFs. 1: use alternative PDF 
suppressBkgSystematic = False # To be used when trying to fit with alternative bkg description, in order to estimate a systematic.If True, it will allow W_pigamma_BR to float negative. 
#Moreover, it will use Signal+Background in the totPDF, so that the fit to the restricted CRs will contain also the POI BR, which will be used to calculate the pull and hence to estimate the systematic

suppressAllSystematics = False


################################################################
#                                                              #
#------------------- Define the observable --------------------#
#                                                              #
################################################################

Wmass = ROOT.RooRealVar("Wmass","m_{#pi#gamma}",50.,100.,"GeV/c^{2}")
Wmass.setRange("LowSideband",50.,65.)
Wmass.setRange("HighSideband",90.,100.)

################################################################
#                                                              #
#--------------------- Retrieve the sample --------------------#
#                                                              #
################################################################

fInput = ROOT.TFile("Tree_input_massfit_Data_3.root")

fInput.cd()

mytree = fInput.Get("minitree")

################################################################
#                                                              #
#-------------------- Define the categories -------------------#
#                                                              #
################################################################

#Define the mu/ele categories
Categorization = ROOT.RooCategory("Categorization","Categorization")
Categorization.defineType("MuonSignal",1)
Categorization.defineType("ElectronSignal",3)

################################################################
#                                                              #
#---------------------- Get the dataset -----------------------#
#                                                              #
################################################################

#Import dataset
data_initial = ROOT.RooDataSet("data","data", ROOT.RooArgSet(Wmass,Categorization), ROOT.RooFit.Import(mytree))
data = data_initial.reduce("(Categorization==Categorization::MuonSignal) || (Categorization==Categorization::ElectronSignal)")
#data = data_initial.reduce("(Categorization==Categorization::MuonSignal)")
#data = data_initial.reduce("(Categorization==Categorization::ElectronSignal)")

print "number of events mu - SR: ", data.sumEntries("Categorization==1")
print "number of events ele - SR: ", data.sumEntries("Categorization==3")
print "Using ", data.numEntries(), " events to fit"

################################################################
#                                                              #
#------------------ Import signal variables -------------------#
#                                                              #
################################################################

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

if selectSigShift == 0: #The nominal case (dCB parameters fixed to their central values from fit on MC signal)
    workspace_sig.var("dCB_pole").setConstant(1)
    workspace_sig.var("dCB_width").setConstant(1)
if selectSigShift == 1: #Both dCB_pole and dCB_width fixed to central value +1sigma
    workspace_sig.var("dCB_pole").setVal(dCB_pole_nominal+dCB_pole_err)
    workspace_sig.var("dCB_width").setVal(dCB_width_nominal+dCB_width_err)
    workspace_sig.var("dCB_pole").setConstant(1)
    workspace_sig.var("dCB_width").setConstant(1)
if selectSigShift == 2: #Both dCB_pole -1sigma and dCB_width +1sigma
    workspace_sig.var("dCB_pole").setVal(dCB_pole_nominal-dCB_pole_err)
    workspace_sig.var("dCB_width").setVal(dCB_width_nominal+dCB_width_err)
    workspace_sig.var("dCB_pole").setConstant(1)
    workspace_sig.var("dCB_width").setConstant(1)
if selectSigShift == 3: #Both dCB_pole +1sigma and dCB_width -1sigma
    workspace_sig.var("dCB_pole").setVal(dCB_pole_nominal+dCB_pole_err)
    workspace_sig.var("dCB_width").setVal(dCB_width_nominal-dCB_width_err)
    workspace_sig.var("dCB_pole").setConstant(1)
    workspace_sig.var("dCB_width").setConstant(1)
if selectSigShift == 4: #Both dCB_pole and dCB_width fixed to central value -1sigma
    workspace_sig.var("dCB_pole").setVal(dCB_pole_nominal-dCB_pole_err)
    workspace_sig.var("dCB_width").setVal(dCB_width_nominal-dCB_width_err)
    workspace_sig.var("dCB_pole").setConstant(1)
    workspace_sig.var("dCB_width").setConstant(1)


################################################################
#                                                              #
#----------------- Import background variables ----------------#
#                                                              #
################################################################

fInput_bkgmodel = ROOT.TFile("fitBackground.root")
fInput_bkgmodel.cd()

workspace_bkg = fInput_bkgmodel.Get("workspace_bkg")

if selectBkgFunction == 0: #Nominal case (Chebychev)

    backPDF = workspace_bkg.pdf("backPDF_cheb")

elif selectBkgFunction == 1: #Use Exponential (estimate systematic on background parametrization)

    backPDF = workspace_bkg.pdf("backPDF_exp")


################################################################
#                                                              #
#------------------ Systematic on ttbar xsec ------------------#
#                                                              #
################################################################
#CMS ttbar measurement/W->lnu BR (it is measured with both W in lnu), in pb
#http://cms-results.web.cern.ch/cms-results/public-results/publications/TOP-16-005/index.html
W_xsec_nominal = (2.+0.3521)*2.*815.*0.1086 #The two factors 2 account for the possible charge signs of the Ws and for the two leptonic decay channels of the tag W. The +0.3521 accounts for the leptonic tau decays
#W_xsec_nominal = (0.3521)*2.*815.*0.1086 #The two factors 2 account for the possible charge signs of the Ws and for the two leptonic decay channels of the tag W. The +0.3521 accounts for the leptonic tau decays
W_xsec_syst    = (43.*2.*(2.+0.3521)*0.1086)/W_xsec_nominal
#W_xsec_syst    = (43.*2.*(0.3521)*0.1086)/W_xsec_nominal

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

totsig_2016 = 80000.*(2.3521/3.) #total number of signal events in 2016
totmu_2016 = 3851.
totel_2016 = 2193.
tot_2016  = totmu_2016 + totel_2016

totsig_2017 = 79985.*(2.3521/3.) #total number of signal events in 2017
totmu_2017 = 3587.
totel_2017 = 2234.
tot_2017  = totmu_2017 + totel_2017

totsig_2018 = 79905.*(2.3521/3.) #total number of signal events in 2018
totmu_2018 = 3719.
totel_2018 = 2170.
tot_2018  = totmu_2018 + totel_2018

#The uncertainties are squared because they will be summed in quadrature
eff_nominal_2016 = tot_2016/totsig_2016
binom_eff_2016 = ((1 - eff_nominal_2016)/totsig_2016)**2 
BDT_syst_2016    = ((0.01*totmu_2016 + 0.02*totel_2016)/tot_2016)**2 
Pythia_syst_pT_2016 = ((0.02*totmu_2016 + 0.04*totel_2016)/tot_2016)**2 
Pythia_syst_angle_2016 = ((0.03*totmu_2016 + 0.06*totel_2016)/tot_2016)**2 
SF_syst_2016 = ((0.014*totmu_2016 + 0.014*totel_2016)/tot_2016)**2
TRK_mischarge_ID_syst_2016 = ((0.01*totmu_2016 + 0.01*totel_2016)/tot_2016)**2
eff_syst_2016    = math.sqrt(binom_eff_2016+BDT_syst_2016+Pythia_syst_pT_2016+Pythia_syst_angle_2016+SF_syst_2016+TRK_mischarge_ID_syst_2016)

eff_nominal_2017 = tot_2017/totsig_2017
binom_eff_2017 = ((1 - eff_nominal_2017)/totsig_2017)**2 
BDT_syst_2017    = ((0.01*totmu_2017 + 0.02*totel_2017)/tot_2017)**2 
Pythia_syst_pT_2017 = ((0.02*totmu_2017 + 0.04*totel_2017)/tot_2017)**2 
Pythia_syst_angle_2017 = ((0.03*totmu_2017 + 0.06*totel_2017)/tot_2017)**2
SF_syst_2017 = ((0.014*totmu_2017 + 0.014*totel_2017)/tot_2017)**2 
TRK_mischarge_ID_syst_2017 = ((0.01*totmu_2017 + 0.01*totel_2017)/tot_2017)**2
eff_syst_2017    = math.sqrt(binom_eff_2017+BDT_syst_2017+Pythia_syst_pT_2017+Pythia_syst_angle_2017+SF_syst_2017+TRK_mischarge_ID_syst_2017)

eff_nominal_2018 = tot_2018/totsig_2018
binom_eff_2018 = ((1 - eff_nominal_2018)/totsig_2018)**2
BDT_syst_2018    = ((0.01*totmu_2018 + 0.02*totel_2018)/tot_2018)**2
Pythia_syst_pT_2018 = ((0.02*totmu_2018 + 0.04*totel_2018)/tot_2018)**2 
Pythia_syst_angle_2018 = ((0.03*totmu_2018 + 0.06*totel_2018)/tot_2018)**2 
SF_syst_2018 = ((0.014*totmu_2018 + 0.014*totel_2018)/tot_2018)**2 
TRK_mischarge_ID_syst_2018 = ((0.01*totmu_2018 + 0.01*totel_2018)/tot_2018)**2
eff_syst_2018    = math.sqrt(binom_eff_2018+BDT_syst_2018+Pythia_syst_pT_2018+Pythia_syst_angle_2018+SF_syst_2018+TRK_mischarge_ID_syst_2018)

efflumi_2016_nominal = eff_nominal_2016*lumi_2016
efflumi_2017_nominal = eff_nominal_2017*lumi_2017
efflumi_2018_nominal = eff_nominal_2018*lumi_2018

efflumi_2016_syst = (eff_syst_2016+lumi_syst_2016)*efflumi_2016_nominal
efflumi_2017_syst = (eff_syst_2017+lumi_syst_2017)*efflumi_2017_nominal
efflumi_2018_syst = (eff_syst_2018+lumi_syst_2018)*efflumi_2018_nominal

efflumi_nominal = efflumi_2016_nominal + efflumi_2017_nominal + efflumi_2018_nominal
efflumi_syst = (efflumi_2016_syst + efflumi_2017_syst + efflumi_2018_syst)/efflumi_nominal

################################################################
#                                                              #
#------------- Systematic on bkg parametrization --------------#
#                                                              #
################################################################

if not suppressSigSystematic:
    sig_syst = 0.069
else:
    sig_syst = 0.0001

if not suppressBkgSystematic:
    bkg_syst = 0.146
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
    W_pigamma_BR = ROOT.RooRealVar("W_pigamma_BR","W_pigamma_BR",0.000005,-0.0001,0.01) # The parameter of interest can go negative when we try the alternative bkg description, to obtain the correct statistic coverage
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
#---------------------------- Fit -----------------------------#
#                                                              #
################################################################

constrained_params = ROOT.RooArgSet()
constrained_params.add(Multi_param_beta)

if not suppressAllSystematics:
    result_dataFit = totPDF.fitTo(data, ROOT.RooFit.Extended(1), ROOT.RooFit.Constrain(constrained_params), ROOT.RooFit.Save() )#For the signal region, I want the fit to be extended (Poisson fluctuation of number of events) to take into account that the total number of events is the sum of signal and background events. Either I do this, or I use a fraction frac*Nbkg+(1-frac)*Nsig, which will become a parameter of the fit and will have a Gaussian behavior (whilst the extended fit preserves the natural Poisson behavior)
else:
    result_dataFit = totPDF.fitTo(data, ROOT.RooFit.Extended(1), ROOT.RooFit.Save() )

################################################################
#                                                              #
#------------------------ Write to file -----------------------#
#                                                              #
################################################################
        
#Save the fit into a root file
fOutput = ROOT.TFile("fitData.root","RECREATE")

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

#xframe = Wmass.frame(55.,95.,15)
xframe = Wmass.frame(50.,100.,20)
xframe.SetTitle(" ")
xframe.SetTitleOffset(1.4,"y")

#################################################

data_reduced = data.reduce("Wmass < 65. || Wmass > 90.")

#data_reduced.plotOn(xframe, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
data.plotOn(xframe, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
totPDF.plotOn(xframe)#, ROOT.RooFit.Range("LowSideband,HighSideband"))#, ROOT.RooFit.LineColor(ROOT.kTeal+10))

chi2 = xframe.chiSquare()
cut_chi2 = "{:.2f}".format(chi2) #Crop the chi2 to 3 decimal digits
label = ROOT.TPaveLabel(0.68,0.4,0.88,0.54,"#chi^{2} = " + cut_chi2,"brNDC")

#fIn_Wmass = ROOT.TFile("../Wmass_4Plotting.root")
#h_Wmass_signal = fIn_Wmass.Get("h_Wmass;1")

#DRAW ON CANVAS
canvas = ROOT.TCanvas()
canvas.cd()
#label.Draw("SAME")
xframe.Draw()

fOut_frame = ROOT.TFile("Wmass_frame.root","RECREATE")
fOut_frame.cd()
xframe.Write("Wmass_frame")

#h_stack.Draw("same,hist")
#h_Wmass_signal.Draw("same,hist")
# Save the plot
canvas.SaveAs("plots/fitData_signalR.pdf")

del workspace_out

raw_input()
