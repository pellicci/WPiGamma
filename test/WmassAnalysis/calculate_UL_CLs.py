import ROOT
import argparse

ROOT.gROOT.ProcessLineSync(".L dCB/RooDoubleCBFast.cc+")

#---------------------------------#
#p = argparse.ArgumentParser(description='Select whether to fit data or MC')
#p.add_argument('runningEra_option', help='Type <<0>> for 2016, <<1>> for 2017, <<2>> for 2016+2017, <<3>> for 2016+2017+2018')
#args = p.parse_args()

#runningEra = int(args.runningEra_option)

suppressAllSystematics = False
#---------------------------------#

#Get the model and the data
fInput = ROOT.TFile("fitData.root")
fInput.cd()

workspace = fInput.Get("workspace")
workspace.Print()
workspace.var("W_pigamma_BR").setRange(0.,0.0001)

#Define the parameter of interest
W_pigamma_BR = workspace.var("W_pigamma_BR")
poi = ROOT.RooArgSet(W_pigamma_BR)

#Define observables
Wmass = workspace.var("Wmass")
Categorization = workspace.cat("Categorization")#All the categories added with defineType in the fit MUST be used for the fit itself (in fitAll_withControlRegions.py), otherwise: SegFault in the UL calculation

observables = ROOT.RooArgSet()
observables.add(Wmass)
observables.add(Categorization)


#Define nuisances
if not suppressAllSystematics:
    constrained_params = ROOT.RooArgSet()
    constrained_params.add(workspace.var("dCB_width"))
    constrained_params.add(workspace.var("dCB_pole"))
    constrained_params.add(workspace.var("W_xsec_beta"))
    constrained_params.add(workspace.var("lumi_beta_2016"))
    constrained_params.add(workspace.var("lumi_beta_2017"))
    constrained_params.add(workspace.var("lumi_beta_2018"))
    constrained_params.add(workspace.var("eff_mu_beta_2016"))
    constrained_params.add(workspace.var("eff_el_beta_2016"))
    constrained_params.add(workspace.var("eff_mu_beta_2017"))
    constrained_params.add(workspace.var("eff_el_beta_2017"))
    constrained_params.add(workspace.var("eff_mu_beta_2018"))
    constrained_params.add(workspace.var("eff_el_beta_2018"))
    constrained_params.add(workspace.var("bkg_param_beta_mu"))
    constrained_params.add(workspace.var("bkg_param_beta_el"))
        
    #Define global observables
    global_params = ROOT.RooArgSet()

    workspace.var("dCB_width_constr").setConstant(1)
    workspace.var("dCB_pole_constr").setConstant(1)
    workspace.var("glb_W_xsec").setConstant(1)
    workspace.var("glb_lumi_2016").setConstant(1)
    workspace.var("glb_eff_mu_2016").setConstant(1)
    workspace.var("glb_eff_el_2016").setConstant(1)
    workspace.var("glb_lumi_2017").setConstant(1)
    workspace.var("glb_eff_mu_2017").setConstant(1)
    workspace.var("glb_eff_el_2017").setConstant(1)
    workspace.var("glb_lumi_2018").setConstant(1)
    workspace.var("glb_eff_mu_2018").setConstant(1)
    workspace.var("glb_eff_el_2018").setConstant(1)
    workspace.var("glb_bkg_param_mu").setConstant(1)
    workspace.var("glb_bkg_param_el").setConstant(1)

    global_params.add(workspace.var("dCB_width_constr"))
    global_params.add(workspace.var("dCB_pole_constr"))
    global_params.add(workspace.var("glb_W_xsec"))
    global_params.add(workspace.var("glb_lumi_2016"))
    global_params.add(workspace.var("glb_eff_mu_2016"))
    global_params.add(workspace.var("glb_eff_el_2016"))
    global_params.add(workspace.var("glb_lumi_2017"))
    global_params.add(workspace.var("glb_eff_mu_2017"))
    global_params.add(workspace.var("glb_eff_el_2017"))
    global_params.add(workspace.var("glb_lumi_2018"))
    global_params.add(workspace.var("glb_eff_mu_2018"))
    global_params.add(workspace.var("glb_eff_el_2018"))
    global_params.add(workspace.var("glb_bkg_param_mu"))
    global_params.add(workspace.var("glb_bkg_param_el"))

else:
    workspace.var("dCB_width").setConstant(1)
    workspace.var("dCB_pole").setConstant(1)
    workspace.var("W_xsec_beta").setConstant(1)
    workspace.var("lumi_beta_2016").setConstant(1)
    workspace.var("lumi_beta_2017").setConstant(1)
    workspace.var("lumi_beta_2018").setConstant(1)
    workspace.var("eff_mu_beta_2016").setConstant(1)
    workspace.var("eff_el_beta_2016").setConstant(1)
    workspace.var("eff_mu_beta_2017").setConstant(1)
    workspace.var("eff_el_beta_2017").setConstant(1)
    workspace.var("eff_mu_beta_2018").setConstant(1)
    workspace.var("eff_el_beta_2018").setConstant(1)
    workspace.var("bkg_param_beta_mu").setConstant(1)
    workspace.var("bkg_param_beta_el").setConstant(1)

#Define the model container
#First the S+B
sbModel = ROOT.RooStats.ModelConfig()
sbModel.SetWorkspace(workspace)
sbModel.SetObservables(observables)
if not suppressAllSystematics:
    sbModel.SetNuisanceParameters(constrained_params)
    sbModel.SetGlobalObservables(global_params)
sbModel.SetPdf("totPDF")
#sbModel.SetPdf("totPDF_el_2016")
sbModel.SetName("S+B Model")
sbModel.SetParametersOfInterest(poi)

bModel = sbModel.Clone()
bModel.SetObservables(observables)
if not suppressAllSystematics:
    bModel.SetNuisanceParameters(constrained_params)
    bModel.SetGlobalObservables(global_params)
bModel.SetPdf("totPDF")
#bModel.SetPdf("totPDF_el_2016")
bModel.SetName("B model")
bModel.SetParametersOfInterest(poi)
oldval = poi.find("W_pigamma_BR").getVal()
poi.find("W_pigamma_BR").setVal(0) #BEWARE that the range of the POI (W_pigamma_BR.setRange) has to contain zero!
bModel.SetSnapshot(poi)
poi.find("W_pigamma_BR").setVal(oldval)

print "Number of events in data = ", workspace.data("data").numEntries()


#use the CLs method with the asymptotic calculator, using Asimov datasets
#See here https://arxiv.org/pdf/1007.1727.pdf

#----------------------------------------------------------------------------------#
fc = ROOT.RooStats.FrequentistCalculator(workspace.data("data"), bModel, sbModel)
fc.SetToys(800,800)
#fc = ROOT.RooStats.AsymptoticCalculator(workspace.data("data"), bModel, sbModel)


#Create hypotest inverter passing desired calculator
calc = ROOT.RooStats.HypoTestInverter(fc)

calc.SetConfidenceLevel(0.95)

#Use CLs
calc.UseCLs(1)

calc.SetVerbose(0)

#Configure ToyMC sampler
#toymc = calc.GetHypoTestCalculator().GetTestStatSampler()
toymc = fc.GetTestStatSampler()

#Set profile likelihood test statistics
profl = ROOT.RooStats.ProfileLikelihoodTestStat(sbModel.GetPdf())
#For CLs (bounded intervals) use one-sided profile likelihood
profl.SetOneSided(1)

#Set the test statistic to use
toymc.SetTestStatistic(profl)

#----------------------------------------------------------------------------------#

#data = workspace.data("data")

# #data = data_initial.reduce("Categorization==Categorization::ElectronSignal_2016")

# fc = ROOT.RooStats.AsymptoticCalculator(data, bModel, sbModel,0)#,15)
# fc.SetOneSided(1)
# #fc.UseSameAltToys()

# #Create hypotest inverter passing the desired calculator 
# calc = ROOT.RooStats.HypoTestInverter(fc)

# #set confidence level (e.g. 95% upper limits)
# calc.SetConfidenceLevel(0.95)
# calc.SetVerbose(1)

# #use CLs
# calc.UseCLs(1)

npoints = 25 #Number of points to scan
# min and max for the scan (better to choose smaller intervals)
poimin = poi.find("W_pigamma_BR").getMin()
poimax = poi.find("W_pigamma_BR").getMax()

min_scan = 0.000004
max_scan = 0.000030
#min_scan = 0.0000005
#max_scan = 0.000003
print "Doing a fixed scan  in interval : ",min_scan, " , ", max_scan
calc.SetFixedScan(npoints,min_scan,max_scan)

# In order to use PROOF, one should install the test statistic on the workers
# pc = ROOT.RooStats.ProofConfig(workspace, 0, "workers=6",0)
# toymc.SetProofConfig(pc)

result = calc.GetInterval() #This is a HypoTestInveter class object

##################################################################
#                                                                #
#--------------------- Only when unblinded ----------------------#
#                                                                #
##################################################################

upperLimit = result.UpperLimit()

#Now let's print the result of the two methods
# print "################"
# print "The observed CLs upper limit is: ", upperLimit

##################################################################

#Compute expected limit
print "Expected upper limits, using the B (alternate) model : "
print " expected limit (median) ", result.GetExpectedUpperLimit(0)
print " expected limit (-1 sig) ", result.GetExpectedUpperLimit(-1)
print " expected limit (+1 sig) ", result.GetExpectedUpperLimit(1)
print "################"

#Plot the results
freq_plot = ROOT.RooStats.HypoTestInverterPlot("HTI_Result_Plot","Frequentist scan result for the W -> pigamma BR",result)

#xPlot in a new canvas with style
canvas = ROOT.TCanvas()
canvas.cd()
#freq_plot.Draw("2CL")
freq_plot.Draw("EXP")
# freq_plot.GetYaxis().SetRangeUser(0.,0.8)
# freq_plot.GetXaxis().SetRange(0.,0.0000107)
canvas.SaveAs("UL_CLs.pdf")

del fc
