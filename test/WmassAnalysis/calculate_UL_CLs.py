import ROOT

ROOT.gROOT.ProcessLineSync(".L dCB/RooDoubleCBFast.cc+")

#Get the model and the data
#fInput = ROOT.TFile("fitMC.root")
fInput = ROOT.TFile("fitData.root")
fInput.cd()

workspace = fInput.Get("workspace")
workspace.Print()
workspace.var("W_pigamma_BR").setRange(0.0000001,0.0001)

#Define the parameter of interest
W_pigamma_BR = workspace.var("W_pigamma_BR")
poi = ROOT.RooArgSet(W_pigamma_BR)

#Define observables
Wmass = workspace.var("Wmass")
Categorization = workspace.cat("Categorization")

observables = ROOT.RooArgSet()
observables.add(Wmass)
observables.add(Categorization)


#Define nuisances
constrained_params = ROOT.RooArgSet()
constrained_params.add(workspace.var("dCB_width"))
constrained_params.add(workspace.var("eta_mu_2016"))
constrained_params.add(workspace.var("eta_el_2016"))
constrained_params.add(workspace.var("eta_mu_2017"))
constrained_params.add(workspace.var("eta_el_2017"))
constrained_params.add(workspace.var("W_xsec_constr"))
constrained_params.add(workspace.var("lumi_constr_2016"))
constrained_params.add(workspace.var("lumi_constr_2017"))
constrained_params.add(workspace.var("eff_mu_constr_2016"))
constrained_params.add(workspace.var("eff_el_constr_2016"))
constrained_params.add(workspace.var("eff_mu_constr_2017"))
constrained_params.add(workspace.var("eff_el_constr_2017"))

#Define global observables
global_params = ROOT.RooArgSet()
global_params.add(workspace.var("dCB_width_constr"))
global_params.add(workspace.var("glb_W_xsec"))
global_params.add(workspace.var("glb_lumi_2016"))
global_params.add(workspace.var("glb_lumi_2017"))
global_params.add(workspace.var("glb_eff_mu_2016"))
global_params.add(workspace.var("glb_eff_el_2016"))
global_params.add(workspace.var("glb_eff_mu_2017"))
global_params.add(workspace.var("glb_eff_el_2017"))

global_params.add(workspace.var("glb_bkg_param_mu_2016"))#Bisogna metterli?
global_params.add(workspace.var("glb_bkg_param_el_2016"))
global_params.add(workspace.var("glb_bkg_param_mu_2017"))
global_params.add(workspace.var("glb_bkg_param_el_2017"))


#Define the model container
#First the S+B
sbModel = ROOT.RooStats.ModelConfig()
sbModel.SetWorkspace(workspace)
sbModel.SetObservables(observables)
sbModel.SetNuisanceParameters(constrained_params)
sbModel.SetGlobalObservables(global_params)
sbModel.SetPdf("totPDF")
sbModel.SetName("S+B Model")
sbModel.SetParametersOfInterest(poi)

bModel = sbModel.Clone()
bModel.SetObservables(observables)
bModel.SetNuisanceParameters(constrained_params)
bModel.SetGlobalObservables(global_params)
bModel.SetPdf("totPDF")
bModel.SetName("B model")
bModel.SetParametersOfInterest(poi)
oldval = poi.find("W_pigamma_BR").getVal()
poi.find("W_pigamma_BR").setVal(0)
bModel.SetSnapshot(poi)
poi.find("W_pigamma_BR").setVal(oldval)

print "Number of events in data = ", workspace.data("data").numEntries()


#use the CLs method with the asymptotic calculator, using Asimov datasets
#See here https://arxiv.org/pdf/1007.1727.pdf

#----------------------------------------------------------------------------------#
# fc = ROOT.RooStats.FrequentistCalculator(workspace.data("data"), bModel, sbModel)
# fc.SetToys(500,500)

# #Create hypotest inverter passing desired calculator
# calc = ROOT.RooStats.HypoTestInverter(fc)

# calc.SetConfidenceLevel(0.95)

# #Use CLs
# calc.UseCLs(1)

# calc.SetVerbose(0)

# #Configure ToyMC sampler
# #toymc = calc.GetHypoTestCalculator().GetTestStatSampler()
# toymc = fc.GetTestStatSampler()

# #Set profile likelihood test statistics
# profl = ROOT.RooStats.ProfileLikelihoodTestStat(sbModel.GetPdf())
# #For CLs (bounded intervals) use one-sided profile likelihood
# profl.SetOneSided(1)

# #Set the test statistic to use
# toymc.SetTestStatistic(profl)

#----------------------------------------------------------------------------------#

data = workspace.data("data")

fc = ROOT.RooStats.AsymptoticCalculator(data, bModel, sbModel,0)
fc.SetOneSided(1)
#fc.UseSameAltToys()

#Create hypotest inverter passing the desired calculator 
calc = ROOT.RooStats.HypoTestInverter(fc)

#set confidence level (e.g. 95% upper limits)
calc.SetConfidenceLevel(0.95)
calc.SetVerbose(0)

#use CLs
calc.UseCLs(1)

npoints = 30 #Number of points to scan
# min and max for the scan (better to choose smaller intervals)
poimin = poi.find("W_pigamma_BR").getMin()
poimax = poi.find("W_pigamma_BR").getMax()

print "Doing a fixed scan  in interval : ", poimin, " , ", poimax
calc.SetFixedScan(npoints,poimin,0.000006)

# In order to use PROOF, one should install the test statistic on the workers
# pc = ROOT.RooStats.ProofConfig(workspace, 0, "workers=6",0)
# toymc.SetProofConfig(pc)

result = calc.GetInterval() #This is a HypoTestInveter class object

##################################################################
#                                                                #
#--------------------- Only when unblinded ----------------------#
#                                                                #
##################################################################

# upperLimit = result.UpperLimit()

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
