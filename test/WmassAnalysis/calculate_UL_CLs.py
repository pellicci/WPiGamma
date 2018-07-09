import ROOT

#Get the model and the data
fInput = ROOT.TFile("fitMC.root")
fInput.cd()

workspace = fInput.Get("workspace")
workspace.Print()
workspace.var("W_pigamma_BR").setRange(0.,0.0001)

#Set a few things constant
#workspace.var("a0_el").setConstant(1)
#workspace.var("a1_el").setConstant(1)
#workspace.var("a0_mu").setConstant(1)
#workspace.var("a1_mu").setConstant(1)
#workspace.var("Nbkg_el").setConstant(1)
#workspace.var("Nbkg_mu").setConstant(1)

#Define the parameter of interest
W_pigamma_BR = workspace.var("W_pigamma_BR")
poi = ROOT.RooArgSet(W_pigamma_BR)

#Define observables
Wmass = workspace.var("Wmass")
isMuon = workspace.cat("isMuon")
observables = ROOT.RooArgSet()
observables.add(Wmass)
observables.add(isMuon)

#Define nuisances
constrained_params = ROOT.RooArgSet()
constrained_params.add(workspace.var("dCB_width"))
constrained_params.add(workspace.var("lumi_constr"))
constrained_params.add(workspace.var("W_xsec_constr"))
constrained_params.add(workspace.var("eff_mu_constr"))
constrained_params.add(workspace.var("eff_el_constr"))

#Define global observables
global_params = ROOT.RooArgSet()
global_params.add(workspace.var("dCB_width_constr"))
global_params.add(workspace.var("glb_W_xsec"))
global_params.add(workspace.var("glb_lumi"))
global_params.add(workspace.var("glb_eff_mu"))
global_params.add(workspace.var("glb_eff_el"))

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

fc = ROOT.RooStats.AsymptoticCalculator(workspace.data("data"), bModel, sbModel,0)
fc.SetOneSided(1)
#fc.UseSameAltToys()

#Create hypotest inverter passing the desired calculator 
calc = ROOT.RooStats.HypoTestInverter(fc)

#set confidence level (e.g. 95% upper limits)
calc.SetConfidenceLevel(0.95)
calc.SetVerbose(0)

#use CLs
calc.UseCLs(1)

npoints = 50 #Number of points to scan
# min and max for the scan (better to choose smaller intervals)
poimin = poi.find("W_pigamma_BR").getMin()
poimax = poi.find("W_pigamma_BR").getMax()

print "Doing a fixed scan  in interval : ", poimin, " , ", poimax
calc.SetFixedScan(npoints,poimin,0.00001)

# In order to use PROOF, one should instal the test statistic on the workers
# pc = ROOT.RooStats.ProofConfig(workspace, 0, "workers=6",0)
# toymc.SetProofConfig(pc)

result = calc.GetInterval() #This is a HypoTestInveter class object
upperLimit = result.UpperLimit()

#Now let's print the result of the two methods
print "################"
print "The observed CLs upper limit is: ", upperLimit

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
freq_plot.Draw("2CL")
canvas.SaveAs("UL_CLs.pdf")

del fc
