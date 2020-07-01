import ROOT

ROOT.gROOT.ProcessLineSync(".L dCB/RooDoubleCBFast.cc+")

suppressAllSystematics = False #Choose if you want to calculate the limit without systematic uncertainties

#Get the model and the data
fInput = ROOT.TFile("fitData.root")
fInput.cd()

workspace = fInput.Get("workspace")
workspace.Print()

data = workspace.data("data")

#Define the parameter of interest
W_pigamma_BR = workspace.var("W_pigamma_BR")
poi = ROOT.RooArgSet(W_pigamma_BR)

#Define observables
Wmass = workspace.var("Wmass")
Categorization = workspace.cat("Categorization")#All the categories added with defineType in the fit MUST be used for the fit itself (in fit_noCats.py), otherwise: SegFault in the UL calculation

observables = ROOT.RooArgSet()
observables.add(Wmass)
observables.add(Categorization)

#Define nuisances
if not suppressAllSystematics:
    constrained_params = ROOT.RooArgSet()
    constrained_params.add(workspace.var("Multi_param_beta"))
    constrained_params.add(workspace.var("Nbkg"))
    constrained_params.add(workspace.var("a0_bkg"))
        
    #Define global observables
    global_params = ROOT.RooArgSet()

    workspace.var("glb_Multi_param").setConstant(1)

    global_params.add(workspace.var("glb_Multi_param"))
else:
    workspace.var("Multi_param_beta").setConstant(1)

#Define the model container
#First the S+B
sbModel = ROOT.RooStats.ModelConfig(workspace)
sbModel.SetObservables(observables)
if not suppressAllSystematics:
    sbModel.SetNuisanceParameters(constrained_params)
    sbModel.SetGlobalObservables(global_params)
sbModel.SetPdf("totPDF")
sbModel.SetName("S+B Model")
sbModel.SetParametersOfInterest(poi)
sbModel.SetSnapshot(poi)

#Then B
bModel = sbModel.Clone()
bModel.SetObservables(observables)
if not suppressAllSystematics:
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

#use the CLs method with the asymptotic calculator
#See here https://arxiv.org/pdf/1007.1727.pdf

fc = ROOT.RooStats.AsymptoticCalculator(data, bModel, sbModel)
fc.SetOneSided(1)

#Create hypotest inverter passing the desired calculator 
calc = ROOT.RooStats.HypoTestInverter(fc)

#Set confidence level (e.g. 95% upper limits)
calc.SetConfidenceLevel(0.95)
calc.SetVerbose(0)

#use CLs
calc.UseCLs(1)

toymc = fc.GetTestStatSampler()

#Set profile likelihood test statistics
profl = ROOT.RooStats.ProfileLikelihoodTestStat(sbModel.GetPdf())
#For CLs (bounded intervals) use one-sided profile likelihood
profl.SetOneSided(1)
#Set the test statistic to use
toymc.SetTestStatistic(profl)

npoints = 100 #Number of points to scan
# min and max for the scan (better to choose smaller intervals)
poimin = poi.find("W_pigamma_BR").getMin()
poimax = poi.find("W_pigamma_BR").getMax()

min_scan = 0.0000001
max_scan = 0.000055
#max_scan = 0.0001

print "Doing a fixed scan  in interval : ",min_scan, " , ", max_scan
calc.SetFixedScan(npoints,min_scan,max_scan)

result = calc.GetInterval() #This is a HypoTestInveter class object

##################################################################
#                                                                #
#--------------------- Only when unblinded ----------------------#
#                                                                #
##################################################################

upperLimit = result.UpperLimit()

#Now print the result of the two methods
print "################"
print "The observed CLs upper limit is: ", upperLimit

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
freq_plot.Draw()
#freq_plot.Draw("EXP")
canvas.SaveAs("UL_CLs.pdf")

del fc
