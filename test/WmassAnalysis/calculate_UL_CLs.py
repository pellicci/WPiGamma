
import ROOT

#Get the model and the data
fInput = ROOT.TFile("fitAllLep.root")
fInput.cd()

workspace = fInput.Get("workspace")
workspace.Print()
workspace.var("W_pigamma_BR").setRange(0.,0.001)

#Set a few things constant
workspace.var("a0_el").setConstant(1)
workspace.var("a1_el").setConstant(1)
workspace.var("a0_mu").setConstant(1)
workspace.var("a1_mu").setConstant(1)
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

#Define the model container
#First the S+B
sbModel = ROOT.RooStats.ModelConfig()
sbModel.SetWorkspace(workspace)
sbModel.SetObservables(observables)
sbModel.SetPdf("totPDF")
sbModel.SetName("S+B Model")
sbModel.SetParametersOfInterest(poi)

bModel = sbModel.Clone()
bModel.SetPdf("totPDF")
sbModel.SetObservables(observables)
bModel.SetParametersOfInterest(poi)
bModel.SetName("B model")
oldval = poi.find("W_pigamma_BR").getVal()
poi.find("W_pigamma_BR").setVal(0)
bModel.SetSnapshot(poi)
poi.find("W_pigamma_BR").setVal(oldval)

print "Number of events in data = ", workspace.data("data").numEntries()

#use the CLs method with the asymptotic calculator, using Asimov datasets
#See here https://arxiv.org/pdf/1007.1727.pdf

fc = ROOT.RooStats.AsymptoticCalculator(workspace.data("data"), bModel, sbModel,0)
fc.SetOneSided(1)
#fc.UseSameAltToys()

#Configure ToyMC Samler
toymcs = fc.GetTestStatSampler()

#Create hypotest inverter passing the desired calculator 
calc = ROOT.RooStats.HypoTestInverter(fc)

#set confidence level (e.g. 95% upper limits)
calc.SetConfidenceLevel(0.95)

#use CLs
calc.UseCLs(1)

#We can use PROOF to speed things along in parallel
pc = ROOT.RooStats.ProofConfig(workspace, 0, "", 0)
toymcs.SetProofConfig(pc)

#Use profile likelihood as test statistics 
profll = ROOT.RooStats.ProfileLikelihoodTestStat(sbModel.GetPdf())

#for CLs (bounded intervals) use one-sided profile likelihood
profll.SetOneSided(1)

#set the test statistic to use for toys
toymcs.SetTestStatistic(profll)

npoints = 10 #8 #Number of points to scan
#x min and max for the scan (better to choose smaller intervals)
poimin = poi.find("W_pigamma_BR").getMin()
poimax = poi.find("W_pigamma_BR").getMax()

print "Doing a fixed scan  in interval : ", poimin, " , ", poimax
calc.SetFixedScan(npoints,poimin,0.0001);

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
#freq_plot = ROOT.RooStats.HypoTestInverterPlot("HTI_Result_Plot","Frequentist scan result for the W -> pigamma BR",result)

#Plot in a new canvas with style
#canvas = ROOT.TCanvas()
#canvas.cd()
#freq_plot.Draw("2CL")
#canvas.SaveAs("calculate_UL_CLs.png")
