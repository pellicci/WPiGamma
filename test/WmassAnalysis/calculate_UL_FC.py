
import ROOT

#Get the model and the data
fInput = ROOT.TFile("fitAllLep.root")
fInput.cd()

workspace = fInput.Get("workspace")

workspace.Print()

workspace.var("W_pigamma_BR").setRange(0.,0.00001)

#Define the parameter of interest
W_pigamma_BR = workspace.var("W_pigamma_BR")
poi = ROOT.RooArgSet(W_pigamma_BR)

#Define the model container
model = ROOT.RooStats.ModelConfig()
model.SetWorkspace(workspace)
model.SetPdf("totPDF")
model.SetParametersOfInterest(poi)

print "Number of events in data = ", workspace.data("data").numEntries()

#Set up the FC calculator
fc = ROOT.RooStats.FeldmanCousins(workspace.data("data"),model)
fc.SetTestSize(0.05)
fc.UseAdaptiveSampling(1)
fc.FluctuateNumDataEntries(1)
fc.SetNBins(10)

#We can use PROOF to speed things along in parallel
pc = ROOT.RooStats.ProofConfig(workspace, 0, "workers=4", 0)
toymcsampler = fc.GetTestStatSampler()
toymcsampler.SetProofConfig(pc)

interval = fc.GetInterval()


print "The interval is [", interval.LowerLimit(W_pigamma_BR), ", ", interval.UpperLimit(W_pigamma_BR), "]"

del fc
