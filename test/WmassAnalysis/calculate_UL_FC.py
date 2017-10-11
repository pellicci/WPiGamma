
import ROOT

#Get the model and the data
fInput = ROOT.TFile("fitData.root")
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
#fc.SetNBins(10)

#Let him know which parameter points to test. This is passed in a RooDataSet with the values of the parameter points
scan_params = poi.snapshot()
points_to_scan = ROOT.RooDataSet("points_to_scan","points_to_scan",scan_params)
scan_params.setRealValue("W_pigamma_BR",0.)
points_to_scan.add(scan_params)
#scan_params.setRealValue("W_pigamma_BR",0.000007)   #7*10-6
#points_to_scan.add(scan_params)
#scan_params.setRealValue("W_pigamma_BR",0.0000075)  #7.5*10-6
#points_to_scan.add(scan_params)
#scan_params.setRealValue("W_pigamma_BR",0.000008)   #8*10-6
#points_to_scan.add(scan_params)
#scan_params.setRealValue("W_pigamma_BR",0.000009)   #9*10-6
#points_to_scan.add(scan_params)
#scan_params.setRealValue("W_pigamma_BR",0.00001)    #1*10-5
#points_to_scan.add(scan_params)
#scan_params.setRealValue("W_pigamma_BR",0.000011)    #1.1*10-5
#points_to_scan.add(scan_params)
#scan_params.setRealValue("W_pigamma_BR",0.000012)    #1.2*10-5
#points_to_scan.add(scan_params)
#scan_params.setRealValue("W_pigamma_BR",0.000013)    #1.3*10-5
#points_to_scan.add(scan_params)
#scan_params.setRealValue("W_pigamma_BR",0.000014)    #1.4*10-5
#points_to_scan.add(scan_params)
#scan_params.setRealValue("W_pigamma_BR",0.000015)    #1.5*10-5
#points_to_scan.add(scan_params)
#scan_params.setRealValue("W_pigamma_BR",0.000016)    #1.6*10-5
#points_to_scan.add(scan_params)
#scan_params.setRealValue("W_pigamma_BR",0.000017)    #1.7*10-5
#points_to_scan.add(scan_params)
#scan_params.setRealValue("W_pigamma_BR",0.000018)    #1.8*10-5
#points_to_scan.add(scan_params)
#scan_params.setRealValue("W_pigamma_BR",0.0000185)    #1.85*10-5
#points_to_scan.add(scan_params)
#scan_params.setRealValue("W_pigamma_BR",0.000019)    #1.9*10-5
#points_to_scan.add(scan_params)
scan_params.setRealValue("W_pigamma_BR",0.00002)     #2.0*10-5
points_to_scan.add(scan_params)
scan_params.setRealValue("W_pigamma_BR",0.00003)     #3.0*10-5
points_to_scan.add(scan_params)
scan_params.setRealValue("W_pigamma_BR",0.000031)     #3.1*10-5
points_to_scan.add(scan_params)
scan_params.setRealValue("W_pigamma_BR",0.000032)     #3.2*10-5
points_to_scan.add(scan_params)
scan_params.setRealValue("W_pigamma_BR",0.000033)     #3.3*10-5
points_to_scan.add(scan_params)
scan_params.setRealValue("W_pigamma_BR",0.000034)     #3.4*10-5
points_to_scan.add(scan_params)
scan_params.setRealValue("W_pigamma_BR",0.000035)    #3.5*10-5
points_to_scan.add(scan_params)

fc.SetPOIPointsToTest(points_to_scan)

#We can use PROOF to speed things along in parallel
pc = ROOT.RooStats.ProofConfig(workspace, 0, "workers=6",0)
toymcsampler = fc.GetTestStatSampler()
toymcsampler.SetProofConfig(pc)

interval = fc.GetInterval()


print "The interval is [", interval.LowerLimit(W_pigamma_BR), ", ", interval.UpperLimit(W_pigamma_BR), "]"

del fc
