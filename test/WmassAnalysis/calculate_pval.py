import ROOT

fInput = ROOT.TFile("fitData.root")
fInput.cd()

workspace = fInput.Get("workspace")

workspace.Print()

#Define the model container
model = ROOT.RooStats.ModelConfig()
model.SetWorkspace(workspace)
model.SetPdf("totPDF")

#Get the parameter of interest and set it for the null hypothesis
W_pigamma_BR = workspace.var("W_pigamma_BR")
poi = ROOT.RooArgSet(W_pigamma_BR)
nullParams = poi.snapshot()
nullParams.setRealValue("W_pigamma_BR",0.)

#Build the profile likelihood calculator
plc = ROOT.RooStats.ProfileLikelihoodCalculator()

plc.SetData(workspace.data("data"))
plc.SetModel(model)
plc.SetParameters(poi)
plc.SetNullParameters(nullParams)

#The the Hypotest result
htr = plc.GetHypoTest()

print "-------------------------------------------------"
print "The p-value for the null is ", htr.NullPValue()
print "Corresponding to a signifcance of ", htr.Significance()
print "-------------------------------------------------"

#PyROOT sometimes fails cleaning memory, this helps
del plc
