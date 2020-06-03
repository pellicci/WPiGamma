import ROOT
ROOT.gROOT.ProcessLineSync(".L dCB/RooDoubleCBFast.cc+")

#Get the model and the data
fInput = ROOT.TFile("fitData.root")
fInput.cd()

workspace = fInput.Get("workspace")
workspace.Print()
#workspace.var("W_pigamma_BR").setRange(-0.0001,0.01)

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

constrained_params = ROOT.RooArgSet()
constrained_params.add(workspace.var("Multi_param_beta"))
constrained_params.add(workspace.var("Nbkg"))
constrained_params.add(workspace.var("a0_bkg"))

#Define global observables
global_params = ROOT.RooArgSet()

workspace.var("glb_Multi_param").setConstant(1)

global_params.add(workspace.var("glb_Multi_param"))

model = ROOT.RooStats.ModelConfig(workspace)
model.SetObservables(observables)

model.SetNuisanceParameters(constrained_params)
model.SetGlobalObservables(global_params)
model.SetPdf("totPDF")
model.SetName("S+B Model")
model.SetParametersOfInterest(poi)
model.SetSnapshot(poi)

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
