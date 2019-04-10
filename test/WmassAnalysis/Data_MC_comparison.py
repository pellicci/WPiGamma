import ROOT
import math

ROOT.gROOT.ProcessLineSync(".L dCB/RooDoubleCBFast.cc+")

fData = ROOT.TFile("fitData.root")

Data = fData.Get("data")

#Define the observable
Wmass = ROOT.RooRealVar("Wmass","m_{#pi#gamma}",50.,100.,"GeV")

#Retrive the sample
fInput_data = ROOT.TFile("Tree_input_massfit_Data.root")

fInput_MC = ROOT.TFile("Tree_input_massfit_MC.root")

#Define the mu/ele category
Categorization = ROOT.RooCategory("Categorization","Categorization")
Categorization.defineType("MuonCR",0)
Categorization.defineType("MuonSignal",1)
Categorization.defineType("ElectronCR",2)
Categorization.defineType("ElectronSignal",3)

isSignal = ROOT.RooCategory("isSignal","isSignal")
isSignal.defineType("Signal",1)
isSignal.defineType("Background",0)

isMuon = ROOT.RooCategory("isMuon","isMuon")
isMuon.defineType("Muon",1)
isMuon.defineType("Electron",0)

minWeight_inMC = -50. # Set the minimum weight an event can have when processing MC. This is useful to remove spikes which cannot be fitted
maxWeight_inMC = 50.  # Set the maximum weight an event can have when processing MC. This is useful to remove spikes which cannot be fitted

#Define the event weight
weight = ROOT.RooRealVar("weight","The event weight",minWeight_inMC,maxWeight_inMC)

#Support variables
BDT_out = ROOT.RooRealVar("BDT_out","Output of BDT",-1.,1.)

#Retrieve data stuff
fInput_data.cd()

mytree_data = fInput_data.Get("minitree")

data_initial_data = ROOT.RooDataSet("data","data", ROOT.RooArgSet(Wmass,Categorization,BDT_out,isMuon), ROOT.RooFit.Import(mytree_data))
#data = ROOT.RooDataSet("data","data", ROOT.RooArgSet(Wmass,Categorization,BDT_out,isMuon), ROOT.RooFit.Import(mytree_data))

#Retrieve MC stuff
fInput_MC.cd()

mytree_MC = fInput_MC.Get("minitree")

data_initial_MC = ROOT.RooDataSet("MC","MC", ROOT.RooArgSet(Wmass,Categorization,weight,BDT_out,isSignal,isMuon), ROOT.RooFit.Import(mytree_MC), ROOT.RooFit.WeightVar("weight"))
#MC = ROOT.RooDataSet("MC","MC", ROOT.RooArgSet(Wmass,Categorization,weight,BDT_out,isSignal,isMuon), ROOT.RooFit.Import(mytree_MC), ROOT.RooFit.WeightVar("weight"))

data = data_initial_data.reduce("(BDT_out > 0.168 && isMuon==isMuon::Muon) || (BDT_out > 0.112 && isMuon==isMuon::Electron)")
MC = data_initial_MC.reduce("(BDT_out > 0.168 && isMuon==isMuon::Muon) || (BDT_out > 0.112 && isMuon==isMuon::Electron)")
#data = data_initial_data.reduce("BDT_out > -3")
#MC = data_initial_MC.reduce("BDT_out > -3")

xframe_data_comp_mu = Wmass.frame(55.,95.,10)
#xframe_data_comp_mu = Wmass.frame(50.,100.,15)
xframe_data_comp_mu.SetTitle(" ")

xframe_data_comp_mu.SetTitleOffset(1.5,"y")
data.plotOn(xframe_data_comp_mu, ROOT.RooFit.Cut("Categorization==0"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson),ROOT.RooFit.MarkerColor(ROOT.kRed))
MC.plotOn(xframe_data_comp_mu, ROOT.RooFit.Cut("Categorization==0"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
xframe_data_comp_mu.SetMaximum(80)
# data.plotOn(xframe_data_comp_mu, ROOT.RooFit.Cut("isMuon==1"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson),ROOT.RooFit.MarkerColor(ROOT.kRed))
# MC.plotOn(xframe_data_comp_mu, ROOT.RooFit.Cut("isMuon==1"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))

xframe_data_comp_ele = Wmass.frame(55.,95.,10)
#xframe_data_comp_ele = Wmass.frame(50.,100.,15)
xframe_data_comp_ele.SetTitle(" ")

xframe_data_comp_ele.SetTitleOffset(1.5,"y")
data.plotOn(xframe_data_comp_ele, ROOT.RooFit.Cut("Categorization==2"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson),ROOT.RooFit.MarkerColor(ROOT.kRed))
MC.plotOn(xframe_data_comp_ele, ROOT.RooFit.Cut("Categorization==2"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
xframe_data_comp_ele.SetMaximum(80)
# data.plotOn(xframe_data_comp_ele, ROOT.RooFit.Cut("isMuon==0"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson),ROOT.RooFit.MarkerColor(ROOT.kRed))
# MC.plotOn(xframe_data_comp_ele, ROOT.RooFit.Cut("isMuon==0"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))

canvas_data_comp = ROOT.TCanvas()
canvas_data_comp.Divide(2,1)
canvas_data_comp.cd(1)
xframe_data_comp_mu.Draw()
canvas_data_comp.cd(2)
xframe_data_comp_ele.Draw()

canvas_data_comp.SaveAs("plots/Data_MC_comp.pdf")

raw_input()
