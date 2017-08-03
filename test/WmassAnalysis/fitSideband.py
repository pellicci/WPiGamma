
#This program only fits a single lepton sample at a time

import ROOT

#Define the observable
Wmass = ROOT.RooRealVar("Wmass","#pi-#gamma invariant mass",50.,100.,"GeV")
Wmass.setRange("LowerRange",50.,65.)
Wmass.setRange("UpperRange",90.,100.)

#Retrive the sample
fInput = ROOT.TFile("Tree_Data.root")
fInput.cd()

mytree = fInput.Get("minitree"))

#Define the mu/ele category
isMuon = ROOT.RooCategory("isMuon","isMuon")
isMuon.defineType("Muon",1)
isMuon.defineType("Electron",0)

#Create the RooDataSet. No need to import weight for signal only analysis
data = ROOT.RooDataSet("data","data", ROOT.RooArgSet(Wmass,isMuon), ROOT.RooFit.Import(mytree))

#Skim one lepton sample only
data_lep = data.reduce("isMuon==1")
#data_lep = data.reduce("isMuon==isMuon::Electron")

print "Using ", data_lep.numEntries(), " events to fit the lepton shape"

#Describe the sidebands
a0 = ROOT.RooRealVar("a0","a0",0.1,-5.,5.)
a1 = ROOT.RooRealVar("a1","a1",0.1,-5.,5.)
b0 = ROOT.RooRealVar("b0","b0",0.1,-5.,5.)
b1 = ROOT.RooRealVar("b1","b1",0.1,-5.,5.)
b2 = ROOT.RooRealVar("b2","b2",0.1,-5.,5.)

sidebandPDF_a = ROOT.RooChebychev("backPDF_a","backPDF_a",Wmass,ROOT.RooArgList(a0,a1))
sidebandPDF_b = ROOT.RooChebychev("backPDF_b","backPDF_b",Wmass,ROOT.RooArgList(b0,b1,b2))

result_a = sidebandPDF_a.fitTo(data_lep, ROOT.RooFit.Range("LowerRange,UpperRange"), ROOT.RooFit.Extended(0), ROOT.RooFit.SumW2Error(1),ROOT.RooFit.Save(1))
result_b = sidebandPDF_b.fitTo(data_lep, ROOT.RooFit.Range("LowerRange,UpperRange"), ROOT.RooFit.Extended(0), ROOT.RooFit.SumW2Error(1),ROOT.RooFit.Save(1))

#Get the minimum of the likelihood from RooFitResult
minNll_a = result_a.minNll()
minNll_b = result_b.minNll()

#Get the number of free parameters
Npars_a = result_a.floatParsFinal().getSize()
Npars_b = result_b.floatParsFinal().getSize()

xframe = Wmass.frame(ROOT.RooFit.Bins(50))
data_lep.plotOn(xframe,ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson),ROOT.RooFit.Name("data"))
sidebandPDF_a.plotOn(xframe,ROOT.RooFit.ProjWData(data),ROOT.RooFit.Name("PDF_a"))
sidebandPDF_b.plotOn(xframe,ROOT.RooFit.ProjWData(data),ROOT.RooFit.Name("PDF_b"),ROOT.RooFit.LineColor(2))

#chi2 = ROOT.RooChi2Var("chi2","chi2",sidebandPDF,data_lep)
chi_2_a = xframe.chiSquare("PDF_a","data",Npars_a)
chi_2_b = xframe.chiSquare("PDF_b","data",Npars_b)

print "chi2_a = ", chi_2_a, "   |||||||   chi2_b", chi_2_b 
print "minNll_a = ", minNll_a, " |||||||   minNll_b = ", minNll_b
print "Deltaloglikelihood = ", -2*ROOT.TMath.Log(minNll_b/minNll_a)


canvas = ROOT.TCanvas()
canvas.cd()
xframe.Draw()
canvas.SaveAs("fitSidebands.png")
