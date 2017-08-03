
#This program only fits a single lepton sample at a time

import ROOT

#Define the observable
Wmass = ROOT.RooRealVar("Wmass","#pi-#gamma invariant mass",50.,100.,"GeV")
Wmass.setRange("LowerRange",50.,65.)
Wmass.setRange("UpperRange",90.,100.)

#Retrive the sample
fInput = ROOT.TFile("Tree_Data.root")
fInput.cd()

mytree = fInput.Get("minitree")

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
j0 = ROOT.RooRealVar("j0","j0",0.1,-5.,5.) # 1 degree polynomial
a0 = ROOT.RooRealVar("a0","a0",0.1,-5.,5.) # 2 degree polynomial
a1 = ROOT.RooRealVar("a1","a1",0.1,-5.,5.)
b0 = ROOT.RooRealVar("b0","b0",0.1,-5.,5.) # 3 degree polynomial
b1 = ROOT.RooRealVar("b1","b1",0.1,-5.,5.)
b2 = ROOT.RooRealVar("b2","b2",0.1,-5.,5.)
c0 = ROOT.RooRealVar("c0","c0",0.1,-5.,5.) # 4 degree polynomial
c1 = ROOT.RooRealVar("c1","c1",0.1,-5.,5.)
c2 = ROOT.RooRealVar("c2","c2",0.1,-5.,5.)
c3 = ROOT.RooRealVar("c3","c3",0.1,-5.,5.)

sidebandPDF_j = ROOT.RooChebychev("backPDF_j","backPDF_j",Wmass,ROOT.RooArgList(j0))
sidebandPDF_a = ROOT.RooChebychev("backPDF_a","backPDF_a",Wmass,ROOT.RooArgList(a0,a1))
sidebandPDF_b = ROOT.RooChebychev("backPDF_b","backPDF_b",Wmass,ROOT.RooArgList(b0,b1,b2))
sidebandPDF_c = ROOT.RooChebychev("backPDF_c","backPDF_c",Wmass,ROOT.RooArgList(c0,c1,c2,c3))

result_j = sidebandPDF_j.fitTo(data_lep, ROOT.RooFit.Range("LowerRange,UpperRange"),ROOT.RooFit.Save(1))
result_a = sidebandPDF_a.fitTo(data_lep, ROOT.RooFit.Range("LowerRange,UpperRange"),ROOT.RooFit.Save(1))
result_b = sidebandPDF_b.fitTo(data_lep, ROOT.RooFit.Range("LowerRange,UpperRange"),ROOT.RooFit.Save(1))
result_c = sidebandPDF_c.fitTo(data_lep, ROOT.RooFit.Range("LowerRange,UpperRange"),ROOT.RooFit.Save(1))

#Get the minimum of the likelihood from RooFitResult
minNll_j = result_j.minNll()
minNll_a = result_a.minNll()
minNll_b = result_b.minNll()
minNll_c = result_c.minNll()

#Get the number of free parameters
Npars_j = result_j.floatParsFinal().getSize()
Npars_a = result_a.floatParsFinal().getSize()
Npars_b = result_b.floatParsFinal().getSize()
Npars_c = result_c.floatParsFinal().getSize()

xframe = Wmass.frame(ROOT.RooFit.Bins(50))
data_lep.plotOn(xframe,ROOT.RooFit.Name("data"))
sidebandPDF_j.plotOn(xframe,ROOT.RooFit.Name("PDF_j"), ROOT.RooFit.Range("Full"),ROOT.RooFit.LineColor(5))
sidebandPDF_a.plotOn(xframe,ROOT.RooFit.Name("PDF_a"), ROOT.RooFit.Range("Full"))
sidebandPDF_b.plotOn(xframe,ROOT.RooFit.Name("PDF_b"), ROOT.RooFit.Range("Full"),ROOT.RooFit.LineColor(2))
sidebandPDF_c.plotOn(xframe,ROOT.RooFit.Name("PDF_c"), ROOT.RooFit.Range("Full"),ROOT.RooFit.LineColor(3))

#chi2 = ROOT.RooChi2Var("chi2","chi2",sidebandPDF,data_lep)
chi_2_j = xframe.chiSquare("PDF_j","data",Npars_j)
chi_2_a = xframe.chiSquare("PDF_a","data",Npars_a)
chi_2_b = xframe.chiSquare("PDF_b","data",Npars_b)
chi_2_c = xframe.chiSquare("PDF_c","data",Npars_c)

print "reduced_chi2_j = ", chi_2_j, "   |||||||   reduced_chi2_a = ", chi_2_a, "   |||||||   reduced_chi2_b = ", chi_2_b , "   |||||||   reduced_chi2_c = ", chi_2_c 
print "minNll_j = ", minNll_j, " |||||||   minNll_a = ", minNll_a, " |||||||   minNll_b = ", minNll_b,  " |||||||   minNll_c = ", minNll_c
print "Deltaloglikelihood_ja = ", -2*ROOT.TMath.Log(minNll_a/minNll_j)
print "Deltaloglikelihood_ab = ", -2*ROOT.TMath.Log(minNll_b/minNll_a)
print "Deltaloglikelihood_bc = ", -2*ROOT.TMath.Log(minNll_c/minNll_b)


canvas = ROOT.TCanvas()
canvas.cd()
xframe.Draw()
canvas.SaveAs("fitSidebands_full_4.png")
