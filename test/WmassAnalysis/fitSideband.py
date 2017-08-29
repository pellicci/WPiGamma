
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
x0 = ROOT.RooRealVar("x0","x0",0.1,0.01,0.5) # exponential

sidebandPDF_j = ROOT.RooBernstein("backPDF_j","backPDF_j",Wmass,ROOT.RooArgList(j0))
sidebandPDF_a = ROOT.RooBernstein("backPDF_a","backPDF_a",Wmass,ROOT.RooArgList(a0,a1))
sidebandPDF_b = ROOT.RooBernstein("backPDF_b","backPDF_b",Wmass,ROOT.RooArgList(b0,b1,b2))
sidebandPDF_c = ROOT.RooBernstein("backPDF_c","backPDF_c",Wmass,ROOT.RooArgList(c0,c1,c2,c3))
sidebandPDF_x = ROOT.RooExponential("backPDF_x","backPDF_x",Wmass,x0)

#Fit the PDFs on data
result_j = sidebandPDF_j.fitTo(data_lep,ROOT.RooFit.Range("LowerRange,UpperRange"),ROOT.RooFit.Save())
result_a = sidebandPDF_a.fitTo(data_lep,ROOT.RooFit.Range("LowerRange,UpperRange"),ROOT.RooFit.Save())
result_b = sidebandPDF_b.fitTo(data_lep,ROOT.RooFit.Range("LowerRange,UpperRange"),ROOT.RooFit.Save())
result_c = sidebandPDF_c.fitTo(data_lep,ROOT.RooFit.Range("LowerRange,UpperRange"),ROOT.RooFit.Save())
result_x = sidebandPDF_x.fitTo(data_lep,ROOT.RooFit.Range("LowerRange,UpperRange"),ROOT.RooFit.Save())

#Get the minimum of the likelihood from RooFitResult
minNll_j = result_j.minNll()
minNll_a = result_a.minNll()
minNll_b = result_b.minNll()
minNll_c = result_c.minNll()
minNll_x = result_x.minNll()

#Get the number of free parameters
Npars_j = result_j.floatParsFinal().getSize()
Npars_a = result_a.floatParsFinal().getSize()
Npars_b = result_b.floatParsFinal().getSize()
Npars_c = result_c.floatParsFinal().getSize()
Npars_x = result_x.floatParsFinal().getSize()

xframe = Wmass.frame(ROOT.RooFit.Bins(30))
data_lep.plotOn(xframe,ROOT.RooFit.Name("data_lep"))
sidebandPDF_j.plotOn(xframe,ROOT.RooFit.Name("PDF_j"),ROOT.RooFit.Range("LowerRange,UpperRange"),ROOT.RooFit.LineColor(5)) # yellow. 1 degree
sidebandPDF_a.plotOn(xframe,ROOT.RooFit.Name("PDF_a"),ROOT.RooFit.Range("LowerRange,UpperRange"),ROOT.RooFit.LineColor(4)) # blue. 2 degree
sidebandPDF_b.plotOn(xframe,ROOT.RooFit.Name("PDF_b"),ROOT.RooFit.Range("LowerRange,UpperRange"),ROOT.RooFit.LineColor(2)) # red. 3 degree
sidebandPDF_c.plotOn(xframe,ROOT.RooFit.Name("PDF_c"),ROOT.RooFit.Range("LowerRange,UpperRange"),ROOT.RooFit.LineColor(3)) # green. 4 degree
sidebandPDF_x.plotOn(xframe,ROOT.RooFit.Name("PDF_x"),ROOT.RooFit.Range("LowerRange,UpperRange"),ROOT.RooFit.LineColor(6)) # purple. exponential

#Redisuals
hresid_j = ROOT.TH1F()
hresid_j = xframe.residHist("data_lep","PDF_j")
hresid_a = xframe.residHist("data_lep","PDF_a")
hresid_b = xframe.residHist("data_lep","PDF_b")
hresid_c = xframe.residHist("data_lep","PDF_c")
"""
res_j = 0.
res_a = 0.
res_b = 0.
res_c = 0.

for x in range (1, 31):
    res_j += hresid_j.getBinContent(x)
    res_a += hresid_a.ROOT.getBinContent(x)
    res_b += hresid_b.ROOT.getBinContent(x)
    res_c += hresid_c.ROOT.getBinContent(x)
"""
#chi2 = ROOT.RooChi2Var("chi2","chi2",sidebandPDF,data_lep)
chi_2_j = xframe.chiSquare("PDF_j","data_lep",Npars_j)
chi_2_a = xframe.chiSquare("PDF_a","data_lep",Npars_a)
chi_2_b = xframe.chiSquare("PDF_b","data_lep",Npars_b)
chi_2_c = xframe.chiSquare("PDF_c","data_lep",Npars_c)
chi_2_x = xframe.chiSquare("PDF_x","data_lep",Npars_x)

Deltaloglikelihood_12 = -2*(minNll_a-minNll_j)
Deltaloglikelihood_23 = -2*(minNll_b-minNll_a)
Deltaloglikelihood_34 = -2*(minNll_c-minNll_b)

prob_12 = ROOT.TMath.Prob(Deltaloglikelihood_12, 1)
prob_23 = ROOT.TMath.Prob(Deltaloglikelihood_23, 1)
prob_34 = ROOT.TMath.Prob(Deltaloglikelihood_34, 1)

print "reduced_chi2_1 = ", chi_2_j, "   |||||||   reduced_chi2_2 = ", chi_2_a, "   |||||||   reduced_chi2_3 = ", chi_2_b , "   |||||||   reduced_chi2_4 = ", chi_2_c, "  |||||||   reduced_chi2_exp = ", chi_2_x
print "minNll_1 = ", minNll_j, " |||||||   minNll_2 = ", minNll_a, " |||||||   minNll_3 = ", minNll_b,  " |||||||   minNll_4 = ", minNll_c, "  |||||||  minNll_exp = ", minNll_x
print "Deltaloglikelihood_12 = ", -2*(minNll_a-minNll_j)
print "Deltaloglikelihood_23 = ", -2*(minNll_b-minNll_a)
print "Deltaloglikelihood_34 = ", -2*(minNll_c-minNll_b)
print "prob_12 = ", prob_12
print "prob_23 = ", prob_23
print "prob_34 = ", prob_34

canvas = ROOT.TCanvas()
canvas.cd()
legend = ROOT.TLegend(0.1,0.7,0.48,0.9)
legend.SetNColumns(2)
legend.AddEntry(0,"yellow: 1 degree","")
legend.AddEntry(0,"blue: 2 degree","")
legend.AddEntry(0,"red: 3 degree","")
legend.AddEntry(0,"green: 4 degree","")
legend.AddEntry(0,"purple: exp","")
xframe.Draw()
legend.Draw()
canvas.SaveAs("fitSidebands_exp_muons.png")
