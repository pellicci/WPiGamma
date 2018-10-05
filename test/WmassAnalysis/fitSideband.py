
#This program only fits a single lepton sample at a time

import ROOT

#Define the observable
Wmass = ROOT.RooRealVar("Wmass","#pi-#gamma invariant mass",50.,100.,"GeV")

fInput = ROOT.TFile("Tree_input_massfit_Data.root")
fInput.cd()

mytree = fInput.Get("minitree")

#Define the mu/ele category
Categorization = ROOT.RooCategory("Categorization","Categorization")
Categorization.defineType("MuonCR",0)
Categorization.defineType("MuonSignal",1)
Categorization.defineType("ElectronCR",2)
Categorization.defineType("ElectronSignal",3)

#Support variables
BDT_out = ROOT.RooRealVar("BDT_out","Output of BDT",-1.,1.)

#Create the RooDataSet. No need to import weight for signal only analysis

data_initial = ROOT.RooDataSet("data_initial","data_initial", ROOT.RooArgSet(Wmass,Categorization,BDT_out), ROOT.RooFit.Import(mytree))

data = data_initial.reduce("BDT_out > -0.1")

print "Using ", data.numEntries(), " events to fit"

#Now describe the background

#First the muon
a0_mu = ROOT.RooRealVar("a0_mu","a0_mu",0.15,-5.,5.)
a1_mu = ROOT.RooRealVar("a1_mu","a1_mu",-0.7,-5.,5.)
a2_mu = ROOT.RooRealVar("a2_mu","a2_mu",-0.2,-5.,5.)
a3_mu = ROOT.RooRealVar("a3_mu","a3_mu",-0.15,-5.,5.)
a4_mu = ROOT.RooRealVar("a4_mu","a4_mu",-0.15,-5.,5.)
a5_mu = ROOT.RooRealVar("a5_mu","a5_mu",0.1,-5.,5.)
a6_mu = ROOT.RooRealVar("a6_mu","a6_mu",0.1,-5.,5.)

a0_alt_mu = ROOT.RooRealVar("a0_alt_mu","a0_alt_mu",0.15,-5.,5.)
a1_alt_mu = ROOT.RooRealVar("a1_alt_mu","a1_alt_mu",-0.7,-5.,5.)
a2_alt_mu = ROOT.RooRealVar("a2_alt_mu","a2_alt_mu",-0.2,-5.,5.)
a3_alt_mu = ROOT.RooRealVar("a3_alt_mu","a3_alt_mu",-0.15,-5.,5.)
a4_alt_mu = ROOT.RooRealVar("a4_alt_mu","a4_alt_mu",-0.15,-1.,1.)
a5_alt_mu = ROOT.RooRealVar("a5_alt_mu","a5_alt_mu",0.1,-5.,5.)
a6_alt_mu = ROOT.RooRealVar("a6_alt_mu","a6_alt_mu",0.1,-5.,5.)

#a0_mu = ROOT.RooRealVar("a0_mu","a0_mu",0.5,0.,1.)
#a1_mu = ROOT.RooRealVar("a1_mu","a1_mu",0.7,0.,1.)
#a2_mu = ROOT.RooRealVar("a2_mu","a2_mu",0.7,0.,1.)
backPDF_mu     = ROOT.RooChebychev("backPDF_mu","backPDF_mu",Wmass,ROOT.RooArgList(a0_mu,a1_mu,a2_mu,a3_mu)) #,a4_mu)) #,a5_mu,a6_mu))
backPDF_alt_mu = ROOT.RooChebychev("backPDF_mu","backPDF_mu",Wmass,ROOT.RooArgList(a0_alt_mu,a1_alt_mu,a2_alt_mu,a3_alt_mu,a4_alt_mu,a5_mu))#,a6_mu))
#backPDF_mu = ROOT.RooBernstein("backPDF_mu","backPDF_mu",Wmass,ROOT.RooArgList(a0_mu,a1_mu,a2_mu)) #,a3_mu)) #,a4_mu)) #,a5_mu,a6_mu))

#Then the electron
a0_el = ROOT.RooRealVar("a0_el","a0_el",0.1,-5.,5.)
a1_el = ROOT.RooRealVar("a1_el","a1_el",-0.6,-5.,5.)
a2_el = ROOT.RooRealVar("a2_el","a2_el",-0.05,-5.,5.)
a3_el = ROOT.RooRealVar("a3_el","a3_el",-0.2,-5.,5.)
a4_el = ROOT.RooRealVar("a4_el","a4_el",0.1,-5.,5.)
a5_el = ROOT.RooRealVar("a5_el","a5_el",0.1,-5.,5.)
a6_el = ROOT.RooRealVar("a6_el","a6_el",0.1,-5.,5.)

a0_alt_el = ROOT.RooRealVar("a0_alt_el","a0_alt_el",0.1,-5.,5.)
a1_alt_el = ROOT.RooRealVar("a1_alt_el","a1_alt_el",-0.6,-5.,5.)
a2_alt_el = ROOT.RooRealVar("a2_alt_el","a2_alt_el",-0.05,-5.,5.)
a3_alt_el = ROOT.RooRealVar("a3_alt_el","a3_alt_el",-0.2,-5.,5.)
a4_alt_el = ROOT.RooRealVar("a4_alt_el","a4_alt_el",-0.15,-1.,1.)
a5_alt_el = ROOT.RooRealVar("a5_alt_el","a5_alt_el",0.1,-5.,5.)
a6_alt_el = ROOT.RooRealVar("a6_alt_el","a6_alt_el",0.1,-5.,5.)

#a0_el = ROOT.RooRealVar("a0_el","a0_el",0.5,0.,1.)
#a1_el = ROOT.RooRealVar("a1_el","a1_el",0.9,0.,1.)
#a2_el = ROOT.RooRealVar("a2_el","a2_el",0.7,0.,1.)
backPDF_el     = ROOT.RooChebychev("backPDF_el","backPDF_el",Wmass,ROOT.RooArgList(a0_el,a1_el,a2_el,a3_el))#a4_el)) #,a5_el,a6_el))
backPDF_alt_el = ROOT.RooChebychev("backPDF_el","backPDF_el",Wmass,ROOT.RooArgList(a0_alt_el,a1_alt_el,a2_alt_el,a3_alt_el,a4_alt_el,a5_el))#,a6_el))
#backPDF_el = ROOT.RooBernstein("backPDF_el","backPDF_el",Wmass,ROOT.RooArgList(a0_el,a1_el,a2_el)) #,a3_el)) #,a4_el)) #,a5_el,a6_el))

###########################################################

#Create the global simultaneous PDF
totPDF = ROOT.RooSimultaneous("totPDF","The total PDF",Categorization)
totPDF.addPdf(backPDF_mu,"MuonCR")
totPDF.addPdf(backPDF_el,"ElectronCR")

totPDF_alt = ROOT.RooSimultaneous("totPDF_alt","The total PDF - alternative",Categorization)
totPDF_alt.addPdf(backPDF_alt_mu,"MuonCR")
totPDF_alt.addPdf(backPDF_alt_el,"ElectronCR")

# constrained_params = ROOT.RooArgSet()
# constrained_params.add(dCB_width)
# constrained_params.add(eta)
# constrained_params.add(W_xsec_constr)
# constrained_params.add(lumi_constr)
# constrained_params.add(eff_mu_constr)
# constrained_params.add(eff_el_constr)

result     = totPDF.fitTo(data, ROOT.RooFit.NumCPU(4),ROOT.RooFit.Save())#,ROOT.RooFit.Constrain(constrained_params) )
result_alt = totPDF_alt.fitTo(data, ROOT.RooFit.NumCPU(4),ROOT.RooFit.Save())#,ROOT.RooFit.Constrain(constrained_params))

#Get the minimum of the likelihood from RooFitResult
# minNll_j = result_j.minNll()
# minNll_a = result_a.minNll()
# minNll_b = result_b.minNll()
# minNll_c = result_c.minNll()
# minNll_x = result_x.minNll()

#Get the number of free parameters
Npars = result.floatParsFinal().getSize()
Npars_alt = result_alt.floatParsFinal().getSize()

Nbins = 30

#Now plot the CR
xframe_mu_CR = Wmass.frame(ROOT.RooFit.Bins(Nbins))
xframe_mu_CR.SetTitle(" ")
#xframe_mu_CR.SetMaximum(30)
xframe_mu_CR.SetTitleOffset(1.4,"y")
data.plotOn(xframe_mu_CR, ROOT.RooFit.Cut("Categorization==0"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson),ROOT.RooFit.Name("data_mu"))
totPDF.plotOn(xframe_mu_CR, ROOT.RooFit.Slice(Categorization,"MuonCR"), ROOT.RooFit.ProjWData(data),ROOT.RooFit.Name("PDF_mu"))
totPDF_alt.plotOn(xframe_mu_CR, ROOT.RooFit.Slice(Categorization,"MuonCR"), ROOT.RooFit.ProjWData(data),ROOT.RooFit.LineColor(2),ROOT.RooFit.Name("PDF_mu_alt"))

xframe_el_CR = Wmass.frame(ROOT.RooFit.Bins(Nbins))
xframe_el_CR.SetTitle(" ")
#xframe_el_CR.SetMaximum(30)
xframe_el_CR.SetTitleOffset(1.4,"y")
data.plotOn(xframe_el_CR, ROOT.RooFit.Cut("Categorization==2"), ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson),ROOT.RooFit.Name("data_el"))
totPDF.plotOn(xframe_el_CR, ROOT.RooFit.Slice(Categorization,"ElectronCR"), ROOT.RooFit.ProjWData(data),ROOT.RooFit.Name("PDF_el"))
totPDF_alt.plotOn(xframe_el_CR, ROOT.RooFit.Slice(Categorization,"ElectronCR"), ROOT.RooFit.ProjWData(data),ROOT.RooFit.LineColor(2),ROOT.RooFit.Name("PDF_el_alt"))

canvas_CR = ROOT.TCanvas()
canvas_CR.Divide(2,1)
canvas_CR.cd(1)
xframe_mu_CR.Draw()
canvas_CR.cd(2)
xframe_el_CR.Draw()


#Redisuals
h_resid_mu     = ROOT.TH1F()
h_resid_mu_alt = ROOT.TH1F()
h_resid_el     = ROOT.TH1F()
h_resid_el_alt = ROOT.TH1F()

h_resid_mu     = xframe_mu_CR.residHist("data_mu","PDF_mu")
h_resid_mu_alt = xframe_mu_CR.residHist("data_mu","PDF_mu_alt")
h_resid_el     = xframe_el_CR.residHist("data_el","PDF_el")
h_resid_el_alt = xframe_el_CR.residHist("data_el","PDF_el_alt")

file_resid = ROOT.TFile("residuals.root","RECREATE")
h_resid_mu.Write("resid_mu")
h_resid_mu_alt.Write("resid_mu_alt")
h_resid_el.Write("resid_el")
h_resid_el_alt.Write("resid_el_alt")


res_mu     = 0.
res_mu_alt = 0.
res_el     = 0.
res_el_alt = 0.

# for x in range (1, Nbins+1):
res_mu     = h_resid_mu.GetPoint(3,ROOT.Double(x),ROOT.Double(y))
res_mu_alt = h_resid_mu_alt.getFitRangeNEvt()
res_el     = h_resid_el.getFitRangeNEvt()
res_el_alt = h_resid_el_alt.getFitRangeNEvt()


print "Npars: ", Npars, "Npars_alt: ", Npars_alt
print "res_mu: ", res_mu, "res_mu_alt: ", res_mu_alt
print "res_el: ", res_el, "res_el_alt: ", res_el_alt

# #chi2 = ROOT.RooChi2Var("chi2","chi2",sidebandPDF,data_lep)
# chi_2_j = xframe.chiSquare("PDF_j","data_lep",Npars_j)
# chi_2_a = xframe.chiSquare("PDF_a","data_lep",Npars_a)
# chi_2_b = xframe.chiSquare("PDF_b","data_lep",Npars_b)
# chi_2_c = xframe.chiSquare("PDF_c","data_lep",Npars_c)
# chi_2_x = xframe.chiSquare("PDF_x","data_lep",Npars_x)

# Deltaloglikelihood_12 = -2*(minNll_a-minNll_j)
# Deltaloglikelihood_23 = -2*(minNll_b-minNll_a)
# Deltaloglikelihood_34 = -2*(minNll_c-minNll_b)

# prob_12 = ROOT.TMath.Prob(Deltaloglikelihood_12, 1)
# prob_23 = ROOT.TMath.Prob(Deltaloglikelihood_23, 1)
# prob_34 = ROOT.TMath.Prob(Deltaloglikelihood_34, 1)

# print "reduced_chi2_1 = ", chi_2_j, "   |||||||   reduced_chi2_2 = ", chi_2_a, "   |||||||   reduced_chi2_3 = ", chi_2_b , "   |||||||   reduced_chi2_4 = ", chi_2_c, "  |||||||   reduced_chi2_exp = ", chi_2_x
# print "minNll_1 = ", minNll_j, " |||||||   minNll_2 = ", minNll_a, " |||||||   minNll_3 = ", minNll_b,  " |||||||   minNll_4 = ", minNll_c, "  |||||||  minNll_exp = ", minNll_x
# print "Deltaloglikelihood_12 = ", -2*(minNll_a-minNll_j)
# print "Deltaloglikelihood_23 = ", -2*(minNll_b-minNll_a)
# print "Deltaloglikelihood_34 = ", -2*(minNll_c-minNll_b)
# print "prob_12 = ", prob_12
# print "prob_23 = ", prob_23
# print "prob_34 = ", prob_34

# canvas = ROOT.TCanvas()
# canvas.cd()
# legend = ROOT.TLegend(0.1,0.7,0.48,0.9)
# legend.SetNColumns(2)
# legend.AddEntry(0,"yellow: 1 degree","")
# legend.AddEntry(0,"blue: 2 degree","")
# legend.AddEntry(0,"red: 3 degree","")
# legend.AddEntry(0,"green: 4 degree","")
# legend.AddEntry(0,"purple: exp","")
# xframe.Draw()
# legend.Draw()
# canvas.SaveAs("plots/fitSidebands_exp_muons.png")

raw_input()
