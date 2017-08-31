import ROOT
import math

fInput = ROOT.TFile("fitAllLep.root")
fInput.cd()
workspace = fInput.Get("workspace")

Wmass = workspace.var("Wmass")
W_xsec_constr = workspace.var("W_xsec_constr")
Nbkg_mu = workspace.var("Nbkg_mu")
Nbkg_el = workspace.var("Nbkg_el")
lumi_constr = workspace.var("lumi_constr")
eff_mu_constr = workspace.var("eff_mu_constr")
eff_el_constr = workspace.var("eff_el_constr")
isMuon = workspace.cat("isMuon")
W_pigamma_BR = workspace.var("W_pigamma_BR")
totPDF = workspace.pdf("totPDF")

constrained_params = ROOT.RooArgSet()
constrained_params.add(W_xsec_constr)
constrained_params.add(lumi_constr)
constrained_params.add(eff_mu_constr)
constrained_params.add(eff_el_constr)

arglist = ROOT.RooArgList(Wmass, isMuon)

mcstudy = ROOT.RooMCStudy(totPDF, ROOT.RooArgSet(arglist), ROOT.RooFit.Silence(), ROOT.RooFit.Extended(1), ROOT.RooFit.FitOptions(ROOT.RooFit.Extended(1),  ROOT.RooFit.Constrain(constrained_params), ROOT.RooFit.Save(1), ROOT.RooFit.PrintEvalErrors(0)))

mcstudy.generateAndFit(3)

W_pigamma_BR_val_frame = mcstudy.plotParam(W_pigamma_BR, ROOT.RooFit.Bins(20))
W_pigamma_BR_err_frame = mcstudy.plotError(W_pigamma_BR, ROOT.RooFit.Bins(20))
W_pigamma_BR_pull_frame = mcstudy.plotPull(W_pigamma_BR, ROOT.RooFit.Bins(20), ROOT.RooFit.FitGauss(1))
Nbkg_mu_frame = mcstudy.plotPull(Nbkg_mu, ROOT.RooFit.Bins(20), ROOT.RooFit.FitGauss(1))
Nbkg_el_frame = mcstudy.plotPull(Nbkg_el, ROOT.RooFit.Bins(20), ROOT.RooFit.FitGauss(1))
W_xsec_constr_frame = mcstudy.plotPull(W_xsec_constr, ROOT.RooFit.Bins(20), ROOT.RooFit.FitGauss(1))

NLLframe = mcstudy.plotNLL(ROOT.RooFit.Bins(20))

#Some settings
W_pigamma_BR_val_frame.SetTitle("")
W_pigamma_BR_val_frame.SetXTitle("BR(W#rightarrow#pi#gamma)")
W_pigamma_BR_err_frame.SetTitle("")
W_pigamma_BR_err_frame.SetXTitle("#sigma_{BR(W#rightarrow#pi#gamma)}")
W_pigamma_BR_pull_frame.SetTitle("")
W_pigamma_BR_pull_frame.SetXTitle("PULL_{BR(W#rightarrow#pi#gamma)}")
NLLframe.SetTitle("")

canvas1 = ROOT.TCanvas("canvas1","canvas1",900,700)
canvas1.Divide(2,2)
canvas1.cd(1)
W_pigamma_BR_val_frame.Draw()
canvas1.cd(2)
W_pigamma_BR_err_frame.Draw()
canvas1.cd(3)
W_pigamma_BR_pull_frame.Draw()
canvas1.cd(4)
NLLframe.Draw()
canvas1.SaveAs("Toy.png")

canvas2 = ROOT.TCanvas()
Nbkg_mu_frame.Draw()
canvas2.SaveAs("Nbkg_mu_pull.png")

canvas3 = ROOT.TCanvas()
Nbkg_el_frame.Draw()
canvas3.SaveAs("Nbkg_el_pull.png")

canvas4 = ROOT.TCanvas()
W_xsec_constr_frame.Draw()
canvas4.SaveAs("W_xsec_constr_pull.png")

