import ROOT
import math

fInput = ROOT.TFile("fitAllLep.root")
fInput.cd()
workspace = fInput.Get("workspace")

Wmass = workspace.var("Wmass")
W_xsec_constr = workspace.var("W_xsec_constr")
Nbkg_mu = workspace.var("Nbkg_mu")
Nbkg_el = workspace.var("Nbkg_el")
isMuon = workspace.cat("isMuon")
W_pigamma_BR = workspace.var("W_pigamma_BR")
totPDF = workspace.pdf("totPDF")

arglist = ROOT.RooArgList(Wmass, isMuon)

mcstudy = ROOT.RooMCStudy(totPDF, ROOT.RooArgSet(arglist), ROOT.RooFit.Silence(), ROOT.RooFit.Extended(1), ROOT.RooFit.FitOptions(ROOT.RooFit.Extended(1),ROOT.RooFit.Save(1), ROOT.RooFit.PrintEvalErrors(0)))

mcstudy.generateAndFit(100)

sigma1val_frame = mcstudy.plotParam(W_pigamma_BR, ROOT.RooFit.Bins(20))
sigma1err_frame = mcstudy.plotError(W_pigamma_BR, ROOT.RooFit.Bins(20))
sigma1pull_frame = mcstudy.plotPull(W_pigamma_BR, ROOT.RooFit.Bins(20), ROOT.RooFit.FitGauss(1))

NLLframe = mcstudy.plotNLL(ROOT.RooFit.Bins(20))

"""
canvas1 = ROOT.TCanvas()
canvas1.Divide(2,2)
canvas1.cd(1)
sigma1val_frame_mu.Draw()
canvas1.cd(2)
sigma1err_frame_mu.Draw()
canvas1.cd(3)
sigma1pull_frame_mu.Draw()
canvas1.cd(4)
NLLframe_mu.Draw()

canvas2 = ROOT.TCanvas()
canvas2.Divide(2,2)
canvas2.cd(1)
sigma1val_frame_mu.Draw()
canvas2.cd(2)
sigma1err_frame_mu.Draw()
canvas2.cd(3)
sigma1pull_frame_mu.Draw()
canvas2.cd(4)
NLLframe_mu.Draw()
"""

canvas1 = ROOT.TCanvas("canvas1","canvas1",900,700)
canvas1.Divide(2,2)
canvas1.cd(1)
sigma1val_frame.Draw()
canvas1.cd(2)
sigma1err_frame.Draw()
canvas1.cd(3)
sigma1pull_frame.Draw()
canvas1.cd(4)
NLLframe.Draw()
canvas1.SaveAs("Toy.png")
