from ROOT import *

def generate_rootfile():

    fIn_ttbar  = TFile("ttbar_pT.root")
    fIn_signal = TFile("signal_pT.root")
    
    h_ttbar  = fIn_ttbar.Get("W_ttbar")
    h_signal = fIn_signal.Get("W_signal")
    
    h_ttbar.Scale(1./h_ttbar.Integral())
    h_signal.Scale(1./h_signal.Integral())
    
    h_ttbar.Divide(h_signal)
    
    canvas = TCanvas()
    gStyle.SetOptStat(0)
    h_ttbar.SetTitle("")
    h_ttbar.SetMarkerStyle(21)
    h_ttbar.SetLineColor(46)
    h_ttbar.GetXaxis().SetTitle("#it{p}_{T}^{W} (GeV)")
    h_ttbar.Draw()
    canvas.SaveAs("ttbar_signal_ratio_2018.pdf")
    
    fOut = TFile.Open("ttbar_signal_ratio_2018.root","RECREATE")
    fOut.cd()
    h_ttbar.Write("ttbar_signal_ratio")
    fOut.Close()

def generate_total_plot():

    fIn_2016 = TFile.Open("ttbar_signal_ratio_2016.root")
    h_ttbar_2016 = fIn_2016.Get("ttbar_signal_ratio")

    fIn_2017 = TFile.Open("ttbar_signal_ratio_2017.root")
    h_ttbar_2017 = fIn_2017.Get("ttbar_signal_ratio")

    fIn_2018 = TFile.Open("ttbar_signal_ratio_2018.root")
    h_ttbar_2018 = fIn_2018.Get("ttbar_signal_ratio")

    canvas_total = TCanvas()
    gStyle.SetOptStat(0)
    h_ttbar_2016.SetTitle("")
    h_ttbar_2016.SetMarkerStyle(21)
    h_ttbar_2016.SetMarkerColor(3)
    h_ttbar_2016.SetLineColor(1)
    h_ttbar_2016.GetXaxis().SetTitle("#it{p}_{T}^{W} (GeV)")
    h_ttbar_2016.GetYaxis().SetTitle("Arbitrary units")

    h_ttbar_2017.SetTitle("")
    h_ttbar_2017.SetMarkerStyle(21)
    h_ttbar_2017.SetMarkerColor(2)
    h_ttbar_2017.SetLineColor(1)
    h_ttbar_2017.GetXaxis().SetTitle("#it{p}_{T}^{W} (GeV)")
    h_ttbar_2017.GetYaxis().SetTitle("Arbitrary units")

    h_ttbar_2018.SetTitle("")
    h_ttbar_2018.SetMarkerStyle(21)
    h_ttbar_2018.SetMarkerColor(1)
    h_ttbar_2018.SetLineColor(1)
    h_ttbar_2018.GetXaxis().SetTitle("#it{p}_{T}^{W} (GeV)")
    h_ttbar_2018.GetYaxis().SetTitle("Arbitrary units")

    legend = TLegend(0.2,0.55,0.35,0.8)
    legend.SetHeader(" ")
    legend.SetFillColor(0)
    legend.SetBorderSize(0)
    legend.SetLineColor(1)
    legend.SetLineStyle(1)
    legend.SetLineWidth(1)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.05)
    legend.AddEntry(h_ttbar_2016,"2016","lep")
    legend.AddEntry(h_ttbar_2017,"2017","lep")
    legend.AddEntry(h_ttbar_2018,"2018","lep")

    h_ttbar_2018.Draw()
    h_ttbar_2017.Draw("SAME")
    h_ttbar_2016.Draw("SAME")
    legend.Draw("SAME")

    canvas_total.SaveAs("ttbar_signal_ratio_total.pdf")

    raw_input()

if __name__ == "__main__":

    #generate_rootfile()
    generate_total_plot()
