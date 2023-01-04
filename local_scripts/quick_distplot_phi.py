#import CombineHarvester.CombineTools.plotting as plot
import ROOT
import argparse

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

parser = argparse.ArgumentParser()
parser.add_argument('--infile',      dest='infile',      default="infile.root",       help="infile")
parser.add_argument('--tree',        dest='tree',        default="shapes_prefit",     help="tree")
parser.add_argument('--regions',     dest='regions',     default=["ptbin0pass"],      help="regions",   nargs='+')
parser.add_argument('--label',       dest='label',       default="prefit_pt0pass",    help="label")
args = parser.parse_args()

#python quick_distplot.py --infile fitDiagnostics.hadmuModel.root --tree shapes_prefit --regions ptbin0fail ptbin0loosepass ptbin0pass ptbin1fail ptbin1loosepass ptbin1pass topCRfail topCRloosepass topCRpass wlnuCRfail wlnuCRloosepass wlnuCRpass --label hadmu_prefit
#python quick_distplot.py --infile fitDiagnostics.hadmuModel.root --tree shapes_fit_b --regions ptbin0fail ptbin0loosepass ptbin0pass ptbin1fail ptbin1loosepass ptbin1pass topCRfail topCRloosepass topCRpass wlnuCRfail wlnuCRloosepass wlnuCRpass --label hadmu_fit_b
#python quick_distplot.py --infile fitDiagnostics.hadmuModel.root --tree shapes_fit_s --regions ptbin0fail ptbin0loosepass ptbin0pass ptbin1fail ptbin1loosepass ptbin1pass topCRfail topCRloosepass topCRpass wlnuCRfail wlnuCRloosepass wlnuCRpass --label hadmu_fit_s

#plot.ModTDRStyle()

canvas = ROOT.TCanvas()
pad1 = ROOT.TPad( 'pad1', 'plot',  0.05, 0.35, 0.95, 0.95 )
pad2 = ROOT.TPad( 'pad2', 'ratio', 0.05, 0.05, 0.95, 0.35 )
pad1.Draw()
pad2.Draw()

fin = ROOT.TFile(args.infile)

first_dir = args.tree
for second_dir in args.regions:
    pad1.cd()
    
    h_err = fin.Get(first_dir + '/' + second_dir + '/total_background')
    h_top = fin.Get(first_dir + '/' + second_dir + '/top')
    h_wlnu = fin.Get(first_dir + '/' + second_dir + '/wlnu')
    h_ztt = fin.Get(first_dir + '/' + second_dir + '/dy')
    h_vvqq = fin.Get(first_dir + '/' + second_dir + '/vvqq')
    h_qcd = fin.Get(first_dir + '/' + second_dir + '/multijet')
    h_h125 = fin.Get(first_dir + '/' + second_dir + '/htt125')
    h_sig = fin.Get(first_dir + '/' + second_dir + '/total_signal')
    h_dat = fin.Get(first_dir + '/' + second_dir + '/data')  # This is a TGraphAsymmErrors, not a TH1F
    
    h_all = ROOT.THStack("h_all","h_all");
    atleast1 = False
    if str(type(h_qcd))=='<class \'ROOT.TH1F\'>':
        h_qcd.SetFillColor(ROOT.kBlue+1);
        h_qcd.SetLineColor(ROOT.kBlack);
        h_all.Add(h_qcd);
        atleast1 = True
    if str(type(h_h125))=='<class \'ROOT.TH1F\'>':
        h_h125.SetFillColor(ROOT.kCyan);
        h_h125.SetLineColor(ROOT.kBlack);
        h_all.Add(h_h125);
        atleast1 = True
    if str(type(h_top))=='<class \'ROOT.TH1F\'>':
        h_top.SetFillColor(ROOT.kRed+1);
        h_top.SetLineColor(ROOT.kBlack);
        h_all.Add(h_top);
        atleast1 = True
    if str(type(h_wlnu))=='<class \'ROOT.TH1F\'>':
        h_wlnu.SetFillColor(ROOT.kGreen+1);
        h_wlnu.SetLineColor(ROOT.kBlack);
        h_all.Add(h_wlnu);
        atleast1 = True
    if str(type(h_ztt))=='<class \'ROOT.TH1F\'>':
        h_ztt.SetFillColor(ROOT.kMagenta+1);
        h_ztt.SetLineColor(ROOT.kBlack);
        h_all.Add(h_ztt);
        atleast1 = True
    if str(type(h_vvqq))=='<class \'ROOT.TH1F\'>':
        h_vvqq.SetFillColor(ROOT.kOrange+1);
        h_vvqq.SetLineColor(ROOT.kBlack);
        h_all.Add(h_vvqq);
        atleast1 = True
    if str(type(h_sig))=='<class \'ROOT.TH1F\'>':
        h_sig.SetLineColor(ROOT.kGray)
        h_sig.SetFillColor(ROOT.kWhite);
        h_sig.SetLineWidth(3)
        h_sig.SetLineStyle(1)
        h_all.Add(h_sig)
        atleast1 = True

    if not atleast1:
        continue
    
    title_str = "";
    if (first_dir=="shapes_prefit"): 
        title_str = "Prefit, ";
    elif (first_dir=="shapes_fit_b"): 
        title_str = "B-only fit, ";
    elif (first_dir=="shapes_fit_s"): 
        title_str = "S+B fit, ";
    h_all.SetTitle(title_str+second_dir);

    h_all.Draw('HIST')
    #h_all.GetXaxis().SetTitle("m_{NN}")
    h_all.GetYaxis().SetTitle("Events / GeV")
    
    h_err.SetFillColorAlpha(12, 0.3)  # Set grey colour (12) and alpha (0.3)
    h_err.SetMarkerSize(0)
    h_err.Draw('E2SAME')
    
    h_dat.SetMarkerStyle(1)
    h_dat.SetMarkerStyle(20)
    h_dat.Draw('PSAME')
    
    h_all.SetMaximum(h_all.GetMaximum() * 1.4)
    
    legend = ROOT.TLegend(0.60, 0.70, 0.90, 0.91, '', 'NBNDC')
    if str(type(h_qcd))=='<class \'ROOT.TH1F\'>':
        legend.AddEntry(h_qcd, 'Multijet', 'F')
    if str(type(h_h125))=='<class \'ROOT.TH1F\'>':
        legend.AddEntry(h_h125, 'H125', 'F')
    if str(type(h_top))=='<class \'ROOT.TH1F\'>':
        legend.AddEntry(h_top, 'Top', 'F')
    if str(type(h_wlnu))=='<class \'ROOT.TH1F\'>':
        legend.AddEntry(h_wlnu, 'W+jets', 'F')
    if str(type(h_ztt))=='<class \'ROOT.TH1F\'>':
        legend.AddEntry(h_ztt, 'Z+jets', 'F')
    if str(type(h_vvqq))=='<class \'ROOT.TH1F\'>':
        legend.AddEntry(h_vvqq, '(V)V(qq)', 'F')
    if str(type(h_sig))=='<class \'ROOT.TH1F\'>':
        legend.AddEntry(h_sig, 'Signal', 'L')
    legend.AddEntry(h_err, 'Background uncertainty', 'F')
    legend.Draw()

    h_div = h_err.Clone()
    h_mc = h_err.Clone()
    h_div.Reset("ICE")
    h_div.Sumw2()
    h_line = h_div.Clone()
    h_line2 = h_div.Clone()
    h_line3 = h_div.Clone()
    h_line4 = h_div.Clone()
    x = ROOT.Double(0.)
    y = ROOT.Double(0.)
    for i in range(h_dat.GetN()):
        h_dat.GetPoint(i, x, y)
        binw = h_div.GetBinWidth(h_div.FindBin(x))
        #print(y,'/',h_err.GetBinContent(h_div.FindBin(x)),'=',y/h_err.GetBinContent(h_div.FindBin(x)) if h_err.GetBinContent(h_div.FindBin(x))>0. else 0.)
        #for j in range(int(y*binw)):
        #    h_div.Fill(x,1./float(binw))
        diff = y-h_mc.GetBinContent(i+1)
        data_err = h_dat.GetErrorYhigh(i) if diff < 0. else h_dat.GetErrorYlow(i)
        err_val = ROOT.TMath.Sqrt(data_err*data_err + h_mc.GetBinError(i+1)*h_mc.GetBinError(i+1))
        err_val = err_val if err_val>0. else 1.
        h_div.SetBinContent(i+1,diff/err_val)
    for i in range(1,h_line.GetNbinsX()+1):
        h_line.Fill(h_line.GetBinCenter(i),1.)
        h_line2.Fill(h_line.GetBinCenter(i),-1.)
        h_line3.Fill(h_line.GetBinCenter(i),2.)
        h_line4.Fill(h_line.GetBinCenter(i),-2.)
    #h_div.Divide(h_err)

    #h_mc.Divide(h_err)

    pad2.cd()
    h_line.SetFillColor(ROOT.kWhite);
    h_line.SetLineColor(ROOT.kBlack);
    h_line.SetLineStyle(2);
    h_line.SetStats(0)
    h_line2.SetFillColor(ROOT.kWhite);
    h_line2.SetLineColor(ROOT.kBlack);
    h_line2.SetLineStyle(2);
    h_line2.SetStats(0)
    h_line3.SetFillColor(ROOT.kWhite);
    h_line3.SetLineColor(ROOT.kBlack);
    h_line3.SetLineStyle(2);
    h_line3.SetStats(0)
    h_line4.SetFillColor(ROOT.kWhite);
    h_line4.SetLineColor(ROOT.kBlack);
    h_line4.SetLineStyle(2);
    h_line4.SetStats(0)
    h_div.SetLineColor(ROOT.kBlack);
    h_div.SetLineWidth(2);
    h_div.SetMarkerStyle(20)
    h_div.SetMarkerSize(1.0)
    h_div.SetMarkerColor(ROOT.kBlack)
    h_div.SetFillColor(ROOT.kRed+1);
    h_div.SetMarkerSize(0.)
    h_mc.SetFillColorAlpha(12, 0.3)  # Set grey colour (12) and alpha (0.3)
    h_mc.SetMarkerSize(0)
    h_div.SetStats(0)
    h_line.Draw("HIST")
    h_line2.Draw("HIST SAME")
    h_line3.Draw("HIST SAME")
    h_line4.Draw("HIST SAME")
    #h_mc.Draw("E2 SAME")
    #h_div.Draw("E1 P0 SAME")
    h_div.Draw("HIST SAME")
    h_line.SetTitle("")
    #h_line.GetYaxis().SetTitle("Data/MC")
    h_line.GetYaxis().SetTitle("#frac{Data-MC}{#sigma}")
    #h_line.GetYaxis().SetRangeUser(0.6,1.4)
    h_line.GetYaxis().SetRangeUser(-3.,3.)
    h_line.GetXaxis().SetTitle("m_{NN}")

    canvas.SaveAs('fitplot_%s%s.pdf'%(second_dir,args.label))
    pad1.SetLogy()
    canvas.SaveAs('fitplot_%s%s_logy.pdf'%(second_dir,args.label))
    pad1.SetLogy(0)
    #canvas.SaveAs('fitplot_%s%s.png'%(second_dir,args.label))
