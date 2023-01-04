
void makePlot(RooWorkspace* model, std::string reg, std::string chan, std::string dirname, std::string sample, std::string syst, std::string namemod, bool mod) {
    TCanvas *can = new TCanvas("can");
    RooPlot *plotmhist = model->var(("massreg"+namemod).c_str())->frame();
    if (model->data((reg+chan+"2017_"+sample).c_str())) {
        model->data((reg+chan+"2017_"+sample).c_str())->plotOn(plotmhist,RooFit::DrawOption("C"),RooFit::LineColor(kBlack),RooFit::LineWidth(2));
    }
    if (model->data((reg+chan+"2017_"+sample+"_"+syst+(mod ? "_"+chan : "")+"Up").c_str())) {
        model->data((reg+chan+"2017_"+sample+"_"+syst+(mod ? "_"+chan : "")+"Up").c_str())->plotOn(plotmhist,RooFit::DrawOption("C"),RooFit::LineColor(kRed),RooFit::LineWidth(2));
    }
    if (model->data((reg+chan+"2017_"+sample+"_"+syst+(mod ? "_"+chan : "")+"Down").c_str())) {
        model->data((reg+chan+"2017_"+sample+"_"+syst+(mod ? "_"+chan : "")+"Down").c_str())->plotOn(plotmhist,RooFit::DrawOption("C"),RooFit::LineColor(kBlue),RooFit::LineWidth(2));
    }
    plotmhist->Draw();
    can->SaveAs((dirname+"/"+reg+chan+"_"+sample+"_"+syst+".pdf").c_str());
    plotmhist->SetMinimum(0.1);
    can->SetLogy();
    can->SaveAs((dirname+"/"+reg+chan+"_"+sample+"_"+syst+"_logy.pdf").c_str());
}

void runPlotLoop(std::string dirname, std::vector<std::string> regions, std::vector<std::string> chans, std::vector<std::string> systs, std::vector<std::string> samples) {
    for (auto & reg : regions) {
        for (auto & chan : chans) {
            TFile *infile = new TFile((dirname+"/"+chan+"Model.root").c_str());
            //TFile *infile = new TFile((dirname+"/"+chan+"Model_m200/"+chan+"Model.root").c_str());
            RooWorkspace *model = (RooWorkspace*)(infile->Get((chan+"Model").c_str()));
            for (auto & syst : systs) {
                for (auto & sample : samples) {
                    std::cout<<reg<<" - "<<chan<<" - "<<syst<<" - "<<sample<<std::endl;
                    std::string namemod = "";
                    const std::string & onechan = "fail";
                    if (chan=="hadhad" && (0 == reg.compare(reg.size()-onechan.size(), onechan.size(), onechan))) {
                    //if ((chan=="hadhad" && (0 == reg.compare(reg.size()-onechan.size(), onechan.size(), onechan))) || strncmp(reg.c_str(), "ptbin", strlen("ptbin"))) {
                    //if ((chan=="hadhad" && (0 == reg.compare(reg.size()-onechan.size(), onechan.size(), onechan))) || (chan!="hadhad" && strncmp(reg.c_str(), "ptbin", strlen("ptbin")))) {
                        std::cout<<"\tonechan"<<std::endl;
                        namemod = "_one";
                    }
                    makePlot(model, reg, chan, dirname, sample, syst, namemod, (syst.find("massscale") != std::string::npos) || !strncmp(syst.c_str(), "qcd", strlen("qcd")) || (syst.find("highmass") != std::string::npos));
                }
            }
        }
    }
}

void shapePlots() {
    std::string dirname = "Dec12_2017_35_m200";
//    runPlotLoop(dirname, {"ptbin0pass"},
//                         {"hadel","hadmu"},
//                         {"qcd_Rpass"},
//                         {"multijet"});
//    runPlotLoop(dirname, {"wlnuCRpass"},
//                         {"hadel","hadmu"},
//                         {"qcd_wlnu_Rpass"},
//                         {"multijet"});
//    runPlotLoop(dirname, {"topCRpass"},
//                         {"hadel","hadmu"},
//                         {"qcd_top_Rpass"},
//                         {"multijet"});
//    runPlotLoop(dirname, {"ptbin0loosepass"},
//                         {"hadel","hadmu"},
//                         {"qcd_Rloosepass"},
//                         {"multijet"});
//    runPlotLoop(dirname, {"wlnuCRloosepass"},
//                         {"hadel","hadmu"},
//                         {"qcd_wlnu_Rloosepass"},
//                         {"multijet"});
//    runPlotLoop(dirname, {"topCRloosepass"},
//                         {"hadel","hadmu"},
//                         {"qcd_top_Rloosepass"},
//                         {"multijet"});
//    runPlotLoop(dirname, {"wlnuCRfail"},
//                         {"hadel","hadmu"},
//                         {"qcd_wlnu_Rfail"},
//                         {"multijet"});
//    runPlotLoop(dirname, {"ptbin0fail"},
//                         {"hadel","hadmu"},
//                         {"qcd_Rfail"},
//                         {"multijet"});
//    runPlotLoop(dirname, {"topCRfail"},
//                         {"hadel","hadmu"},
//                         {"qcd_top_Rfail"},
//                         {"multijet"});

    runPlotLoop(dirname, {"ptbin0pass"},
                         {"hadhad"},
                         {"qcd_Rpass"},
                         {"multijet"});
    runPlotLoop(dirname, {"wlnuCRpass"},
                         {"hadhad"},
                         {"qcd_wlnu_Rpass"},
                         {"multijet"});
    runPlotLoop(dirname, {"topCRpass"},
                         {"hadhad"},
                         {"qcd_top_Rpass"},
                         {"multijet"});
    runPlotLoop(dirname, {"ptbin0loosepass"},
                         {"hadhad"},
                         {"qcd_Rloosepass"},
                         {"multijet"});
    runPlotLoop(dirname, {"wlnuCRloosepass"},
                         {"hadhad"},
                         {"qcd_wlnu_Rloosepass"},
                         {"multijet"});
    runPlotLoop(dirname, {"topCRloosepass"},
                         {"hadhad"},
                         {"qcd_top_Rloosepass"},
                         {"multijet"});
    runPlotLoop(dirname, {"wlnuCRfail"},
                         {"hadhad"},
                         {"qcd_wlnu_Rfail"},
                         {"multijet"});
    runPlotLoop(dirname, {"ptbin0fail"},
                         {"hadhad"},
                         {"qcd_Rfail"},
                         {"multijet"});
    runPlotLoop(dirname, {"topCRfail"},
                         {"hadhad"},
                         {"qcd_top_Rfail"},
                         {"multijet"});

/*
    runPlotLoop(dirname, {"ptbin0fail","topCRfail","wlnuCRfail","ptbin0loosepass","topCRloosepass","wlnuCRloosepass","ptbin0pass","topCRpass","wlnuCRpass"},
                         {"hadhad","hadmu","hadel"},
                         {"jescale","uescale","l1prefire","jeresol"},
                         //{"top","htt125","dy","wlnu","vvqq"});
                         {"top","htt125","dy","wlnu"});
    runPlotLoop(dirname, {"ptbin0pass","topCRpass","wlnuCRpass","ptbin0loosepass","topCRloosepass","wlnuCRloosepass"},
                         {"hadhad","hadmu","hadel"},
                         {"massscale"},
                         {"htt125","dy","phitt200"});
    runPlotLoop(dirname, {"ptbin0pass","topCRpass","wlnuCRpass","ptbin0loosepass","topCRloosepass","wlnuCRloosepass","ptbin0fail","topCRfail","wlnuCRfail"},
                         {"hadhad","hadmu","hadel"},
                         {"massscale_bkg"},
                         {"wlnu","top"});

    runPlotLoop(dirname, {"ptbin0fail","topCRfail","wlnuCRfail","ptbin0loosepass","topCRloosepass","wlnuCRloosepass","ptbin0pass","topCRpass","wlnuCRpass"},
                         {"hadhad","hadmu","hadel"},
                         {"toppt"},
                         {"top"});

    runPlotLoop(dirname, {"ptbin0pass","topCRpass","wlnuCRpass","ptbin0loosepass","topCRloosepass","wlnuCRloosepass"},
                         {"hadhad","hadmu","hadel"},
                         {"wlnu_highmass_bin14","wlnu_highmass_bin15","wlnu_highmass_bin16","wlnu_highmass_bin17"},
                         {"wlnu"});
    runPlotLoop(dirname, {"ptbin0pass","topCRpass","wlnuCRpass","ptbin0loosepass","topCRloosepass","wlnuCRloosepass"},
                         {"hadhad","hadmu","hadel"},
                         {"top_highmass_bin14","top_highmass_bin15","top_highmass_bin16","top_highmass_bin17"},
                         {"top"});
*/

}
