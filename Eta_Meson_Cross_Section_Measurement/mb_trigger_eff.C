#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include <math.h>
#include <vector>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TObject.h>
#include <TString.h>
#include <map>
#include <sstream>

using namespace std;

void mb_trigger_eff(){
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    read_input_files();
    cout << "Program ends" << endl;
}

void read_input_files() {
    map<string, double> liveRateMap;
    ifstream liveRatesFile("liverates.txt");

    if (!liveRatesFile.is_open()) {
        cerr << "Error: Could not open liverates.txt" << endl;
        return;
    }

    string line;
    while (getline(liveRatesFile, line)) {
        if (line.empty()) continue;  // Skip empty lines
        istringstream iss(line);
        string runnumber;
        double liveRates11, liveRates1;

        if (!(iss >> runnumber >> liveRates11 >> liveRates1)) {
            cerr << "Skipping invalid line: " << line << endl;
            continue;
        }

        if (liveRates1 == 0) {
            liveRateMap[runnumber] = -1.0; // divide by zero
        }else{
             double ratio = liveRates1 / liveRates11;
             liveRateMap[runnumber] = ratio;
        }
    }
    liveRatesFile.close();

    cout << "runnumber -> liverates_1/liverates_11 Ratio" << endl;

    if (liveRateMap.empty()) {
        cerr << "Error: No valid data found in the file!" << endl;
    }
    //printing map values
    for (std::map<string, double>::iterator it = liveRateMap.begin(); it != liveRateMap.end(); ++it) {
        std::cout << it->first << ": " << it->second << std::endl;
    }
    

    // Create TChain objects for the trees you want to read
    TChain* chain_pair = new TChain("Pair");

    // Get a list of all files in the directory
    TString directoryPath = "";
    TSystemDirectory dir("data", directoryPath);
    TList *files = dir.GetListOfFiles();

    if (!files) {
        std::cout << "No files found in the directory." << std::endl;
        return;
    }

    // Loop through the list of files
   
    TIter next(files);
    TObject *file;
    int j = 0;
    while ((file = next())) {
        // if (j > 10) break;
        TString fileName = file->GetName();
        if (!fileName.Contains("BadRuns/") && 
            !fileName.Contains("se/") && 
            !fileName.Contains("erttrigefficiencies") &&
            !fileName.Contains("mixedEventPairs") &&
            fileName.EndsWith(".root")) {
            cout << fileName << endl;
            TString fullPath = directoryPath + fileName;
            chain_pair->Add(fullPath);
            j++;
        }
    }
    std::cout << "Number of files in Pair chain: " << chain_pair->GetNtrees() << std::endl;
    std::cout << "Number of entries in Pair chain: " << chain_pair->GetEntries() << std::endl;
    
    Float_t pt, mass;
    Int_t sector_i, evtfire_ertb, evtfire_bbcll1narrowvtx_live, evtfire_bbcll1novtx_live;
    // Set the branch address to the 'pt' branch
    chain_pair->SetBranchAddress("pt", &pt);
    chain_pair->SetBranchAddress("mass", &mass);
    chain_pair->SetBranchAddress("sector_i", &sector_i);
    chain_pair->SetBranchAddress("evtfire_ertb", &evtfire_ertb);
    chain_pair->SetBranchAddress("evtfire_bbcll1narrowvtx_live", &evtfire_bbcll1narrowvtx_live);
    chain_pair->SetBranchAddress("evtfire_bbcll1novtx_live", &evtfire_bbcll1novtx_live);
    const int nBinsPt = 11;  // Total number of bins
    Float_t ptBins[nBinsPt + 1] = {2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 
                                   7.5, 8.0};

    // First row: sector_i = 0-5 are PbSc
    // Second row: sector_i = 6-7 are PbGl
    TH1D* PbSc_mass[nBinsPt];
    TH1D* PbSc_mass_tot[nBinsPt];
    TH1D* PbGl_mass[nBinsPt];
    TH1D* PbGl_mass_tot[nBinsPt];
    for (j = 0; j < nBinsPt; j++) {
        TString histNameSc = Form("Triggered_Mass_PbSc_ptBin%d", j);
        TString histTitleSc = Form("Triggered Mass distribution for PbSc: %.1f <= pt < %.1f", 
                                   ptBins[j], ptBins[j + 1]);
        PbSc_mass[j] = new TH1D(histNameSc, histTitleSc, 50, 0, 1.5);
        PbSc_mass[j]->GetXaxis()->SetTitle("Invariant Mass");
        PbSc_mass[j]->GetYaxis()->SetTitle("Entries");

        TString histNameGl = Form("Triggered_Mass_PbGl_ptBin%d", j);
        TString histTitleGl = Form("Triggered Mass distribution for PbGl: %.1f <= pt < %.1f", 
                                   ptBins[j], ptBins[j + 1]);
        PbGl_mass[j] = new TH1D(histNameGl, histTitleGl, 50, 0, 1.5);
        PbGl_mass[j]->GetXaxis()->SetTitle("Invariant Mass");
        PbGl_mass[j]->GetYaxis()->SetTitle("Entries");

        TString histNameSc_tot = Form("Total_Mass_PbSc_ptBin%d", j);
        TString histTitleSc_tot = Form("Total Mass distribution for PbSc: %.1f <= pt < %.1f", 
                                   ptBins[j], ptBins[j + 1]);
        PbSc_mass_tot[j] = new TH1D(histNameSc_tot, histTitleSc_tot, 50, 0, 1.5);
        PbSc_mass_tot[j]->GetXaxis()->SetTitle("Invariant Mass");
        PbSc_mass_tot[j]->GetYaxis()->SetTitle("Entries");

        TString histNameGl_tot = Form("Total_Mass_PbGl_tot_ptBin%d", j);
        TString histTitleGl_tot = Form("Total Mass distribution for PbGl: %.1f <= pt < %.1f", 
                                   ptBins[j], ptBins[j + 1]);
        PbGl_mass_tot[j] = new TH1D(histNameGl_tot, histTitleGl_tot, 50, 0, 1.5);
        PbGl_mass_tot[j]->GetXaxis()->SetTitle("Invariant Mass");
        PbGl_mass_tot[j]->GetYaxis()->SetTitle("Entries");

        // enable error propagation
        PbSc_mass[j]->Sumw2();
        PbSc_mass_tot[j]->Sumw2();
        PbGl_mass[j]->Sumw2();
        PbGl_mass_tot[j]->Sumw2();
    }

    // Loop through the entries and fill the appropriate histogram based on pt
    Long64_t nEntries = chain_pair->GetEntries();
    for (int i = 0; i < nEntries; i++) {
        chain_pair->GetEntry(i);
        if (i % 100000 == 0) {
            std::cout << "Processing entry " << i << " out of " << nEntries << "..." << std::endl;
        }
        if (evtfire_ertb){
             for (int bin = 0; bin < nBinsPt; bin++) {
                 if (pt >= ptBins[bin] && pt < ptBins[bin + 1]) {
                     if (sector_i <= 5 && sector_i >= 0){ //PbSc
                         PbSc_mass_tot[bin]->Fill(mass);
                         // if (evtfire_bbcll1narrowvtx_live){
                         if (evtfire_bbcll1novtx_live){
                             PbSc_mass[bin]->Fill(mass);
                         }
                     }else if(sector_i == 6 || sector_i == 7){ //PbGl at west arm
                         PbGl_mass_tot[bin]->Fill(mass);
                         // if (evtfire_bbcll1narrowvtx_live){
                         if (evtfire_bbcll1novtx_live){
                             PbGl_mass[bin]->Fill(mass);
                         }
                     }
                     break;
                 } 
             }
        }
    }
    Float_t eff[2][nBinsPt];
    Float_t eff_tot[2][nBinsPt];
    for (int bin = 0; bin < nBinsPt; bin++) {
        eff[0][bin] = PbSc_mass[bin]->GetEntries();
        eff_tot[0][bin] = PbSc_mass_tot[bin]->GetEntries();
        eff[1][bin] = PbGl_mass[bin]->GetEntries();
        eff_tot[1][bin] = PbGl_mass_tot[bin]->GetEntries();   
    }
    Float_t centerPt[nBinsPt];
    Float_t binwidth[nBinsPt];
    for (int bin = 0; bin < nBinsPt; bin++) {
        binwidth[bin] = ptBins[bin+1] - ptBins[bin];
        centerPt[bin] = (ptBins[bin] + ptBins[bin + 1]) / 2.0;
    }
    TFile *outputFile = new TFile("mb_trig_eff_ERT_result_no_frange.root", "RECREATE");
     TCanvas *ceff = new TCanvas("ceff", "Efficiencies", 800, 600);

    TGraphErrors *graphEffSc = new TGraphErrors(nBinsPt);
    TGraphErrors *graphEffGl = new TGraphErrors(nBinsPt);

    for (int bin = 0; bin < nBinsPt; bin++) {
        graphEffSc->SetPoint(bin, centerPt[bin], eff[0][bin] / eff_tot[0][bin]);
        graphEffGl->SetPoint(bin, centerPt[bin], eff[1][bin] / eff_tot[1][bin]);
        graphEffSc->SetPointError(bin, binwidth[bin], eff[0][bin] / eff_tot[0][bin] * sqrt(1./ eff[0][bin] + 1./eff_tot[0][bin]));
        graphEffGl->SetPointError(bin, binwidth[bin], eff[1][bin] / eff_tot[1][bin] * sqrt(1./ eff[1][bin] + 1./eff_tot[1][bin]));
    }
    graphEffSc->SetName("graphEffSc");
    graphEffGl->SetName("graphEffGl");
    graphEffSc->SetTitle("");
    graphEffGl->SetTitle("");
    graphEffSc->GetYaxis()->SetRangeUser(0, 1.5);
    graphEffSc->GetXaxis()->SetRangeUser(0, 45);
    graphEffSc->SetLineColor(kRed);
    graphEffGl->SetLineColor(kGreen+2);

    graphEffSc->SetLineWidth(2);
    graphEffGl->SetLineWidth(2);

    graphEffSc->GetXaxis()->SetTitle("pT");
    graphEffSc->GetYaxis()->SetTitle("Efficiency");

    graphEffSc->Write();
    graphEffGl->Write();

    graphEffSc->Draw("APE");
    graphEffGl->Draw("PE same");
    
    TLegend *legend = new TLegend(0.7, 0.1, 0.9, 0.3);
    legend->AddEntry(graphEffSc, "PbSc Eff", "l");
    legend->AddEntry(graphEffGl, "PbGl Eff", "l");
    legend->Draw();
    ceff->Write();

    outputFile->Close();




    // Invariant Mass Fit
    /*
    TF1* drawboth[2][nBinsPt];
    TF1* drawboth_tot[2][nBinsPt];
    
    for (int bin = 0; bin < nBinsPt; bin++) {
        double drawmin = 0.25;
        double drawmax = 0.9;
        if (pt <= 4){
            drawmax = 0.75;
        }
        if (bin % 1 == 0) {
            std::cout << "Fitting pT bin = " << bin  << " out of " << nBinsPt << "..." << std::endl;
        }
        if(bin >= 7){
            // drawmin = 0.3;
            drawmax = 1.1;
        }else if (bin > 4){
            drawmax = 0.9;
        }
        // Define the fit function
        drawboth[0][bin] = new TF1(Form("fit_PbSc_ptBin%d", bin), "[0]*TMath::GammaDist(x,[1],[2],[3]) + [4]/[6]*exp(-0.5*pow((x-[5])/[6], 2))/sqrt(2 * TMath::Pi())", drawmin, drawmax);
        drawboth[1][bin] = new TF1(Form("fit_PbGl_ptBin%d", bin), "[0]*TMath::GammaDist(x,[1],[2],[3]) + [4]/[6]*exp(-0.5*pow((x-[5])/[6], 2))/sqrt(2 * TMath::Pi())", drawmin, drawmax);
        drawboth_tot[0][bin] = new TF1(Form("fit_PbSc_tot_ptBin%d", bin), "[0]*TMath::GammaDist(x,[1],[2],[3]) + [4]/[6]*exp(-0.5*pow((x-[5])/[6], 2))/sqrt(2 * TMath::Pi())", drawmin, drawmax);
        drawboth_tot[1][bin] = new TF1(Form("fit_PbGl_tot_ptBin%d", bin), "[0]*TMath::GammaDist(x,[1],[2],[3]) + [4]/[6]*exp(-0.5*pow((x-[5])/[6], 2))/sqrt(2 * TMath::Pi())", drawmin, drawmax);

        // Set parameters
        if (bin == 0){
            drawboth[0][bin]->SetParameters(1e4, 2.0, 0.15, 0.15, 5000, 0.55, 0.04);
            drawboth[1][bin]->SetParameters(1e4, 2.0, 0.15, 0.15, 5000, 0.55, 0.04);
            drawboth_tot[0][bin]->SetParameters(1e4, 2.0, 0.15, 0.15, 5000, 0.55, 0.04);
            drawboth_tot[1][bin]->SetParameters(1e4, 2.0, 0.15, 0.15, 5000, 0.55, 0.04);
            std::cout << "TEST"  << std::endl;
        }else{
            drawboth[0][bin]->SetParameters(drawboth[0][bin-1]->GetParameter(0), drawboth[0][bin-1]->GetParameter(1), drawboth[0][bin-1]->GetParameter(2), drawboth[0][bin-1]->GetParameter(3), drawboth[0][bin-1]->GetParameter(4), drawboth[0][bin-1]->GetParameter(5), drawboth[0][bin-1]->GetParameter(6));
            drawboth[1][bin]->SetParameters(drawboth[1][bin-1]->GetParameter(0), drawboth[1][bin-1]->GetParameter(1), drawboth[1][bin-1]->GetParameter(2), drawboth[1][bin-1]->GetParameter(3), drawboth[1][bin-1]->GetParameter(4), drawboth[1][bin-1]->GetParameter(5), drawboth[1][bin-1]->GetParameter(6));
            drawboth_tot[0][bin]->SetParameters(drawboth_tot[0][bin-1]->GetParameter(0), drawboth_tot[0][bin-1]->GetParameter(1), drawboth_tot[0][bin-1]->GetParameter(2), drawboth_tot[0][bin-1]->GetParameter(3), drawboth_tot[0][bin-1]->GetParameter(4), drawboth_tot[0][bin-1]->GetParameter(5), drawboth_tot[0][bin-1]->GetParameter(6));
            drawboth_tot[1][bin]->SetParameters(drawboth_tot[1][bin-1]->GetParameter(0), drawboth_tot[1][bin-1]->GetParameter(1), drawboth_tot[1][bin-1]->GetParameter(2), drawboth_tot[1][bin-1]->GetParameter(3), drawboth_tot[1][bin-1]->GetParameter(4), drawboth_tot[1][bin-1]->GetParameter(5), drawboth_tot[1][bin-1]->GetParameter(6));
        }
        
        drawboth[0][bin]->SetParLimits(0, 0, 1000000);
        drawboth[0][bin]->SetParLimits(1, 0.6, 50.5);
        drawboth[0][bin]->SetParLimits(2, 0, 0.24);
        drawboth[0][bin]->SetParLimits(3, 0.01, 0.25);
        drawboth[0][bin]->SetParLimits(4, 0, 10000000);
        drawboth[0][bin]->SetParLimits(5, 0.53, 0.57);
        drawboth[0][bin]->SetParLimits(6, 0.02, 0.05);

        drawboth[1][bin]->SetParLimits(0, 0, 1000000);
        drawboth[1][bin]->SetParLimits(1, 0.6, 50.5);
        drawboth[1][bin]->SetParLimits(2, 0, 0.24);
        drawboth[1][bin]->SetParLimits(3, 0.01, 0.25);
        drawboth[1][bin]->SetParLimits(4, 0, 10000000);
        drawboth[1][bin]->SetParLimits(5, 0.53, 0.57);
        drawboth[1][bin]->SetParLimits(6, 0.02, 0.05);

        drawboth_tot[0][bin]->SetParLimits(0, 0, 1000000);
        drawboth_tot[0][bin]->SetParLimits(1, 0.6, 50.5);
        drawboth_tot[0][bin]->SetParLimits(2, 0, 0.24);
        drawboth_tot[0][bin]->SetParLimits(3, 0.01, 0.25);
        drawboth_tot[0][bin]->SetParLimits(4, 0, 10000000);
        drawboth_tot[0][bin]->SetParLimits(5, 0.53, 0.57);
        drawboth_tot[0][bin]->SetParLimits(6, 0.02, 0.05);

        drawboth_tot[1][bin]->SetParLimits(0, 0, 1000000);
        drawboth_tot[1][bin]->SetParLimits(1, 0.6, 50.5);
        drawboth_tot[1][bin]->SetParLimits(2, 0, 0.24);
        drawboth_tot[1][bin]->SetParLimits(3, 0.01, 0.25);
        drawboth_tot[1][bin]->SetParLimits(4, 0, 10000000);
        drawboth_tot[1][bin]->SetParLimits(5, 0.53, 0.57);
        drawboth_tot[1][bin]->SetParLimits(6, 0.02, 0.05);


        PbSc_mass[bin]->Fit(drawboth[0][bin], "RMQ");
        PbSc_mass_tot[bin]->Fit(drawboth_tot[0][bin], "RMQ");
        PbGl_mass[bin]->Fit(drawboth[1][bin], "RMQ");
        PbGl_mass_tot[bin]->Fit(drawboth_tot[1][bin], "RMQ");
       

        outputFile->cd();

        TCanvas *c_PbSc = new TCanvas(Form("c_PbSc_bin%d", bin), "PbSc Mass Plot", 800, 600);
        PbSc_mass[bin]->Draw();
        drawboth[0][bin]->Draw("same");
        c_PbSc->Write();

        TCanvas *c_PbSc_tot = new TCanvas(Form("c_PbSc_tot_bin%d", bin), "PbSc Mass Tot Plot", 800, 600);
        PbSc_mass_tot[bin]->Draw();
        drawboth_tot[0][bin]->Draw("same");
        c_PbSc_tot->Write();

        TCanvas *c_PbGl = new TCanvas(Form("c_PbGl_bin%d", bin), "PbGl Mass Plot", 800, 600);
        PbGl_mass[bin]->Draw();
        drawboth[1][bin]->Draw("same");
        c_PbGl->Write();

        TCanvas *c_PbGl_tot = new TCanvas(Form("c_PbGl_tot_bin%d", bin), "PbGl Mass Tot Plot", 800, 600);
        PbGl_mass_tot[bin]->Draw();
        drawboth_tot[1][bin]->Draw("same");
        c_PbGl_tot->Write();

        delete c_PbSc;
        delete c_PbSc_tot;
        delete c_PbGl;
        delete c_PbGl_tot;
    }

    // plot with widths
    Float_t width[2][nBinsPt];
    Float_t widthError[2][nBinsPt];
    Float_t width_tot[2][nBinsPt];
    Float_t widthError_tot[2][nBinsPt];

    Float_t peak[2][nBinsPt];
    Float_t peak_tot[2][nBinsPt];

    for (int i = 0; i < 2; i++){
        for (int j = 0; j < nBinsPt; j++){
            width[i][j] = drawboth[i][j]->GetParameter(6);
            widthError[i][j] = drawboth[i][j]->GetParError(6);
            width_tot[i][j] = drawboth_tot[i][j]->GetParameter(6);
            widthError_tot[i][j] = drawboth_tot[i][j]->GetParError(6);
            peak[i][j] = drawboth[i][j]->GetParameter(5);
            peak_tot[i][j] = drawboth_tot[i][j]->GetParameter(5);
        }
    }
    TCanvas *cWidth = new TCanvas("cWidth", "Mass Width vs Pt", 1200, 800);
    cWidth->cd();

    TGraphErrors *graph1_width = new TGraphErrors(nBinsPt, centerPt, width[0], 0, widthError[0]);
    graph1_width->SetName("PeakWidth_vs_Pt_PbSc");
    graph1_width->SetLineColor(kRed);
    graph1_width->SetMarkerColor(kRed);
    graph1_width->SetMarkerStyle(20);
    graph1_width->SetTitle(""); 
    graph1_width->GetXaxis()->SetTitle("Pt");
    graph1_width->GetYaxis()->SetTitle("Mass Width");
    graph1_width->GetYaxis()->SetRangeUser(0.0, 0.1);
    graph1_width->Draw("AP");

    // Graph for PbGl (individual)
    TGraphErrors *graph2_width = new TGraphErrors(nBinsPt, centerPt, width[1], 0, widthError[1]);
    graph2_width->SetName("PeakWidth_vs_Pt_PbGl");
    graph2_width->SetLineColor(kBlue);
    graph2_width->SetMarkerColor(kBlue);
    graph2_width->SetMarkerStyle(21);
    graph2_width->Draw("P same");

    // Graph for PbSc (total)
    TGraphErrors *graph1_width_tot = new TGraphErrors(nBinsPt, centerPt, width_tot[0], 0, widthError_tot[0]);
    graph1_width_tot->SetName("PeakWidth_vs_Pt_PbSc_Tot");
    graph1_width_tot->SetLineColor(kGreen+2);
    graph1_width_tot->SetMarkerColor(kGreen+2);
    graph1_width_tot->SetMarkerStyle(22);
    graph1_width_tot->Draw("P same");

    // Graph for PbGl (total)
    TGraphErrors *graph2_width_tot = new TGraphErrors(nBinsPt, centerPt, width_tot[1], 0, widthError_tot[1]);
    graph2_width_tot->SetName("PeakWidth_vs_Pt_PbGl_Tot");
    graph2_width_tot->SetLineColor(kMagenta);
    graph2_width_tot->SetMarkerColor(kMagenta);
    graph2_width_tot->SetMarkerStyle(23);
    graph2_width_tot->Draw("P same");

    // Add a legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(graph1_width, "PbSc (Individual)", "lp");
    legend->AddEntry(graph2_width, "PbGl (Individual)", "lp");
    legend->AddEntry(graph1_width_tot, "PbSc (Total)", "lp");
    legend->AddEntry(graph2_width_tot, "PbGl (Total)", "lp");
    legend->SetBorderSize(0);
    legend->SetTextSize(0.03);
    legend->Draw();

    cWidth->Write();

    // integrate the histograms and calculate the efficiencies   
    
    for (int bin = 0; bin < nBinsPt; bin++) {
        eff[0][bin] = PbSc_mass[bin]->Integral(PbSc_mass[bin]->FindBin(peak[0][bin] - 3*width[0][bin]), 
                                               PbSc_mass[bin]->FindBin(peak[0][bin] + 3*width[0][bin]));
        eff_tot[0][bin] = PbSc_mass_tot[bin]->Integral(PbSc_mass_tot[bin]->FindBin(peak[0][bin] - 3*width[0][bin]),
                                                       PbSc_mass_tot[bin]->FindBin(peak[0][bin] + 3*width[0][bin]));
        eff[1][bin] = PbGl_mass[bin]->Integral(PbGl_mass[bin]->FindBin(peak[1][bin] - 3*width[1][bin]), 
                                               PbGl_mass[bin]->FindBin(peak[1][bin] + 3*width[1][bin]));
        eff_tot[1][bin] = PbGl_mass_tot[bin]->Integral(PbGl_mass_tot[bin]->FindBin(peak[1][bin] - 3*width[1][bin]),
                                                       PbGl_mass_tot[bin]->FindBin(peak[1][bin] + 3*width[1][bin]));
    }
*/
}
