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

using namespace std;

void mid_rapidity_eta_meson_mix(){
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    read_mix_file();
    cout << "Program ends" << endl;
}


void read_mix_file() {
    // Get a list of all files in the directory
    TString directoryPath = "";
    TSystemDirectory dir("data", directoryPath);
    TList *files = dir.GetListOfFiles();

    if (!files) {
        std::cout << "No files found in the directory." << std::endl;
        return;
    }

    TString fileName = "mixedEventPairs.root";
    TFile* file;
    
    // Loop through the files and open the desired file
    TSystemFile* fileItem;
    TIter next(files);
    while ((fileItem = (TSystemFile*)next())) {
        if (fileItem->IsDirectory()) {
            continue; // Skip directories
        }
        
        if (fileItem->GetName() == fileName) {
            file = TFile::Open(directoryPath + fileName);
            if (file) {
                std::cout << "Opened file: " << directoryPath + fileName << std::endl;
                break; // Exit loop once file is found and opened
            } else {
                std::cerr << "Error: Could not open file " << directoryPath + fileName << std::endl;
            }
        }
    }

    // If file is successfully opened, you can proceed to access the objects inside
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file " << directoryPath + fileName << std::endl;
        return;
    }
    
    TTree* mix_tree = (TTree*)file->Get("Mix");  // Replace "myTree" with the actual tree name
    if (!mix_tree) {
        std::cerr << "Error" << std::endl;
        file->Close();
        return;
    }

    Float_t pt, mass;
    Int_t sector;
    // Set the branch address to the 'pt' branch
    mix_tree->SetBranchAddress("pt", &pt);
    mix_tree->SetBranchAddress("mass", &mass);
    mix_tree->SetBranchAddress("sector", &sector);

    const int nBinsPt = 28;  // Total number of bins
    Float_t ptBins[nBinsPt + 1] = {2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 
                                   7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 12.0, 14.0, 
                                   16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 
                                   35.0, 40.0, 45.0};

    // First row: sector_i = 0-5 are PbSc
    // Second row: sector_i = 6-7 are PbGl
    TH1F* hMassByPt_PbSc[nBinsPt];
    TH1F* hMassByPt_PbGl[nBinsPt];

    for (int j = 0; j < nBinsPt; j++) {
        TString histNameSc = Form("Invariant_Mass_PbSc_ptBin%d", j);  // Unique name for PbSc
        TString histTitleSc = Form("Invariant Mass distribution for PbSc: %.1f <= pt < %.1f", 
                                   ptBins[j], ptBins[j + 1]);
        hMassByPt_PbSc[j] = new TH1F(histNameSc, histTitleSc, 50, 0, 1.5);
        hMassByPt_PbSc[j]->GetXaxis()->SetTitle("Invariant Mass");
        hMassByPt_PbSc[j]->GetYaxis()->SetTitle("Entries");

        TString histNameGl = Form("Invariant_Mass_PbGl_ptBin%d", j);  // Unique name for PbGl
        TString histTitleGl = Form("Invariant Mass distribution for PbGl: %.1f <= pt < %.1f", 
                                   ptBins[j], ptBins[j + 1]);
        hMassByPt_PbGl[j] = new TH1F(histNameGl, histTitleGl, 50, 0, 1.5);
        hMassByPt_PbGl[j]->GetXaxis()->SetTitle("Invariant Mass");
        hMassByPt_PbGl[j]->GetYaxis()->SetTitle("Entries");
    }


    // Loop through the entries and fill the appropriate histogram based on pt
    Long64_t nEntries = mix_tree->GetEntries();
    for (j = 0; j < nEntries; j++) {
        mix_tree->GetEntry(j);

        // Find which pt bin this entry belongs to and fill the corresponding mass histogram
        for (int bin = 0; bin < nBinsPt; bin++) {
            if (pt >= ptBins[bin] && pt < ptBins[bin + 1]) {
                if (sector <= 5) {
                    hMassByPt_PbSc[bin]->Fill(mass);  // Fill PbSc histogram
                } else {
                    hMassByPt_PbGl[bin]->Fill(mass);  // Fill PbGl histogram
                }
                break;
            }
        }
    }

    // mix-event background subtraction
    TFile* filePbSc = TFile::Open("hMassByPt_PbSc_v2.root");
    TFile* filePbGl = TFile::Open("hMassByPt_PbGl_v2.root");

    TH1F* PbSc_arr[nBinsPt];
    TH1F* PbGl_arr[nBinsPt];
    TFile *outputFile_PbSc = new TFile("hMassByPt_PbSc_sub.root", "RECREATE");
    TFile *outputFile_PbGl = new TFile("hMassByPt_PbGl_sub.root", "RECREATE");
    int upper_bin_bound = 10; // up to pT = 7.0 inclusive
    TFile *outputFile_ratio = new TFile("output_ratios.root", "RECREATE");
    for (j = 0; j < upper_bin_bound; j++) {
        TH1F* hPbSc = (TH1F*)filePbSc->Get(Form("Invariant_Mass_PbSc_ptBin%d", j));
        TH1F* hPbGl = (TH1F*)filePbGl->Get(Form("Invariant_Mass_PbGl_ptBin%d", j));
        hPbSc->SetDirectory(0);
        hPbGl->SetDirectory(0);
        /*
        hPbSc->Divide(hMassByPt_PbSc[j]);
        hPbGl->Divide(hMassByPt_PbGl[j]);
        c_ratio->cd(j+1);
        hPbSc->Draw();
        */
        // Clone histograms for division
        TH1F* hPbSc_ratio = (TH1F*)hMassByPt_PbSc[j]->Clone(Form("hPbSc_ratio_ptBin%d", j));
        TH1F* hPbGl_ratio = (TH1F*)hMassByPt_PbGl[j]->Clone(Form("hPbGl_ratio_ptBin%d", j));

        // Perform division
        hPbSc_ratio->Divide(hPbSc);
        hPbGl_ratio->Divide(hPbGl);

        hPbSc_ratio->Write();
        hPbGl_ratio->Write();

        /*
        // normalize by pT bin by pT bin using mass between 0.65, 0.7, 0.75, 0.8
        double PbSc_ratio = hPbSc->Integral(hPbSc->FindBin(0.65), hPbSc->FindBin(0.8)) / hMassByPt_PbSc[j]->Integral(hMassByPt_PbSc[j]->FindBin(0.65), hMassByPt_PbSc[j]->FindBin(0.8));
        double PbGl_ratio = hPbGl->Integral(hPbGl->FindBin(0.65), hPbGl->FindBin(0.8)) / hMassByPt_PbGl[j]->Integral(hMassByPt_PbGl[j]->FindBin(0.65), hMassByPt_PbGl[j]->FindBin(0.8));
        std::cout << "ptbin = " << j << ", PbSc ratio = " << PbSc_ratio << ", PbGl_ratio = " << PbGl_ratio << std::endl;
        // Scale the background histograms by the calculated ratios
        hMassByPt_PbSc[j]->Scale(PbSc_ratio);
        hMassByPt_PbGl[j]->Scale(PbGl_ratio);

        // Subtract the scaled background from the original histograms
        TH1F* hPbSc_sub = (TH1F*)hPbSc->Clone(Form("hPbSc_sub_ptBin%d", j));
        hPbSc_sub->Add(hMassByPt_PbSc[j], -1); // Subtract background
        TH1F* hPbGl_sub = (TH1F*)hPbGl->Clone(Form("hPbGl_sub_ptBin%d", j));
        hPbGl_sub->Add(hMassByPt_PbGl[j], -1); // Subtract background
        
        // Zero out any negative values in the subtracted histograms
        for (int bin = 1; bin <= hPbSc_sub->GetNbinsX(); bin++) {
            if (hPbSc_sub->GetBinContent(bin) < 0) {
                hPbSc_sub->SetBinContent(bin, 0);
            }
        }

        for (int bin = 1; bin <= hPbGl_sub->GetNbinsX(); bin++) {
            if (hPbGl_sub->GetBinContent(bin) < 0) {
                hPbGl_sub->SetBinContent(bin, 0);
            }
        }

        PbSc_arr[j] = hPbSc_sub;
        PbGl_arr[j] = hPbGl_sub;

        outputFile_PbSc->cd();
        PbSc_arr[j]->Write();
        outputFile_PbGl->cd();
        PbGl_arr[j]->Write();
        */
    }
    for (j = upper_bin_bound; j < nBinsPt; j++) {
        TH1F* hPbSc = (TH1F*)filePbSc->Get(Form("Invariant_Mass_PbSc_ptBin%d", j));
        TH1F* hPbGl = (TH1F*)filePbGl->Get(Form("Invariant_Mass_PbGl_ptBin%d", j));
        hPbSc->SetName(Form("hPbSc_sub_ptBin%d", j));  // Renaming the histogram
        hPbGl->SetName(Form("hPbGl_sub_ptBin%d", j));
        outputFile_PbSc->cd();
        hPbSc->Write();
        outputFile_PbGl->cd();
        hPbGl->Write();
    }

    outputFile_PbSc->Close();
    outputFile_PbGl->Close();
    outputFile_ratio->Close();
    filePbSc->Close();
    filePbGl->Close();
}
