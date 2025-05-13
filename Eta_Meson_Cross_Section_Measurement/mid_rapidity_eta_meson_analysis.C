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

void mid_rapidity_eta_meson_analysis(){
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    read_input_files();
    cout << "Program ends" << endl;
}


void read_input_files() {
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
    Long64_t j = 0;
    while ((file = next())) {
        TString fileName = file->GetName();
        if (!fileName.Contains("BadRuns/") && 
            !fileName.Contains("se/") && 
            !fileName.Contains("erttrigefficiencies") &&
            !fileName.Contains("mixedEventPairs") &&
            fileName.EndsWith(".root")) {
            cout << fileName << endl;
            TString fullPath = directoryPath + fileName;
            chain_pair->Add(fullPath);
        }
        j++;
    }
    std::cout << "Number of files in Pair chain: " << chain_pair->GetNtrees() << std::endl;
    std::cout << "Number of entries in Pair chain: " << chain_pair->GetEntries() << std::endl;
    
    Float_t pt, mass;
    Int_t sector_i, clustfire_ertgamma3_i, clustfire_ertgamma3_j, evtfire_bbcll1narrowvtx_live;
    // Set the branch address to the 'pt' branch
    chain_pair->SetBranchAddress("pt", &pt);
    chain_pair->SetBranchAddress("mass", &mass);
    chain_pair->SetBranchAddress("sector_i", &sector_i);
    chain_pair->SetBranchAddress("clustfire_ertgamma3_i", &clustfire_ertgamma3_i);
    chain_pair->SetBranchAddress("clustfire_ertgamma3_j", &clustfire_ertgamma3_j);
    chain_pair->SetBranchAddress("evtfire_bbcll1narrowvtx_live", &evtfire_bbcll1narrowvtx_live);
    const int nBinsPt = 28;  // Total number of bins
    Float_t ptBins[nBinsPt + 1] = {2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 
                                   7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 12.0, 14.0, 
                                   16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 
                                   35.0, 40.0, 45.0};

    // First row: sector_i = 0-5 are PbSc
    // Second row: sector_i = 6-7 are PbGl
    TH1F* hMassByPt_PbSc[nBinsPt];
    TH1F* hMassByPt_PbGl[nBinsPt];

    for (j = 0; j < nBinsPt; j++) {
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
    Long64_t nEntries = chain_pair->GetEntries();
    for (j = 0; j < nEntries; j++) {
        chain_pair->GetEntry(j);
        // Cuts
        if (clustfire_ertgamma3_i == 0 && clustfire_ertgamma3_j == 0){
            continue;
        }
        if (evtfire_bbcll1narrowvtx_live == 0){
            continue;
        }
        // Find which pt bin this entry belongs to and fill the corresponding mass histogram
        for (int bin = 0; bin < nBinsPt; bin++) {
            if (pt >= ptBins[bin] && pt < ptBins[bin + 1]) {
                if (sector_i <= 5) {
                    hMassByPt_PbSc[bin]->Fill(mass);  // Fill PbSc histogram
                } else {
                    hMassByPt_PbGl[bin]->Fill(mass);  // Fill PbGl histogram
                }
                break;
            }
        }
    }

    TFile *outputFile_PbSc = new TFile("hMassByPt_PbSc_v2.root", "RECREATE");
    outputFile_PbSc->cd();  // Set PbSc output file as the current active file

    // Write the PbSc histograms to their file
    for (int i = 0; i < nBinsPt; i++) {
        hMassByPt_PbSc[i]->Write();
    }
    outputFile_PbSc->Close();

    TFile *outputFile_PbGl = new TFile("hMassByPt_PbGl_v2.root", "RECREATE");
    outputFile_PbGl->cd();  // Set PbGl output file as the current active file

    // Write the PbGl histograms to their file
    for (int i = 0; i < nBinsPt; i++) {
        hMassByPt_PbGl[i]->Write();
    }
    outputFile_PbGl->Close();

}
