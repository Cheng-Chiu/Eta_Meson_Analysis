#include <vector>
#include <TString.h>
#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>

void comp_ERTB_GAMMA3(){

    TFile *ERTB_file = TFile::Open("invmass_tofshapecut.root");
    TFile *GAMMA3_file = TFile::Open("invmass_tofshapecut.root");

    if (!ERTB_file || !GAMMA3_file) {
        cout << "Error: Could not open ERTB or GAMMA3 file!" << endl;
        return;
    }
    TFile *outputFile = new TFile("comparison_plots.root", "RECREATE");
    vector<TString> histTypes;
    histTypes.push_back("PbScWest");
    histTypes.push_back("PbScEast");
    histTypes.push_back("PbGl");

    for (int i = 1; i <= 27; i++) {
        for (int j = 0; j < 3; j++) {
            TString type = histTypes[j];
            TString histName = Form("Invariant_Mass_%s_ptBin%d", type.Data(), i);

            // Retrieve histograms from both files
            TH1F *ERTB_hist = (TH1F*)ERTB_file->Get(histName);
            TH1F *GAMMA3_hist = (TH1F*)GAMMA3_file->Get(histName);

            if (!ERTB_hist || !GAMMA3_hist) {
                cout << "Warning: Histogram " << histName << " not found in one of the files!" << endl;
                continue;
            }
      
            // Set styles
            ERTB_hist->SetLineColor(kRed);
            ERTB_hist->SetLineWidth(2);

            GAMMA3_hist->SetLineColor(kBlue);
            GAMMA3_hist->SetLineWidth(2);
            GAMMA3_hist->SetLineStyle(2); // Dashed line

            TCanvas *c = new TCanvas(histName, histName, 800, 600);
            ERTB_hist->Draw();
            GAMMA3_hist->Draw("same");

            TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
            leg->AddEntry(ERTB_hist, "ERTB", "l");
            leg->AddEntry(GAMMA3_hist, "GAMMA3", "l");
            leg->Draw();
            
            outputFile->cd();
            c->Write();
        }
    }

    ERTB_file->Close();
    GAMMA3_file->Close();
}
