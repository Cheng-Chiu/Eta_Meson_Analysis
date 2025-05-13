#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>

void yield_plots() {
    
    
    
    const char* fileNames[4] = {
        "eta_invmass_fits_together_v2_minbias.root",
        "eta_invmass_fits_separate_v2_minbias.root",
        "eta_invmass_fits_together_v2_gamma_minbias.root",
        "eta_invmass_fits_separate_v2_gamma_minbias.root"
    };

    /*
    const char* fileNames[4] = {
        "eta_invmass_fits_together_v2.root",
        "eta_invmass_fits_separate_v2.root",
        "eta_invmass_fits_together_v2_gamma.root",
        "eta_invmass_fits_separate_v2_gamma.root"
    };
    */
 
    const char* fileTitles[4] = {
    "together, bg = P_{2}",
    "separate, bg = P_{2}",
    "together, bg = \\Gamma",
    "separate, bg = \\Gamma"
};
    const char* graphNames[2] = {
        "Yield_vs_Pt_PbSc",
        "Yield_vs_Pt_PbGl"
    };
    
    TLatex latex;
    latex.SetNDC();  // Set NDC (Normalized Device Coordinates)
    latex.SetTextSize(0.05);

    // Colors for the histograms
    Color_t colors[4] = {kRed, kBlue, kGreen+2, kMagenta};

    // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "PbSc", 1200, 600);
    c1->Divide(2, 1);
    // Create a legend for each sub-canvas
    TLegend *legend1 = new TLegend(0.1, 0.1, 0.6, 0.35);
    // TLegend *legend1 = new TLegend(0.4, 0.65, 0.9, 0.9);
    legend1->SetTextSize(0.05);
    TLegend *legend2 = new TLegend(0.1, 0.1, 0.6, 0.35);
    // TLegend *legend2 = new TLegend(0.4, 0.65, 0.9, 0.9);
    legend2->SetTextSize(0.05);

    // Loop over files
    for (int i = 0; i < 4; i++) {
        // Open the ROOT file
        TFile *file = new TFile(fileNames[i]);
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Could not open file " << fileNames[i] << std::endl;
            continue;
        }
        
        // Loop over the TGraph names within each file
        for (int j = 0; j < 2; j++) {
            TGraphErrors *graph = (TGraphErrors*)file->Get(graphNames[j]);
            if (!graph) {
                std::cerr << "Error: Could not find TGraph " << graphNames[j] 
                          << " in file " << fileNames[i] << std::endl;
                continue;
            }

            TGraphErrors *filteredGraph = new TGraphErrors();
            int numPoints = graph->GetN();
            for (int k = 0; k < numPoints; k++) {
                double x, y;
                graph->GetPoint(k, x, y);
                if (x >= 0) {
                // if (x >= 5) {
                    int n = filteredGraph->GetN();
                    double ex = graph->GetErrorX(k);
                    double ey = graph->GetErrorY(k);
                    filteredGraph->SetPoint(n, x, y);
                    filteredGraph->SetPointError(n, ex, ey);
                        
                }
            }

            filteredGraph->SetLineColor(colors[i]);
            filteredGraph->SetLineWidth(2);
            filteredGraph->SetMarkerStyle(20 + j);
            filteredGraph->SetMarkerColor(colors[i]);
            filteredGraph->SetMarkerSize(0.7);
            filteredGraph->SetTitle("");
            filteredGraph->GetYaxis()->SetLabelOffset(0.001);
            filteredGraph->GetYaxis()->SetTitleOffset(0.9);
            filteredGraph->GetXaxis()->SetTitleOffset(0.7);
            filteredGraph->GetXaxis()->SetTitle("pt");
            filteredGraph->GetYaxis()->SetTitle("Yield");
            filteredGraph->GetXaxis()->SetTitleSize(0.05);
            filteredGraph->GetYaxis()->SetTitleSize(0.05);
            
            if (j == 1){
                filteredGraph->GetYaxis()->SetRangeUser(0.1, 120);
            }
            

            c1->cd(j + 1);
            gPad->SetLogy();
            gPad->SetGrid();
            
            if (i == 0) {
                filteredGraph->Draw("APE");
            } else {
                filteredGraph->Draw("PE SAME");
            }
           
            if (j == 0) {
                legend1->AddEntry(filteredGraph, fileTitles[i], "lp");
                latex.DrawLatex(0.75, 0.55, Form("PbSc"));
            } else {
                legend2->AddEntry(filteredGraph, fileTitles[i], "lp");
                latex.DrawLatex(0.75, 0.55, Form("PbGl"));
            }

            
        }
        file->Close();
    }  
    c1->cd(1); 
    legend1->Draw();
    c1->SetTitle("MB Trigger Yields");
    // c1->SetTitle("gamma3 Trigger Yields");
    
    
    c1->cd(2); 
    legend2->Draw();

    c1->SaveAs("yield_dif_fits_MB.pdf");
    // c1->SaveAs("yield_dif_fits_gamma3.pdf");

}
