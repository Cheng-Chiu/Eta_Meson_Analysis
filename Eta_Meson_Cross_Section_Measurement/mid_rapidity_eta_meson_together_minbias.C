#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <algorithm>



void mid_rapidity_eta_meson_together_minbias() {
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1112);

    total_fit();
}


void total_fit(){
    const int nBinsPt = 11;  // Total number of bins
    TF1* drawboth[2][nBinsPt];
    TF1* background[2][nBinsPt];
    Float_t ptBins[nBinsPt + 1] = {2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 
                                   7.5, 8.0};
    Float_t centerPt[nBinsPt];
    Float_t binwidth[nBinsPt];
    for (int i = 0; i < nBinsPt; i++) {
        binwidth[i] = ptBins[i+1] - ptBins[i];
        centerPt[i] = (ptBins[i] + ptBins[i + 1]) / 2.0;
    }
    

    TFile* filePbSc = TFile::Open("hMassByPt_PbSc_minbias.root");
    TFile* filePbGl = TFile::Open("hMassByPt_PbGl_minbias.root");
    TFile *outfile = new TFile("eta_invmass_fits_together_v2_minbias.root","RECREATE");
    if (!filePbSc || !filePbGl) {
        std::cerr << "Error opening files!" << std::endl;
        return;
    }
    TCanvas* cPbSc = new TCanvas("cFit_PbSc", "Fit for PbSc ptBin", 3600, 2400);
    cPbSc->Divide(4,3);
    TCanvas* cPbGl = new TCanvas("cFit_PbGl", "Fit for PbGl ptBin", 3600, 2400);
    cPbGl->Divide(4,3);

    // Loop over the bins and fit histograms
    for (int pt = 0; pt < nBinsPt; pt++) {
        double drawmin = 0.25;
        double drawmax = 0.75;
        if(pt >= 7){
            // drawmin = 0.3;
            drawmax = 1.1;
        }else if (pt > 4){
            drawmax = 0.9;
        }

        TH1F* hPbSc = (TH1F*)filePbSc->Get(Form("Invariant_Mass_PbSc_ptBin%d", pt));
        TH1F* hPbGl = (TH1F*)filePbGl->Get(Form("Invariant_Mass_PbGl_ptBin%d", pt));
        
        if (!hPbSc || !hPbGl) {
            std::cerr << "Error retrieving histograms!" << std::endl;
            continue;
        }
        // Define the fit function
        drawboth[0][pt] = new TF1(Form("fit_PbSc_ptBin%d", pt), "pol2 + [4]/[6]*exp(-0.5*pow((x-[5])/[6], 2))/sqrt(2 * TMath::Pi())", drawmin, drawmax);
        drawboth[1][pt] = new TF1(Form("fit_PbGl_ptBin%d", pt), "pol2 + [4]/[6]*exp(-0.5*pow((x-[5])/[6], 2))/sqrt(2 * TMath::Pi())", drawmin, drawmax);
        
        // Set parameters
        drawboth[0][pt]->FixParameter(3,0);
        drawboth[1][pt]->FixParameter(3,0);
        if (pt == 0){
            drawboth[0][pt]->SetParameters(1e6, -1e6, 1e6, 0, 5000, 0.55, 0.04);
            drawboth[1][pt]->SetParameters(1e6, -1e6, 1e6, 0, 5000, 0.55, 0.04);
        }else{
            drawboth[0][pt]->SetParameters(drawboth[0][pt-1]->GetParameter(0), drawboth[0][pt-1]->GetParameter(1), drawboth[0][pt-1]->GetParameter(2),
                                           drawboth[0][pt-1]->GetParameter(3), drawboth[0][pt-1]->GetParameter(4), drawboth[0][pt-1]->GetParameter(5), 
                                           drawboth[0][pt-1]->GetParameter(6));
            drawboth[1][pt]->SetParameters(drawboth[1][pt-1]->GetParameter(0), drawboth[1][pt-1]->GetParameter(1), drawboth[1][pt-1]->GetParameter(2),
                                           drawboth[1][pt-1]->GetParameter(3), drawboth[1][pt-1]->GetParameter(4), drawboth[1][pt-1]->GetParameter(5), 
                                           drawboth[1][pt-1]->GetParameter(6));
        }

        drawboth[0][pt]->SetParLimits(4, 0, 10000000);
        drawboth[0][pt]->SetParLimits(5, 0.53, 0.57);
        drawboth[0][pt]->SetParLimits(6, 0.02, 0.05);

        drawboth[1][pt]->SetParLimits(4, 0, 10000000);
        drawboth[1][pt]->SetParLimits(5, 0.53, 0.57);
        drawboth[1][pt]->SetParLimits(6, 0.02, 0.05);

        // Fit the histograms
        hPbSc->Fit(drawboth[0][pt], "RMQ");
        hPbGl->Fit(drawboth[1][pt], "RMQ");

        background[0][pt] = new TF1(Form("background_PbSc_ptBin%d", pt), "pol2", drawmin, drawmax);
        background[0][pt]->SetParameters(drawboth[0][pt]->GetParameter(0), drawboth[0][pt]->GetParameter(1), drawboth[0][pt]->GetParameter(2), drawboth[0][pt]->GetParameter(3));
        background[1][pt] = new TF1(Form("background_PbGl_ptBin%d", pt), "pol2", drawmin, drawmax);
        background[1][pt]->SetParameters(drawboth[1][pt]->GetParameter(0), drawboth[1][pt]->GetParameter(1), drawboth[1][pt]->GetParameter(2), drawboth[1][pt]->GetParameter(3));

        
        // PbGl
        hPbGl->Draw();
        drawboth[1][pt]->Draw("same");
        background[1][pt]->SetLineColor(kGreen);
        background[1][pt]->Draw("same");

        TLegend* legendGl = new TLegend(0.7, 0.15, 0.9, 0.3);
        legendGl->AddEntry(hPbGl, "Data", "l");
        legendGl->AddEntry(drawboth[1][pt], "Fit", "l");
        legendGl->AddEntry(background[1][pt], "B", "l");
        legendGl->Draw("same");        

        // PbSc
        cPbSc->cd(pt+1);
        hPbSc->Draw();
        drawboth[0][pt]->Draw("same");
        background[0][pt]->SetLineColor(kGreen);
        background[0][pt]->Draw("same");

        TLegend* legendSc = new TLegend(0.7, 0.15, 0.9, 0.3);
        legendSc->AddEntry(hPbSc, "Data", "l");
        legendSc->AddEntry(drawboth[0][pt], "Fit", "l");
        legendSc->AddEntry(background[0][pt], "B", "l");
        legendSc->Draw("same");
        
        hPbSc->Write();
        hPbGl->Write();
    }

    cPbGl->SaveAs(Form("fit_PbGl_together_minbias.png"));
    cPbGl->SaveAs(Form("fit_PbGl_together_minbias.pdf"));
    cPbSc->SaveAs(Form("fit_PbSc_together_minbias.png"));
    cPbSc->SaveAs(Form("fit_PbSc_together_minbias.pdf"));
    
    // Plot Yield
    Float_t yield[2][nBinsPt];
    Float_t yieldError[2][nBinsPt];
    Float_t rapd = 0.35;
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < nBinsPt; j++){
            yield[i][j] = drawboth[i][j]->GetParameter(4)/rapd/binwidth[j];
            yieldError[i][j] = drawboth[i][j]->GetParError(4)/rapd/binwidth[j];
        }
    }

    TCanvas *cYield = new TCanvas("cYield", "Yield vs Pt", 1200, 800);
    cYield->cd();
    cYield->SetLogy();

    TGraphErrors *graph1 = new TGraphErrors(nBinsPt, centerPt, yield[0], 0, yieldError[0]);
    graph1->SetName("Yield_vs_Pt_PbSc");
    graph1->SetTitle("Yield_vs_Pt_PbSc");
    graph1->SetLineColor(kRed);       // Set line color to red
    graph1->SetMarkerColor(kRed);     // Set marker color to red
    graph1->SetMarkerStyle(20);       // Set marker style to circle (20)
    graph1->SetTitle(""); 
    graph1->GetYaxis()->SetRangeUser(1e-2, 1e5);
    graph1->GetXaxis()->SetTitle("Pt");
    graph1->GetYaxis()->SetTitle("Yield");
    graph1->Draw("AP");               // Draw with axis and points

    TGraphErrors *graph2 = new TGraphErrors(nBinsPt, centerPt, yield[1], 0, yieldError[1]);
    graph2->SetName("Yield_vs_Pt_PbGl");
    graph2->SetTitle("Yield_vs_Pt_PbGl");
    graph1->GetXaxis()->SetTitle("Pt");
    graph1->GetYaxis()->SetTitle("Yield");
    graph2->SetLineColor(kBlue);      // Set line color to blue
    graph2->SetMarkerColor(kBlue);    // Set marker color to blue
    graph2->SetMarkerStyle(20);       // Set marker style to circle (20)
    graph2->GetYaxis()->SetRangeUser(1e-2, 1e5);
    graph2->Draw("P same");           // Draw the second graph on the same plot

    // Add a legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(graph1, "Yield for PbSc", "p");
    legend->AddEntry(graph2, "Yield for PbGl", "p");
    legend->Draw();
    
    graph1->Write();
    graph2->Write();
    // cYield->SaveAs("Yield_vs_Pt.png");
    // cYield->SaveAs("Yield_vs_Pt.pdf");


    // Plot Peak Position
    Float_t peak[2][nBinsPt];
    Float_t peakError[2][nBinsPt];
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < nBinsPt; j++){
            peak[i][j] = drawboth[i][j]->GetParameter(5);
            peakError[i][j] = drawboth[i][j]->GetParError(5);
        }
    }

    TCanvas *cPeak = new TCanvas("cPeak", "Peak Position vs Pt", 1200, 800);
    cPeak->cd();

    TGraphErrors *graph1_peak = new TGraphErrors(nBinsPt, centerPt, peak[0], 0, peakError[0]);
    graph1_peak->SetName("PeakPosition_vs_Pt_PbSc");
    graph1_peak->SetTitle("PeakPosition_vs_Pt_PbSc");
    graph1_peak->SetLineColor(kRed);       // Set line color to red
    graph1_peak->SetMarkerColor(kRed);     // Set marker color to red
    graph1_peak->SetMarkerStyle(20);       // Set marker style to circle (20)
    graph1_peak->SetTitle(""); 
    graph1_peak->GetXaxis()->SetTitle("Pt");
    graph1_peak->GetYaxis()->SetTitle("Peak Position");
    graph1_peak->GetYaxis()->SetRangeUser(0.5, 0.6);
    graph1_peak->Draw("AP");               // Draw with axis and points

    TGraphErrors *graph2_peak = new TGraphErrors(nBinsPt, centerPt, peak[1], 0, peakError[1]);
    graph2_peak->SetName("PeakPosition_vs_Pt_PbGl");
    graph2_peak->SetTitle("PeakPosition_vs_Pt_PbGl");
    graph2_peak->GetXaxis()->SetTitle("Pt");
    graph2_peak->GetYaxis()->SetTitle("Peak Position");
    graph2_peak->SetLineColor(kBlue);      // Set line color to blue
    graph2_peak->SetMarkerColor(kBlue);    // Set marker color to blue
    graph2_peak->SetMarkerStyle(20);       // Set marker style to circle (20)
    graph2_peak->GetYaxis()->SetRangeUser(0.5, 0.6);
    graph2_peak->Draw("P same");           // Draw the second graph on the same plot

    // Add a legend
    TLegend *legend_peak = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_peak->AddEntry(graph1_peak, "Peak Position for PbSc", "p");
    legend_peak->AddEntry(graph2_peak, "Peak Position for PbGl", "p");
    legend_peak->Draw();
    
    graph1_peak->Write();
    graph2_peak->Write();
    // cPeak->SaveAs("Peak_vs_Pt.png");
    // cPeak->SaveAs("Peak_vs_Pt.pdf");


    // Plot Mass Width
    Float_t width[2][nBinsPt];
    Float_t widthError[2][nBinsPt];
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < nBinsPt; j++){
            width[i][j] = drawboth[i][j]->GetParameter(6);
            widthError[i][j] = drawboth[i][j]->GetParError(6);
        }
    }

    TCanvas *cWidth = new TCanvas("cWidth", "Mass Width vs Pt", 1200, 800);
    cWidth->cd();

    TGraphErrors *graph1_width = new TGraphErrors(nBinsPt, centerPt, width[0], 0, widthError[0]);
    graph1_width->SetName("PeakWidth_vs_Pt_PbSc");
    graph1_width->SetTitle("PeakWidth_vs_Pt_PbSc");
    graph1_width->SetLineColor(kRed);       // Set line color to red
    graph1_width->SetMarkerColor(kRed);     // Set marker color to red
    graph1_width->SetMarkerStyle(20);       // Set marker style to circle (20)
    graph1_width->SetTitle(""); 
    graph1_width->GetXaxis()->SetTitle("Pt");
    graph1_width->GetYaxis()->SetTitle("Mass Width");
    graph1_width->GetYaxis()->SetRangeUser(0.0, 0.1);
    graph1_width->Draw("AP");               // Draw with axis and points

    TGraphErrors *graph2_width = new TGraphErrors(nBinsPt, centerPt, width[1], 0, widthError[1]);
    graph2_width->SetName("PeakWidth_vs_Pt_PbGl");
    graph2_width->SetTitle("PeakWidth_vs_Pt_PbGl");
    graph2_width->GetXaxis()->SetTitle("Pt");
    graph2_width->GetYaxis()->SetTitle("Mass Width");
    graph2_width->SetLineColor(kBlue);      // Set line color to blue
    graph2_width->SetMarkerColor(kBlue);    // Set marker color to blue
    graph2_width->SetMarkerStyle(20);       // Set marker style to circle (20)
    graph2_width->GetYaxis()->SetRangeUser(0.0, 0.1);
    graph2_width->Draw("P same");           // Draw the second graph on the same plot

    // Add a legend
    TLegend *legend_width = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_width->AddEntry(graph1_width, "Mass Width for PbSc", "p");
    legend_width->AddEntry(graph2_width, "Mass Width for PbGl", "p");
    legend_width->Draw();
    
    graph1_width->Write();
    graph2_width->Write();
    // cWidth->SaveAs("Width_vs_Pt.png");
    // cWidth->SaveAs("Width_vs_Pt.pdf");


    // Plot Reduced Chisquared
    Float_t chi2[2][nBinsPt];
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < nBinsPt; j++){
            chi2[i][j] = drawboth[i][j]->GetChisquare() / drawboth[i][j]->GetNDF();
        }
    }

    TCanvas *cChi2 = new TCanvas("cChi2", "Mass Chi2 vs Pt", 1200, 800);
    cChi2->cd();

    TGraph *graph1_chi2 = new TGraph(nBinsPt, centerPt, chi2[0]);
    graph1_chi2->SetName("ReducedChi2_vs_Pt_PbSc");
    graph1_chi2->SetTitle("ReducedChi2_vs_Pt_PbSc");
    graph1_chi2->SetLineColor(kRed);       // Set line color to red
    graph1_chi2->SetMarkerColor(kRed);     // Set marker color to red
    graph1_chi2->SetMarkerStyle(20);       // Set marker style to circle (20)
    graph1_chi2->SetTitle(""); 
    graph1_chi2->GetXaxis()->SetTitle("Pt");
    graph1_chi2->GetYaxis()->SetTitle("Chi2");
    graph1_chi2->GetYaxis()->SetRangeUser(0.0, 80);
    graph1_chi2->Draw("AP");               // Draw with axis and points

    TGraph *graph2_chi2 = new TGraph(nBinsPt, centerPt, chi2[1]);
    graph2_chi2->SetName("ReducedChi2_vs_Pt_PbGl");
    graph2_chi2->SetTitle("ReducedChi2_vs_Pt_PbGl");
    graph2_chi2->SetLineColor(kBlue);      // Set line color to blue
    graph2_chi2->SetMarkerColor(kBlue);    // Set marker color to blue
    graph2_chi2->SetMarkerStyle(20);       // Set marker style to circle (20)
    graph2_chi2->GetXaxis()->SetTitle("Pt");
    graph2_chi2->GetYaxis()->SetTitle("Chi2");
    graph2_chi2->GetYaxis()->SetRangeUser(0.0, 80);
    graph2_chi2->Draw("P same");           // Draw the second graph on the same plot

    // Add a legend
    TLegend *legend_chi2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_chi2->AddEntry(graph1_chi2, "Chi2 for PbSc", "p");
    legend_chi2->AddEntry(graph2_chi2, "Chi2 for PbGl", "p");
    legend_chi2->Draw();
    
    graph1_chi2->Write();
    graph2_chi2->Write();
    // cChi2->SaveAs("Chi2_vs_Pt.png");
    // cChi2->SaveAs("Chi2_vs_Pt.pdf");

    filePbSc->Close();
    filePbGl->Close();

    outfile->Write();
    outfile->Close(); 
}

