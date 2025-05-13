void gamma3_trig_eff_summary() {
    TFile *file1 = TFile::Open("trig_eff_minbias_result.root");
    TFile *file2 = TFile::Open("trig_eff_result_mfit_hclust_corrected.root");
    // 1 : MB
    // 2 : EAST & WEST ARM
    TGraphErrors *graphEffSc1 = (TGraphErrors*)file1->Get("graphEffSc");
    TGraphErrors *graphEffGl1 = (TGraphErrors*)file1->Get("graphEffGl");
    TGraphErrors *graphEffSc2 = (TGraphErrors*)file2->Get("graphEffSc");
    TGraphErrors *graphEffGl2 = (TGraphErrors*)file2->Get("graphEffGl");

    TCanvas *ceff = new TCanvas("ceff", "Efficiencies", 800, 600);

    graphEffSc1->SetLineColor(kRed);
    graphEffGl1->SetLineColor(kBlue);
    graphEffSc2->SetLineColor(kRed);
    graphEffGl2->SetLineColor(kBlue);

    graphEffSc1->SetMarkerStyle(24);
    graphEffGl1->SetMarkerStyle(24);
    graphEffSc2->SetMarkerStyle(21);
    graphEffGl2->SetMarkerStyle(21);

    graphEffSc1->SetMarkerColor(kRed);
    graphEffGl1->SetMarkerColor(kBlue);
    graphEffSc2->SetMarkerColor(kRed);
    graphEffGl2->SetMarkerColor(kBlue);

    graphEffSc1->SetLineWidth(2);
    graphEffGl1->SetLineWidth(2);
    graphEffSc2->SetLineWidth(2);
    graphEffGl2->SetLineWidth(2);
    

    TF1 *doubleErfFit_Sc1 = new TF1("sigmoidFit_Sc1", "[0] / (1 + exp(-(x - [1]) / [2]))", 2, 8);
    TF1 *doubleErfFit_Sc2 = new TF1("sigmoidFit_Sc2", "[0] / (1 + exp(-(x - [1]) / [2]))", 2, 15);
    TF1 *doubleErfFit_Gl1 = new TF1("sigmoidFit_Gl1", "[0] / (1 + exp(-(x - [1]) / [2]))", 2, 8);
    TF1 *doubleErfFit_Gl2 = new TF1("sigmoidFit_Gl2", "[0] / (1 + exp(-(x - [1]) / [2]))", 2, 15);

    doubleErfFit_Sc1->SetParameters(1.0, 3.0, 2.0);
    doubleErfFit_Sc2->SetParameters(1.0, 3.0, 2.0);
    doubleErfFit_Gl1->SetParameters(1.0, 3.0, 2.0);
    doubleErfFit_Gl2->SetParameters(1.0, 3.0, 2.0);
    
    graphEffSc1->Fit(doubleErfFit_Sc1, "R+");
    graphEffSc2->Fit(doubleErfFit_Sc2, "R+");
    graphEffGl1->Fit(doubleErfFit_Gl1, "R+");
    graphEffGl2->Fit(doubleErfFit_Gl2, "R+");

    graphEffSc2->SetTitle("");
    graphEffSc2->GetXaxis()->SetTitle("pT");
    graphEffSc2->GetYaxis()->SetTitle("Trigger Efficiency");
    graphEffSc2->GetYaxis()->SetRangeUser(0, 1.35);
    graphEffSc2->GetXaxis()->SetRangeUser(0, 45);

    graphEffSc2->Draw("AP");
    graphEffGl1->Draw("P same");
    graphEffSc1->Draw("P same");
    graphEffGl2->Draw("P same");

    doubleErfFit_Sc1->SetLineColor(kGreen+2);
    doubleErfFit_Sc2->SetLineColor(kGreen);
    doubleErfFit_Gl1->SetLineColor(kOrange+7);
    doubleErfFit_Gl2->SetLineColor(kMagenta);


    doubleErfFit_Sc1->SetLineWidth(3);
    doubleErfFit_Sc2->SetLineWidth(3);
    doubleErfFit_Gl1->SetLineWidth(3);
    doubleErfFit_Gl2->SetLineWidth(3);

    doubleErfFit_Sc1->Draw("same");
    doubleErfFit_Sc2->Draw("same");
    doubleErfFit_Gl1->Draw("same");
    doubleErfFit_Gl2->Draw("same");

    TLegend *legend = new TLegend(0.35, 0.13, 0.65, 0.35);
    legend->AddEntry(graphEffSc1, "MB PbSc", "p");
    legend->AddEntry(graphEffGl1, "MB PbGl", "p");
    legend->AddEntry(graphEffSc2, "E&W Arm PbSc", "p");
    legend->AddEntry(graphEffGl2, "E&W Arm PbGl", "p");
    legend->SetTextSize(0.04);
    legend->Draw();

    ceff->SaveAs("gamma3_Comb_TrigEff.pdf");

    file1->Close();
    file2->Close();
}
