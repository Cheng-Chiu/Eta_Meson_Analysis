void plot_together(){
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    
    TFile* file1 = TFile::Open("trig_eff_result.root"); 
    TFile* file2 = TFile::Open("trig_eff_minbias_result.root");

    TGraphErrors* graph1 = (TGraphErrors*)file1->Get("ceff");
    TGraphErrors* graph2 = (TGraphErrors*)file2->Get("ceff");

    TCanvas* c1 = new TCanvas("c1", "Two Efficiencies", 800, 600); 
    graph1->Draw("APE");
    graph2->Draw("PE SAME");

    cout << "Program ends" << endl;
}
