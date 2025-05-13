#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include <math.h>
void analyze_eta_initial_subfrac_wpT(){	
	gStyle->SetOptStat(0);
	TFile* f = new TFile("hardqcd_for_15.root");
	double eta_pT, eta_part1_pdg, eta_part2_pdg;
	TTree* ETA = (TTree*) f->Get("eta_tree");
	ETA->SetBranchAddress("eta_pT",&eta_pT);
	ETA->SetBranchAddress("eta_part1_pdg",&eta_part1_pdg);
	ETA->SetBranchAddress("eta_part2_pdg",&eta_part2_pdg);

	TH1F *qq_pT = new TH1F("qq", "", 50, 0., 4.);
	TH1F *gg_pT = new TH1F("gg", "", 50, 0., 4.);
	TH1F *qg_pT = new TH1F("qg", "", 50, 0., 4.);
	TH1F *Total = new TH1F("Total", "Total", 50, 0., 4.);
	qq_pT->Sumw2();
	gg_pT->Sumw2();
	qg_pT->Sumw2();
	Total->Sumw2();


	int Eta_Entries = ETA->GetEntries();
	for (Int_t i = 0; i < Eta_Entries; ++i){
		ETA->GetEntry(i);
		//excluding diffractive events
		if (eta_part1_pdg == 2212 || eta_part1_pdg == 9902210) continue;
		if (eta_part2_pdg == 2212 || eta_part2_pdg == 9902210) continue;

		Total->Fill(eta_pT);
		if (fabs(eta_part1_pdg) <= 6 && fabs(eta_part2_pdg) <= 6){ //qq
			qq_pT->Fill(eta_pT);
		}else if (fabs(eta_part1_pdg) > 6 && fabs(eta_part2_pdg) > 6){ //gg
			gg_pT->Fill(eta_pT);
		}else{ //qg
			qg_pT->Fill(eta_pT);
		}
	}
	
	qq_pT->Divide(Total);
	gg_pT->Divide(Total);
	qg_pT->Divide(Total);
	

	TCanvas *c1 = new TCanvas("c1","");

	qq_pT->SetLineColor(4);
	qq_pT->SetMarkerStyle(21);
	qq_pT->SetMarkerColor(4);
	qq_pT->SetMarkerSize(1.1);
	gg_pT->SetLineColor(2);
	gg_pT->SetMarkerStyle(21);
	gg_pT->SetMarkerColor(2);
	gg_pT->SetMarkerSize(1.1);
	qg_pT->SetLineColor(3);
	qg_pT->SetMarkerStyle(21);
	qg_pT->SetMarkerColor(3);
	qg_pT->SetMarkerSize(1.1);
	
	qq_pT->GetYaxis()->SetTitle("Subprocess Fraction");
	qq_pT->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	qq_pT->GetYaxis()->SetRangeUser(0, 1);
	qq_pT->GetXaxis()->SetRangeUser(0, 4);
	qq_pT->GetXaxis()->SetTitleSize(0.047);
	qq_pT->GetXaxis()->SetTitleOffset(0.9);
	qq_pT->GetYaxis()->SetTitleSize(0.047);
	qq_pT->GetYaxis()->SetTitleOffset(0.9);

	qq_pT->Draw("P");
	gg_pT->Draw("PSAME");
	qg_pT->Draw("PSAME");
	
	TLatex latex;
	latex.SetTextSize(0.045);
	latex.SetTextAlign(13);
	latex.DrawLatex(.2,.95,"p + p #rightarrow #eta + X (#sqrt{s} = 500 GeV)");
	latex.DrawLatex(.2,.85,"3.0 < |#eta| < 3.8");
	//latex.DrawLatex(.2,.85,"|#eta| < 1");

	latex.DrawLatex(.5,.13,"#color[4]{qq}");
	latex.DrawLatex(.5,.63,"#color[2]{gg}");
	latex.DrawLatex(.5,.45,"#color[3]{qg}");
	c1->SaveAs("eta_initial_wpT_colitype_h_f_15.pdf");
	
}