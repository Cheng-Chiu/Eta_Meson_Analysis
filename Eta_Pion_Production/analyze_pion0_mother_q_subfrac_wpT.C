#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include <math.h>
void analyze_pion0_mother_q_subfrac_wpT(){
	gStyle->SetOptStat(0);
	TFile* f = new TFile("hardqcd_for_15.root");
	//TFile* f = new TFile("softqcd_for_15.root");


	TTree* PION0 = (TTree*) f->Get("pion0_tree");

	double pion0_mother_pdg, pion0_pT, pion0_part1_pdg, pion0_part2_pdg;
	PION0->SetBranchAddress("pion0_mother_pdg",&pion0_mother_pdg);
	PION0->SetBranchAddress("pion0_pT",&pion0_pT);
	PION0->SetBranchAddress("pion0_part1_pdg",&pion0_part1_pdg);
	PION0->SetBranchAddress("pion0_part2_pdg",&pion0_part2_pdg);

	TH1F *G = new TH1F("g-fragmentation", "", 25, 0., 6.5);
	TH1F *D = new TH1F("d-fragmentation", "", 25, 0., 6.5);
	TH1F *U = new TH1F("u-fragmentation", "", 25, 0., 6.5);
	TH1F *S = new TH1F("s-fragmentation", "", 25, 0., 6.5);

	TH1F *aD = new TH1F("#bar{d}-fragmentation", "", 25, 0., 6.5);
	TH1F *aU = new TH1F("#bar{u}-fragmentation", "", 25, 0., 6.5);
	TH1F *aS = new TH1F("#bar{s}-fragmentation", "", 25, 0., 6.5);

	TH1F *Total = new TH1F("Total_fraction", "Total_fraction", 25, 0., 6.5);
	G->Sumw2();
	D->Sumw2();
	U->Sumw2();
	S->Sumw2();
	Total->Sumw2();

	
	int PION0_Entries = PION0->GetEntries();
	for (int i = 0; i < PION0_Entries; i++){
		PION0->GetEntry(i);
		//excluding diffractive events
		if (pion0_part1_pdg == 2212 || pion0_part1_pdg == 9902210) continue;
		if (pion0_part2_pdg == 2212 || pion0_part2_pdg == 9902210) continue;

		Total->Fill(pion0_pT);
		if (pion0_mother_pdg == 21){
			G->Fill(pion0_pT);
		}else if(pion0_mother_pdg == 1){
			D->Fill(pion0_pT);
		}else if(pion0_mother_pdg == 2){
			U->Fill(pion0_pT);
		}else if(pion0_mother_pdg == 3){
			S->Fill(pion0_pT);
		}else if (pion0_mother_pdg == 21){
			G->Fill(pion0_pT);
		}
	}
	
	G->Divide(Total);
	D->Divide(Total);
	U->Divide(Total);
	S->Divide(Total);
	

	TCanvas *c1 = new TCanvas("c1","");

	G->GetYaxis()->SetTitle("Subprocess Fraction");
	G->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	G->GetYaxis()->SetRangeUser(0, 1);
	G->GetXaxis()->SetRangeUser(0, 6.5);
	G->GetXaxis()->SetTitleSize(0.047);
	G->GetXaxis()->SetTitleOffset(0.9);
	G->GetYaxis()->SetTitleSize(0.047);
	G->GetYaxis()->SetTitleOffset(0.9);

	G->SetLineColor(1);
	G->SetMarkerStyle(21);
	G->SetMarkerColor(1);
	G->SetMarkerSize(1.1);

	D->SetLineColor(2);
	D->SetMarkerStyle(21);
	D->SetMarkerColor(2);
	D->SetMarkerSize(1.1);

	U->SetLineColor(3);
	U->SetMarkerStyle(21);
	U->SetMarkerColor(3);
	U->SetMarkerSize(1.1);

	S->SetLineColor(880);
	S->SetMarkerStyle(21);
	S->SetMarkerColor(880);
	S->SetMarkerSize(1.1);
	

	G->Draw("P");
	D->Draw("PSAME");
	U->Draw("PSAME");
	S->Draw("PSAME");


	TLatex latex;
	latex.SetTextSize(0.045);
	latex.SetTextAlign(13);
	latex.DrawLatex(.3,.95,"p + p #rightarrow #pi^{0} + X (#sqrt{s} = 500 GeV)");
	latex.DrawLatex(.3,.85,"3.0 < |#eta| < 3.8");
	//latex.DrawLatex(.3,.85,"|#eta| < 1");

	latex.DrawLatex(.5,.43,"#color[2]{d-fragmentation}");
	latex.DrawLatex(.5,.65,"#color[3]{u-fragmentation}");
	latex.DrawLatex(.5,.18,"#color[1]{g-fragmentation}");
	latex.DrawLatex(5.25,.05,"#color[880]{s-fragmentation}");
	//latex.DrawLatex(2.3,.05,"#color[880]{s-fragmentation}");
	c1->SaveAs("pion0_mother_wpT_h_f_15.pdf");
}