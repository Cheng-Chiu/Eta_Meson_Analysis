#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include <math.h>
void analyze_pion0_initial_subfrac_wxBj(){
	gStyle->SetOptStat(0);
	TFile* f = new TFile("hardqcd_for.root");
	//TFile* f = new TFile("hardqcd_mid.root");
	//TFile* f = new TFile("hardsoftqcd_for.root");
	//TFile* f = new TFile("hardsoftqcd_mid.root");
	double pion0_part1_xBj, pion0_part2_xBj, pion0_part1_pdg, pion0_part2_pdg;
	TTree* PION0 = (TTree*) f->Get("pion0_tree");
	PION0->SetBranchAddress("pion0_part1_xBj",&pion0_part1_xBj);
	PION0->SetBranchAddress("pion0_part2_xBj",&pion0_part2_xBj);
	PION0->SetBranchAddress("pion0_part1_pdg",&pion0_part1_pdg);
	PION0->SetBranchAddress("pion0_part2_pdg",&pion0_part2_pdg);

	TH1F *qq_xbj = new TH1F("qq", "qq", 100, 0., 1.);
	TH1F *gg_xbj = new TH1F("gg", "gg", 100, 0., 1.);
	TH1F *qg_xbj = new TH1F("qg", "qg", 100, 0., 1.);
	TH1F *Total = new TH1F("Total", "Total", 100, 0., 1.);


	int Pion0_Entries = PION0->GetEntries();
	for (Int_t i = 0; i < Pion0_Entries; ++i){
		PION0->GetEntry(i);
		//excluding diffractive events
		if (pion0_part1_pdg == 2212 || pion0_part1_pdg == 9902210) continue;
		if (pion0_part2_pdg == 2212 || pion0_part2_pdg == 9902210) continue;

		Total->Fill(pion0_part1_xBj);
		Total->Fill(pion0_part2_xBj);
		if (fabs(pion0_part1_pdg) <= 6 && abs(pion0_part2_pdg) <= 6){ //qq
			qq_xbj->Fill(pion0_part1_xBj);
			qq_xbj->Fill(pion0_part2_xBj);
		}else if (fabs(pion0_part1_pdg) > 6 && abs(pion0_part2_pdg) > 6){ //gg
			gg_xbj->Fill(pion0_part1_xBj);
			gg_xbj->Fill(pion0_part2_xBj);
		}else{
			qg_xbj->Fill(pion0_part1_xBj);
			qg_xbj->Fill(pion0_part2_xBj);
		}
	}

	qq_xbj->Divide(Total);
	gg_xbj->Divide(Total);
	qg_xbj->Divide(Total);

	qq_xbj->SetTitle("x_{Bj} sorted with collision type, Hardqcd, Forward, Diffractive removed");
	qq_xbj->GetYaxis()->SetTitleOffset(1.);
	qq_xbj->GetXaxis()->SetTitleOffset(1.);
	qq_xbj->GetYaxis()->SetTitle("Subprocess Fraction");
	qq_xbj->GetXaxis()->SetTitle("x_{Bj}");
	qq_xbj->SetLineColor(kBlue);
	gg_xbj->SetLineColor(kRed);
	qg_xbj->SetLineColor(kGreen);
	qq_xbj->Draw();
	gg_xbj->Draw("SAME");
	qg_xbj->Draw("SAME");
	qq_xbj->GetYaxis()->SetRangeUser(0, 1);
	TLegend *leg = new TLegend(0.2,0.65,0.3,0.9);
	leg->AddEntry(qq_xbj, "qq","l");
	leg->AddEntry(gg_xbj, "gg","l");
	leg->AddEntry(qg_xbj, "qg","l");
	leg->Draw();
	//gPad->SetLogx();
}