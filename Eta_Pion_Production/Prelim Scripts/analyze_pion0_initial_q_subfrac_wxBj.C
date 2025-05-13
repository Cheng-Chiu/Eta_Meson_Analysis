#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include <math.h>
void analyze_pion0_initial_q_subfrac_wxBj(){
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

	TH1F *D = new TH1F("Down", "Down", 100, 0., 1.);
	TH1F *U = new TH1F("Up", "Up", 100, 0., 1.);
	TH1F *S = new TH1F("Strange", "Strange", 100, 0., 1.);
	TH1F *C = new TH1F("Charm", "Charm", 100, 0., 1.);
	TH1F *B = new TH1F("Bottom", "Bottom", 100, 0., 1.);
	TH1F *T = new TH1F("Top", "Top", 100, 0., 1.);

	TH1F *_D = new TH1F("Anti-Down", "Anti-Down", 100, 0., 1.);
	TH1F *_U = new TH1F("Anti-Up", "Anti-Up", 100, 0., 1.);
	TH1F *_S = new TH1F("Anti-Strange", "Anti-Strange", 100, 0., 1.);
	TH1F *_C = new TH1F("Anti-Charm", "Anti-Charm", 100, 0., 1.);
	TH1F *_B = new TH1F("Anti-Bottom", "Anti-Bottom", 100, 0., 1.);
	TH1F *_T = new TH1F("Anti-Top", "Anti-Top", 100, 0., 1.);

	TH1F *Total = new TH1F("Total", "Total", 100, 0., 1.);


	int Pion0_Entries = PION0->GetEntries();
	for (Int_t i = 0; i < Pion0_Entries; i++){
		PION0->GetEntry(i);
		//excluding diffractive events
		if (pion0_part1_pdg == 2212 || pion0_part1_pdg == 9902210) continue;
		if (pion0_part2_pdg == 2212 || pion0_part2_pdg == 9902210) continue;

		if (fabs(pion0_part1_pdg) <= 6){
			Total->Fill(pion0_part1_xBj);
			switch((int)pion0_part1_pdg){
				case 1: 
					D->Fill(pion0_part1_xBj);
					break;
				case 2: 
					U->Fill(pion0_part1_xBj);
					break;
				case 3: 
					S->Fill(pion0_part1_xBj);
					break;
				case 4: 
					C->Fill(pion0_part1_xBj);
					break;
				case 5: 
					B->Fill(pion0_part1_xBj);
					break;
				case 6: 
					T->Fill(pion0_part1_xBj);
					break;
				case -1: 
					_D->Fill(pion0_part1_xBj);
					break;
				case -2: 
					_U->Fill(pion0_part1_xBj);
					break;
				case -3: 
					_S->Fill(pion0_part1_xBj);
					break;
				case -4: 
					_C->Fill(pion0_part1_xBj);
					break;
				case -5: 
					_B->Fill(pion0_part1_xBj);
					break;
				case -6: 
					_T->Fill(pion0_part1_xBj);
					break;
			}
		}

		if (fabs(pion0_part2_pdg) <= 6){
			Total->Fill(pion0_part2_xBj);
			switch((int)pion0_part2_pdg){
				case 1: 
					D->Fill(pion0_part2_xBj);
					break;
				case 2: 
					U->Fill(pion0_part2_xBj);
					break;
				case 3: 
					S->Fill(pion0_part2_xBj);
					break;
				case 4: 
					C->Fill(pion0_part2_xBj);
					break;
				case 5: 
					B->Fill(pion0_part2_xBj);
					break;
				case 6: 
					T->Fill(pion0_part2_xBj);
					break;
				case -1: 
					_D->Fill(pion0_part2_xBj);
					break;
				case -2: 
					_U->Fill(pion0_part2_xBj);
					break;
				case -3: 
					_S->Fill(pion0_part2_xBj);
					break;
				case -4: 
					_C->Fill(pion0_part2_xBj);
					break;
				case -5: 
					_B->Fill(pion0_part2_xBj);
					break;
				case -6: 
					_T->Fill(pion0_part2_xBj);
					break;
			}
		}
	}

	D->Divide(Total);
	U->Divide(Total);
	S->Divide(Total);
	C->Divide(Total);
	B->Divide(Total);
	T->Divide(Total);

	_D->Divide(Total);
	_U->Divide(Total);
	_S->Divide(Total);
	_C->Divide(Total);
	_B->Divide(Total);
	_T->Divide(Total);


	T->SetTitle("x_{Bj} sorted with (only) flavors of Quarks, Hardqcd, Forward, Diffractive removed");
	T->GetYaxis()->SetTitleOffset(1.);
	T->GetXaxis()->SetTitleOffset(1.);
	T->GetYaxis()->SetTitle("Subprocess Fraction of Quarks");
	T->GetXaxis()->SetTitle("x_{Bj}");
	D->SetLineColor(1);
	U->SetLineColor(2);
	S->SetLineColor(3);
	C->SetLineColor(4);
	B->SetLineColor(40);
	T->SetLineColor(7);
	_D->SetLineColor(1);
	_U->SetLineColor(2);
	_S->SetLineColor(3);
	_C->SetLineColor(4);
	_B->SetLineColor(40);
	_T->SetLineColor(7);
	_D->SetLineStyle(2);
	_U->SetLineStyle(2);
	_S->SetLineStyle(2);
	_C->SetLineStyle(2);
	_B->SetLineStyle(2);
	_T->SetLineStyle(2);
	T->Draw();
	D->Draw("SAME");
	U->Draw("SAME");
	S->Draw("SAME");
	C->Draw("SAME");
	B->Draw("SAME");
	_D->Draw("SAME");
	_U->Draw("SAME");
	_S->Draw("SAME");
	_C->Draw("SAME");
	_B->Draw("SAME");
	_T->Draw("SAME");

	T->GetYaxis()->SetRangeUser(0, 1);
	
	TLegend *leg = new TLegend(0.2,0.65,0.3,0.9);
	leg->AddEntry(D, "Down", "l");
	leg->AddEntry(U, "Up", "l");
	leg->AddEntry(S, "Strange", "l");
	leg->AddEntry(C, "Charm", "l");
	leg->AddEntry(B, "Bottom", "l");
	leg->AddEntry(T, "Top", "l");

	leg->AddEntry(_D, "Anti-Down", "l");
	leg->AddEntry(_U, "Anti-Up", "l");
	leg->AddEntry(_S, "Anti-Strange", "l");
	leg->AddEntry(_C, "Anti-Charm", "l");
	leg->AddEntry(_B, "Anti-Bottom", "l");
	leg->AddEntry(_T, "Anti-Top", "l");
	leg->Draw();
}