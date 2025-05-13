#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include <math.h>
void analyze_pion0_initial_q_hist_wxBj(){
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

	//set log bins
	const int nbins = 100;
	Double_t xbins[nbins + 1];
	double dx = 1./nbins;
	double L10 = TMath::Log(1E5);
	double xmin = 1E-5;
	for (int i = 0; i <= nbins; i++) {
		xbins[i] = xmin*TMath::Exp(L10 * i * dx);
	}

	TH1F *G_xbj = new TH1F("Gluon_xBj", "Gluon_xBj", nbins, xbins);
	TH1F *Q_xbj = new TH1F("Quark_xBj", "Quark_xBj", nbins, xbins);
	TH1F *D = new TH1F("Down", "Down", nbins, xbins);
	TH1F *U = new TH1F("Up", "Up", nbins, xbins);
	TH1F *S = new TH1F("Strange", "Strange", nbins, xbins);
	TH1F *C = new TH1F("Charm", "Charm", nbins, xbins);
	TH1F *B = new TH1F("Bottom", "Bottom", nbins, xbins);
	TH1F *T = new TH1F("Top", "Top", nbins, xbins);

	TH1F *_D = new TH1F("Anti-Down", "Anti-Down", nbins, xbins);
	TH1F *_U = new TH1F("Anti-Up", "Anti-Up", nbins, xbins);
	TH1F *_S = new TH1F("Anti-Strange", "Anti-Strange", nbins, xbins);
	TH1F *_C = new TH1F("Anti-Charm", "Anti-Charm", nbins, xbins);
	TH1F *_B = new TH1F("Anti-Bottom", "Anti-Bottom", nbins, xbins);
	TH1F *_T = new TH1F("Anti-Top", "Anti-Top", nbins, xbins);


	int Pion0_Entries = PION0->GetEntries();
	for (Int_t i = 0; i < Pion0_Entries; i++){
		PION0->GetEntry(i);
		//excluding diffractive events
		if (pion0_part1_pdg == 2212 || pion0_part1_pdg == 9902210) continue;
		if (pion0_part2_pdg == 2212 || pion0_part2_pdg == 9902210) continue;

		if (fabs(pion0_part1_pdg) <= 6){
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
	T->SetTitle("x_{Bj} from Quarks, Hardqcd, Forward, Diffractive removed");
	T->GetYaxis()->SetTitleOffset(1.);
	T->GetXaxis()->SetTitleOffset(1.);
	T->GetYaxis()->SetTitle("Magnitude");
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
	
	TLegend *leg = new TLegend(0.2,0.45,0.5,0.9);
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

	T->GetYaxis()->SetRangeUser(0., 5.5E4);
	T->GetXaxis()->SetRangeUser(1E-2, 1);
	T->GetYaxis()->SetTitleOffset(1.4);
	gPad->SetLogx();
}