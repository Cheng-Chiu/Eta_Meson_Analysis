#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include <math.h>
void analyze_pion0_initial_hist_wxBj(){
	gStyle->SetOptStat(0);
	//TFile* f = new TFile("hardqcd_for.root");
	TFile* f = new TFile("hardqcd_mid.root");
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

	int Pion0_Entries = PION0->GetEntries();
	for (Int_t i = 0; i < Pion0_Entries; i++){
		PION0->GetEntry(i);
		//excluding diffractive events
		if (pion0_part1_pdg == 2212 || pion0_part1_pdg == 9902210) continue;
		if (pion0_part2_pdg == 2212 || pion0_part2_pdg == 9902210) continue;

		if (fabs(pion0_part1_pdg) <= 6){
			Q_xbj->Fill(pion0_part1_xBj);
		}else{
			G_xbj->Fill(pion0_part1_xBj);
		}
		if (fabs(pion0_part2_pdg) <= 6){
			Q_xbj->Fill(pion0_part2_xBj);
		}else{
			G_xbj->Fill(pion0_part2_xBj);
		}
	}
	
	Q_xbj->SetTitle("x_{Bj} of Initial State Parton, Hard+softqcd, Mid");
	Q_xbj->GetYaxis()->SetTitleOffset(1.);
	Q_xbj->GetXaxis()->SetTitleOffset(1.);
	Q_xbj->GetYaxis()->SetTitle("Magnitude");
	Q_xbj->GetXaxis()->SetTitle("x_{Bj}");
	Q_xbj->SetLineColor(kBlue);
	G_xbj->SetLineColor(kRed);
	Q_xbj->Draw();
	G_xbj->Draw("SAME");
	Q_xbj->GetYaxis()->SetRangeUser(0., 220000.);
	Q_xbj->GetYaxis()->SetTitleOffset(1.4);
	TLegend *leg = new TLegend(0.2,0.5,0.4,0.75);
	leg->AddEntry(Q_xbj, "Quark","l");
	leg->AddEntry(G_xbj, "Gluon","l");
	leg->Draw();
	gPad->SetLogx();
}