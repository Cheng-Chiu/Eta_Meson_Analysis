#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include <math.h>
void analyze_eta_initial_hist_wxBj(){
	gStyle->SetOptStat(0);
	TFile* f = new TFile("hardqcd_for_15.root");
	//TFile* f = new TFile("hardqcd_for_15.root");
	//TFile* f = new TFile("softqcd_for_15..root");
	//TFile* f = new TFile("softqcd_for_15..root");
	double eta_part1_xBj, eta_part2_xBj, eta_part1_pdg, eta_part2_pdg;
	TTree* ETA = (TTree*) f->Get("eta_tree");
	ETA->SetBranchAddress("eta_part1_xBj",&eta_part1_xBj);
	ETA->SetBranchAddress("eta_part2_xBj",&eta_part2_xBj);
	ETA->SetBranchAddress("eta_part1_pdg",&eta_part1_pdg);
	ETA->SetBranchAddress("eta_part2_pdg",&eta_part2_pdg);

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

	int Eta_Entries = ETA->GetEntries();
	for (Int_t i = 0; i < Eta_Entries; i++){
		ETA->GetEntry(i);
		//excluding diffractive events
		if (eta_part1_pdg == 2212 || eta_part1_pdg == 9902210) continue;
		if (eta_part2_pdg == 2212 || eta_part2_pdg == 9902210) continue;

		if (fabs(eta_part1_pdg) <= 6){
			Q_xbj->Fill(eta_part1_xBj);
		}else{
			G_xbj->Fill(eta_part1_xBj);
		}
		if (fabs(eta_part2_pdg) <= 6){
			Q_xbj->Fill(eta_part2_xBj);
		}else{
			G_xbj->Fill(eta_part2_xBj);
		}
	}
	Q_xbj->SetTitle("x_{Bj} of Initial State Parton, Hard+softqcd, Forward");
	Q_xbj->GetYaxis()->SetTitleOffset(1.);
	Q_xbj->GetXaxis()->SetTitleOffset(1.);
	Q_xbj->GetYaxis()->SetTitle("Magnitude");
	Q_xbj->GetXaxis()->SetTitle("x_{Bj}");
	Q_xbj->SetLineColor(kBlue);
	G_xbj->SetLineColor(kRed);
	Q_xbj->Draw();
	G_xbj->Draw("SAME");
	Q_xbj->GetYaxis()->SetRangeUser(0., 5.5E4);
	Q_xbj->GetYaxis()->SetTitleOffset(1.4);
	TLegend *leg = new TLegend(0.2,0.5,0.4,0.75);
	leg->AddEntry(Q_xbj, "Quark","l");
	leg->AddEntry(G_xbj, "Gluon","l");
	leg->Draw();
	gPad->SetLogx();
}