#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include <math.h>
void analyze_eta_mother_wxBj_cuts_soft(){
	gStyle->SetOptStat(0);
	
	/*
	TFile* f = new TFile("hardqcd_mid.root");
	TFile* f2 = new TFile("hardqcd_mid_cut2.root");
	TFile* f3 = new TFile("hardqcd_mid_cut3.root");
	TFile* f5 = new TFile("hardqcd_mid_cut5.root");
	*/
	
	
	
	TFile* f = new TFile("hardqcd_for.root");
	TFile* f2 = new TFile("hardqcd_for_cut2.root");
	TFile* f3 = new TFile("hardqcd_for_cut3.root");
	TFile* f5 = new TFile("hardqcd_for_cut5.root");
	
	
	
	/*
	TFile* f = new TFile("hardsoftqcd_for.root");
	TFile* f2 = new TFile("hardsoftqcd_mid.root");
	*/

	TTree* ETA = (TTree*) f->Get("eta_tree");
	TTree* ETA2 = (TTree*) f2->Get("eta_tree");
	TTree* ETA3 = (TTree*) f3->Get("eta_tree");
	TTree* ETA5 = (TTree*) f5->Get("eta_tree");
	
	double eta_mother_pdg, eta_part1_xBj, eta_part2_xBj, eta_part1_pdg, eta_part2_pdg;
	ETA->SetBranchAddress("eta_mother_pdg",&eta_mother_pdg);
	ETA->SetBranchAddress("eta_part1_xBj",&eta_part1_xBj);
	ETA->SetBranchAddress("eta_part2_xBj",&eta_part2_xBj);
	ETA->SetBranchAddress("eta_part1_pdg",&eta_part1_pdg);
	ETA->SetBranchAddress("eta_part2_pdg",&eta_part2_pdg);

	double eta2_mother_pdg, eta2_part1_xBj, eta2_part2_xBj, eta2_part1_pdg, eta2_part2_pdg;
	ETA2->SetBranchAddress("eta_mother_pdg",&eta2_mother_pdg);
	ETA2->SetBranchAddress("eta_part1_xBj",&eta2_part1_xBj);
	ETA2->SetBranchAddress("eta_part2_xBj",&eta2_part2_xBj);
	ETA2->SetBranchAddress("eta_part1_pdg",&eta2_part1_pdg);
	ETA2->SetBranchAddress("eta_part2_pdg",&eta2_part2_pdg);
	
	double eta3_mother_pdg, eta3_part1_xBj, eta3_part2_xBj, eta3_part1_pdg, eta3_part2_pdg;
	ETA3->SetBranchAddress("eta_mother_pdg",&eta3_mother_pdg);
	ETA3->SetBranchAddress("eta_part1_xBj",&eta3_part1_xBj);
	ETA3->SetBranchAddress("eta_part2_xBj",&eta3_part2_xBj);
	ETA3->SetBranchAddress("eta_part1_pdg",&eta3_part1_pdg);
	ETA3->SetBranchAddress("eta_part2_pdg",&eta3_part2_pdg);

	double eta5_mother_pdg, eta5_part1_xBj, eta5_part2_xBj, eta5_part1_pdg, eta5_part2_pdg;
	ETA5->SetBranchAddress("eta_mother_pdg",&eta5_mother_pdg);
	ETA5->SetBranchAddress("eta_part1_xBj",&eta5_part1_xBj);
	ETA5->SetBranchAddress("eta_part2_xBj",&eta5_part2_xBj);
	ETA5->SetBranchAddress("eta_part1_pdg",&eta5_part1_pdg);
	ETA5->SetBranchAddress("eta_part2_pdg",&eta5_part2_pdg);
	
	TH1F *G = new TH1F("g-fragmentation", "g-fragmentation", 100, 0., 1.);
	TH1F *Q = new TH1F("q-fragmentation", "q-fragmentation", 100, 0., 1.);
	TH1F *Total = new TH1F("Total_fraction", "Total_fraction", 100, 0., 1.);

	TH1F *G2 = new TH1F("g-fragmentation2", "g-fragmentation2", 100, 0., 1.);
	TH1F *Q2 = new TH1F("q-fragmentation2", "q-fragmentation2", 100, 0., 1.);
	TH1F *Total2 = new TH1F("Total_fraction2", "Total_fraction2", 100, 0., 1.);

	TH1F *G3 = new TH1F("g-fragmentation3", "g-fragmentation3", 100, 0., 1.);
	TH1F *Q3 = new TH1F("q-fragmentation3", "q-fragmentation3", 100, 0., 1.);
	TH1F *Total3 = new TH1F("Total_fraction3", "Total_fraction3", 100, 0., 1.);

	TH1F *G5 = new TH1F("g-fragmentation5", "g-fragmentation5", 100, 0., 1.);
	TH1F *Q5 = new TH1F("q-fragmentation5", "q-fragmentation5", 100, 0., 1.);
	TH1F *Total5 = new TH1F("Total_fraction5", "Total_fraction5", 100, 0., 1.);

	int Eta_Entries = ETA->GetEntries();
	for (int i = 0; i < Eta_Entries; i++){
		if (eta_part1_pdg == 2212 || eta_part1_pdg == 9902210) continue;
		if (eta_part2_pdg == 2212 || eta_part2_pdg == 9902210) continue;
		ETA->GetEntry(i);
		Total->Fill(eta_part1_xBj);
		Total->Fill(eta_part2_xBj);
		if (fabs(eta_mother_pdg)<=6){
			Q->Fill(eta_part1_xBj);
			Q->Fill(eta_part2_xBj);
		}else{
			G->Fill(eta_part1_xBj);
			G->Fill(eta_part2_xBj);
		}
	}
	Q->Divide(Total);
	G->Divide(Total);
	
	int Eta2_Entries = ETA2->GetEntries();
	for (int i = 0; i < Eta2_Entries; i++){
		if (eta2_part1_pdg == 2212 || eta2_part1_pdg == 9902210) continue;
		if (eta2_part2_pdg == 2212 || eta2_part2_pdg == 9902210) continue;
		ETA2->GetEntry(i);
		Total2->Fill(eta2_part1_xBj);
		Total2->Fill(eta2_part2_xBj);
		if (fabs(eta2_mother_pdg)<=6){
			Q2->Fill(eta2_part1_xBj);
			Q2->Fill(eta2_part2_xBj);
		}else{
			G2->Fill(eta2_part1_xBj);
			G2->Fill(eta2_part2_xBj);
		}
	}
	Q2->Divide(Total2);
	G2->Divide(Total2);

	int Eta3_Entries = ETA3->GetEntries();
	for (int i = 0; i < Eta3_Entries; i++){
		if (eta3_part1_pdg == 2212 || eta3_part1_pdg == 9902210) continue;
		if (eta3_part2_pdg == 2212 || eta3_part2_pdg == 9902210) continue;
		ETA3->GetEntry(i);
		Total3->Fill(eta3_part1_xBj);
		Total3->Fill(eta3_part2_xBj);
		if (fabs(eta3_mother_pdg)<=6){
			Q3->Fill(eta3_part1_xBj);
			Q3->Fill(eta3_part2_xBj);
		}else{
			G3->Fill(eta3_part1_xBj);
			G3->Fill(eta3_part2_xBj);
		}
	}
	Q3->Divide(Total3);
	G3->Divide(Total3);

	int Eta5_Entries = ETA5->GetEntries();
	for (int i = 0; i < Eta5_Entries; i++){
		if (eta5_part1_pdg == 2212 || eta5_part1_pdg == 9902210) continue;
		if (eta5_part2_pdg == 2212 || eta5_part2_pdg == 9902210) continue;
		ETA5->GetEntry(i);
		Total5->Fill(eta5_part1_xBj);
		Total5->Fill(eta5_part2_xBj);
		if (fabs(eta5_mother_pdg)<=6){
			Q5->Fill(eta5_part1_xBj);
			Q5->Fill(eta5_part2_xBj);
		}else{
			G5->Fill(eta5_part1_xBj);
			G5->Fill(eta5_part2_xBj);
		}
	}
	Q5->Divide(Total5);
	G5->Divide(Total5);


	G->SetTitle("Final State Mother Particle of #eta Meson at Forward Rapidity");
	G->GetYaxis()->SetTitleOffset(1.);
	G->GetXaxis()->SetTitleOffset(1.);
	G->GetYaxis()->SetTitle("Subprocess Fraction");
	G->GetXaxis()->SetTitle("x_{Bj}");
	G->SetLineColor(kBlue);
	Q->SetLineColor(kBlue);
	G2->SetLineColor(kRed);
	Q2->SetLineColor(kRed);
	G3->SetLineColor(kGreen);
	Q3->SetLineColor(kGreen);
	G5->SetLineColor(kBlack);
	Q5->SetLineColor(kBlack);

	Q->SetLineStyle(2);
	Q2->SetLineStyle(2);
	Q3->SetLineStyle(2);
	Q5->SetLineStyle(2);

	G->Draw();
	Q->Draw("SAME");
	G2->Draw("SAME");
	Q2->Draw("SAME");
	G3->Draw("SAME");
	Q3->Draw("SAME");
	G5->Draw("SAME");
	Q5->Draw("SAME");

	G->GetXaxis()->SetRangeUser(0, 1);
	G->GetYaxis()->SetRangeUser(0, 1);

	
	TLegend *leg = new TLegend(0.2,0.3,0.5,0.7);
	leg->AddEntry(Q, "q","l");
	leg->AddEntry(G, "g","l");
	leg->AddEntry(Q2, "q #hatp_{T} #geq 2 GeV","l");
	leg->AddEntry(G2, "g #hatp_{T} #geq 2 GeV","l");
	leg->AddEntry(Q3, "q #hatp_{T} #geq 3 GeV","l");
	leg->AddEntry(G3, "g #hatp_{T} #geq 3 GeV","l");
	leg->AddEntry(Q5, "q #hatp_{T} #geq 5 GeV","l");
	leg->AddEntry(G5, "g #hatp_{T} #geq 5 GeV","l");

	leg->Draw();
	
}