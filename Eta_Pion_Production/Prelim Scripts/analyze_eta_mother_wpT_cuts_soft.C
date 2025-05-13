#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include <math.h>
void analyze_eta_mother_wpT_cuts_soft(){
	gStyle->SetOptStat(0);
	
	/*
	TFile* f = new TFile("hardqcd_mid.root");
	TFile* f2 = new TFile("hardqcd_mid_cut2.root");
	TFile* f3 = new TFile("hardqcd_mid_cut3.root");
	TFile* f5 = new TFile("hardqcd_mid_cut5.root");
	*/
	
	
	/*
	TFile* f = new TFile("hardqcd_for.root");
	TFile* f2 = new TFile("hardqcd_for_cut2.root");
	TFile* f3 = new TFile("hardqcd_for_cut3.root");
	TFile* f5 = new TFile("hardqcd_for_cut5.root");
	*/
	

	TFile* f = new TFile("hardsoftqcd_for.root");
	TFile* f2 = new TFile("hardsoftqcd_mid.root");

	TTree* ETA = (TTree*) f->Get("eta_tree");
	TTree* ETA2 = (TTree*) f2->Get("eta_tree");
	/*
	TTree* ETA3 = (TTree*) f3->Get("eta_tree");
	TTree* ETA5 = (TTree*) f5->Get("eta_tree");
	*/
	double eta_mother_pdg, eta_pT;
	ETA->SetBranchAddress("eta_mother_pdg",&eta_mother_pdg);
	ETA->SetBranchAddress("eta_mother_pT",&eta_pT);

	double eta2_mother_pdg, eta2_pT;
	ETA2->SetBranchAddress("eta_mother_pdg",&eta2_mother_pdg);
	ETA2->SetBranchAddress("eta_mother_pT",&eta2_pT);
	/*
	double eta3_mother_pdg, eta3_pT;
	ETA3->SetBranchAddress("eta_mother_pdg",&eta3_mother_pdg);
	ETA3->SetBranchAddress("eta_mother_pT",&eta3_pT);

	double eta5_mother_pdg, eta5_pT;
	ETA5->SetBranchAddress("eta_mother_pdg",&eta5_mother_pdg);
	ETA5->SetBranchAddress("eta_mother_pT",&eta5_pT);
	*/
	TH1F *G = new TH1F("g-fragmentation", "g-fragmentation", 50, 0., 7.);
	TH1F *Q = new TH1F("q-fragmentation", "q-fragmentation", 50, 0., 7.);
	TH1F *Total = new TH1F("Total_fraction", "Total_fraction", 50, 0., 7.);

	TH1F *G2 = new TH1F("g-fragmentation2", "g-fragmentation2", 50, 0., 7.);
	TH1F *Q2 = new TH1F("q-fragmentation2", "q-fragmentation2", 50, 0., 7.);
	TH1F *Total2 = new TH1F("Total_fraction2", "Total_fraction2", 50, 0., 7.);
	/*
	TH1F *G3 = new TH1F("g-fragmentation3", "g-fragmentation3", 100, 0., 10.);
	TH1F *Q3 = new TH1F("q-fragmentation3", "q-fragmentation3", 100, 0., 10.);
	TH1F *Total3 = new TH1F("Total_fraction3", "Total_fraction3", 100, 0., 10.);

	TH1F *G5 = new TH1F("g-fragmentation5", "g-fragmentation5", 100, 0., 10.);
	TH1F *Q5 = new TH1F("q-fragmentation5", "q-fragmentation5", 100, 0., 10.);
	TH1F *Total5 = new TH1F("Total_fraction5", "Total_fraction5", 100, 0., 10.);
	*/
	int Eta_Entries = ETA->GetEntries();
	for (int i = 0; i < Eta_Entries; i++){
		ETA->GetEntry(i);
		Total->Fill(eta_pT);
		if (fabs(eta_mother_pdg)<=6){
			Q->Fill(eta_pT);
		}else{
			G->Fill(eta_pT);
		}
	}
	Q->Divide(Total);
	G->Divide(Total);
	
	int Eta2_Entries = ETA2->GetEntries();
	for (int i = 0; i < Eta2_Entries; i++){
		ETA2->GetEntry(i);
		Total2->Fill(eta2_pT);
		if (fabs(eta2_mother_pdg)<=6){
			Q2->Fill(eta2_pT);
		}else{
			G2->Fill(eta2_pT);
		}
	}
	Q2->Divide(Total2);
	G2->Divide(Total2);
	/*
	int Eta3_Entries = ETA3->GetEntries();
	for (int i = 0; i < Eta3_Entries; i++){
		ETA3->GetEntry(i);
		Total3->Fill(eta3_pT);
		if (fabs(eta3_mother_pdg)<=6){
			Q3->Fill(eta3_pT);
		}else{
			G3->Fill(eta3_pT);
		}
	}
	Q3->Divide(Total3);
	G3->Divide(Total3);

	int Eta5_Entries = ETA5->GetEntries();
	for (int i = 0; i < Eta5_Entries; i++){
		ETA5->GetEntry(i);
		Total5->Fill(eta5_pT);
		if (fabs(eta5_mother_pdg)<=6){
			Q5->Fill(eta5_pT);
		}else{
			G5->Fill(eta5_pT);
		}
	}
	Q5->Divide(Total5);
	G5->Divide(Total5);
	*/
	G->SetTitle("Hard+softqcd");
	G->GetYaxis()->SetTitleOffset(1.);
	G->GetXaxis()->SetTitleOffset(1.);
	G->GetYaxis()->SetTitle("Subprocess Fraction");
	G->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	G->SetLineColor(kBlue);
	Q->SetLineColor(kBlue);
	G2->SetLineColor(kRed);
	Q2->SetLineColor(kRed);
	/*
	G3->SetLineColor(kGreen);
	Q3->SetLineColor(kGreen);
	G5->SetLineColor(kBlack);
	Q5->SetLineColor(kBlack);
	*/
	G->Draw();
	Q->Draw("SAME");
	G2->Draw("SAME");
	Q2->Draw("SAME");
	/*
	G3->Draw("SAME");
	Q3->Draw("SAME");
	G5->Draw("SAME");
	Q5->Draw("SAME");
	*/
	G->GetYaxis()->SetRangeUser(0, 1);

	
	TLegend *leg = new TLegend(0.2,0.4,0.5,0.6);
	leg->AddEntry(Q, "q-forward","l");
	leg->AddEntry(G, "g-forward","l");
	leg->AddEntry(Q2, "q-mid","l");
	leg->AddEntry(G2, "g-mid","l");
	/*
	leg->AddEntry(Q3, "q #hatp_{T} #geq 3 GeV","l");
	leg->AddEntry(G3, "g #hatp_{T} #geq 3 GeV","l");
	leg->AddEntry(Q5, "q #hatp_{T} #geq 5 GeV","l");
	leg->AddEntry(G5, "g #hatp_{T} #geq 5 GeV","l");
	*/
	leg->Draw();
	
}