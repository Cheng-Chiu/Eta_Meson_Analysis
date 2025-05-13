#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include <math.h>
void analyze_pion0_mother_wpT(){
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
	
	TTree* PION0 = (TTree*) f->Get("pion0_tree");
	TTree* PION02 = (TTree*) f2->Get("pion0_tree");
	/*
	TTree* PION03 = (TTree*) f3->Get("pion0_tree");
	TTree* PION05 = (TTree*) f5->Get("pion0_tree");
	*/

	double pion0_mother_pdg, pion0_pT;
	PION0->SetBranchAddress("pion0_mother_pdg",&pion0_mother_pdg);
	PION0->SetBranchAddress("pion0_mother_pT",&pion0_pT);
	
	double pion02_mother_pdg, pion02_pT;
	PION02->SetBranchAddress("pion0_mother_pdg",&pion02_mother_pdg);
	PION02->SetBranchAddress("pion0_mother_pT",&pion02_pT);
	/*
	double pion03_mother_pdg, pion03_pT;
	PION03->SetBranchAddress("pion0_mother_pdg",&pion03_mother_pdg);
	PION03->SetBranchAddress("pion0_mother_pT",&pion03_pT);

	double pion05_mother_pdg, pion05_pT;
	PION05->SetBranchAddress("pion0_mother_pdg",&pion05_mother_pdg);
	PION05->SetBranchAddress("pion0_mother_pT",&pion05_pT);
	*/
	TH1F *G = new TH1F("g-fragmentation", "g-fragmentation", 100, 0., 10.);
	TH1F *Q = new TH1F("q-fragmentation", "q-fragmentation", 100, 0., 10.);
	TH1F *Total = new TH1F("Total_fraction", "Total_fraction", 100, 0., 10.);
	
	TH1F *G2 = new TH1F("g-fragmentation2", "g-fragmentation2", 100, 0., 10.);
	TH1F *Q2 = new TH1F("q-fragmentation2", "q-fragmentation2", 100, 0., 10.);
	TH1F *Total2 = new TH1F("Total_fraction2", "Total_fraction2", 100, 0., 10.);
	/*
	TH1F *G3 = new TH1F("g-fragmentation3", "g-fragmentation3", 100, 0., 10.);
	TH1F *Q3 = new TH1F("q-fragmentation3", "q-fragmentation3", 100, 0., 10.);
	TH1F *Total3 = new TH1F("Total_fraction3", "Total_fraction3", 100, 0., 10.);

	TH1F *G5 = new TH1F("g-fragmentation5", "g-fragmentation5", 100, 0., 10.);
	TH1F *Q5 = new TH1F("q-fragmentation5", "q-fragmentation5", 100, 0., 10.);
	TH1F *Total5 = new TH1F("Total_fraction5", "Total_fraction5", 100, 0., 10.);
	*/
	int Pion0_Entries = PION0->GetEntries();
	for (int i = 0; i < Pion0_Entries; i++){
		PION0->GetEntry(i);
		Total->Fill(pion0_pT);
		if (fabs(pion0_mother_pdg)<=6){
			Q->Fill(pion0_pT);
		}else{
			G->Fill(pion0_pT);
		}
	}
	Q->Divide(Total);
	G->Divide(Total);
	
	int Pion02_Entries = PION02->GetEntries();
	for (int i = 0; i < Pion02_Entries; i++){
		PION02->GetEntry(i);
		Total2->Fill(pion02_pT);
		if (fabs(pion02_mother_pdg)<=6){
			Q2->Fill(pion02_pT);
		}else{
			G2->Fill(pion02_pT);
		}
	}
	Q2->Divide(Total2);
	G2->Divide(Total2);
	/*
	int Pion03_Entries = PION03->GetEntries();
	for (int i = 0; i < Pion03_Entries; i++){
		PION03->GetEntry(i);
		Total3->Fill(pion03_pT);
		if (fabs(pion03_mother_pdg)<=6){
			Q3->Fill(pion03_pT);
		}else{
			G3->Fill(pion03_pT);
		}
	}
	Q3->Divide(Total3);
	G3->Divide(Total3);

	int Pion05_Entries = PION05->GetEntries();
	for (int i = 0; i < Pion05_Entries; i++){
		PION05->GetEntry(i);
		Total5->Fill(pion05_pT);
		if (fabs(pion05_mother_pdg)<=6){
			Q5->Fill(pion05_pT);
		}else{
			G5->Fill(pion05_pT);
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