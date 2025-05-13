#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include <math.h>
void analyze_eta_initial_subfrac_wxT(){
	gStyle->SetOptStat(0);
	TFile* f = new TFile("hardqcd_mid34.root");
	//TFile* f = new TFile("hardqcd_for.root");
	//TFile* f = new TFile("hardqcd_mid.root");
	//TFile* f = new TFile("hardsoftqcd_for.root");
	//TFile* f = new TFile("hardsoftqcd_mid.root");
	double eta_pT, eta_part3_pdg, eta_part4_pdg;
	TTree* ETA = (TTree*) f->Get("eta_tree");
	ETA->SetBranchAddress("eta_pT",&eta_pT);
	ETA->SetBranchAddress("eta_part3_pdg",&eta_part3_pdg);
	ETA->SetBranchAddress("eta_part4_pdg",&eta_part4_pdg);

	TH1F *qq_xT = new TH1F("qq", "qq", 100, 0., 6.5);
	TH1F *gg_xT = new TH1F("gg", "gg", 100, 0., 6.5);
	TH1F *qg_xT = new TH1F("qg", "qg", 100, 0., 6.5);
	TH1F *Total = new TH1F("Total", "Total", 100, 0., 6.5);


	int Eta_Entries = ETA->GetEntries();
	for (Int_t i = 0; i < Eta_Entries; ++i){
		ETA->GetEntry(i);
		Total->Fill(eta_pT);
		if (fabs(eta_part3_pdg) <= 6 && abs(eta_part4_pdg) <= 6){ //qq
			qq_xT->Fill(eta_pT);
		}else if (fabs(eta_part3_pdg) > 6 && abs(eta_part4_pdg) > 6){ //gg
			gg_xT->Fill(eta_pT);
		}else{
			qg_xT->Fill(eta_pT);
		}
	}

	qq_xT->Divide(Total);
	gg_xT->Divide(Total);
	qg_xT->Divide(Total);

	qq_xT->SetTitle("p_{T} sorted with collision type, Hardqcd, Mid");
	qq_xT->GetYaxis()->SetTitleOffset(1.);
	qq_xT->GetXaxis()->SetTitleOffset(1.);
	qq_xT->GetYaxis()->SetTitle("Subprocess Fraction");
	qq_xT->GetXaxis()->SetTitle("p_{T}");
	qq_xT->SetLineColor(kBlue);
	gg_xT->SetLineColor(kRed);
	qg_xT->SetLineColor(kGreen);
	qq_xT->Draw();
	gg_xT->Draw("SAME");
	qg_xT->Draw("SAME");
	qq_xT->GetYaxis()->SetRangeUser(0, 1);
	TLegend *leg = new TLegend(0.85,0.65,1,0.9);
	leg->AddEntry(qq_xT, "qq","l");
	leg->AddEntry(gg_xT, "gg","l");
	leg->AddEntry(qg_xT, "qg","l");
	leg->Draw();
	//gPad->SetLogx();
}