#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include <math.h>
void analyze_eta_momentum_fraction(){
	gStyle->SetOptStat(0);
	TFile* f = new TFile("hardqcd_mid.root");
	double eta_pT, eta_mother_pdg, eta_part1_xBj, eta_part2_xBj;
	TTree* ETA = (TTree*) f->Get("eta_tree");
	ETA->SetBranchAddress("eta_part1_xBj",&eta_part1_xBj);
	ETA->SetBranchAddress("eta_part2_xBj",&eta_part2_xBj);
	ETA->SetBranchAddress("eta_mother_pdg",&eta_mother_pdg);
	ETA->SetBranchAddress("eta_pT",&eta_pT);

	TH1F *Gluon_fraction = new TH1F("Gluon_fraction", "Gluon_fraction", 100, 0., 7.);
	TH1F *Quark_fraction = new TH1F("Quark_fraction", "Quark_fraction", 100, 0., 7.);
	TH1F *Total_fraction = new TH1F("Total_fraction", "Total_fraction", 100, 0., 7.);

	Int_t Eta_Entries = ETA->GetEntries();
	for (Int_t i = 0; i < Eta_Entries; i++){
		ETA->GetEntry(i);
		Total_fraction->Fill(eta_pT);
		if (fabs(eta_mother_pdg) <= 6){
			Quark_fraction->Fill(eta_pT);
		}else{
			Gluon_fraction->Fill(eta_pT);
		}
	}
	Quark_fraction->Divide(Total_fraction);
	Gluon_fraction->Divide(Total_fraction);

	Quark_fraction->Draw();
	Gluon_fraction->Draw("SAME");
}