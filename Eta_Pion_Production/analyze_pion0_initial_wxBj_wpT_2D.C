#include "TSystem.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include <math.h>
void analyze_pion0_initial_wxBj_wpT_2D(){
	gStyle->SetOptStat(0);
	TFile* f = new TFile("hardqcd_mid_15.root");
	//TFile* f = new TFile("softqcd_mid_15.root");
	double pion0_part1_xBj, pion0_part2_xBj, pion0_pT;
	double pion0_part1_pdg, pion0_part2_pdg;
	TTree* PION0 = (TTree*) f->Get("pion0_tree");
	PION0->SetBranchAddress("pion0_part1_xBj",&pion0_part1_xBj);
	PION0->SetBranchAddress("pion0_part2_xBj",&pion0_part2_xBj);
	PION0->SetBranchAddress("pion0_part1_pdg",&pion0_part1_pdg);
	PION0->SetBranchAddress("pion0_part2_pdg",&pion0_part2_pdg);
	const int nbins = 50;
	Double_t xbins[nbins + 1];
	double dx = 1./nbins;
	double L10 = TMath::Log(1E4);
	double xmin = 1E-4;
	for (int i = 0; i <= nbins; i++) {
		xbins[i] = xmin*TMath::Exp(L10 * i * dx);
	}
	TH2F *xBj_contour = new TH2F("X1X2", "", nbins, xbins, nbins, xbins);	


	int PION0_Entries = PION0->GetEntries();
	for (Int_t i = 0; i < PION0_Entries; i++){
		PION0->GetEntry(i);
		//excluding diffractive events
		if (pion0_part1_pdg == 2212 || pion0_part1_pdg == 9902210) continue;
		if (pion0_part2_pdg == 2212 || pion0_part2_pdg == 9902210) continue;
		xBj_contour->Fill(pion0_part1_xBj, pion0_part2_xBj);
	}

	TCanvas *c1 = new TCanvas("c1","",800,600);
	c1->SetLogx();
	c1->SetLogy();
	xBj_contour->Draw("COLZ");
	xBj_contour->SetContour(0);
	xBj_contour->GetXaxis()->SetTitle("x_{Bj, 1}");
	xBj_contour->GetXaxis()->SetTitleSize(0.05);
	xBj_contour->GetXaxis()->SetTitleOffset(0.9);
	xBj_contour->GetYaxis()->SetTitle("x_{Bj, 2}");
	xBj_contour->GetYaxis()->SetTitleSize(0.05);
	xBj_contour->GetYaxis()->SetTitleOffset(0.9);

	TLatex latex;
	latex.SetTextSize(0.045);
	latex.SetTextAlign(13);
	//hard
	latex.DrawLatex(.0002,.001,"p + p #rightarrow #pi^{0} + X (#sqrt{s} = 500 GeV)");
	//latex.DrawLatex(.0002,.00053, "3.0 < |#eta| < 3.8");
	latex.DrawLatex(.0002,.00053,"|#eta| < 1");
	latex.DrawLatex(.0002,.00028, "HardQCD:all = on");

	//soft
	//latex.DrawLatex(.000002,.00004,"p + p #rightarrow #pi^{0} + X (#sqrt{s} = 500 GeV)");
	//latex.DrawLatex(.000002,.000016, "3.0 < |#eta| < 3.8");
	//latex.DrawLatex(.000002,.000016,"|#eta| < 1");
	//latex.DrawLatex(.000002,.000006, "SoftQCD:nonDiffractive = on");

	c1->SetRightMargin(0.13);

    c1->SaveAs("pion0_initial_x1x2_h_m_15.pdf");

	

}
