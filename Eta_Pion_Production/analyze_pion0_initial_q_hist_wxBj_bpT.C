#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include <math.h>
#include <vector>
void analyze_pion0_initial_q_hist_wxBj_bpT(){
	TFile* f = new TFile("hardqcd_for_15.root");
	//TFile* f = new TFile("softqcd_mid_15.root");

	
	double pion0_part1_xBj, pion0_part2_xBj, pion0_part1_pdg, pion0_part2_pdg, pion0_pT;
	TTree* PION0 = (TTree*) f->Get("pion0_tree");
	PION0->SetBranchAddress("pion0_part1_xBj",&pion0_part1_xBj);
	PION0->SetBranchAddress("pion0_part2_xBj",&pion0_part2_xBj);
	PION0->SetBranchAddress("pion0_part1_pdg",&pion0_part1_pdg);
	PION0->SetBranchAddress("pion0_part2_pdg",&pion0_part2_pdg);
	PION0->SetBranchAddress("pion0_pT",&pion0_pT);

	//set log bins
	const int nbins = 50;
	Double_t xbins[nbins + 1];
	double dx = 1./nbins;
	double L10 = TMath::Log(1E5);
	double xmin = 1E-5;
	for (int i = 0; i <= nbins; i++) {
		xbins[i] = xmin*TMath::Exp(L10 * i * dx);
	}

	vector<vector<TH1F*>> hist_pT(4, vector<TH1F*>(4));
	for (size_t i = 0; i < 4; ++i){
		vector<TH1F*> temp;
		string g = "Gluon" + to_string(i);
		string d = "Down" + to_string(i);
		string u = "Up" + to_string(i);
		string s = "Sea Quark" + to_string(i);
		temp.push_back(new TH1F(g.c_str(), "", nbins, xbins));
		temp.push_back(new TH1F(d.c_str(), "", nbins, xbins));
		temp.push_back(new TH1F(u.c_str(), "", nbins, xbins));
		temp.push_back(new TH1F(s.c_str(), "", nbins, xbins));
		hist_pT[i] = temp;
	}
	std::vector<size_t> entries_pT(4);
	int PION0_Entries = PION0->GetEntries();
	for (Int_t i = 0; i < PION0_Entries; i++){
		PION0->GetEntry(i);
		//excluding diffractive events
		if (pion0_part1_pdg == 2212 || pion0_part1_pdg == 9902210) continue;
		if (pion0_part2_pdg == 2212 || pion0_part2_pdg == 9902210) continue;
		uint32_t pT_region = 4;
		if (pion0_pT >= 1 && pion0_pT < 2){
			pT_region = 0;
			entries_pT[0]++;
		}else if (pion0_pT >=2 && pion0_pT < 3){
			pT_region = 1;
			entries_pT[1]++;
		}else if (pion0_pT >= 3 && pion0_pT < 4){
			pT_region = 2;
			entries_pT[2]++;
		}else if (pion0_pT >= 4 && pion0_pT < 6){
			pT_region = 3;
			entries_pT[3]++;
		}else{
			continue;
		}
		if (fabs(pion0_part1_pdg) <= 6 || pion0_part1_pdg == 21){
			switch(static_cast<int>(pion0_part1_pdg)){
				case 21:
					hist_pT[pT_region][0]->Fill(pion0_part1_xBj);
					break;
				case 1: 
					hist_pT[pT_region][1]->Fill(pion0_part1_xBj);
					break;
				case 2:
					hist_pT[pT_region][2]->Fill(pion0_part1_xBj);
					break;
				default:
					hist_pT[pT_region][3]->Fill(pion0_part1_xBj);
			}
		}
		
		if (fabs(pion0_part2_pdg) <= 6 || pion0_part2_pdg == 21){
			switch(static_cast<int>(pion0_part2_pdg)){
				case 21:
					hist_pT[pT_region][0]->Fill(pion0_part2_xBj);
					break;
				case 1: 
					hist_pT[pT_region][1]->Fill(pion0_part2_xBj);
					break;
				case 2:
					hist_pT[pT_region][2]->Fill(pion0_part2_xBj);
					break;
				default:
					hist_pT[pT_region][3]->Fill(pion0_part2_xBj);
			}
		}
	}
	for (size_t i = 0; i < 4; i++){
		for (size_t j = 0; j < 4; j++){
			hist_pT[i][j]->Scale(1./entries_pT[i]);
		}
	}

	//plot
	TCanvas *c1 = new TCanvas("c1","multipads",1000,600);
   	gStyle->SetOptStat(0);
   	c1->Divide(2,2,1E-4,1E-8);
	for (size_t i = 0; i < 4; ++i){
		c1->cd(i+1);
   		gPad->SetTickx(2);
   		gPad->SetLogx();
   		hist_pT[i][0]->GetYaxis()->SetRangeUser(0., 0.15);
	   	hist_pT[i][0]->GetXaxis()->SetRangeUser(1E-4, 1);
	   	for (size_t j = 0; j < 4; ++j){
	   		hist_pT[i][j]->Draw("SAME");
	   		hist_pT[i][j]->SetLineColor(j+1);
	   		hist_pT[i][j]->SetMarkerStyle(21);
			hist_pT[i][j]->SetMarkerColor(j+1);
			hist_pT[i][j]->SetMarkerSize(0.6);
	   	}
	}

	c1->cd(1);
	hist_pT[0][0]->GetXaxis()->SetTitle("x_{Bj}");
	hist_pT[0][0]->GetXaxis()->SetTitleSize(0.06);
	hist_pT[0][0]->GetXaxis()->SetTitleOffset(0.7);

	TLegend *leg = new TLegend(0.15,0.4,0.35,0.7);
	leg->AddEntry(hist_pT[0][0], "g", "pe");
	leg->AddEntry(hist_pT[0][2], "u", "pe");
	leg->AddEntry(hist_pT[0][1], "d", "pe");
	leg->AddEntry(hist_pT[0][3], "sea", "pe");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.05);
	leg->Draw();
	
	TLatex latex;
	latex.SetTextSize(0.045);
	latex.SetTextAlign(13);
	latex.DrawLatex(.00015,.143,"p + p #rightarrow #pi^{0} + X (#sqrt{s} = 500 GeV)");
	latex.DrawLatex(.00015,.133,"3.0 < |#eta| < 3.8");
	//latex.DrawLatex(.00015,.133,"|#eta| < 1");
   	latex.DrawLatex(.00015,.123,"1 #leq p_{T} < 2 (GeV/c)");

   	c1->cd(2);
   	latex.DrawLatex(.00015,.123,"2 #leq p_{T} < 3 (GeV/c)");
   	c1->cd(3);
   	latex.DrawLatex(.00015,.123,"3 #leq p_{T} < 4 (GeV/c)");
   	c1->cd(4);
   	latex.DrawLatex(.00015,.123,"4 #leq p_{T} < 6 (GeV/c)");
   	c1->SaveAs("pion0_initial_wxBj_hist_h_f_4pTwindow_15.pdf");
}