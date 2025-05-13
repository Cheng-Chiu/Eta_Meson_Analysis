#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include <math.h>
void pythia8_part_pTcut_ratio(){
  double arr_x[8] = {0.24903574202822917, 0.35749586925379195, 0.5593241243560779, 0.9294462918156647, 2.479662209584589, 4.3389353336827705, 8.491256550216814, 0.14858077415017457};
  double arr_y[8] = {0.04194491333525463, 0.06471406809666924, 0.11388200419943184, 0.17183745728519084, 0.36000823418008154, 0.4778558195067727, 0.5005064020750135, 0.028045617357651542};

  TGraph *graph = new TGraph(8, arr_x, arr_y);
  graph->SetMarkerStyle(21);
  
  TFile *f_eta = new TFile("eta_mid_rap_min_pThat10.root");
  TFile *f_pi0 = new TFile("pi0_mid_rap_min_pThat10.root");
  TH1D *h_eta = (TH1D*)f_eta->Get("ETA");
  TH1D *h_pi0 = (TH1D*)f_pi0->Get("ETA");

 
  TH1D *h_ratio = new TH1D(*h_eta);
  h_ratio->Divide(h_pi0);
  h_ratio->SetNameTitle("", "ETA / PION0 at Mid Rapidity, pT Min Cut 10 GeV");
  h_ratio->GetYaxis()->SetTitleOffset(1.);
  h_ratio->GetXaxis()->SetTitleOffset(1.);
  h_ratio->GetYaxis()->SetTitle("Ratio of # of ETA over PION0");
  h_ratio->GetXaxis()->SetTitle("Pseudorapidity");
  auto c1 = new TCanvas("c1","c1");
  h_ratio->Draw();
  graph->Draw("PSAME");
  TLegend *le = new TLegend(0.1,0.7,0.4,0.9);
  le->AddEntry(h_ratio, "Pythia Simulated", "l");
  le->AddEntry(graph, "Measured", "p");
  le->Draw();
  c1->SaveAs("ETA_PION0_mid_10.jpg");
  
}
