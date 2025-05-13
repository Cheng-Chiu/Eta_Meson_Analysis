#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include <math.h>
void pythia8_part_ratio(){
  double arr_x[8] = {0.24903574202822917, 0.35749586925379195, 0.5593241243560779, 0.9294462918156647, 2.479662209584589, 4.3389353336827705, 8.491256550216814, 0.14858077415017457};
  double arr_y[8] = {0.04194491333525463, 0.06471406809666924, 0.11388200419943184, 0.17183745728519084, 0.36000823418008154, 0.4778558195067727, 0.5005064020750135, 0.028045617357651542};

  TGraph *graph = new TGraph(8, arr_x, arr_y);
  graph->SetMarkerStyle(21);
  
  TFile *f_mid_eta = new TFile("eta_mid_rap.root");
  TFile *f_mid_pi0 = new TFile("pi0_mid_rap.root");
  TFile *f_for_eta = new TFile("eta_for_rap.root");
  TFile *f_for_pi0 = new TFile("pi0_for_rap.root");
  
  TH1D *h_mid_eta = (TH1D*)f_mid_eta->Get("ETA");
  TH1D *h_mid_pi0 = (TH1D*)f_mid_pi0->Get("ETA");
  TH1D *h_for_eta = (TH1D*)f_for_eta->Get("ETA");
  TH1D *h_for_pi0 = (TH1D*)f_for_pi0->Get("ETA");
  
  TH1D *h_mid_ratio = new TH1D(*h_mid_eta); //
  h_mid_ratio->Divide(h_mid_pi0);
  h_mid_ratio->SetNameTitle("Mid_Rapid_Ratio", "eta / pi0 at Mid Rapidity");
  
  TH1D *h_for_ratio = new TH1D(*h_for_eta);
  h_for_ratio->Divide(h_for_pi0);
  h_for_ratio->SetNameTitle("For_Rapid_Ratio", "eta / pi0 at Forward Rapidity");
  

  h_for_ratio->SetFillColor(kRed);
  auto c1 = new TCanvas("c1","c1");
  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry("Mid_Rapid_Ratio","Mid","f");
  legend->AddEntry("For_Rapid_Ratio","Forward","l");
  
  h_mid_ratio->Draw();
  h_for_ratio->Draw("SAME");
  graph->Draw("PSAME");
  
  
  
}
