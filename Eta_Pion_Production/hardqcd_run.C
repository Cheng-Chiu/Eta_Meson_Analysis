#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include <math.h>
void hardqcd_run(Int_t nev  = 15000000, Int_t ndeb = 1){

//Load libraries
	gSystem->Load("libEG");
	gSystem->Load("libEGPythia8");
	gStyle->SetOptFit(1111);

//Array of particles
	TClonesArray* particles = new TClonesArray("TParticle", 1000);
//Create pythia8 object
	TPythia8* pythia8 = new TPythia8();

#if PYTHIA_VERSION_INTEGER == 8235
//Pythia 8.235 is known to cause crashes:
	printf("ABORTING PYTHIA8 TUTORIAL!\n");
	printf("The version of Pythia you use is known to case crashes due to memory errors.\n");
	printf("They have been reported to the authors; the Pythia versions 8.1... are known to work.\n");
	return;
#endif

//Configure
	pythia8->ReadString("HardQCD:all = on");
	//pythia8->ReadString("PhaseSpace:pTHatMin = 2.");
	//pythia8->ReadString("Random:setSeed = on");
	//pythia8->ReadString("Random:seed = 42");


//Initialize
	double root_s = 500.;
	pythia8->Initialize(2212 /* p */, 2212 /* p */, root_s /* GeV */);

	TFile *output = new TFile("hardqcd_for_15.root", "RECREATE");
   
//Book trees
	TTree *eta_tree = new TTree("eta_tree", "ETA");
	double eta_pT, eta_prap, eta_engy, eta_xF, eta_mother_pT, eta_mother_pdg, eta_mother_engy, 
	eta_part1_pdg, eta_part2_pdg, eta_part1_xBj, eta_part2_xBj,
	eta_part3_pdg, eta_part4_pdg, eta_part3_xBj, eta_part4_xBj;
	eta_tree->Branch("eta_pT",&eta_pT,"eta_pT/D");
	eta_tree->Branch("eta_prap",&eta_prap,"eta_prap/D");
	eta_tree->Branch("eta_engy",&eta_engy,"eta_engy/D");
	eta_tree->Branch("eta_xF",&eta_xF,"eta_xF/D");
	eta_tree->Branch("eta_mother_pT",&eta_mother_pT,"eta_mother_pT/D");
	eta_tree->Branch("eta_mother_pdg",&eta_mother_pdg,"eta_mother_pdg/D");
	eta_tree->Branch("eta_mother_engy",&eta_mother_engy,"eta_mother_energy/D");
	eta_tree->Branch("eta_part1_pdg",&eta_part1_pdg,"eta_part1_pdg/D");
	eta_tree->Branch("eta_part2_pdg",&eta_part2_pdg,"eta_part2_pdg/D");
	eta_tree->Branch("eta_part1_xBj",&eta_part1_xBj,"eta_part1_xBj/D");
	eta_tree->Branch("eta_part2_xBj",&eta_part2_xBj,"eta_part2_xBj/D");
	eta_tree->Branch("eta_part3_pdg",&eta_part3_pdg,"eta_part3_pdg/D");
	eta_tree->Branch("eta_part4_pdg",&eta_part4_pdg,"eta_part4_pdg/D");
	eta_tree->Branch("eta_part3_xBj",&eta_part3_xBj,"eta_part3_xBj/D");
	eta_tree->Branch("eta_part4_xBj",&eta_part4_xBj,"eta_part4_xBj/D");

	
	TTree *pion0_tree = new TTree("pion0_tree", "PION0");
	double pion0_pT, pion0_prap, pion0_engy, pion0_xF, pion0_mother_pT, pion0_mother_pdg, pion0_mother_engy,
	pion0_part1_pdg, pion0_part2_pdg, pion0_part1_xBj, pion0_part2_xBj,
	pion0_part3_pdg, pion0_part4_pdg, pion0_part3_xBj, pion0_part4_xBj;
	pion0_tree->Branch("pion0_pT",&pion0_pT,"pion0_pT/D");
	pion0_tree->Branch("pion0_prap",&pion0_prap,"pion0_prap/D");
	pion0_tree->Branch("pion0_engy",&pion0_engy,"pion0_engy/D");
	pion0_tree->Branch("pion0_xF",&pion0_xF,"pion0_xF/D");
	pion0_tree->Branch("pion0_mother_pT",&pion0_mother_pT,"pion0_mother_pT/D");
	pion0_tree->Branch("pion0_mother_pdg",&pion0_mother_pdg,"pion0_mother_pdg/D");
	pion0_tree->Branch("pion0_mother_engy",&pion0_mother_engy,"pion0_mother_energy/D");
	pion0_tree->Branch("pion0_part1_pdg",&pion0_part1_pdg,"pion0_part1_pdg/D");
	pion0_tree->Branch("pion0_part2_pdg",&pion0_part2_pdg,"pion0_part2_pdg/D");
	pion0_tree->Branch("pion0_part1_xBj",&pion0_part1_xBj,"pion0_part1_xBj/D");
	pion0_tree->Branch("pion0_part2_xBj",&pion0_part2_xBj,"pion0_part2_xBj/D");
	pion0_tree->Branch("pion0_part3_pdg",&pion0_part3_pdg,"pion0_part3_pdg/D");
	pion0_tree->Branch("pion0_part4_pdg",&pion0_part4_pdg,"pion0_part4_pdg/D");
	pion0_tree->Branch("pion0_part3_xBj",&pion0_part3_xBj,"pion0_part3_xBj/D");
	pion0_tree->Branch("pion0_part4_xBj",&pion0_part4_xBj,"pion0_part4_xBj/D");	

	
//Event loop
	for (Int_t iev = 0; iev < nev; iev++) {
		pythia8->GenerateEvent();
		if (iev < ndeb) pythia8->EventListing();
		pythia8->ImportParticles(particles,"All");
		Int_t np = particles->GetEntriesFast();

//Particle loop
		for (Int_t ip = 0; ip < np; ip++) {
			TParticle* part = (TParticle*) particles->At(ip);
			int pdg = part->GetPdgCode();
			if (fabs(part->Eta()) > 3.8 || fabs(part->Eta()) < 3.) continue; //forward
	 		//if (fabs(part->Eta()) > 1) continue; //mid
	 		if (pdg != 221 && pdg != 111) continue; //pi0 = 111, eta = 221
	 		Int_t motherindex = part->GetFirstMother();
	 		TParticle* mother = (TParticle*) particles->At(motherindex);
	 		int mother_pdg = mother->GetPdgCode();
	 		if (mother_pdg != 21 && abs(mother_pdg) > 6) continue; //select only quarks and gluons
	 		TParticle* proton1 = (TParticle*) particles->At(0);
			TParticle* proton2 = (TParticle*) particles->At(1);
			TParticle* part1 = (TParticle*) particles->At(2);
			TParticle* part2 = (TParticle*) particles->At(3);
			TParticle* part3 = (TParticle*) particles->At(4);
			TParticle* part4 = (TParticle*) particles->At(5);
	 		if (pdg == 221){
	 			eta_pT = part->Pt();
		 		eta_prap = part->Eta();
		 		eta_xF = part->Pz()/(0.5*root_s);
		 		eta_mother_pT = mother->Pt();
		 		eta_mother_pdg = mother_pdg;
		 		eta_mother_engy = mother->Energy();
		 		eta_part1_pdg = part1->GetPdgCode();
		 		eta_part2_pdg = part2->GetPdgCode();
		 		eta_part3_pdg = part3->GetPdgCode();
		 		eta_part4_pdg = part4->GetPdgCode();
				eta_part1_xBj = part1->Pz()/proton1->Pz();
		 		eta_part2_xBj = part2->Pz()/proton2->Pz();
		 		eta_part3_xBj = part3->Pz()/proton1->Pz();
		 		eta_part4_xBj = part4->Pz()/proton2->Pz();
		 		eta_tree->Fill();
		 	}else{
		 		pion0_pT = part->Pt();
		 		pion0_prap = part->Eta();
		 		pion0_engy = part->Energy();
		 		pion0_xF = part->Pz()/(0.5*root_s);
		 		pion0_mother_pT = mother->Pt();
		 		pion0_mother_pdg = mother_pdg;
		 		pion0_mother_engy = mother->Energy();
		 		pion0_part1_pdg = part1->GetPdgCode();
		 		pion0_part2_pdg = part2->GetPdgCode();
		 		pion0_part3_pdg = part3->GetPdgCode();
		 		pion0_part4_pdg = part4->GetPdgCode();
		 		pion0_part1_xBj = part1->Pz()/proton1->Pz();
		 		pion0_part2_xBj = part2->Pz()/proton2->Pz();
		 		pion0_part3_xBj = part3->Pz()/proton1->Pz();
		 		pion0_part4_xBj = part4->Pz()/proton2->Pz();
		 		pion0_tree->Fill();
		 	}
		 }
	 }
	eta_tree->Write();
	pion0_tree->Write();
	output->Write();
	output->Close();
 }
