#include "HistMaker.h"
#include "Event.h"

ClassImp(HistMaker)

//=================HistMaker part========================
HistMaker::HistMaker(){
	//track wise
	h_FS_pt = new TH1D("h_FS_pt","pt for the final state particles;pt;counts",100, 0, 25);
	h_FS_eta = new TH1D("h_FS_eta","#eta for final state particles;#eta;counts",150, -15, 15);
	h_FS_rapidity = new TH1D("h_FS_rapidity","y for final state particles;y;counts",150, -15, 15);
	h_FS_pt_Pi = new TH1D("h_FS_pt_Pi","pt for final charged pions;pt;counts",100, 0, 25);
	h_FS_pt_K = new TH1D("h_FS_pt_K","pt for final charged kaons;pt;counts",100, 0, 25);
	h_FS_pt_P = new TH1D("h_FS_pt_P","pt for final p/pbar;pt;counts",100, 0, 25);
	h_FS_PID = new TH1D("h_FS_PID","PID for final state particles;PID;counts",6000, -3000, 3000);

	//event wise
	h_FS_num = new TH1D("h_FS_num","Multiplicity;N_{FS};counts",100, 0, 200);
	h_process = new TH1D("h_process","process;process id;counts",100, 1, 100);
	h_pt2hat = new TH1D("h_pt2hat","pt2hat;pt2hat;counts",200, 0, 100);

}

HistMaker::~HistMaker(){

}

void HistMaker::Hfill(const Event *event){
	//Fill histograms
	h_FS_num->Fill(event->GetFSNTrack());
	h_process->Fill(event->GetProcess());
	h_pt2hat->Fill(event->GetPt2_hat());

	for(Int_t itrack = 0; itrack<event->GetNTrack(); itrack++){
		ParticleMC* particle = event->GetTrack(itrack);
		if(particle->GetKS()==1){
			h_FS_pt->Fill(particle->GetPt());
			h_FS_eta->Fill(particle->GetEta());
			h_FS_rapidity->Fill(particle->GetRapidity());
			h_FS_PID->Fill(particle->GetPid());

			if(abs(particle->GetPid())==211)
				h_FS_pt_Pi->Fill(particle->GetPt());
			if(abs(particle->GetPid())==321)
				h_FS_pt_K->Fill(particle->GetPt());
			if(abs(particle->GetPid())==2212)
				h_FS_pt_P->Fill(particle->GetPt());
		}//if
	}//for
}

void HistMaker::Hwrite(TFile *file){
	TDirectory *save = gDirectory;
	file->cd();

	//track wise
	h_FS_pt->Write();
	h_FS_eta->Write();
	h_FS_rapidity->Write();
	h_FS_pt_Pi->Write();
	h_FS_pt_K->Write();
	h_FS_pt_P->Write();
	h_FS_PID->Write();

	//event wise
	h_FS_num->Write();
	h_process->Write();
	h_pt2hat->Write();

	//cd back to original dir
	save->cd();
}

