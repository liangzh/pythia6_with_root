#include "Analysis.h"
#include "math.h"
#include "TLorentzVector.h"
#include "TVector3.h"

ClassImp(Analysis)

TH1D *getLogXBin(TString histname, TString histaxis, 
		int bin_num, double bin_min, double bin_max){

		double x_edge[bin_num+1];
		for(int i(0); i<=bin_num; i++){
			double step = (log10(bin_max)-log10(bin_min))/bin_num;
			x_edge[i]=pow(10, log10(bin_min)+i*step);
		}

		TH1D *htemp = new TH1D(histname, histaxis, bin_num, x_edge);
		return htemp;
}

//=================Initializer========================
Analysis::Analysis(){
	TH1::SetDefaultSumw2();
	//track wise
	h_theta_scattered = new TH1D("h_theta_scattered","theta of scattered electron;#theta;counts",100, 0, 3.14159);
	h_D0_eta = new TH1D("h_D0_eta", "#eta;#eta;d#sigma/d#eta [nb]", 50, -4.5, 4.5);
	h_D0_z = new TH1D("h_D0_z", "z;z;d#sigma/dz [nb]", 10, 0, 1);
	h_D0_pt = getLogXBin("h_D0_pt", "pt;pt;d#sigma/dpt [nb/GeV]", 40, 0.1, 20);
	//pair wise
	h_phi_pair = new TH1D("h_phi_pair", ";#phi;N(#phi)", 10, 0, 2*3.14);
	h_phi_pair_PGF = new TH1D("h_phi_pair_PGF", ";#phi;N(#phi)", 10, 0, 2*3.14);
	h_deltaphi = new TH1D("h_deltaphi", ";#Delta#phi;N(#Delta#phi)", 50, -0.5*3.142, 1.5*3.142);
	h_deltaphi_PGF = new TH1D("h_deltaphi_PGF", ";#Delta#phi;N(#Delta#phi)", 50, -0.5*3.142, 1.5*3.142);


	//event wise
	h_QSquared = new TH1D("h_QSquared", "Q^{2};log_{10}(Q^{2});counts", 50, 0, 4);
	h_xBj = new TH1D("h_xBj", "x_{Bj};log_{10}(x_{Bj});counts", 50, -5, 0);
	h_W2 = new TH1D("h_W2", "W^{2};log_{10}(W^{2});counts", 100, 0, 6);
	h_y = new TH1D("h_y", "y;log_{10}(y);counts", 50, -5, 0);
	h_Q2VsxBj = new TH2D("h_Q2VsxBj", "log_{10}(Q^{2}) vs log_{10}(x_{Bj});log_{10}(x_{Bj});log_{10}(Q^{2})", 50, -5, 0, 50, 0, 4);
	h_Q2VsxBj_D0 = new TH2D("h_Q2VsxBj_D0", "log_{10}(Q^{2}) vs log_{10}(x_{Bj});log_{10}(x_{Bj});log_{10}(Q^{2})", 50, -5, 0, 50, 0, 4);

	h_alpha = new TH1D("h_alpha", "alpha_{em};#alpha_{em};counts", 50, 0, 0.01);

	status = false;
}

//destructor
Analysis::~Analysis(){

}


//analyse current event and fill the corresponding histograms
//if return 1, this event is not analysed, throw another event
//if return 0, this event is accepted and do analysis
void Analysis::FillHist(const Event *event){
	status = false;
	double Q2 = pyint1.vint[306];
	double y = pyint1.vint[308];
	double sqs = pypars.pari[100];

	double xBj = Q2/sqs/sqs/y;
	double W = sqrt(sqs*sqs*y);

	h_xBj->Fill(log10(xBj));
	h_QSquared->Fill(log10(Q2));
	h_W2->Fill(log10(y*sqs*sqs));
	h_y->Fill(log10(y));
	h_Q2VsxBj->Fill(log10(xBj), log10(Q2));

	double alpha_em = pydat1.paru[107];
	h_alpha->Fill(alpha_em);

	//find the scattered electron and fill the theta
	TLorentzVector v_ebeam;
	TLorentzVector v_e_out;
	TLorentzVector v_pbeam;
	TLorentzVector vPh, vPhTrue;

	double Ee_beam = -1;
	//fill e/p beam and out e, exchanged photon info
	for(int itrack=0; itrack<event->GetNTrack(); itrack++){
		ParticleMC *particle = event->GetTrack(itrack);
		if(itrack==0) {
			Ee_beam=particle->GetE();
			v_ebeam.SetPxPyPzE(particle->GetPx(), particle->GetPy(), particle->GetPz(), particle->GetE());
		}
		if(itrack==1) {
			v_pbeam.SetPxPyPzE(particle->GetPx(), particle->GetPy(), particle->GetPz(), particle->GetE());
		}
		if(particle->GetKS()==1&&particle->GetPid()==11&&particle->GetParent()==3){
			v_e_out.SetPxPyPzE(particle->GetPx(), particle->GetPy(), particle->GetPz(), particle->GetE());
			h_theta_scattered->Fill(particle->GetTheta());
			break;
		}
	}//for

	vPh = v_ebeam-v_e_out;
	vPhTrue = vPh;

	//get lorentz boost vector to gamma p center of mass frame
	TVector3 b=-(v_pbeam.Vect()+vPhTrue.Vect())*(1./(v_pbeam.E()+vPhTrue.E()));
	vPh.Boost(b);

	TLorentzVector vhadron, vTrig, vTotal;
	int trigId = 0;
	int assoId = 0;

	for(int itrack=0; itrack<event->GetNTrack(); itrack++){
		ParticleMC *particle = event->GetTrack(itrack);
		if(particle->GetKS()==1&&abs(particle->GetPid())==421&&
				fabs(particle->GetEta())<4.5){
			vhadron.SetPxPyPzE(particle->GetPx(), particle->GetPy(), particle->GetPz(), particle->GetE());
			double z1 = (particle->GetE()-particle->GetPz())/(2*y*Ee_beam);
			double z = (vhadron*v_pbeam)/(v_pbeam*vPhTrue);
			vhadron.Boost(b);
			vhadron.RotateZ(-vPh.Phi());
			vhadron.RotateY(-vPh.Theta());
			h_D0_eta->Fill(particle->GetEta());
			h_D0_pt->Fill(particle->GetPt());
			h_D0_z->Fill(z);

			//z cut for the particles
			if(z>0.75||z<0.25) continue;

			//look for trig
			if( trigId==0 && vhadron.Pt()>3 ){
				trigId = particle->GetPid();
				vTrig=vhadron;
			}//if trig

			//look for associate
			if( trigId!=0 && trigId == - particle->GetPid() ){
				vTotal = vTrig+vhadron;
				double kperp = vTotal.Pt();
				double Pperp = (vTrig-vhadron).Pt()/2.;
				//correlation limit
				if(kperp/Pperp<0.5){
					assoId=particle->GetPid();
					double kt_phi = vTotal.Phi();
					double delta_phi = vTrig.DeltaPhi(vhadron);

					if(delta_phi < -3.14159/2.) delta_phi+=2*3.14159;
					if(delta_phi > 3*3.14159/2.) delta_phi-=2*3.14159;
					if(kt_phi<0) kt_phi+=2*3.14159;

					h_phi_pair->Fill(kt_phi);
					h_deltaphi->Fill(delta_phi);
					if(event->GetProcess()>132 && abs(event->GetOutParton1())==4){
						h_phi_pair_PGF->Fill(kt_phi);
						h_deltaphi_PGF->Fill(delta_phi);
					}
				}
			}//if assoc
		}//if
	}//for
	
//	std::cout<<num_D0<<std::endl;

	if(assoId!=0) {
		status = true;
		h_Q2VsxBj_D0->Fill(log10(xBj), log10(Q2));
	}

}


void Analysis::WriteResults(TFile *file){
	TDirectory *save = gDirectory;
	file->cd();

	double total_cs = pypars.pari[0]; //in mb
//	double total_nevt = pypars.msti[4]; //number of gen events
	double total_ntrials = pyint5.ngen[2][0]; //number of gen events
	double lumi = total_ntrials/total_cs; //in mb^-1
	//1E6 added to convert mb to nb
	h_D0_eta->Scale(1E6/lumi,"width"); //dsigma/dX in nb
	h_D0_pt->Scale(1E6/lumi,"width"); //dsigma/dX in nb
	h_D0_z->Scale(1E6/lumi,"width"); //dsigma/dX in nb

	//track wise
	h_theta_scattered->Write();
	h_D0_eta->Write();
	h_D0_pt->Write();
	h_D0_z->Write();
	//pair wise
	h_phi_pair->Write();
	h_phi_pair_PGF->Write();
	h_deltaphi->Write();
	h_deltaphi_PGF->Write();


	//event wise
	h_QSquared->Write();
	h_xBj->Write();
	h_W2->Write();
	h_y->Write();
	h_Q2VsxBj->Write();
	h_Q2VsxBj_D0->Write();
	h_alpha->Write();

	//D0: 421, D*: 413
	std::cout<<"421 comp code:"<<call_pycomp(421)<<std::endl;

	//cd back to original dir
	save->cd();
}

