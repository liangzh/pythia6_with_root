#include "Analysis.h"
#include "math.h"

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
	//track wise
	h_theta_scattered = new TH1D("h_theta_scattered","theta of scattered electron;#theta;counts",100, 0, 3.14159);
	h_Dstar_eta = new TH1D("h_Dstar_eta", "#eta;#eta;d#sigma/d#eta [nb]", 50, -1.8, 1.8);
	h_Dstar_z = new TH1D("h_Dstar_z", "z;z;d#sigma/dz [nb]", 10, 0, 1);
	h_Dstar_pt = getLogXBin("h_Dstar_pt", "pt;pt;d#sigma/dpt [nb/GeV]", 40, 1.25, 20);

	//event wise
	h_QSquared = new TH1D("h_QSquared", "Q^{2};log_{10}(Q^{2});counts", 50, 0, 4);
	h_xBj = new TH1D("h_xBj", "x_{Bj};log_{10}(x_{Bj});counts", 50, -5, 0);
	h_W2 = new TH1D("h_W2", "W^{2};log_{10}(W^{2});counts", 100, 0, 6);
	h_y = new TH1D("h_y", "y;log_{10}(y);counts", 50, -5, 0);
	h_Q2VsxBj = new TH2D("h_Q2VsxBj", "log_{10}(Q^{2}) vs log_{10}(x_{Bj});log_{10}(x_{Bj});log_{10}(Q^{2})", 50, -5, 0, 50, 0, 4);
	h_Q2VsxBj_Dstar = new TH2D("h_Q2VsxBj_Dstar", "log_{10}(Q^{2}) vs log_{10}(x_{Bj});log_{10}(x_{Bj});log_{10}(Q^{2})", 50, -5, 0, 50, 0, 4);

	h_alpha = new TH1D("h_alpha", "alpha_{em};#alpha_{em};counts", 50, 0, 0.01);

	for(int i(0); i<7; i++){
		h_sig_reduced_cc[i] = getLogXBin(Form("h_sig_reduced_cc%d",i), ";x;#sigma_{ccbar}^{red}", 50, 1E-5, 1);
		h_dsig[i] = getLogXBin(Form("h_dsig%d",i), ";x;d#sigma_{ccbar}/dx/dQ^{2}", 50, 1E-5, 1);
	}

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

	double prefactor = xBj*Q2*Q2/2/3.14159/alpha_em/alpha_em/(1+(1-y)*(1-y));

	//fill reduced cross section for ccbar process
	if(event->GetOutParton1()==-event->GetOutParton2() &&
			abs(event->GetOutParton1())==4) {
		for(int i(0); i<7; i++){
			int icut = 2*i;
			if(Q2<Q2_cut[icut+1] && Q2>Q2_cut[icut]){
				h_dsig[i]->Fill(xBj);
				h_sig_reduced_cc[i]->Fill(xBj, prefactor);
			}
		}
	}

	//if no charm production, skip analysis
//	if(abs(event->GetOutParton1())!=4 && abs(event->GetOutParton1())!=4)
//		return 1;

	int num_Dstar=0;
	double Ee_beam = -1;
	//find the scattered electron and fill the theta
	if(Q2>5&&Q2<100&&y>0.02&&y<0.7){
//	if(Q2<2&&W>100&&W<285){
		for(int itrack=0; itrack<event->GetNTrack(); itrack++){
			ParticleMC *particle = event->GetTrack(itrack);
			if(itrack==0) Ee_beam=particle->GetE();
			if(particle->GetKS()==1&&particle->GetPid()==11&&particle->GetParent()==3)
				h_theta_scattered->Fill(particle->GetTheta());
			if(particle->GetKS()==1&&abs(particle->GetPid())==413
//					&&particle->GetPt()>1.8&&fabs(particle->GetEta())<1.5){
					&&particle->GetPt()>1.25&&fabs(particle->GetEta())<1.8){
				h_Dstar_eta->Fill(particle->GetEta());
				h_Dstar_pt->Fill(particle->GetPt());
				double z = (particle->GetE()-particle->GetPz())/(2*y*Ee_beam);
				h_Dstar_z->Fill(z);
				num_Dstar++;
//				std::cout<<"found 1 Dstar"<<std::endl;
			}
		}
	}
	
//	std::cout<<num_Dstar<<std::endl;

	if(num_Dstar>0) {
		status = true;
		h_Q2VsxBj_Dstar->Fill(log10(xBj), log10(Q2));
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
	h_Dstar_eta->Scale(1E6/lumi,"width"); //dsigma/dX in nb
	h_Dstar_pt->Scale(1E6/lumi,"width"); //dsigma/dX in nb
	h_Dstar_z->Scale(1E6/lumi,"width"); //dsigma/dX in nb

	//track wise
	h_theta_scattered->Write();
	h_Dstar_eta->Write();
	h_Dstar_pt->Write();
	h_Dstar_z->Write();

	//event wise
	h_QSquared->Write();
	h_xBj->Write();
	h_W2->Write();
	h_y->Write();
	h_Q2VsxBj->Write();
	h_Q2VsxBj_Dstar->Write();
	h_alpha->Write();

	std::cout<<"421 comp code:"<<call_pycomp(-421)<<std::endl;

	for(int i(0); i<7; i++){
		int icut = 2*i;
		double Q2width = Q2_cut[icut+1]-Q2_cut[icut];
		double mb_to_GeV2 = 1./0.3894; //mb to GeV^-2
		h_dsig[i]->Scale(1./(lumi*Q2width),"width");
		h_sig_reduced_cc[i]->Scale(1./(lumi*Q2width)*mb_to_GeV2,"width");
		h_dsig[i]->Write();
		h_sig_reduced_cc[i]->Write();
	}

	//cd back to original dir
	save->cd();
}

