#include "math.h"

#include "Analysis.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "Event.h"

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

TH2D *getLogXYBin(TString histname, TString histaxis, 
		int bin_num_X, double bin_min_X, double bin_max_X,
		int bin_num_Y, double bin_min_Y, double bin_max_Y){

		double x_edge[bin_num_X+1];
		for(int i(0); i<=bin_num_X; i++){
			double step = (log10(bin_max_X)-log10(bin_min_X))/bin_num_X;
			x_edge[i]=pow(10, log10(bin_min_X)+i*step);
		}

		double y_edge[bin_num_Y+1];
		for(int i(0); i<=bin_num_Y; i++){
			double step = (log10(bin_max_Y)-log10(bin_min_Y))/bin_num_Y;
			y_edge[i]=pow(10, log10(bin_min_Y)+i*step);
		}


		TH2D *htemp = new TH2D(histname, histaxis, bin_num_X, x_edge, bin_num_Y, y_edge);
		return htemp;
}

TH2D *getLinXLogYBin(TString histname, TString histaxis, 
		int bin_num_X, double bin_min_X, double bin_max_X,
		int bin_num_Y, double bin_min_Y, double bin_max_Y){

		double x_edge[bin_num_X+1];
		for(int i(0); i<=bin_num_X; i++){
			double step = (bin_max_X-bin_min_X)/bin_num_X;
			x_edge[i]=(bin_min_X+i*step);
		}

		double y_edge[bin_num_Y+1];
		for(int i(0); i<=bin_num_Y; i++){
			double step = (log10(bin_max_Y)-log10(bin_min_Y))/bin_num_Y;
			y_edge[i]=pow(10, log10(bin_min_Y)+i*step);
		}

		TH2D *htemp = new TH2D(histname, histaxis, bin_num_X, x_edge, bin_num_Y, y_edge);
		return htemp;
}


//=================Initializer========================
Analysis::Analysis(){
	TH1::SetDefaultSumw2();
	//track wise
	h_theta_scattered = new TH1D("h_theta_scattered","theta of scattered electron;#theta;counts",100, 0, 3.14159);
	h_FS_PID_PGF = new TH1D("h_FS_PID_PGF","PID for final state particles;PID;counts",6000, -3000, 3000);
	h_zVsRapidity = new TH2D("h_zVsRapidity", ";y;z", 50, -3, 3, 50, 0, 1);
	h_zVsEta = new TH2D("h_zVsEta", ";#eta;z", 50, -5, 5, 50, 0, 1);

	//tune test
	h_Prim_dNdEta_1 = new TH1D("h_Prim_dNdEta_1", ";#eta;dN^{ch}/(N_{e}d#eta)", 20, 0, 6);
	h_Prim_dNdEta_2 = new TH1D("h_Prim_dNdEta_2", ";#eta;dN^{ch}/(N_{e}d#eta)", 20, 0, 6);
	h_Prim_dNdPt_1 = getLogXBin("h_Prim_dNdPt_1", ";p_{T};dN^{ch}/(N_{e}dp_{T})", 20, 0.1, 10);
	h_Prim_dNdPt_2 = getLogXBin("h_Prim_dNdPt_2", ";p_{T};dN^{ch}/(N_{e}dp_{T})", 20, 0.1, 10);
	h_Q2_tune = new TH1D("h_Q2_tune", ";Q2;counts", 20, 0, 20);

	//all
	h_Kaon_eta = new TH1D("h_Kaon_eta", "#eta;#eta;d#sigma/d#eta [nb]", 50, -4.5, 4.5);
	h_Kaon_z = new TH1D("h_Kaon_z", "z;z;d#sigma/dz [nb]", 10, 0, 1);
	h_Kaon_pt = getLogXBin("h_Kaon_pt", "pt;pt;d#sigma/dpt [nb/GeV]", 40, 0.1, 20);
	h_Pion_eta = new TH1D("h_Pion_eta", "#eta;#eta;d#sigma/d#eta [nb]", 50, -4.5, 4.5);
	h_Pion_z = new TH1D("h_Pion_z", "z;z;d#sigma/dz [nb]", 10, 0, 1);
	h_Pion_pt = getLogXBin("h_Pion_pt", "pt;pt;d#sigma/dpt [nb/GeV]", 40, 0.1, 20);
	//LODIS
	h_Kaon_eta_LODIS = new TH1D("h_Kaon_eta_LODIS", "#eta;#eta;d#sigma/d#eta [nb]", 50, -4.5, 4.5);
	h_Kaon_z_LODIS = new TH1D("h_Kaon_z_LODIS", "z;z;d#sigma/dz [nb]", 10, 0, 1);
	h_Kaon_pt_LODIS = getLogXBin("h_Kaon_pt_LODIS", "pt;pt;d#sigma/dpt [nb/GeV]", 40, 0.1, 20);
	h_Pion_eta_LODIS = new TH1D("h_Pion_eta_LODIS", "#eta;#eta;d#sigma/d#eta [nb]", 50, -4.5, 4.5);
	h_Pion_z_LODIS = new TH1D("h_Pion_z_LODIS", "z;z;d#sigma/dz [nb]", 10, 0, 1);
	h_Pion_pt_LODIS = getLogXBin("h_Pion_pt_LODIS", "pt;pt;d#sigma/dpt [nb/GeV]", 40, 0.1, 20);
	//Direct
	h_Kaon_eta_Direct = new TH1D("h_Kaon_eta_Direct", "#eta;#eta;d#sigma/d#eta [nb]", 50, -4.5, 4.5);
	h_Kaon_z_Direct = new TH1D("h_Kaon_z_Direct", "z;z;d#sigma/dz [nb]", 10, 0, 1);
	h_Kaon_pt_Direct = getLogXBin("h_Kaon_pt_Direct", "pt;pt;d#sigma/dpt [nb/GeV]", 40, 0.1, 20);
	h_Kaon_pt_PGF = getLogXBin("h_Kaon_pt_PGF", "pt;pt;d#sigma/dpt [nb/GeV]", 40, 0.1, 20);
	h_Kaon_pt_QCDC = getLogXBin("h_Kaon_pt_QCDC", "pt;pt;d#sigma/dpt [nb/GeV]", 40, 0.1, 20);
	h_Pion_eta_Direct = new TH1D("h_Pion_eta_Direct", "#eta;#eta;d#sigma/d#eta [nb]", 50, -4.5, 4.5);
	h_Pion_z_Direct = new TH1D("h_Pion_z_Direct", "z;z;d#sigma/dz [nb]", 10, 0, 1);
	h_Pion_pt_Direct = getLogXBin("h_Pion_pt_Direct", "pt;pt;d#sigma/dpt [nb/GeV]", 40, 0.1, 20);
	//Resolved
	h_Kaon_eta_Resolved = new TH1D("h_Kaon_eta_Resolved", "#eta;#eta;d#sigma/d#eta [nb]", 50, -4.5, 4.5);
	h_Kaon_z_Resolved = new TH1D("h_Kaon_z_Resolved", "z;z;d#sigma/dz [nb]", 10, 0, 1);
	h_Kaon_pt_Resolved = getLogXBin("h_Kaon_pt_Resolved", "pt;pt;d#sigma/dpt [nb/GeV]", 40, 0.1, 20);
	h_Pion_eta_Resolved = new TH1D("h_Pion_eta_Resolved", "#eta;#eta;d#sigma/d#eta [nb]", 50, -4.5, 4.5);
	h_Pion_z_Resolved = new TH1D("h_Pion_z_Resolved", "z;z;d#sigma/dz [nb]", 10, 0, 1);
	h_Pion_pt_Resolved = getLogXBin("h_Pion_pt_Resolved", "pt;pt;d#sigma/dpt [nb/GeV]", 40, 0.1, 20);
	//PGF ssbar
	h_Kaon_eta_PGFss = new TH1D("h_Kaon_eta_PGFss", "#eta;#eta;d#sigma/d#eta [nb]", 50, -4.5, 4.5);
	h_Kaon_z_PGFss = new TH1D("h_Kaon_z_PGFss", "z;z;d#sigma/dz [nb]", 10, 0, 1);
	h_Kaon_pt_PGFss = getLogXBin("h_Kaon_pt_PGFss", "pt;pt;d#sigma/dpt [nb/GeV]", 40, 0.1, 20);
	h_Pion_eta_PGFss = new TH1D("h_Pion_eta_PGFss", "#eta;#eta;d#sigma/d#eta [nb]", 50, -4.5, 4.5);
	h_Pion_z_PGFss = new TH1D("h_Pion_z_PGFss", "z;z;d#sigma/dz [nb]", 10, 0, 1);
	h_Pion_pt_PGFss = getLogXBin("h_Pion_pt_PGFss", "pt;pt;d#sigma/dpt [nb/GeV]", 40, 0.1, 20);


	//D0 in PGF
	h_D0_pt = getLogXBin("h_D0_pt", "pt;pt;d#sigma/dpt [nb/GeV]", 40, 0.1, 20);
	h_D0_eta = new TH1D("h_D0_eta", "#eta;#eta;d#sigma/d#eta [nb]", 50, -4.5, 4.5);
	h_D0_z = new TH1D("h_D0_z", "z;z;d#sigma/dz [nb]", 10, 0, 1);
	h_D0_pt_Dstar_feed = getLogXBin("h_D0_pt_Dstar_feed", "pt from D* feed down;pt;d#sigma/dpt [nb/GeV]", 40, 0.1, 20);
	h_D0_eta_Dstar_feed = new TH1D("h_D0_eta_Dstar_feed", "#eta from D* feed down;#eta;d#sigma/d#eta [nb]", 50, -4.5, 4.5);
	h_D0_z_Dstar_feed = new TH1D("h_D0_z_Dstar_feed", "z from D* feed down;z;d#sigma/dz [nb]", 10, 0, 1);
	h_D0_decay_pi = getLinXLogYBin("h_D0_decay_pi", "p vs #eta for pi from D0; #eta; p [GeV]",
			50, -5, 5, 50, 0.1, 100);
	h_D0_decay_kaon = getLinXLogYBin("h_D0_decay_kaon", "p vs #eta for pi from D0; #eta; p [GeV]",
			50, -5, 5, 50, 0.1, 100);


	//pair wise
	h_Pair_Eta = new TH1D("h_Pair_Eta", ";#eta_{trig}; #eta_{assoc}", 50, -5, 5);
	h_RapidityCorre = new TH2D("h_RapidityCorre", ";y_trig; y_assoc", 50, -5, 5, 50, -5, 5);
	h_RapidityCorre_PGF = new TH2D("h_RapidityCorre_PGF", ";y_trig; y_assoc", 50, -5, 5, 50, -5, 5);
	h_RapidityCorre_QCDC = new TH2D("h_RapidityCorre_QCDC", ";y_trig; y_assoc", 50, -5, 5, 50, -5, 5);
	h_RapidityPair_PGF = new TH1D("h_RapidityPair_PGF", ";y_trig or y_assoc", 50, -5, 5);
	h_RapidityPair_QCDC = new TH1D("h_RapidityPair_QCDC", ";y_trig or y_assoc", 50, -5, 5);
	h_kperp_Pperp_all = new TH2D("h_kperp_Pperp_all", "no kt or Pt cut;k_{#perp} [GeV];P_{#perp} [GeV]", 50, 0, 5, 50, 0, 10);
	h_dphi_kpOverPp = new TH2D("h_dphi_kpOverPp", ";#Delta#phi;k_{#perp}/P_{#perp}", 50, -0.5*3.142, 1.5*3.142, 50, 0, 2);
	h_kperp_Pperp_all_PGF = new TH2D("h_kperp_Pperp_all_PGF", "no kt or Pt cut;k_{#perp} [GeV];P_{#perp} [GeV]", 50, 0, 5, 50, 0, 10);
	h_kperp_Pperp = new TH2D("h_kperp_Pperp", ";k_{#perp} [GeV];P_{#perp} [GeV]", 50, 0, 5, 50, 0, 10);
	h_kperp_Pperp_PGF = new TH2D("h_kperp_Pperp_PGF", ";k_{#perp} [GeV];P_{#perp} [GeV]", 50, 0, 5, 50, 0, 10);
	h_kperp_PGF = new TH1D("h_kperp_PGF", ";k_{#perp} [GeV];", 50, 0, 5);
	h_kperp_QCDC = new TH1D("h_kperp_QCDC", ";k_{#perp} [GeV];", 50, 0, 5);
	h_kperp_q_chan = new TH1D("h_kperp_q_chan", ";k_{#perp} [GeV];", 50, 0, 5);
	h_kperp_g_chan = new TH1D("h_kperp_g_chan", ";k_{#perp} [GeV];", 50, 0, 5);
	h_deltaphi_Pperp = new TH2D("h_deltaphi_Pperp", ";#Delta#phi;P_{#perp}", 50, -0.5*3.142, 1.5*3.142, 50, 0, 10);
	h_deltaphi_Pperp_PGF = new TH2D("h_deltaphi_Pperp_PGF", ";#Delta#phi;P_{#perp}", 50, -0.5*3.142, 1.5*3.142, 50, 0, 10);
	h_phi_pair_1 = new TH1D("h_phi_pair_1", ";#phi;N(#phi)", 10, 0, 2*3.14);
	h_phi_pair_QCDC_1 = new TH1D("h_phi_pair_QCDC_1", ";#phi;N(#phi)", 10, 0, 2*3.14);
	h_phi_pair_PGF_1 = new TH1D("h_phi_pair_PGF_1", ";#phi;N(#phi)", 10, 0, 2*3.14);
	h_phi_pair_q_chan_1 = new TH1D("h_phi_pair_q_chan_1", ";#phi;N(#phi)", 10, 0, 2*3.14);
	h_phi_pair_g_chan_1 = new TH1D("h_phi_pair_g_chan_1", ";#phi;N(#phi)", 10, 0, 2*3.14);
	h_deltaphi_1 = new TH1D("h_deltaphi_1", ";#Delta#phi;N(#Delta#phi)", 50, -0.5*3.142, 1.5*3.142);
	h_deltaphi_Pperp_1 = new TH2D("h_deltaphi_Pperp_1", ";#Delta#phi;P_{#perp}", 50, -0.5*3.142, 1.5*3.142, 50, 0, 10);
	h_deltaphi_PGF_1 = new TH1D("h_deltaphi_PGF_1", ";#Delta#phi;N(#Delta#phi)", 50, -0.5*3.142, 1.5*3.142);
	h_deltaphi_Pperp_PGF_1 = new TH2D("h_deltaphi_Pperp_PGF_1", ";#Delta#phi;P_{#perp}", 50, -0.5*3.142, 1.5*3.142, 50, 0, 10);
	h_phi_pair_2 = new TH1D("h_phi_pair_2", ";#phi;N(#phi)", 10, 0, 2*3.14);
	h_phi_pair_QCDC_2 = new TH1D("h_phi_pair_QCDC_2", ";#phi;N(#phi)", 10, 0, 2*3.14);
	h_phi_pair_PGF_2 = new TH1D("h_phi_pair_PGF_2", ";#phi;N(#phi)", 10, 0, 2*3.14);
	h_phi_pair_q_chan_2 = new TH1D("h_phi_pair_q_chan_2", ";#phi;N(#phi)", 10, 0, 2*3.14);
	h_phi_pair_g_chan_2 = new TH1D("h_phi_pair_g_chan_2", ";#phi;N(#phi)", 10, 0, 2*3.14);
	h_deltaphi_2 = new TH1D("h_deltaphi_2", ";#Delta#phi;N(#Delta#phi)", 50, -0.5*3.142, 1.5*3.142);
	h_deltaphi_Pperp_2 = new TH2D("h_deltaphi_Pperp_2", ";#Delta#phi;P_{#perp}", 50, -0.5*3.142, 1.5*3.142, 50, 0, 10);
	h_deltaphi_PGF_2 = new TH1D("h_deltaphi_PGF_2", ";#Delta#phi;N(#Delta#phi)", 50, -0.5*3.142, 1.5*3.142);
	h_deltaphi_Pperp_PGF_2 = new TH2D("h_deltaphi_Pperp_PGF_2", ";#Delta#phi;P_{#perp}", 50, -0.5*3.142, 1.5*3.142, 50, 0, 10);
	h_bkgVsPperp = new TH2D("h_bkgVsPperp", ";P_{#perp}; bkg num", 50, 0, 10, 50, 0, 40);
	h_bkg_process = new TH1D("h_bkg_process",";process id; counts", 400, 0, 200);

	h_pt_D0_jet = new TH2D("h_pt_D0_jet", ";p_{T}^{D0};p_{T}^{jet}", 50, 0, 10, 50, 0, 10);
	h_dphi_D0_jet = new TH1D("h_dphi_D0_jet", ";#phi_{jet}-#phi_{D0}", 50, -3.142, 3.142);

	//event wise
	h_QSquared = new TH1D("h_QSquared", "Q^{2};log_{10}(Q^{2});counts", 50, 0, 4);
	h_xBj = new TH1D("h_xBj", "x_{Bj};log_{10}(x_{Bj});counts", 50, -5, 0);
	h_W2 = new TH1D("h_W2", "W^{2};log_{10}(W^{2});counts", 100, 0, 6);
	h_y = new TH1D("h_y", "y;log_{10}(y);counts", 50, -5, 0);
	h_Q2VsxBj = new TH2D("h_Q2VsxBj", "log_{10}(Q^{2}) vs log_{10}(x_{Bj});log_{10}(x_{Bj});log_{10}(Q^{2})", 50, -5, 0, 50, 0, 4);
//	h_Q2VsxBj_Pair = new TH2D("h_Q2VsxBj_Pair", "log_{10}(Q^{2}) vs log_{10}(x_{Bj});log_{10}(x_{Bj});log_{10}(Q^{2})", 50, -5, 0, 50, 0, 4);

	h_Q2VsxBj_D0 = getLogXYBin("h_Q2VsxBj_D0", "Q^{2} vs x_{Bj};x_{Bj};Q^{2}", 50, 1E-5, 1, 50, 0.1, 1E3);
	h_Q2VsxBj_Pair = getLogXYBin("h_Q2VsxBj_Pair", "Q^{2} vs x_{Bj};x_{Bj};Q^{2}", 40, 1E-5, 1, 40, 0.1, 1E3);
	h_Q2VsxBj_Pair_PGF = getLogXYBin("h_Q2VsxBj_Pair_PGF", "Q^{2} vs x_{Bj};x_{Bj};Q^{2}", 40, 1E-5, 1, 40, 0.1, 1E3);
	h_Q2VsxBj_Pair_1 = getLogXYBin("h_Q2VsxBj_Pair_1", "Q^{2} vs x_{Bj};x_{Bj};Q^{2}", 40, 1E-5, 1, 40, 0.1, 1E3);
	h_Q2VsxBj_Pair_2 = getLogXYBin("h_Q2VsxBj_Pair_2", "Q^{2} vs x_{Bj};x_{Bj};Q^{2}", 40, 1E-5, 1, 40, 0.1, 1E3);
	h_Q2VsxBj_Pair_xg = getLogXYBin("h_Q2VsxBj_Pair_xg", "Q^{2} vs x_{Bj};x_{Bj};Q^{2}", 40, 1E-5, 1, 40, 0.1, 1E3);
	h_Q2VsxBj_qqbar = getLogXYBin("h_Q2VsxBj_qqbar", "Q^{2} vs x_{Bj};x_{Bj};Q^{2}", 50, 1E-5, 1, 50, 0.1, 1E3);
	h_Parton_PGF = new TH1D("h_Parton_PGF", ";parton id;counts", 100, -23, 23);
	h_xg_Pair = getLogXBin("h_xg_Pair", "x_{g}; x_{g}; counts", 50, 1E-5, 1);

	h_alpha = new TH1D("h_alpha", "alpha_{em};#alpha_{em};counts", 50, 0, 0.01);

	this->status = false;
}

//destructor
Analysis::~Analysis(){

}


//analyse current event and fill the corresponding histograms
//status gives the whether the pair is found or not
void Analysis::FillHist(const Event *event){
	this->status = false;
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
	int trigIndex = -1;

//	int selectPid=421;
//	int selectQuark=4;
	int selectPid=321;
	int selectQuark=3;

	if((event->GetProcess()>132 || event->GetProcess()==84)){
		h_Parton_PGF->Fill(event->GetOutParton1());
		h_Parton_PGF->Fill(event->GetOutParton2());
		if(abs(event->GetOutParton1())==selectQuark)
			h_Q2VsxBj_qqbar->Fill(xBj, Q2);
	}

	//this part is used for the test the tuning parameters with H1 data 1302.1321~~~~~~~~~~~
	if(Q2>5&&Q2<10&&xBj>5E-4&&xBj<2E-3&&(event->GetProcess()>=99||event->GetProcess()<90)) h_Q2_tune->Fill(Q2);

	double foundIndex[100]={0};
	int foundNum = 0;
	int bkg_num = 0;
	//first loop to find all the wanted particles
	for(int itrack=0; itrack<event->GetNTrack(); itrack++){
		ParticleMC *particle = event->GetTrack(itrack);
		//do acceptance cut
//		if(fabs(particle->GetEta())<4.5){
		if(fabs(particle->GetEta())<15){ //no acceptance cut at this stage
			vhadron.SetPxPyPzE(particle->GetPx(), particle->GetPy(), particle->GetPz(), particle->GetE());
//			double z1 = (particle->GetE()-particle->GetPz())/(2*y*Ee_beam);
			double z = (vhadron*v_pbeam)/(v_pbeam*vPhTrue);
			vhadron.Boost(b);
			vhadron.RotateZ(-vPh.Phi());
			vhadron.RotateY(-vPh.Theta());

			//this part is used for the test the tuning parameters with H1 data 1302.1321~~~~~~~~~~~
			if(particle->GetKS()==1&&(abs(particle->GetPid())==211||abs(particle->GetPid())==321||abs(particle->GetPid())==2212||abs(particle->GetPid())==11)){
				if(Q2>5&&Q2<10&&xBj>5E-4&&xBj<2E-3&&(event->GetProcess()>=99||event->GetProcess()<90)){
					if(vhadron.Pt()>0&&vhadron.Pt()<1) h_Prim_dNdEta_1->Fill(vhadron.Eta());
					if(vhadron.Pt()>1&&vhadron.Pt()<10) h_Prim_dNdEta_2->Fill(vhadron.Eta());
					if(vhadron.Eta()>0&&vhadron.Eta()<1.5) h_Prim_dNdPt_1->Fill(vhadron.Pt());
					if(vhadron.Eta()>1.5&&vhadron.Eta()<5) h_Prim_dNdPt_2->Fill(vhadron.Pt());
				}
			}



			if(abs(particle->GetPid())==211||abs(particle->GetPid())==321){
				bkg_num++;
			}

			if(particle->GetKS()==1&&abs(particle->GetPid())==321){
				h_Kaon_eta->Fill(particle->GetEta());
				h_Kaon_pt->Fill(vhadron.Pt());
				h_Kaon_z->Fill(z);
			}
			else if(particle->GetKS()==1&&abs(particle->GetPid())==211){
				h_Pion_eta->Fill(particle->GetEta());
				h_Pion_pt->Fill(vhadron.Pt());
				h_Pion_z->Fill(z);
			}

			if(event->GetProcess()==99){
				if(abs(particle->GetPid())==321){
					h_Kaon_eta_LODIS->Fill(particle->GetEta());
					h_Kaon_pt_LODIS->Fill(vhadron.Pt());
					h_Kaon_z_LODIS->Fill(z);
				}
				else if(abs(particle->GetPid())==211){
					h_Pion_eta_LODIS->Fill(particle->GetEta());
					h_Pion_pt_LODIS->Fill(vhadron.Pt());
					h_Pion_z_LODIS->Fill(z);
				}
			}
			else if(event->GetProcess()>130){
				if(abs(particle->GetPid())==321){
					h_Kaon_eta_Direct->Fill(particle->GetEta());
					h_Kaon_pt_Direct->Fill(vhadron.Pt());
					if(event->GetProcess()>132)
						h_Kaon_pt_PGF->Fill(vhadron.Pt());
					else
						h_Kaon_pt_QCDC->Fill(vhadron.Pt());

					h_Kaon_z_Direct->Fill(z);
				}
				else if(abs(particle->GetPid())==211){
					h_Pion_eta_Direct->Fill(particle->GetEta());
					h_Pion_pt_Direct->Fill(vhadron.Pt());
					h_Pion_z_Direct->Fill(z);
				}
			}
			else if(event->GetProcess()<90){
				if(abs(particle->GetPid())==321){
					h_Kaon_eta_Resolved->Fill(particle->GetEta());
					h_Kaon_pt_Resolved->Fill(vhadron.Pt());
					h_Kaon_z_Resolved->Fill(z);
				}
				else if(abs(particle->GetPid())==211){
					h_Pion_eta_Resolved->Fill(particle->GetEta());
					h_Pion_pt_Resolved->Fill(vhadron.Pt());
					h_Pion_z_Resolved->Fill(z);
				}
			}

			if(event->GetProcess()>132 || event->GetProcess()==84){

				//check the final state particle types in ssbar PGF channel
				if(abs(event->GetOutParton1())==3&&particle->GetKS()==1)
					h_FS_PID_PGF->Fill(particle->GetPid());

				if(abs(particle->GetPid())==321&&abs(event->GetOutParton1())==3){
					h_Kaon_eta_PGFss->Fill(particle->GetEta());
					h_Kaon_pt_PGFss->Fill(vhadron.Pt());
					h_Kaon_z_PGFss->Fill(z);
				}//kaon
				else if(abs(particle->GetPid())==211&&abs(event->GetOutParton1())==3){
					h_Pion_eta_PGFss->Fill(particle->GetEta());
					h_Pion_pt_PGFss->Fill(vhadron.Pt());
					h_Pion_z_PGFss->Fill(z);
				}//pion
				else if(abs(particle->GetPid())==421){
					h_D0_pt->Fill(vhadron.Pt());
					h_D0_eta->Fill(particle->GetEta());
					h_D0_z->Fill(z);
					if(abs(event->GetTrack(particle->GetParent()-1)->GetPid())==413){
						h_D0_pt_Dstar_feed->Fill(vhadron.Pt());
						h_D0_eta_Dstar_feed->Fill(particle->GetEta());
						h_D0_z_Dstar_feed->Fill(z);
					}//D0 from D*

					h_Q2VsxBj_D0->Fill(xBj, Q2);
				}//D0

			}

			//z cut for the particles
//			if(z>0.75||z<0.25) continue;
//			if(z<0.25) continue; //used for D0 pair
//			if(z<0.2) continue; //used for D0 pair
//			if(fabs(particle->GetEta())>1) continue; //used for charged K pair
//			if(z<0.1||particle->GetKS()!=1||abs(particle->GetPid())!=selectPid||vhadron.Pt()<1.7) continue; //used for charged K pair
			if(particle->GetKS()!=1||vhadron.Pt()<1.7) continue; //all primary charge

			int Pidtemp = particle->GetPid(); //all primary charge
			if(abs(Pidtemp)!=211&&abs(Pidtemp)!=321&&abs(Pidtemp)!=2212) continue; //all primary charge
			h_zVsRapidity->Fill(vhadron.Rapidity(), z);
			if(z<0.1) continue; //all primary charge
			h_zVsEta->Fill(particle->GetEta(), z);

//			if(z<0.2||z>0.4) continue; //match with Bowen comparison

			if(fabs(particle->GetEta())>4.5) continue; //all primary charge
			
			//pt cut for the D0 particles
//			if(vhadron.Pt()<0.5||abs(particle->GetPid())!=selectPid) continue;
//			if(abs(particle->GetPid())!=selectPid) continue;

/*
			//make sure only 2 decay products for D0
			if(particle->GetDaughter2()<=particle->GetDaughter1()||particle->GetDaughter2()>particle->GetDaughter1()+1) continue;

			//===========decay products cut==================================//
			int id1 = event->GetTrack(particle->GetDaughter1()-1)->GetPid();
			int id2 = event->GetTrack(particle->GetDaughter2()-1)->GetPid();
			double eta1 = event->GetTrack(particle->GetDaughter1()-1)->GetEta();
			double eta2 = event->GetTrack(particle->GetDaughter2()-1)->GetEta();
			double mom1 = event->GetTrack(particle->GetDaughter1()-1)->GetP();
			double mom2 = event->GetTrack(particle->GetDaughter2()-1)->GetP();
			double pt1 = event->GetTrack(particle->GetDaughter1()-1)->GetPt();
			double pt2 = event->GetTrack(particle->GetDaughter2()-1)->GetPt();

			if(abs(id1)!=211&&abs(id1)!=321)
				continue;
			else if(abs(id1)==211&&abs(id2)!=321)
				continue;
			else if(abs(id1)==321&&abs(id2)!=211)
				continue;

			bool kaon_mom_cut=false;

			//1st pion, 2nd kaon
			if(abs(id1)==211){
				h_D0_decay_pi->Fill(eta1, mom1);
				h_D0_decay_kaon->Fill(eta2, mom2);
				//if(mom2<0.2) kaon_mom_cut=true;
				if(pt1<0.2||pt2<0.2) kaon_mom_cut=true; //do mom cut for both pi and k
			}
			//1st kaon, 2nd pion
			else{
				h_D0_decay_pi->Fill(eta2, mom2);
				h_D0_decay_kaon->Fill(eta1, mom1);
				//if(mom1<0.2) kaon_mom_cut=true;
				if(pt1<0.2||pt2<0.2) kaon_mom_cut=true; //do mom cut for both pi and k
			}

			//do kaon decay product cut
			if(kaon_mom_cut) continue;

			if(fabs(eta1)>4.5||fabs(eta2)>4.5) continue;
			//===========decay products cut end===============================//
*/

//			std::cout<<"id1:"<<id1<<" id2:"<<id2<<std::endl;

			foundIndex[foundNum]=itrack;
			foundNum++;
			//look for trig
/*			if( trigId==0 && vhadron.Pt()>1 ){
				trigId = particle->GetPid();
				vTrig=vhadron;
				trigIndex = itrack;
			}//if trig
*/
		}//if accepted
	}//for

	//make pairs according to the particle list accepted above
	for(int itrig=0; itrig<foundNum-1; itrig++){
		ParticleMC *particle1 = event->GetTrack(foundIndex[itrig]);
		vTrig.SetPxPyPzE(particle1->GetPx(), particle1->GetPy(), particle1->GetPz(), particle1->GetE());

		vTrig.Boost(b);
		vTrig.RotateZ(-vPh.Phi());
		vTrig.RotateY(-vPh.Theta());

		for(int iasso=itrig+1; iasso<foundNum; iasso++){
			ParticleMC *particle2 = event->GetTrack(foundIndex[iasso]);
			vhadron.SetPxPyPzE(particle2->GetPx(), particle2->GetPy(), particle2->GetPz(), particle2->GetE());
//			if(particle1->GetPid()==-particle2->GetPid()){
			if(true){ //no antiparticle-particle correlation: for all primary charge
				vhadron.Boost(b);
				vhadron.RotateZ(-vPh.Phi());
				vhadron.RotateY(-vPh.Theta());

				h_RapidityCorre->Fill(vTrig.Rapidity(), vhadron.Rapidity());

				vTotal = vTrig+vhadron;
				double kt_phi = vTotal.Phi();
				double delta_phi = vTrig.DeltaPhi(vhadron);
				double kperp = vTotal.Pt();
				double Pperp = (vTrig-vhadron).Pt()/2.;

				//test with Bowen's calculation
//				if(Pperp<4.) continue;

				if(delta_phi < -3.14159/2.) delta_phi+=2*3.14159;
				if(delta_phi > 3*3.14159/2.) delta_phi-=2*3.14159;
				if(kt_phi<0) kt_phi+=2*3.14159;

				h_kperp_Pperp_all->Fill(kperp, Pperp);
				h_dphi_kpOverPp->Fill(delta_phi, kperp/Pperp);
				if(event->GetProcess()>132) h_kperp_Pperp_all_PGF->Fill(kperp, Pperp);
//				std::cout<<"found 1 pair, kperp:"<<kperp<<" Pperp:"<<Pperp<<std::endl;
				//correlation limit
//				if(kperp/Pperp<0.4){
				if(kperp/Pperp<0.7){
					std::cout<<"accept 1 pair"<<std::endl;

					h_Pair_Eta->Fill(particle1->GetEta());
					h_Pair_Eta->Fill(particle2->GetEta());

					if(event->GetProcess()>132){
						h_RapidityCorre_PGF->Fill(vTrig.Rapidity(), vhadron.Rapidity());
						h_RapidityPair_PGF->Fill(vTrig.Rapidity());
						h_RapidityPair_PGF->Fill(vhadron.Rapidity());
					}
					else if(event->GetProcess()>130){
						h_RapidityCorre_QCDC->Fill(vTrig.Rapidity(), vhadron.Rapidity());
						h_RapidityPair_QCDC->Fill(vTrig.Rapidity());
						h_RapidityPair_QCDC->Fill(vhadron.Rapidity());
					}

					int parentId1 = event->GetTrack(particle1->GetParent()-1)->GetPid();
					int parentId2 = event->GetTrack(particle2->GetParent()-1)->GetPid();
					if(parentId1==91||parentId1==92){
						ParticleMC *jet = event->GetTrack(event->GetTrack(particle1->GetParent()-1)->GetParent()-1);
						if(abs(jet->GetPid())<10){
							//found 1 mother jet
							h_pt_D0_jet->Fill(jet->GetPt(), particle1->GetPt());
							h_dphi_D0_jet->Fill(jet->GetPhi() - particle1->GetPhi());
						}
					}
					if(parentId2==91||parentId2==92){
						ParticleMC *jet = event->GetTrack(event->GetTrack(particle2->GetParent()-1)->GetParent()-1);
						if(abs(jet->GetPid())<10){
							h_pt_D0_jet->Fill(jet->GetPt(), particle2->GetPt());
							h_dphi_D0_jet->Fill(jet->GetPhi() - particle2->GetPhi());
						}
					}

					h_bkgVsPperp->Fill(Pperp, bkg_num);

					h_kperp_Pperp->Fill(kperp, Pperp);
					h_deltaphi_Pperp->Fill(delta_phi, Pperp);

					if(event->GetProcess()>130&&event->GetProcess()<133) h_kperp_QCDC->Fill(kperp);
					if(event->GetProcess()>133) h_kperp_PGF->Fill(kperp);
					if(event->GetTargetParton()==21&&(event->GetProcess()<90||event->GetProcess()>98)) h_kperp_g_chan->Fill(kperp);
					if(abs(event->GetTargetParton())<10&&(event->GetProcess()<90||event->GetProcess()>98)) h_kperp_q_chan->Fill(kperp);


					if(kperp<1&&kperp>0.5){
						h_phi_pair_1->Fill(kt_phi);
						if(event->GetProcess()>130&&event->GetProcess()<133) h_phi_pair_QCDC_1->Fill(kt_phi);
						if(event->GetTargetParton()==21&&(event->GetProcess()<90||event->GetProcess()>98)) h_phi_pair_g_chan_1->Fill(kt_phi);
						if(abs(event->GetTargetParton())<10&&(event->GetProcess()<90||event->GetProcess()>98)) h_phi_pair_q_chan_1->Fill(kt_phi);
						h_deltaphi_1->Fill(delta_phi);
						h_Q2VsxBj_Pair_1->Fill(xBj, Q2);
						h_deltaphi_Pperp_1->Fill(delta_phi, Pperp);
					}
					else if(kperp>1&&kperp<2){
						if(event->GetProcess()>130&&event->GetProcess()<133) h_phi_pair_QCDC_2->Fill(kt_phi);
						if(event->GetTargetParton()==21&&(event->GetProcess()<90||event->GetProcess()>98)) h_phi_pair_g_chan_2->Fill(kt_phi);
						if(abs(event->GetTargetParton())<10&&(event->GetProcess()<90||event->GetProcess()>98)) h_phi_pair_q_chan_2->Fill(kt_phi);
	
						h_phi_pair_2->Fill(kt_phi);
						h_deltaphi_2->Fill(delta_phi);
						h_Q2VsxBj_Pair_2->Fill(xBj, Q2);
						h_deltaphi_Pperp_2->Fill(delta_phi, Pperp);
					}

			//		if((event->GetProcess()>132||event->GetProcess()==84) && abs(event->GetOutParton1())==selectQuark){
					if((event->GetProcess()>132||event->GetProcess()==84)){
						h_kperp_Pperp_PGF->Fill(kperp, Pperp);
						h_deltaphi_Pperp_PGF->Fill(delta_phi, Pperp);
						h_Q2VsxBj_Pair_PGF->Fill(xBj, Q2);
						if(kperp<1&&kperp>0.5){
							h_phi_pair_PGF_1->Fill(kt_phi);
							h_deltaphi_PGF_1->Fill(delta_phi);
							h_deltaphi_Pperp_PGF_1->Fill(delta_phi, Pperp);
						}
						else if(kperp>1&&kperp<2){
							h_phi_pair_PGF_2->Fill(kt_phi);
							h_deltaphi_PGF_2->Fill(delta_phi);
							h_deltaphi_Pperp_PGF_2->Fill(delta_phi, Pperp);
						}
					}
					else
						h_bkg_process->Fill(event->GetProcess());

					h_Q2VsxBj_Pair->Fill(xBj, Q2);
					h_Q2VsxBj_Pair_xg->Fill(xBj, Q2, event->GetTargetPartonX());
					h_xg_Pair->Fill(event->GetTargetPartonX());
					this->status = true;//this event found one pair
				}//kperp/Pperp cut
			}//if assoc id== trig id
		}//for asso loop
	}//for trig loop
	
//	std::cout<<num_Kaon<<std::endl;

}


void Analysis::WriteResults(TFile *file){
	TDirectory *save = gDirectory;
	file->cd();

	double total_cs = pypars.pari[0]; //in mb
//	double total_nevt = pypars.msti[4]; //number of gen events
	double total_ntrials = pyint5.ngen[2][0]; //number of gen events
	double lumi = total_ntrials/total_cs; //in mb^-1
	//1E6 added to convert mb to nb
//	h_Kaon_eta->Scale(1E6/lumi,"width"); //dsigma/dX in nb
//	h_Kaon_pt->Scale(1E6/lumi,"width"); //dsigma/dX in nb
//	h_Kaon_z->Scale(1E6/lumi,"width"); //dsigma/dX in nb
//	h_Pion_eta->Scale(1E6/lumi,"width"); //dsigma/dX in nb
//	h_Pion_pt->Scale(1E6/lumi,"width"); //dsigma/dX in nb
//	h_Pion_z->Scale(1E6/lumi,"width"); //dsigma/dX in nb
	h_D0_eta->Scale(1E6/lumi,"width"); //dsigma/dX in nb
	h_D0_pt->Scale(1E6/lumi,"width"); //dsigma/dX in nb
	h_D0_z->Scale(1E6/lumi,"width"); //dsigma/dX in nb
	h_D0_eta_Dstar_feed->Scale(1E6/lumi,"width"); //dsigma/dX in nb
	h_D0_pt_Dstar_feed->Scale(1E6/lumi,"width"); //dsigma/dX in nb
	h_D0_z_Dstar_feed->Scale(1E6/lumi,"width"); //dsigma/dX in nb

	//this part is used for the test the tuning parameters with H1 data 1302.1321~~~~~~~~~~~
	h_Prim_dNdEta_1->Scale(1./h_Q2_tune->GetEntries(), "width");
	h_Prim_dNdEta_2->Scale(1./h_Q2_tune->GetEntries(), "width");
	h_Prim_dNdPt_1->Scale(1./h_Q2_tune->GetEntries(), "width");
	h_Prim_dNdPt_2->Scale(1./h_Q2_tune->GetEntries(), "width");

	//tune test
	h_Prim_dNdEta_1->Write();
	h_Prim_dNdEta_2->Write();
	h_Prim_dNdPt_1->Write();
	h_Prim_dNdPt_2->Write();

	//track wise
	h_theta_scattered->Write();
	h_FS_PID_PGF->Write();

	h_zVsRapidity->Write();
	h_zVsEta->Write();

	h_Kaon_eta->Write();
	h_Kaon_pt->Write();
	h_Kaon_z->Write();
	h_Pion_eta->Write();
	h_Pion_pt->Write();
	h_Pion_z->Write();

	h_Kaon_eta_LODIS->Write();
	h_Kaon_pt_LODIS->Write();
	h_Kaon_z_LODIS->Write();
	h_Pion_eta_LODIS->Write();
	h_Pion_pt_LODIS->Write();
	h_Pion_z_LODIS->Write();

	h_Kaon_eta_Direct->Write();
	h_Kaon_pt_Direct->Write();
	h_Kaon_pt_PGF->Write();
	h_Kaon_pt_QCDC->Write();
	h_Kaon_z_Direct->Write();
	h_Pion_eta_Direct->Write();
	h_Pion_pt_Direct->Write();
	h_Pion_z_Direct->Write();

	h_Kaon_eta_Resolved->Write();
	h_Kaon_pt_Resolved->Write();
	h_Kaon_z_Resolved->Write();
	h_Pion_eta_Resolved->Write();
	h_Pion_pt_Resolved->Write();
	h_Pion_z_Resolved->Write();

	h_Kaon_eta_PGFss->Write();
	h_Kaon_pt_PGFss->Write();
	h_Kaon_z_PGFss->Write();
	h_Pion_eta_PGFss->Write();
	h_Pion_pt_PGFss->Write();
	h_Pion_z_PGFss->Write();

	h_D0_eta->Write();
	h_D0_pt->Write();
	h_D0_z->Write();
	h_D0_eta_Dstar_feed->Write();
	h_D0_pt_Dstar_feed->Write();
	h_D0_z_Dstar_feed->Write();

	h_D0_decay_pi->Write();
	h_D0_decay_kaon->Write();

	//pair wise
	h_Pair_Eta->Write();
	h_RapidityCorre->Write();
	h_RapidityCorre_PGF->Write();
	h_RapidityCorre_QCDC->Write();
	h_RapidityPair_PGF->Write();
	h_RapidityPair_QCDC->Write();
	h_dphi_kpOverPp->Write();
	h_kperp_Pperp_all->Write();
	h_kperp_Pperp_all_PGF->Write();
	h_kperp_Pperp->Write();
	h_kperp_Pperp_PGF->Write();
	h_kperp_PGF->Write();
	h_kperp_QCDC->Write();
	h_kperp_q_chan->Write();
	h_kperp_g_chan->Write();
	h_deltaphi_Pperp->Write();
	h_deltaphi_Pperp_PGF->Write();
	h_phi_pair_1->Write();
	h_phi_pair_QCDC_1->Write();
	h_phi_pair_PGF_1->Write();
	h_phi_pair_g_chan_1->Write();
	h_phi_pair_q_chan_1->Write();
	h_deltaphi_1->Write();
	h_deltaphi_Pperp_1->Write();
	h_deltaphi_PGF_1->Write();
	h_deltaphi_Pperp_PGF_1->Write();
	h_phi_pair_2->Write();
	h_phi_pair_QCDC_2->Write();
	h_phi_pair_PGF_2->Write();
	h_phi_pair_g_chan_2->Write();
	h_phi_pair_q_chan_2->Write();
	h_deltaphi_2->Write();
	h_deltaphi_Pperp_2->Write();
	h_deltaphi_PGF_2->Write();
	h_deltaphi_Pperp_PGF_2->Write();
	h_bkgVsPperp->Write();
	h_bkg_process->Write();

	//event wise
	h_QSquared->Write();
	h_xBj->Write();
	h_W2->Write();
	h_y->Write();
	h_Q2VsxBj->Write();
	h_Q2VsxBj_Pair->Write();
	h_Q2VsxBj_Pair_PGF->Write();
	h_Q2VsxBj_Pair_1->Write();
	h_Q2VsxBj_Pair_2->Write();
	h_Q2VsxBj_Pair_xg->Write();
	h_Q2VsxBj_D0->Write();
	h_Q2VsxBj_qqbar->Write();
	h_Parton_PGF->Write();
	h_alpha->Write();
	h_xg_Pair->Write();
	h_pt_D0_jet->Write();
	h_dphi_D0_jet->Write();

	//D0: 421, D*: 413
	std::cout<<"421 comp code:"<<call_pycomp(421)<<std::endl;

	//cd back to original dir
	save->cd();
}
