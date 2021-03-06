#ifndef Analysis_H
#define Analysis_H

#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <Rtypes.h>
#include <TObject.h>

//#include "Event.h"

class Event;

class Analysis {

private:
	//tune test
	TH1D *h_Prim_dNdEta_1;
	TH1D *h_Prim_dNdEta_2;
	TH1D *h_Prim_dNdPt_1;
	TH1D *h_Prim_dNdPt_2;
	TH1D *h_Q2_tune;

	//track wise
	TH1D *h_theta_scattered;
	TH1D *h_FS_PID_PGF;
	TH1D *h_Kaon_pt;
	TH1D *h_Kaon_eta;
	TH1D *h_Kaon_z;
	TH1D *h_Pion_pt;
	TH1D *h_Pion_eta;
	TH1D *h_Pion_z;
	TH1D *h_charge_pt;
	TH1D *h_charge_eta;
	TH1D *h_charge_z;


	TH2D *h_zVsRapidity;
	TH2D *h_zVsEta;

	TH1D *h_Kaon_pt_LODIS;
	TH1D *h_Kaon_eta_LODIS;
	TH1D *h_Kaon_z_LODIS;
	TH1D *h_Pion_pt_LODIS;
	TH1D *h_Pion_eta_LODIS;
	TH1D *h_Pion_z_LODIS;
	TH1D *h_charge_pt_LODIS;
	TH1D *h_charge_eta_LODIS;
	TH1D *h_charge_z_LODIS;


	TH1D *h_Kaon_pt_Direct;
	TH1D *h_Kaon_pt_PGF;
	TH1D *h_Kaon_pt_QCDC;
	TH1D *h_Kaon_eta_Direct;
	TH1D *h_Kaon_z_Direct;
	TH1D *h_Pion_pt_Direct;
	TH1D *h_Pion_eta_Direct;
	TH1D *h_Pion_z_Direct;
	TH1D *h_charge_pt_Direct;
	TH1D *h_charge_pt_PGF;
	TH1D *h_charge_pt_QCDC;
	TH1D *h_charge_eta_Direct;
	TH1D *h_charge_z_Direct;


	TH1D *h_Kaon_pt_Resolved;
	TH1D *h_Kaon_eta_Resolved;
	TH1D *h_Kaon_z_Resolved;
	TH1D *h_Pion_pt_Resolved;
	TH1D *h_Pion_eta_Resolved;
	TH1D *h_Pion_z_Resolved;
	TH1D *h_charge_pt_Resolved;
	TH1D *h_charge_eta_Resolved;
	TH1D *h_charge_z_Resolved;


	TH1D *h_Kaon_pt_PGFss;
	TH1D *h_Kaon_eta_PGFss;
	TH1D *h_Kaon_z_PGFss;
	TH1D *h_Pion_pt_PGFss;
	TH1D *h_Pion_eta_PGFss;
	TH1D *h_Pion_z_PGFss;

	TH1D *h_D0_pt;
	TH1D *h_D0_eta;
	TH1D *h_D0_z;
  TH1D *h_D0_pt_Dstar_feed;
	TH1D *h_D0_eta_Dstar_feed;
	TH1D *h_D0_z_Dstar_feed;
	TH2D *h_D0_decay_pi;
	TH2D *h_D0_decay_kaon;
	TH2D *h_D0_pipt_kaonpt;
	//pair wise
	TH1D *h_Pair_Eta;
	TH2D *h_RapidityCorre;
	TH2D *h_EtaCorre;
	TH2D *h_RapidityCorre_PGF;
	TH2D *h_RapidityCorre_QCDC;
	TH1D *h_RapidityPair_PGF;
	TH1D *h_RapidityPair_QCDC;
	TH2D *h_kperp_Pperp_all;
	TH2D *h_dphi_kpOverPp;
	TH2D *h_kperp_Pperp_all_PGF;
	TH2D *h_kperp_Pperp;
	TH2D *h_kperp_Pperp_PGF;
	TH1D *h_kperp_PGF;
	TH1D *h_kperp_QCDC;
	TH1D *h_kperp_q_chan;
	TH1D *h_kperp_g_chan;
	TH2D *h_deltaphi_Pperp;
	TH2D *h_deltaphi_Pperp_PGF;
	TH1D *h_phi_pair_1;
	TH1D *h_phi_pair_QCDC_1;
	TH1D *h_phi_pair_PGF_1;
	TH1D *h_phi_pair_q_chan_1;
	TH1D *h_phi_pair_g_chan_1;
	TH1D *h_deltaphi_1;
	TH2D *h_deltaphi_Pperp_1;
	TH1D *h_deltaphi_PGF_1;
	TH2D *h_deltaphi_Pperp_PGF_1;
	TH1D *h_phi_pair_2;
	TH1D *h_phi_pair_QCDC_2;
	TH1D *h_phi_pair_PGF_2;
	TH1D *h_phi_pair_q_chan_2;
	TH1D *h_phi_pair_g_chan_2;
	TH1D *h_deltaphi_2;
	TH2D *h_deltaphi_Pperp_2;
	TH1D *h_deltaphi_PGF_2;
	TH2D *h_deltaphi_Pperp_PGF_2;
	TH2D *h_bkgVsPperp;
	TH1D *h_bkg_process;

	TH2D *h_pt_D0_jet;
	TH1D *h_dphi_D0_jet;
		
	//event wise
	TH2D *h_oute_pVsEta;
	TH1D *h_QSquared;
	TH1D *h_xBj;
	TH1D *h_W2;
	TH1D *h_W;
	TH1D *h_y;
	TH2D *h_Q2VsxBj;
	TH1D *h_alpha;
	TH2D *h_Q2VsxBj_Pair;
	TH2D *h_Q2VsxBj_Pair_PGF;
	TH2D *h_Q2VsxBj_Pair_1;
	TH2D *h_Q2VsxBj_Pair_2;
	TH2D *h_Q2VsxBj_Pair_xg;
	TH2D *h_Q2VsxBj_D0;
	TH2D *h_Q2VsxBj_qqbar;
	TH1D *h_Parton_PGF;
	TH1D *h_xg_Pair;

	// if this analysis is accepted
	bool status;


public:
	Analysis();
	virtual ~Analysis();


	/**
	 * fill the current event info to some QA hist
	 */
	void FillHist(const Event *event);

	/**
	 * write QA hist to the output root file
	 */
	void WriteResults(TFile *file);

	/**
	 * get status for the current analysis
	 */
	bool GetStatus() { return status; }


	ClassDef(Analysis,1)
};

#endif
