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
	//track wise
	TH1D *h_theta_scattered;
	TH1D *h_FS_PID_PGF;
	TH1D *h_Kaon_pt;
	TH1D *h_Kaon_eta;
	TH1D *h_Kaon_z;
	TH1D *h_Pion_pt;
	TH1D *h_Pion_eta;
	TH1D *h_Pion_z;
	TH1D *h_D0_pt;
	TH1D *h_D0_eta;
	TH1D *h_D0_z;
	//pair wise
	TH2D *h_kperp_Pperp;
	TH2D *h_kperp_Pperp_PGF;
	TH2D *h_deltaphi_Pperp;
	TH2D *h_deltaphi_Pperp_PGF;
	TH1D *h_phi_pair_1;
	TH1D *h_phi_pair_PGF_1;
	TH1D *h_deltaphi_1;
	TH2D *h_deltaphi_Pperp_1;
	TH1D *h_deltaphi_PGF_1;
	TH2D *h_deltaphi_Pperp_PGF_1;
	TH1D *h_phi_pair_2;
	TH1D *h_phi_pair_PGF_2;
	TH1D *h_deltaphi_2;
	TH2D *h_deltaphi_Pperp_2;
	TH1D *h_deltaphi_PGF_2;
	TH2D *h_deltaphi_Pperp_PGF_2;
	TH2D *h_bkgVsPperp;
	TH1D *h_bkg_process;
		
	//event wise
	TH1D *h_QSquared;
	TH1D *h_xBj;
	TH1D *h_W2;
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
