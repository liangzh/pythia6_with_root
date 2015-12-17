#ifndef Analysis_H
#define Analysis_H

#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <Rtypes.h>
#include <TObject.h>

#include "Event.h"

class Analysis {

private:
	//track wise
	TH1D *h_theta_scattered;
	TH1D *h_Dstar_pt;
	TH1D *h_Dstar_eta;
	TH1D *h_Dstar_z;
	
	//event wise
	TH1D *h_QSquared;
	TH1D *h_xBj;
	TH1D *h_W2;
	TH1D *h_y;
	TH2D *h_Q2VsxBj;
	TH2D *h_Q2VsxBj_Dstar;
	TH1D *h_alpha;

	//inclusive cross section
	TH1D *h_dsig[7];
	//reduced cross section
	TH1D *h_sig_reduced_cc[7];
	double Q2_cut[14]={2.4, 2.6, 
										4.5, 5.5, 
										5, 9,
										10, 14,
										15, 21,
										30, 34,
										180, 220};

public:
	Analysis();
	virtual ~Analysis();

	// if this analysis is accepted
	bool status;

	/**
	 * fill the current event info to some QA hist
	 */
	void FillHist(const Event *event);

	/**
	 * write QA hist to the output root file
	 */
	void WriteResults(TFile *file);

	ClassDef(Analysis,1);
};

#endif
