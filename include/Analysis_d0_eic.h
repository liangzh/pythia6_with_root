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
	TH1D *h_D0_pt;
	TH1D *h_D0_eta;
	TH1D *h_D0_z;
	//pair wise
	TH1D *h_phi_pair;
	TH1D *h_phi_pair_PGF;
	TH1D *h_deltaphi;
	TH1D *h_deltaphi_PGF;
	
	//event wise
	TH1D *h_QSquared;
	TH1D *h_xBj;
	TH1D *h_W2;
	TH1D *h_y;
	TH2D *h_Q2VsxBj;
	TH2D *h_Q2VsxBj_D0;
	TH1D *h_alpha;


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
