#ifndef HistMaker_H
#define HistMaker_H

#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <Rtypes.h>
#include <TObject.h>

//#include "Event.h"

class Event;

class HistMaker {

private:
	//track wise
	TH1D *h_FS_pt;
	TH1D *h_FS_eta;
	TH1D *h_FS_rapidity;
	TH1D *h_FS_pt_Pi;
	TH1D *h_FS_pt_K;
	TH1D *h_FS_pt_P;
	TH1D *h_FS_PID;
	//event wise
	TH1D *h_FS_num;
	TH1D *h_process;
	TH1D *h_pt2hat;

public:
	HistMaker();
	virtual ~HistMaker();

	/**
	 * fill the current event info to some QA hist
	 */
	void Hfill(const Event *event);

	/**
	 * write QA hist to the output root file
	 */
	void Hwrite(TFile *file);

	ClassDef(HistMaker,1)
};

#endif
