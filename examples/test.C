#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include "/home/adminuser/Software/pythia_with_root/pythia/include/Event.h"

int test() {

	cout<<"This is in the test module"<<endl;
	gSystem->Load("/home/adminuser/Software/pythia_with_root/pythia/lib/libEvent.so");

	TFile *file = new TFile("../pythia.root");
	TTree *tree = (TTree*)file->Get("mcTree");

	Event *event = new Event();
	tree->SetBranchAddress("event", &event);

	if(tree->GetEntries()<=0) {
		cout<< "file not opened!!" <<endl;
		return 1;
	}

	tree->Draw("fPt");

	return 0;

}
