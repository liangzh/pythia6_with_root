/* This is a example showing how to read the
 * generated root files. Accessible functions
 * or data variables can be found in the
 * include/Event.h decalaration.
 */


//wrap up the header files for compiler
//but hiding from the CINT

#ifndef __CINT__

#include <iostream>

#include <TSystem.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

//used to read a event tree
#include "/home/adminuser/Software/pythia_with_root/pythia/include/Event.h"
#include "/home/adminuser/Software/pythia_with_root/pythia/include/HistMaker.h"



#endif

using std::cout;
using std::endl;


void read(){

	//make the chain
	TChain *tree = new TChain("mcTree");
	tree->Add("../pythia.root");

	if(tree->GetEntries()<=0) cout<<" File not opened!"<<endl;

	//assign the branch address
	Event *event = new Event();
	tree->SetBranchAddress("event", &event);

	int iLoop = 0;

	//temp particle looper
	ParticleMC *particle;
	while( tree->GetEntry(iLoop)>0 ){
		iLoop++;
		if(iLoop%100==0) cout<<iLoop<<" done"<<endl;
		for(int itrack(0); itrack<event->GetNTrack(); itrack++){
			particle = event->GetTrack(itrack);
			//analyse the final state particles here
			if(particle->GetKS()==1){
				cout << "id=" << particle->GetPid() << " KS=" << particle->GetKS() << " orig=" << particle->GetParent() << " daughter1=" << particle->GetDaughter1()
					<< " daughter2=" << particle->GetDaughter2() << " E=" << particle->GetE() << endl;
			}
		}//for
	}//while
}
