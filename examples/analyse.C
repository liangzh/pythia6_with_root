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
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>

//used to read a event tree
#include "/home/adminuser/Software/pythia_with_root/pythia/include/Event.h"
#include "/home/adminuser/Software/pythia_with_root/pythia/include/HistMaker.h"

#endif

using std::cout;
using std::endl;

void getFmoments(const TH1D*, const int&, TH1D*);
double getProduct(const int&, const double&);
double getFq(const TH1D*, const TH1D*, const double&);



void analyse(){

	//make the chain
	TChain *tree = new TChain("mcTree");
	//read the data file
//	tree->Add("../pythia_BES.root");
	tree->Add("../pythia_ee_200GeV.root");

	if(tree->GetEntries()<=0) cout<<" File not opened!"<<endl;

	//assign the branch address
	Event *event = new Event();
	tree->SetBranchAddress("event", &event);


	TH1::SetDefaultSumw2();//to get the bin error saved
	//book hist
	int Npoints=20;
	TH1D *h_final_eta_per_event[Npoints];
	TH1D *h_final_eta[Npoints];
	TH1D *h_product_per_event[Npoints];
	TH1D *h_product[Npoints];

	TH1D *h_moments[Npoints];

	int power_start=0;
	for(int i(0); i<Npoints; i++){
//		int Nbins=pow(2, i+power_start);
		int Nbins=i+1;
		h_final_eta_per_event[i] = new TH1D(Form("h_final_eta_per_event%d",i), "final #eta; #eta; counts", Nbins, -5, 5);
		h_product_per_event[i] = new TH1D(Form("h_product_per_event%d",i), "final #eta; #eta; counts", Nbins, -5, 5);

		h_final_eta[i] = new TH1D(Form("h_final_eta%d",i), "final #eta; #eta; counts", Nbins, -5, 5);
		h_product[i] = new TH1D(Form("h_product%d",i), "final #eta; #eta; counts", Nbins, -5, 5);
	}

	//kinematic cut
	double etaCut=5; //acceptance cut

	//moment index
	int q_index = 2; //F moments index q

	long iLoop = 0; //event looper
	ParticleMC *particle;//temp particle looper

	while( tree->GetEntry(iLoop)>0 ){
		iLoop++;
		if(iLoop%10000==0) cout<<iLoop<<" done"<<endl;
		
		//particle level loop
		for(int itrack(0); itrack<event->GetNTrack(); itrack++){
			particle = event->GetTrack(itrack);
			//analyse the final state particles
			if(particle->GetKS()==1&&fabs(particle->GetEta())<etaCut&&particle->GetPt()<3){
				//exclude n because it can not be observed in the detector
//				if(abs(particle->GetPid())!=2112){
			
					for(int iset(0); iset<Npoints; iset++)
						h_final_eta_per_event[iset]->Fill(particle->GetRapidity());
//				}
			}
		}//for

		//get the F moments for different bin size
		for(int ipoint(0); ipoint<Npoints; ipoint++){
			//now get factorial moments according to the one dimensional distribution
			getFmoments(h_final_eta_per_event[ipoint], q_index, h_product_per_event[ipoint]);

			h_product[ipoint]->Add(h_product_per_event[ipoint]);
			h_final_eta[ipoint]->Add(h_final_eta_per_event[ipoint]);

			h_product_per_event[ipoint]->Reset();
			h_final_eta_per_event[ipoint]->Reset();
		}

	}//while

	//save the final result in a graph
	TGraphErrors *gr = new TGraphErrors();
	//average the hist
	for(int ipoint(0); ipoint<Npoints; ipoint++){
		//average for the events
		h_product[ipoint]->Scale(1./tree->GetEntries());
		h_final_eta[ipoint]->Scale(1./tree->GetEntries());

		//clone the product content to the h_moments
		h_moments[ipoint] = (TH1D*)h_product[ipoint]->Clone(Form("h_moment%d",ipoint));

		//now get the Fq for test(by hand)
		double Fq = getFq(h_product[ipoint], h_final_eta[ipoint], q_index);
		for(int i_index(0); i_index<q_index; i_index++){
			h_moments[ipoint]->Divide(h_final_eta[ipoint]);
		}
		cout<<"index="<<pow(2, ipoint+power_start)<<" Fq="<<Fq<<endl;
		
		//sum over the kinematic bin to get Fq and its error
		double error = -1;
		double value = h_moments[ipoint]->IntegralAndError(-1,-1,error)/h_moments[ipoint]->GetNbinsX();

		//set the value and error to the graph
		gr->SetPoint(ipoint, h_moments[ipoint]->GetNbinsX(), value);
		gr->SetPointError(ipoint, 0, error);
	}

	TCanvas *c1 = new TCanvas();
	gr->Draw("apl");


	//save the hist to an root file
/*	TFile *outFile = new TFile("analyse_hist.root","recreate");
	gr->Write();
	for(int ipoint(0); ipoint<Npoints; ipoint++){
		h_final_eta[ipoint]->Write();
		h_product[ipoint]->Write();
		h_moments[ipoint]->Write();
	}
*/
}


//get the hout with bin content set by the product of the moments
//h_in: the input hist with multiplicity distribution
//q: the q index
//h_out: the output hist with the content set by the products
void getFmoments(const TH1D *h_in, const int &q, TH1D *h_out){
	for(int ibin(0); ibin<h_in->GetSize(); ibin++){
		double content = getProduct(q, h_in->GetBinContent(ibin));
		h_out->SetBinContent(ibin, content);
	}
}

//get the q-th product of a count
//q: the q index
//content: the bin content of a certain multplicity
double getProduct(const int &q, const double &content){
	double prod = content;
	for(int i_q(1); i_q<q; i_q++){
		prod*=(content - i_q);
	}

	return prod;
}


//the test routine to calculate Fq
double getFq(const TH1D *h_moment, const TH1D *h_aver, const double &q){

	double Fq = 0;
	for(int ibin(0); ibin<h_moment->GetSize(); ibin++){
		if(h_aver->GetBinContent(ibin)>0)
			Fq += h_moment->GetBinContent(ibin)/pow(h_aver->GetBinContent(ibin), q);
	}

	return Fq/h_aver->GetSize();
}
