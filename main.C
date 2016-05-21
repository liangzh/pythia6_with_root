/** 
 * This is the main program coordinating
 * I/O routines, event generation process
 * and preselector functions for the saved
 * events.
 */

#include <iostream>
#include <sstream>
#include <string>

// GNU long option variation of getopt, for command line option parsing.
#include <getopt.h>

#include <TSystem.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>

//used to wrap up fortran routines used
//in PYTHIA
#include "pythiaWrapper.h"

//used to create a event tree
#include "Event.h"
//used to create some QA plots
#include "HistMaker.h"
//used to do analysis
#include "Analysis.h"

using std::cout;
using std::endl;


// Forward declarations
std::string initialise(int argc, char* argv[]); //process command line arguments
void write_run_data(TFile&); //write cross section to output root file

//preselector to choose which events
//will be written to the output root file
bool selector(Event *e);
static int doSelect = 0; //<selector switch---1:on; 0:off


int main(int argc, char** argv){

	//switch to select those events to be saved on disk 
	//if 1, do preselection; if 0, write all the
	//events to the disk
	doSelect = 0; //disable preselection by default

	std::string outfile = initialise(argc, argv);

	int nEvent = -1;
	int nPrint = -1;//step number to print when running
	runpyinit(&nEvent, &nPrint);

	cout<<"outfile:"<<outfile<<endl;

	cout<<"============================================="<<endl;


	if(nEvent<0) {
		cout << "Invalid event number!" << endl;
		return 1;
	}

	//make the output tree file
	TFile file(outfile.c_str(),"recreate");
	TTree tree("mcTree","a pythia MC tree");

	//make hist manager, save some QA plots
	HistMaker hist;
	Analysis worker;

	Event *event = new Event();
	tree.Branch("event", event->ClassName(), &event, 32000, 99);
	//Auto-save every 500 MB
	tree.SetAutoSave(500LL * 1024LL * 1024LL);
	
	int iLoop = 0;//event accept loop
	long long iTrial = 0;//event generation loop, must be consistent with pythia internal loop
	int iLoopLast = -1;
	while(iLoop<nEvent){

		//clear the current event container
		if(iTrial>0) event->Clear();

		if(iLoop%nPrint==0&&iLoop!=iLoopLast) cout << iLoop << " events generated" << endl;
		iLoopLast = iLoop;
		//generate 1 event
		call_pyevnt();
		call_setphispin();
		//if bad event skip this one
		if(pypars.msti[60]==1) continue;

		event->Build(iTrial);//dump the pythia event data to a Event object

		if(iTrial<5) {
			call_pylist(2);
			event->PrintParticles();
			cout<<"iLoop="<<iLoop<<" iTrial:"<<iTrial<<endl;
		}


		iTrial++;
		//dump the current event data to the QA hist
		hist.Hfill(event);

		worker.FillHist(event);
		//test if the current event is what we want
		if(worker.GetStatus()==false) continue; //if analysis not accepted generate another

//		if( doSelect == 1 && !selector(event) ) continue;

		//fill the tree
		tree.Fill();

		if(iLoop<20)
			cout<<"accept 1 event!: iLoop="<<iLoop<<endl;
		//update loop flag
		iLoop++;

	}//event loop

	//write tree data and hist data
	file.cd();
	tree.Write();
	hist.Hwrite(&file);
	worker.WriteResults(&file);
	//write cross section data
	write_run_data(file);

	file.Close();

	//print cross sections
	call_pystat(1);
	call_pystat(4);

	cout<<"============================================="<<endl;
	cout << "Pythia total cross section normalisation: "
		 << pypars.pari[0] << " mb" << endl;
	cout << "Total Number of generated events: " << pypars.msti[4] << endl;
	cout << "Total Number of trials: " << pyint5.ngen[2][0] << endl;
	cout << "Integrated Lumi: " << pyint5.ngen[2][0]/pypars.pari[0]*1E-12 << " fb^-1" << endl;

	call_pdfsta();


	return 0;
}


//selector can be used to trigger any events you
//want with the selection rules you defined
bool selector(Event *e){
	bool isPass = false;

	ParticleMC *particle;
	//check if a jpsi exists in the current event list
	for(Int_t itrack(0); itrack<e->GetNTrack(); itrack++){
		particle = e->GetTrack(itrack);
		//if found a Jpsi in the particle list, trigger this event
		if(abs(particle->GetPid()) == 413) {
		//if found a final state pi+, trigger this event
//		if(particle->GetPt()>2&&particle->GetPid()==211&&particle->GetKS()==1) {
			isPass = true;
			break;
		}//if
	}//for

	return isPass;
}

/**
 Print the help message.
 */
void print_help(std::ostream& out) {
   std::cout << "PP implementation of PYTHIA 6" << std::endl;
   std::cout << "Usage: main [options] < input-file" << std::endl;
   std::cout << "Options:" << std::endl;
   std::cout << " --help    print this message and exit" << std::endl;
   std::cout << " --filter  enable selector on the output events" << std::endl;
   std::cout << " --out     set the name of the output ROOT file (default: pythia.root)" << std::endl;
}

/**
 Initialise the program, processing command line arguments.
 Returns the name of the ROOT file to generate.
 */
std::string initialise(int argc, char* argv[]) {
   std::string filename("pythia.root");
   static struct option long_options[] = {
	  // These set a flag
      {"filter", no_argument, &doSelect, 1},
      // These don't set a flag.
      // The arguments are processed below.
      {"out", required_argument, NULL, 'o'},
      {"help", no_argument, NULL, 'h'},
      {NULL, 0, NULL, 0}
    };
   // Loop through options.
   int option_index = 0;
   int code(0);
   while((code = getopt_long(argc, argv, "o:r:",
                             long_options, &option_index)) not_eq -1) {
      switch(code) {
         case 0:
            if(long_options[option_index].flag not_eq 0) {
               break;
            } // if
            printf("option %s", long_options[option_index].name);
            if(optarg) {
               printf (" with arg %s", optarg);
            } // if
            printf("\n");
            break;
        case 'o':
            filename = optarg;
            break;
         case 'h':
            print_help(std::cout);
            exit(0);
         default:
            abort();
      } // switch
   } // while
   // Now that we have processed all command line arguments
   // input/output filename needs to be specified in the command line
	return filename;
}

///write cross section to the root file
void write_run_data(TFile& file) {
    TObjString text;
    std::stringstream stream;
    // Total cross section is stored in PYTHIA in pari(1) in millibarns.
    stream << pypars.pari[0];
    text.SetString(stream.str().c_str());
    file.WriteObject(&text, "crossSection");
    // Total number of generated events is stored in PYTHIA in msti(5).
    stream.str("");
    stream.clear();
    stream << pypars.msti[4];
    text.SetString(stream.str().c_str());
    file.WriteObject(&text, "nEvents");
    // Total number of generateds is stored in PYTHIA in ngen(0, 3).
    stream.str("");
    stream.clear();
    stream << pyint5.ngen[2][0];
    text.SetString(stream.str().c_str());
    file.WriteObject(&text, "nTrials");
}

