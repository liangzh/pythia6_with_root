#include <math.h>

#include <TVector2.h>

#include "Event.h"

ClassImp(Event)
ClassImp(ParticleMC)

//=================ParticleMC part========================
///Particle initialize
ParticleMC::ParticleMC(){
	fPx = 0;
	fPy = 0;
	fPz = 0;
	fE = 0;
	fMass = 0;
	fVx = 0;
	fVy = 0;
	fVz = 0;

	fPid = -999;
	fKS = -999;
	fIndex = -999;
	fOrig = -999;
	fDaughter1 = -999;
	fDaughter2 = -999;

   //derived quantities
    fPt = 0; 
    fP = 0;          
	fRapidity = -19;	 
	fEta = -19;
	fTheta = 0;
	fPhi = 0;	  


}

///particle copier
ParticleMC::ParticleMC(const ParticleMC &track){
	fPx = track.GetPx();
	fPy = track.GetPy();
	fPz = track.GetPz();
	fE = track.GetE();
	fMass = track.GetMass();
	fVx = track.GetVx();
	fVy = track.GetVy();
	fVz = track.GetVz();

	fIndex = track.GetIndex();
	fKS = track.GetKS();
	fPid = track.GetPid();
	fOrig = track.GetParent();
	fDaughter1 = track.GetDaughter1();
	fDaughter2 = track.GetDaughter2();
	
	ComputeDerived(); //compute derived

}

void ParticleMC::Clear(Option_t *){
	TObject::Clear();
}

/** 
 * set pythia variable to particle
 * hard coded with pythia data structure
 * better use the model independent
 * data assignment method
 */
void ParticleMC::Set(const Int_t &index){
	fPx = pyjets.p[0][index];
	fPy = pyjets.p[1][index];
	fPz = pyjets.p[2][index];
	fE = pyjets.p[3][index];
	fMass = pyjets.p[4][index];
	fVx = pyjets.v[0][index];
	fVy = pyjets.v[1][index];
	fVz = pyjets.v[2][index];
	fIndex = index;
	fPid = pyjets.k[1][index];
	fKS = pyjets.k[0][index];
	fOrig = pyjets.k[2][index];
	fDaughter1 = pyjets.k[3][index];
	fDaughter2 = pyjets.k[4][index];

	ComputeDerived(); //compute derived
}

void ParticleMC::ComputeDerived(){
	fPt = TMath::Sqrt(fPx*fPx + fPy*fPy);
	fP = TMath::Sqrt(fPx*fPx + fPy*fPy + fPz*fPz);
	// Rapidity and pseudorapidity
	Double_t Epluspz = fE + fPz; 
	Double_t Eminuspz = fE - fPz; 
	Double_t Ppluspz = fP + fPz; 
	Double_t Pminuspz = fP - fPz; 
	if (Eminuspz <= 0.0 || Pminuspz == 0.0 ||
	 Ppluspz == 0.0 || Epluspz <= 0.0) {
	// Dummy values to avoid zero or infinite arguments in calculations
		fRapidity = -19.;
		fEta = -19.;
	} 
	else {
		fRapidity = 0.5 * log(Epluspz / Eminuspz);
		fEta = 0.5 * log(Ppluspz / Pminuspz);
	}  // if
	fTheta = atan2(fPt, fPz);
	fPhi = TVector2::Phi_0_2pi(atan2(fPy, fPx));
	
}

//=================Event part========================
///Event initialize
///assign value to event wise variables and instantiate fParticles
Event::Event()
: fParticles("ParticleMC", 500) {
	fNTracks = 0;
	fFSNTracks = 0;
	fProcess = -999;
	fTargetParton = -999;
	fTargetPartonX = 0;
	fProjParton = -999;
	fProjPartonX = 0;
	fPt2_hat = 0;
}

void Event::Clear(Option_t *){
	fNTracks = 0;
	fFSNTracks = 0;
	fProcess = -999;
	fTargetParton = -999;
	fTargetPartonX = 0;
	fProjParton = -999;
	fProjPartonX = 0;
	fPt2_hat = 0;
	fParticles.Clear();
}

///fill the content of this event
void Event::Build(const Int_t &index){
	fNEvent = index;
	fNTracks = pyjets.n;
	fProcess = pypars.msti[0];
	fTargetParton = pypars.msti[15];
	fTargetPartonX = pypars.pari[33];
	fProjParton = pypars.msti[14];
	fProjPartonX = pypars.pari[32];
	fPt2_hat = pypars.pari[17];
	
	//add tracks to the current event
	for(int i=0; i<fNTracks; i++){
		ParticleMC particle;
		particle.Set(i);
		AddTrack(particle);
		fLastTrack = &particle;
		//count number of final state particles
		if(particle.GetKS()==1)
			fFSNTracks++;
	}
}

///add one track to the particle list
void Event::AddTrack(const ParticleMC &particle ){
   // To avoid calling the very time consuming operator new for each track,
   // the standard but not well know C++ operator "new with placement"
   // is called. If tracks[i] is 0, a new Track object will be created
   // otherwise the previous Track[i] will be overwritten.
	new(fParticles[fParticles.GetEntries()]) ParticleMC(particle);
}

///get one track from the particle list
ParticleMC* Event::GetTrack(const Int_t &indx) const{
	if(indx<fParticles.GetEntries())
		return static_cast<ParticleMC*>(fParticles.At(indx));
	else
		return NULL;
}

///get the particle list
std::vector<const ParticleMC*> Event::GetTracks() const{
	std::vector<const ParticleMC*> tracks;
	TObject *obj(NULL);
	TIter next(&fParticles);
	while((obj = next())){
		tracks.push_back(static_cast<ParticleMC*>(obj));
	}
	return tracks;
}

