#ifndef Event_H
#define Event_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Event for the pythia pp simulations                                  //
//                                                                      //
// Description of the event and track parameters                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>

#include <TObject.h>
#include <TClonesArray.h>
#include <TRefArray.h>
#include <TRef.h>
#include <TH1.h>
#include <TBits.h>
#include <TMath.h>
#include <Rtypes.h> //for class Def

#include "pythiaWrapper.h" //for pythia data struct


class ParticleMC : public TObject {

private:
	///read in identities 
	Double_t    fPx;           ///<X component of the momentum P(I,1)
	Double_t    fPy;           ///<Y component of the momentum P(I,2)
	Double_t    fPz;           ///<Z component of the momentum P(I,3)
	Double_t    fE;           ///<E f the particle P(I,4)
	Double_t	fMass;        ///<The mass of this particle P(I,5)
	Double_t    fVx;           ///<X component of the position V(I,1)
	Double_t    fVy;           ///<Y component of the position V(I,2)
	Double_t    fVz;           ///<Z component of the position V(I,3)

	Int_t	   fPid;		///<MC Pid of this particle K(I,2)
	Int_t	   fKS;		///<MC status of this particle K(I,1)
	Int_t	   fIndex;		///<index of this particle in the MC file list, I
	Int_t	   fOrig;		///<mother index of this particle K(I,3)
	Int_t	   fDaughter1;	///<first daughter index of this particle K(I,4)
	Int_t	   fDaughter2;	///<last daughter index of this particle K(I,5)


    ///derived quantities
	Double_t    fPt;           ///<T component of the momentum
	Double_t    fP;           ///<magnitude of the momentum
	Double_t	fEta;		 ///<Eta of the particle
	Double_t	fRapidity;	 ///<Rapidity of the particle
	Double_t	fTheta;		 ///<Theta of the particle
	Double_t	fPhi;	 ///<Phi of the particle

public:
	///constructuor
	ParticleMC();

	///allow copy another particle
	ParticleMC(const ParticleMC &track);

	///destructor
	virtual ~ParticleMC() { this->Clear();}

   //operator= definition, assign the variables for another particle
//   ParticleMC &operator=(const ParticleMC &orig);

	//inline functions
	///read data
	Double_t       GetPx() const { return fPx; }
	Double_t       GetPy() const { return fPy; }
	Double_t       GetPz() const { return fPz; }
	Double_t	   GetE() const { return fE; }
	Double_t       GetPt() const { return fPt; }
	Double_t       GetP() const { return  fP; }
	Double_t	   GetTheta() const { return fTheta; }
	Double_t	   GetPhi() const { return fPhi; }
	Double_t	   GetEta() const { return fEta; }
	Double_t	   GetRapidity() const { return fRapidity; }
	Double_t       GetMass() const { return fMass; }
	Double_t       GetVx() const { return fPx; }
	Double_t       GetVy() const { return fPy; }
	Double_t       GetVz() const { return fPz; }

	Int_t		GetIndex() const { return fIndex; }
	Int_t		GetParent() const { return fOrig; }
	Int_t		GetKS() const { return fKS; }
	Int_t		GetPid() const { return fPid; }
	Int_t		GetDaughter1() const { return fDaughter1; }
	Int_t		GetDaughter2() const { return fDaughter2; }

	///set data
	void SetPx(const Double_t &px) { fPx = px; }
	void SetPy(const Double_t &py) { fPy = py; }
	void SetPz(const Double_t &pz) { fPz = pz; }
	void SetE(const Double_t &E) { fE = E; }
	void SetMass(const Double_t &M) { fMass = M; }
	void SetVx(const Double_t &vx) { fVx = vx; }
	void SetVy(const Double_t &vy) { fVy = vy; }
	void SetVz(const Double_t &vz) { fVz = vz; }
	void SetIndex(const Double_t &index) { fIndex = index; }
	void SetPid(const Int_t &pid) { fPid = pid; }
	void SetKS(const Int_t &KS) { fKS = KS; }
	void SetParent(const Int_t &orig) { fOrig =  orig; }
	void SetDaughter1(const Int_t &daugther1) { fDaughter1 = daugther1; }
	void SetDaughter2(const Int_t &daugther2) { fDaughter2 = daugther2; }
	void Set(const Int_t &index);
	void ComputeDerived();

	void Clear(Option_t *option="");

	ClassDef(ParticleMC,1);  //A MC particle
};


class Event : public TObject {

private:
	Int_t fNEvent; ///<event index
	Int_t fProcess; ///<PYTHIA process code for the current event MSTI(1)
	Int_t fNTracks; ///<number of total particles in the event (Final + intermediate)
	Int_t fFSNTracks; ///<number of FS particles in the event
	Int_t fTargetParton; ///<parton pdg from target MSTI(16)
	Int_t fProjParton; ///<parton pdg from projectile MSTI(15)
	Double_t fTargetPartonX; ///<parton x from target PARI(34)
	Double_t fProjPartonX; ///<parton x from projectile PARI(33)

	Double_t fPt2_hat; ///<pt2_hat at hard interaction PARI(18)

	TClonesArray fParticles; ///<particle list
	TRef fLastTrack; ///<index to the last track

public:
	Event();
	virtual ~Event() { this->Clear(); };

	void Clear(Option_t *option="");

	void Build(const Int_t &index);

	void AddTrack(const ParticleMC &particle);

	/**
	 * Get a pointer to one track from the current Event
	 * returns NULL if the track can not be located in the Event
	 */
	ParticleMC* GetTrack(const Int_t &indx) const;

	Int_t GetNTrack() const { return fNTracks; }

	Int_t GetFSNTrack() const { return fFSNTracks; }

	//   TClonesArray *GetTracks() const { return fParticles; }
	std::vector<const ParticleMC*> GetTracks() const;

	Int_t GetProcess() const { return fProcess; }

	Int_t GetEventIndex() const { return fNEvent; }

	Int_t GetTargetParton() const { return fTargetParton; }

	Int_t GetProjParton() const { return fProjParton; }

	Double_t GetTargetPartonX() const { return fTargetPartonX; }

	Double_t GetProjPartonX() const { return fProjPartonX; }

	Double_t GetPt2_hat() const { return fPt2_hat; }

	ClassDef(Event,1);
};

#endif
