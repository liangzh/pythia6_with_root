/**
 * pythia c++ wrapper which converts pythia6 data and functions
 * to accessible c++ data types and functions
 */

//--------------------------------------------------------------------------
#ifndef PYTHIA_WRAPPER_H
#define PYTHIA_WRAPPER_H

#include <ctype.h>
#include <cstring>

#include <algorithm>
#include <iostream>


//--------------------------------------------------------------------------
// Initialization routine

extern "C" {
    void initpydata(void);
}
#define initpydata initpydata_

//--------------------------------------------------------------------------
// PYTHIA Common Block Declarations

//index in PYHIA common blocks starts from 0 as a data array
//while the fortran block starts from 1

const int pyjets_maxn =4000;
// COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
extern "C" {
    extern struct {
        int n, npad, k[5][pyjets_maxn];
        double p[5][pyjets_maxn], v[5][pyjets_maxn];
    } pyjets_;
}
#define pyjets pyjets_

// COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
extern "C" {
    extern struct {
        int mstu[200];
        double paru[200];
        int mstj[200];
        double parj[200];
    } pydat1_;
}
#define pydat1 pydat1_

// COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
extern "C" {
    extern struct {
        int kchg[4][500];
        double pmas[4][500], parf[2000], vckm[4][4];  
    } pydat2_;
}
#define pydat2 pydat2_

// COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
extern "C" {
    extern struct {
        int mdcy[3][500], mdme[2][8000];
        double brat[8000];
        int kfdp[5][8000];
    } pydat3_;
}
#define pydat3 pydat3_

// COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
extern "C" {
    extern struct {
        int msel, mselpd, msub[500], kfin[81][2];
        double ckin[200];
    } pysubs_;
}
#define pysubs pysubs_

// COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
extern "C" {
    extern struct {
        int mstp[200];
        double parp[200];
        int msti[200];
        double pari[200];
    } pypars_;
}
#define pypars pypars_

// COMMON/PYINT1/MINT(400),VINT(400)
extern "C" {
    extern struct {
        int mint[400];
        double vint[400];
    } pyint1_;
}
#define pyint1 pyint1_

// COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
extern "C" {
    extern struct {
        int iset[500], kfpr[2][500];
        double coef[20][500];
        int icol[2][4][40];       // was [320] was [40][4][2]
    } pyint2_;
}
#define pyint2 pyint2_

// COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
extern "C" {
    extern struct {
        int ngenpd, ngen[3][501];
        double xsec[3][501];
    } pyint5_;
}
#define pyint5 pyint5_

//--------------------------------------------------------------------------
// PYTHIA routines declaration

#define pycomp pycomp_
#define pyhepc pyhepc_ 
#define pyinit pyinit_
#define pylist pylist_
#define pyedit pyedit_
#define pystat pystat_
#define pyevnt pyevnt_
#define pyeevt pyeevt_
#define upinit upinit_
#define upevnt upevnt_
    extern "C" {
				int pycomp(int*);
        void pyhepc(int*);
        void pyinit(const char*,const char*,const char*,double*,int,int,int);
        void pylist(int*);
        void pyedit(int*);
        void pystat(int*);
        void pyeevt(int*, double*);
        void pyevnt();
        void upinit();
        void upevnt();
    }

// define methods to hide the subtle syntax necessary to call fortran from C++
inline int call_pycomp( int pid ) { return pycomp( &pid ); }
inline void call_pyhepc( int mode ){ pyhepc( &mode ); }
inline void call_pyinit( const char* frame, const char* beam, const char* target,
                  double win ) 
{ pyinit( frame,beam,target,&win,strlen(frame),strlen(beam),strlen(target) ); }
inline void call_pylist( int mode ){ pylist( &mode ); }
inline void call_pystat( int mode ){ pystat( &mode ); }
inline void call_pyedit( int mode ){ pyedit( &mode ); }
inline void call_pyevnt(){ pyevnt(); }
inline void call_pyeevt(int flavor, double ecm){ pyeevt(&flavor, &ecm); }

#define pysphe pysphe_
extern "C" {
	void pysphe(double*, double*);
}

inline void call_pysphe(double sph, double apl){ pysphe(&sph, &apl); }


#define pdfsta pdfsta_
extern "C" {
	void pdfsta();
}
inline void call_pdfsta(){ pdfsta(); }

//--------------------------------------------------------------------------
// PYTHIA block data
// ( with gcc it works to initialize the block data by calling 
//   "pydata();" at beginning, but this fails for f77, so the fortran routine
//   initpydata.f is supplied ... call it instead for platform independent
//   behaviour )

#define pydata pydata_
extern "C" {
    void pydata(void);
}

// originally strip from the fortran pythia main function
// read in the input file parameters setup and call initialization
// routines for pythia6
// 1st par: number of events to run
// 2nd par: every number of events to print
#define runpyinit runpyinit_
extern "C" {
	void runpyinit(int*, int*);
}


#endif  // PYTHIA_WRAPPER_H
//--------------------------------------------------------------------------
