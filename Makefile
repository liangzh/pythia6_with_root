########################################################################
#                                                                      #
#                      makefile for pythia6.4                          #
#                                                                      #
#                      (zheng.l  04-05-2014)                           #
#                                                                      #
########################################################################

MAIN = main
#PYUSR = pyMainRHIC
PYUSR = pyMain
PYTHIA = pythia-6.4.25
EVENT = Event
HIST = HistMaker
ANA = Analysis

TOP = $(PWD)
EXE = $(MAIN)

######################### compiler - options ###########################
FOPT = -c -Wall -fPIC
COPT = -g -fPIC
GFC = gfortran -g 
GCC = g++ 

#ROOTCFLAGS =  $(COPT) $(shell root-config --cflags)
ROOTCOPT = -pthread -m64 -std=c++11 -fPIC -g

############################# libraries ################################

#LHAPDF=/afs/rhic.bnl.gov/eic/PACKAGES/LHAPDF-5.8.6/lib
#LHAPDF=/direct/eic+data/zhengl/LHAPDF-5.9.1/
LHAPDF=/home/adminuser/lib/LHAPDF-5.9/lib/
FASTJET=/home/adminuser/lib/fastjet3/
LIB = -L$(LHAPDF)/ -lLHAPDF -L$(FASTJET)/lib/ -lfastjet

ROOTGLIBS = $(shell root-config --glibs)

############################# include ################################
#it is very important to set include path to be absolute, otherwise
#the later usage will cause some trouble when dealing the load of
#shared library
INC = -I$(TOP)/include/
FASTJETINC = -I$(FASTJET)/include/
ROOTINC = -I$(shell root-config --incdir) 

############################# linker - options #########################
LDOPT = -lgfortran -g $(LIB) $(ROOTLDFLAG) -fPIC
ROOTLDFLAG = $(ROOTGLIBS) $(shell root-config --ldflags)

LD = g++ 

ROOTCFLAGS = $(ROOTCOPT) $(ROOTINC)

########################################################################
all:  $(EXE) 

clean:
	rm -f $(EXE) $(MAIN).o $(PYTHIA).o $(PYUSR).o src/$(EVENT).o src/$(HIST).o src/$(ANA).o lib/libEvent.so cint/*Dict*


$(EXE): $(MAIN).o $(PYUSR).o src/$(EVENT).o src/$(HIST).o src/$(ANA).o lib/libEvent.so
	$(LD) -o $(EXE)  $(MAIN).o $(PYUSR).o src/$(EVENT).o src/$(HIST).o src/$(ANA).o $(TOP)/lib/libEvent.so $(LDOPT)

$(PYTHIA).o: $(PYTHIA).f
	$(GFC) -c $(FOPT) $(INC) $<

$(PYUSR).o: $(PYUSR).f
	$(GFC) -c $(FOPT) $(INC) $<

src/$(EVENT).o: src/$(EVENT).cxx include/$(EVENT).h
	$(GCC) -c $(ROOTCFLAGS) $(INC) $<  -o src/$(EVENT).o

src/$(HIST).o: src/$(HIST).cxx include/$(HIST).h
	$(GCC) -c $(ROOTCFLAGS) $(INC) $<  -o src/$(HIST).o

src/$(ANA).o: src/$(ANA).cxx include/$(ANA).h
	$(GCC) -c $(ROOTCFLAGS) $(INC) $(FASTJETINC) $<  -o src/$(ANA).o


#make dictionary for the Event class
#should not use rootcint, otherwise can not Branch a tree when libEvent.so loaded
cint/MyDict.cxx: include/HistMaker.h include/Analysis.h include/pythiaWrapper.h include/Event.h cint/Linkdef.h
	rootcling -f $@ -c $(ROOTCFLAGS)  $(INC) -p HistMaker.h Analysis.h Event.h cint/Linkdef.h
#cint/MyDict.cxx: include/HistMaker.h include/Analysis.h include/pythiaWrapper.h include/Event.h cint/Linkdef.h
#	rootcint -f $@ -c $(ROOTCFLAGS)  $(INC) -p HistMaker.h Analysis.h Event.h cint/Linkdef.h


cint/MyDict.o: cint/MyDict.cxx
	$(GCC) -c $(ROOTCFLAGS) $(INC) $< -o cint/MyDict.o $(LDOPT)
	 
lib/libEvent.so: cint/MyDict.o src/$(EVENT).o src/$(HIST).o src/$(ANA).o $(PYTHIA).o
	mkdir -p lib; g++  $(ROOTCFLAGS) $(INC) $^ -shared -o $@  $(LDOPT); cp cint/*pcm lib/

$(MAIN).o: $(MAIN).C include/pythiaWrapper.h
	$(GCC) -c $(ROOTCFLAGS) $(INC) $<

