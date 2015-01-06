########################################################################
#                                                                      #
#                      makefile for pythia6.4                          #
#                                                                      #
#                      (zheng.l  04-05-2014)                           #
#                                                                      #
########################################################################

MAIN = main
#PYUSR = pyMainRHIC
PYUSR = pyMainLEP
PYTHIA = pythia-6.4.25
EVENT = Event
HIST = HistMaker

TOP = $(PWD)
EXE = $(MAIN)

######################### compiler - options ###########################
FOPT = -c -Wall -fPIC
COPT = -g -fPIC
GFC = gfortran -g 
GCC = g++ 

#ROOTCFLAGS =  $(COPT) $(shell root-config --cflags)
ROOTCOPT = -pthread -m64 -std=c++11 -fPIC

############################# libraries ################################

#LHAPDF=/afs/rhic.bnl.gov/eic/PACKAGES/LHAPDF-5.8.6/lib
#LHAPDF=/direct/eic+data/zhengl/LHAPDF-5.9.1/
LHAPDF=/home/adminuser/lib/LHAPDF-5.9/lib/
LIB = -L$(LHAPDF)/ -lLHAPDF

ROOTGLIBS = $(shell root-config --glibs)

############################# include ################################
INC = -Iinclude/
ROOTINC = -I$(shell root-config --incdir) 

############################# linker - options #########################
LDOPT = -lgfortran -g $(LIB) $(ROOTLDFLAG) -fPIC
ROOTLDFLAG = $(ROOTGLIBS) $(shell root-config --ldflags)

LD = g++ 

ROOTCFLAGS = $(ROOTCOPT) $(ROOTINC)

########################################################################
all:  $(EXE) 

clean:
	rm -f $(EXE) $(MAIN).o $(PYTHIA).o $(PYUSR).o src/$(EVENT).o src/$(HIST).o lib/libEvent.so cint/*Dict*


$(EXE): $(MAIN).o $(PYUSR).o src/$(EVENT).o src/$(HIST).o lib/libEvent.so
	$(LD) -o $(EXE)  $(MAIN).o $(PYUSR).o src/$(EVENT).o src/$(HIST).o lib/libEvent.so $(LDOPT)

$(PYTHIA).o: $(PYTHIA).f
	$(GFC) -c $(FOPT) $(INC) $<

$(PYUSR).o: $(PYUSR).f
	$(GFC) -c $(FOPT) $(INC) $<

src/$(EVENT).o: src/$(EVENT).cxx include/$(EVENT).h
	$(GCC) -c $(ROOTCFLAGS) $(INC) $<  -o src/$(EVENT).o

src/$(HIST).o: src/$(HIST).cxx include/$(HIST).h
	$(GCC) -c $(ROOTCFLAGS) $(INC) $<  -o src/$(HIST).o


#make dictionary for the Event class
cint/MyDict.cxx: include/HistMaker.h include/pythiaWrapper.h include/Event.h cint/Linkdef.h
	rootcint -f $@ -c $(ROOTCFLAGS)  $(INC) -p HistMaker.h Event.h cint/Linkdef.h

cint/MyDict.o: cint/MyDict.cxx
	$(GCC) -c $(ROOTCFLAGS) $(INC) $< -o cint/MyDict.o
	 
lib/libEvent.so: cint/MyDict.o src/$(EVENT).o src/$(HIST).o $(PYTHIA).o
	mkdir -p lib; g++ -shared -o $@  $(ROOTCFLAGS) $(INC) $^

$(MAIN).o: $(MAIN).C include/pythiaWrapper.h
	$(GCC) -c $(ROOTCFLAGS) $(INC) $<

