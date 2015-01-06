C...This routine initialize pythia6 with parameters read from 
C...the input unit
*======================================================================
      subroutine runpyinit(NEV, NPRT)

      include 'pythia.inc'              ! All PYTHIA commons blocks

C...Added by liang 1/3/12
C...Switches for nuclear correction
      COMMON/PYNUCL/INUMOD,CHANUM,ORDER
      SAVE /PYNUCL/
      !switch to define A(INUMOD), Z(CHANUM)
      DOUBLE PRECISION INUMOD,CHANUM
      !switch to define LO/NLO; ORDER=x*100+y, x=1/2:LO/NLO, y:1~30
      !error sets.
      INTEGER ORDER 

      integer NEV, NPRT, ievent, genevent, I, tracknr
      integer lastgenevent
      REAL pbeam1, pbeam2, emass, pmass
      DOUBLE PRECISION sqrts
      CHARACTER PARAM*100
      CHARACTER OutPutFile*100

c ---------------------------------------------------------------------
c     Run parameter
c ---------------------------------------------------------------------
      integer*4 today(3), now(3)
c---------------------------------------------------------------------
c     ASCII output file
c ---------------------------------------------------------------------
      integer inUnit
      parameter (inUnit=29)
      CHARACTER*256 inname

       pbeam1=100. 
       pbeam2=100.
       pmass1=PYMASS(11)
       pmass2=PYMASS(11)
       ievent=0
       genevent=0
       lastgenevent=0

      

C...Read switch for output file yes or no
c       READ(*,*) outname
C...Read parameters for PYINIT call ( beam and target particle energy).
c       READ(inUnit,*) pbeam1, pbeam2
       READ(*,*) pbeam1, pbeam2
C...Read number of events to generate, and to print.
c       READ(inUnit,*) NEV,NPRT
       READ(*,*) NEV,NPRT
C...Read number of the order used in EPS09, x*100+y:x=1/2,LO/NLO;y=0-30
C...y=1 central err, 2-31 different error sets
c       READ(inUnit,*) ORDER
       READ(*,*) ORDER
C...Read nPDF parameters for A and Z (eg. A=1,Z=1 for p; A=197, Z=79 for
C...Au)
c       READ(inUnit,*) INUMOD,CHANUM
       READ(*,*) INUMOD,CHANUM
C...Read information for cross section used in radgen
c  100  READ(inUnit,'(A)',END=200) PARAM
  100  READ(*,'(A)',END=200) PARAM
       CALL PYGIVE(PARAM)
       GOTO 100
c ---------------------------------------------------------------------
C...Initialize PYTHIA.      
c ---------------------------------------------------------------------
  200  write(*,*) '*********************************************'
       write(*,*) 'NOW all parameters are read by PYTHIA'
       write(*,*) '*********************************************'
C       call PYLIST(11)
C       call PYLIST(12)
       
c...initialie done, close input unit
c       close(inUnit)

C!     
C!     Getting the date and time of the event generation
C!
C        
C      call idate(today)   ! today(1)=day, (2)=month, (3)=year
C      call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
C        
C!      
C!     Take date as the SEED for the random number generation
C!      
C      initseed = today(1) + 10*today(2) + today(3) + now(1) + 5*now(3)
C      write(6,*) 'SEED = ', initseed
C      call rndmq (idum1,idum2,initseed,' ')
C
      sqrts=sqrt(4*pbeam1*pbeam2)
      write(*,*) '*********************************************'
      write(*,*) '1st proton beam energy:', pbeam1, 'GeV'
      write(*,*) '2nd proton beam energy:', pbeam2, 'GeV'
      write(*,*) 'resulting sqrt(s):', sqrts, 'GeV'
      write(*,*) '*********************************************'
C     beam1 is defined in positive z and as target
      P(2,1)=0.0  
      P(2,2)=0.0  
      P(2,3)=pbeam1
C     beam2 is defined in negative z and as beam
      P(1,1)=0.0  
      P(1,2)=0.0  
      P(1,3)=-pbeam2
      call pyinit ('3MOM','e+','e-',WIN)

c ---------------------------------------------------------------------
C...Event generation loop
c ---------------------------------------------------------------------

C   This is what we write in the ascii-file

c        write(29,*)' PYTHIA EVENT FILE '
c        write(29,*)'============================================'
c        write(29,30) 
c30      format('I, ievent, genevent, subprocess, nucleon,
c     &  targetparton, xtargparton, beamparton, xbeamparton,
c     &  thetabeamprtn, pt2_hat, Q2_hat, nrTracks')
c        write(29,*)'============================================'
c
c        write(29,*)' I  K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5)
c     &  P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5)  V(I,1)  V9I,2)  V(I,3)'
c        write(29,*)'============================================'
c
c       DO 300 IEV=1,NEV
c999      CALL PYEVNT
c         if (MSTI(61).eq.1) then
c            write(*,*) 'go back to PYEVNT call'
c            goto 999
c         endif
cC         CALL PYLIST(2)
c
c         ievent=IEV       
c         genevent=NGEN(0,3)-lastgenevent
c
c         tracknr=N
c
c         write(29,32) 0,ievent,genevent, msti(1), msti(12), 
c     &        msti(16), pari(34), msti(15), pari(33), pari(53), 
c     &        pari(18),  pari(22), tracknr
c 32      format((I4,1x,$),(I10,1x,$),4(I4,1x,$),f9.6,1x,$,I12,1x,$,
c     &           4(f12.6,1x,$),I12,/)
c         write(29,*)'============================================'
c
c         DO I=1,tracknr
c         if (K(I,3).le.tracknr) then
c         write(29,34) I,K(I,1),K(I,2),K(I,3),K(I,4),K(I,5),
c     &        P(I,1),P(I,2),P(I,3),P(I,4),P(I,5),
c     &        V(I,1),V(I,2),V(I,3)
c         endif
c         ENDDO
c 34      format(2(I6,1x,$),I10,1x,$,3(I6,1x,$),8(f15.6,1x,$),/)
c         write(29,*)'=============== Event finished ==============='
c         lastgenevent=NGEN(0,3)
c
c  300  CONTINUE
c      
cC...Print cross sections.
c       CALL PYSTAT(1)
c       CALL PYSTAT(4)
cC...Print the Pythia cross section which is needed to get an absolut 
cC   normalisation the number is in microbarns
c       write(*,*)'==================================================='
c       write(*,*)'Pythia total cross section normalisation:',
c     +            pari(1)*1000, ' microbarn'
c       write(*,*)'Total Number of generated events', MSTI(5)
c       write(*,*)'Total Number of trials', NGEN(0,3)
c       write(*,*)'==================================================='
c       close(29)
c

      RETURN
      END
