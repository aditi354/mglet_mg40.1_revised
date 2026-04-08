










CCCCC        DEFINITIONEN FUER C-PREPROZESSOR
C            MASCHINE

C              MESSAGE PASSING INTERFACE

C            RAEUMLICHE DISKRETISIERUNG
C
***********************************************************************
C                       KOMPAKT-UPWIND in X-RICHTUNG (3-TER ORDNUNG) fuer
C                       die U-Komponente (V und W bleiben Kompakt 4ter ordnung)
C                       falls im Stroemungsfeld eine Koerper vorhanden ist
C
***********************************************************************
C                       KOMPAKTVERF. in XYZ-RICHTUNG (4-TER ORDNUNG)


**********************************************************************
C                      Preprocessing with ADM for 2nd Order Central
***********************************************************************
CCC                     Bei periodischen Randbedingungen in X-Richtung
CCC                     ist eine hoehere O
C
CCC                     Bei periodischen Randbedingungen in Y-Richtung
CCC                     ist eine hoehere Ordnung moeglich
C
CCC                     Bei periodischen Randbedingungen in Y-Richtung:
CCC                     Zentraldiff. 4-ter Ordnung (Parallelisieung moeglich)
C
***********************************************************************
C
C
C            INTERPOLATION AN GITTER-GRENZEN

C            BEHANDLUNG DER TOPAR-RANDBEDINGUNG

C            ZEITLICHE DISKRETISIERUNG

C            FEINSTRUKTURMODELL

C            NICHT-NEWTONSCHE SPANNUNGEN


C            TRANSPORT UND ORIENTIERUNG VON PARTIKELN


C            SCALAR HEAT/TEMPERATURE TRANSPORT
C            APROXIMATE DECONVOLUTION FOR SCALAR
C            TURBULENT PRANDTL NUMBER 
C            PLOT LOCAL SCALAR CONVECTION DIFFUSION EXTREMES 
C            ANALYZE AND PLOT NEAR WALL GRID RESLOLUTION 
C            WRITE SPECIAL 1D LINE FOR TMIX FOR SPECTRA AND PDF
C            SCALAR TIME ADVANCEMENT/DISCRETISATION METHOD

C            STROEMUNG NACH OBEN?


C            BEHANDLUNG DER FLUKTUATIONS-RANDBEDINGUNG

C            BEHANDLUNG VON PPHYS ALS QUELLTERM IN TSTLE2

C            POSITIONIERUNG DER EINSTROEMPROFILE

C           STATISTIK


C                 GEMISCHTES

C                GRDFMI WIRD NICHT VERWENDET, DAHER EINSPARUNG DER FELDER
C                IGRHF UND GRHF

C                VERGROESSERN DER KJI-RICHTUNG BEI FORTSETZUNGSLAUF ERLAUBT!

C                RAUSSCHREIBEN VON ZEITRECORDS

C                FUER KORRELATIONEN UND HAEFIGKEITSVERTEILUNGEN
      SUBROUTINE OPTFREQ(VAROPT,NDIMVAROPT
     +                   ,IDIM2D,NTO1,TIMEPH,ITTOT)
C
C
C*MGLET***************************************************************
C        O P T F R E Q     OPTIMIERUNG DER KONTROLLFREQUENZ 
C*MGLET***************************************************************
C
C PARAM: NDIMVAROPT     - ANZAHL DER STICHPROBEN FUER ANALYSE VON XR
C        VAROPT(2,NDIMVAROPT)    - WERTE DER VARIABLEN, DIE AUSGEWERTET WIRD
C        IDIM2D         - DIMENSIONIERUNG DES FELDES VAROPT
C        NTO1           - ZEITSCHRITT, AB DEM OPTIMIERT WERDEN SOLL
C        TIMEPH         - AKTUELLE PHYSIKALISCHE ZEIT
C	
C ALGOTRITHMUS:   NACH JEWEILS .PER. PERIODEN DER EINGESTELLTEN FREQUENZ
C                 WIRD EINE AUSWERTUNG VON VAROPT DURCHGEFUEHRT.
C                 DABEI WIRD DIE FREQUENZ VERAENDERT, WENN SICH AUFGRUND
C                 DIESER VERAENDERUNG EINE ERHOEHUNG/ERNIEDRIGUNG DER
C                 AUGENBLICKLICHEN REZIRKULATIONSLAENGE ERGIBT.
C
C*MGLET***************************************************************
C
      REAL VAROPT(2,IDIM2D*8)


      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


      INTEGER MAXLVLMAX,MINLVLMIN,MAXGRDSOFLVL

      PARAMETER ( MAXLVLMAX = 4,  MINLVLMIN = -4 )
      PARAMETER ( MAXGRDSOFLVL =260)

      COMMON /COLEVEL/ MAXLEVEL,MINLEVEL
      COMMON /COLEVEL/ NOFLEVEL,IGRDOFLEVEL,NOFPSDIR,IGRDOFPSDIR
      COMMON /COLEVEL/ NOFVPIT ,IGRDOFVPIT ,NOFTST,  IGRDOFTST
      COMMON /COLEVEL/ NOFSLCED,IGRDOFSLCED,NOFSLICE,IGRDOFSLICE
      COMMON /COLEVEL/ NOFSCAI, IGRDOFSCAI,NGRIDPERLEVEL,NOFTHISGRID

      INTEGER MAXLEVEL,MINLEVEL

      INTEGER NOFLEVEL   (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFLEVEL(MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFVPIT    (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFVPIT (MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFTST     (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFTST  (MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFSCAI    (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFSCAI (MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFSLCED   (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFSLCED(MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFSLICE   (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFSLICE(MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFPSDIR   (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFPSDIR(MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NGRIDPERLEVEL    
      INTEGER NOFTHISGRID(MAXGRDSOFLVL)

      INTEGER MAXPROCS,NUMPROCS,MYID,IDIM_MG_BUFS,IDIM_MG_BIG
      PARAMETER (MAXPROCS = 64)
      PARAMETER ( IDIM_MG_BUFS=  9000000 , IDIM_MG_BIG= 19500000 )
      COMMON /COIDTHR/ NUMPROCS,MYID,IDPROCOFGRD
      COMMON /COIDTHR/ NVPOFPROC,IVPOFPROC,NTSOFPROC,ITSOFPROC
      COMMON /COIDTHR/ NSIOFPROC,ISIOFPROC
      INTEGER IDPROCOFGRD(MAXGRIDS)
      INTEGER IVPOFPROC  (MAXGRDSOFLVL,0:MAXPROCS)
      INTEGER NVPOFPROC  (             0:MAXPROCS,MINLVLMIN:MAXLVLMAX)
      INTEGER ITSOFPROC  (MAXGRDSOFLVL,0:MAXPROCS)
      INTEGER NTSOFPROC  (             0:MAXPROCS,MINLVLMIN:MAXLVLMAX)
      INTEGER ISIOFPROC  (MAXGRDSOFLVL,0:MAXPROCS)
      INTEGER NSIOFPROC  (             0:MAXPROCS,MINLVLMIN:MAXLVLMAX)
      COMMON /MG_BUFS/ SENDBUF,RECVBUF,BIGBUF
      REAL SENDBUF(IDIM_MG_BUFS),RECVBUF(IDIM_MG_BUFS)
      REAL BIGBUF(IDIM_MG_BIG)

      COMMON /COSTRLES/
     $                  VERS,   NRRUN,
     $                  DREAD,  DWRITE, DCONT,  DCUB,   
     $                  LBODYINSERT,
     $                  LGRIDNEW,    LGRIDREMOVE,
     $                  NPRNEU, FPRNEU, MTSTEP,
     $                  DT,     ITPRIN, IPINF,  ITINT,
     $                  IPRINT_WSS, IPRINT_WNS, ITFLUC, ITMIT,  
     $                  MPCORR, EPCORR, EPFAK,  MPCVOR, MPCNACH,
     $                  IVPINF, 
     $                  OMG,    LDIMLO, ISETRE,
     $                  LINPRN, IWRB,   LREC

      INTEGER 
     $                  NRRUN,
     $                  MTURB,  NPRNEU, MTSTEP,
     $                          ITPRIN, IPINF,  ITINT,
     $                  ITFLUC, ITMIT,  MPCORR,        
     $                                          ISETRE,
     $                          MSLIN,  IWRB          

      LOGICAL
     $                  DREAD,  DWRITE, DCONT,  DCUB,   
     $                  LBODYINSERT,
     $                  LGRIDNEW,       LGRIDREMOVE,
     $                  LDIMLO, LINPRN, LREC

      REAL
     $                  DT,     EPCORR,  OMG, FPRNEU

      CHARACTER (LEN=8)       VERS

C
C     VERS:   28.09.95 (AO) JET VARIABLES INTRODUCED
C             12.12.02 (TB) SCALAR PARAMETERS (FIX VALUE, GRADIENT) ADDED
C
C     XMPOS1,XMPOS2,YMPOS1,YMPOS2,ZMPOS1,ZMPOS2: BEREICH DER EFFEKTMESSUNG 
C                                                AUFGRUND DER MANIPULATION

      COMMON /COBOUND/
     &                 NBOCD,
     &                NBOCONDS,     LARBOCONDS,     ITYPBOCONDS,
     &                LBOGRIDS,    LPOSBOGRIDS,
     &                 FRONT,    BACK,     RIGHT,     LEFT,
     &                 BOTTOM,   TOP,      CUBE,
     &                 IBPOS,   JBPOS,    KBPOS,
     &                 IBANF,   JBANF,    KBANF,
     &                 IBEND,   JBEND,    KBEND,
     &                 XBANF,   YBANF,    ZBANF,
     &                 XBEND,   YBEND,    ZBEND,
     &                 ANIVEAU,
     &                 AUB, AVB, AWB,
     &                 FREQB,  FLOWTYP, WAVENUMBER,
     &                 LOPT,NTOPT1,NTOPT2,
     &                 XMPOS1,YMPOS1,ZMPOS1,
     &                 XMPOS2,YMPOS2,ZMPOS2,
     &                 TIMEALT,XRTALT1,XRTALT2,FREQOPT,PERIODE,ITALT,
     &                 FREQALT1,FREQALT2,XRMIN,FXRMIN,NXRMIN,
     &                 RANNUM,PHASE


      INTEGER
     &       NBOCD   (                9, MAXGRIDS),
     &       NBOCONDS(                9, MAXGRIDS),
     &     LARBOCONDS( 6, MAXBOCONDS, 9, MAXGRIDS),
     &    ITYPBOCONDS(    MAXBOCONDS, 9, MAXGRIDS),
     &       LBOGRIDS(    MAXBOCONDS, 9, MAXGRIDS),
     &    LPOSBOGRIDS( 3, MAXBOCONDS, 9, MAXGRIDS),
     &   IBPOS(MAXBOCONDS,9,MAXGRIDS),  JBPOS(MAXBOCONDS,9,MAXGRIDS),
     &   KBPOS(MAXBOCONDS,9,MAXGRIDS),
     &   IBANF(MAXBOCONDS,9,MAXGRIDS),  IBEND(MAXBOCONDS,9,MAXGRIDS),
     &   JBANF(MAXBOCONDS,9,MAXGRIDS),  JBEND(MAXBOCONDS,9,MAXGRIDS),
     &   KBANF(MAXBOCONDS,9,MAXGRIDS),  KBEND(MAXBOCONDS,9,MAXGRIDS),
     &   LOPT(MAXBOCONDS,MAXGRIDS),
     &   NTOPT1(MAXBOCONDS,MAXGRIDS),NTOPT2(MAXBOCONDS,MAXGRIDS)



      CHARACTER (LEN=16)
     &      FRONT(MAXBOCONDS,MAXGRIDS),   BACK(MAXBOCONDS,MAXGRIDS),
     &      RIGHT(MAXBOCONDS,MAXGRIDS),   LEFT(MAXBOCONDS,MAXGRIDS),
     &     BOTTOM(MAXBOCONDS,MAXGRIDS),    TOP(MAXBOCONDS,MAXGRIDS),
     &       CUBE(MAXBOCONDS,MAXGRIDS),
     &      FLOWTYP(MAXBOCONDS,9,MAXGRIDS)

      REAL  ANIVEAU(MAXBOCONDS,9,MAXGRIDS), AUB(MAXBOCONDS,9,MAXGRIDS),
     &      AVB(MAXBOCONDS,9,MAXGRIDS), AWB(MAXBOCONDS,9,MAXGRIDS),
     &      FREQB(MAXBOCONDS,9,MAXGRIDS),
     &      XBANF(MAXBOCONDS,9,MAXGRIDS),XBEND(MAXBOCONDS,9,MAXGRIDS),
     &      YBANF(MAXBOCONDS,9,MAXGRIDS),YBEND(MAXBOCONDS,9,MAXGRIDS),
     &      ZBANF(MAXBOCONDS,9,MAXGRIDS),ZBEND(MAXBOCONDS,9,MAXGRIDS),
     &      XM1(MAXBOCONDS,9,MAXGRIDS),XM2(MAXBOCONDS,9,MAXGRIDS),
     &      XM3(MAXBOCONDS,9,MAXGRIDS),
     &      YM1(MAXBOCONDS,9,MAXGRIDS),YM2(MAXBOCONDS,9,MAXGRIDS),
     &      YM3(MAXBOCONDS,9,MAXGRIDS),
     &      ZM1(MAXBOCONDS,9,MAXGRIDS),ZM2(MAXBOCONDS,9,MAXGRIDS),
     &      ZM3(MAXBOCONDS,9,MAXGRIDS),
     &      WAVENUMBER(MAXBOCONDS,9,MAXGRIDS),
     &      XMPOS1(MAXBOCONDS,MAXGRIDS),
     &      XMPOS2(MAXBOCONDS,MAXGRIDS),
     &      YMPOS1(MAXBOCONDS,MAXGRIDS),
     &      YMPOS2(MAXBOCONDS,MAXGRIDS),
     &      ZMPOS1(MAXBOCONDS,MAXGRIDS),
     &      ZMPOS2(MAXBOCONDS,MAXGRIDS),
     &      RANNUM,
     &      PHASE(MAXBOCONDS,9,MAXGRIDS)

C*******************************************************************
C*******************************************************************
C          1. PHASE: WENN NTO1 DANN DEFINIERE ANFANGSPARAMETER
C*******************************************************************
C*******************************************************************
      IF(ITTOT.EQ.NTO1) THEN
         FREQOPT = 1.00
         TIMEALT = TIMEPH
         ITALT   = ITTOT+1
         PERIODE = MIN( MAX(5.0,20.0*FREQOPT) , 80.0*FREQOPT)
         NDIMVAROPT = INT(PERIODE/(FREQOPT*DT)+1.)
         TRAN = TIMEPH+NDIMVAROPT*DT
         WRITE(6,*)'OPTF: 1.PHASE >>>> INITIALISIERUNG <<<<!'
         WRITE(6,*)'OPTF: BIS           ',TRAN 
         WRITE(6,*)'OPTF: TIMEPH,OPTFREQ',TIMEALT,FREQOPT

         IF (FREQOPT.GT.0.4) THEN
           XRTALT1 = 2.0
           XRTALT2 = 1.0
         ELSE
           XRTALT1 = 100.0
           XRTALT2 = 150.0
         ENDIF
         FREQALT1 = 0.0
	 XRMIN    = 100.0
	 NXRMIN   = 0
	 GOTO 10
      ENDIF

C*******************************************************************
C*******************************************************************
C           2. PHASE: WENN VAROPT VOLLSTAENDIG, DANN ANALYSIERE
C                     VAROPT
C*******************************************************************
C*******************************************************************
      IF(ITTOT.GE.ITALT+NDIMVAROPT) THEN
        FREQRATIO1 = 0.75
        FREQRATIO11= 0.55
        FREQRATIO2 = 1.36
        FREQRATIO22= 1.80

        CALL WRIVAROPT(VAROPT,IDIM2D,NDIMVAROPT)
        WRITE(6,*)'OPTF: 2.PHASE >>>> XR ANALYSE <<<<!' 
        FREQALT2=FREQALT1
        FREQALT1=FREQOPT

	XRSUM=0.0
        DO I=1,NDIMVAROPT
          XRSUM = XRSUM + VAROPT(2,I)
        ENDDO
        XRT = XRSUM / (FLOAT(NDIMVAROPT))

        NXRMIN=NXRMIN+1
        IF (XRT.LE.XRMIN) THEN
           XRMIN  = XRT
           FXRMIN = FREQOPT
           NXRMIN = 0
        ENDIF

        WRITE(6,*)'OPTF: TIMEPH,OPTFREQ',TIMEPH,FREQOPT 
        WRITE(6,*)'OPTF: XRT',TIMEPH,XRT,XRTALT1,FREQOPT
        WRITE(6,*)'OPTF: XRMIN',TIMEPH,XRMIN,FXRMIN,NXRMIN

C*******************************************************************
C                 DECISION:
C                     CASE X: NXRMIN > 5 => JUMP BACK TO FXRMIN
C                     CASE A: XRT < XRTALT1
C                     CASE B: XRT > XRTALT1
C                     CASE C: XRT ~ XRTALT1  &&  XRTALT1 ~ XRTALT2
C                     CASE D: XRT ~ XRTALT1  &&  XRTALT1 <> XRTALT2 
C*******************************************************************

        IF(NXRMIN.GT.5) THEN
          FREQOPT=FXRMIN
	  NXRMIN =0
          WRITE(6,*)'OPTF: ST:=ST(XRMIN) (CASE X)'
     +               ,TIMEPH,FREQOPT 
        ELSEIF(XRT*0.99.GE.XRTALT1) THEN 
          IF(FREQALT2.GT.FREQOPT) THEN
            FREQOPT=FREQOPT*FREQRATIO2
            WRITE(6,*)'OPTF: ST UP (CASE A1)'
     +                 ,TIMEPH,FREQOPT 
          ELSE
            FREQOPT=FREQOPT*FREQRATIO1
            WRITE(6,*)'OPTF: ST DOWN (CASE A2)'
     +                 ,TIMEPH,FREQOPT 
          ENDIF
        ELSEIF(XRT*1.01.LT.XRTALT1) THEN 
          IF(FREQALT2.GT.FREQOPT) THEN       
            FREQOPT=FREQOPT*FREQRATIO1      
            WRITE(6,*)'OPTF: ST DOWN (CASE B1)'
     +                 ,TIMEPH,FREQOPT 
          ELSE
            FREQOPT=FREQOPT*FREQRATIO2
            WRITE(6,*)'OPTF: ST UP (CASE B2)'
     +                 ,TIMEPH,FREQOPT 
          ENDIF
        ELSEIF( (XRTALT1*0.99.LT.XRTALT2) .AND. 
     +          (XRTALT1*1.01.GE.XRTALT2) ) THEN
          IF(FREQALT2.GT.FREQOPT) THEN
            FREQOPT=FREQOPT*FREQRATIO11
            WRITE(6,*)'OPTF: XR GLEICH GEBLIEBEN'
            WRITE(6,*)'OPTF: ST DOWN_DOWN (CASE C1)'
     +                 ,TIMEPH,FREQOPT 
          ELSE
            FREQOPT=FREQOPT*FREQRATIO22
            WRITE(6,*)'OPTF: XR GLEICH GEBLIEBEN'
            WRITE(6,*)'OPTF: ST UP_UP (CASE C2)'
     +                 ,TIMEPH,FREQOPT 
          ENDIF
        ELSE
          IF(FREQALT2.GT.FREQOPT) THEN       
            FREQOPT=FREQOPT*FREQRATIO1      
            WRITE(6,*)'OPTF: ST DOWN (CASE D1)'
     +                 ,TIMEPH,FREQOPT 
          ELSE
            FREQOPT=FREQOPT*FREQRATIO2
            WRITE(6,*)'OPTF: ST UP (CASE D2)'
     +                 ,TIMEPH,FREQOPT 
          ENDIF
        ENDIF

C*******************************************************************
C                                            ENDE DER IF-ABFRAGEN
C*******************************************************************
                                         
        DO I=1,IDIM2D*8
          VAROPT(1,I)  = 0.0
          VAROPT(2,I)  = 0.0
        ENDDO
        PERIODE = MIN( MAX(5.0,20.0*FREQOPT) , 80.0*FREQOPT)
        NDIMVAROPT = INT(PERIODE/(FREQOPT*DT)+1.)
        TIMEALT = TIMEPH
        ITALT   = ITTOT
        XRTALT2 = XRTALT1
        XRTALT1 = XRT


      ENDIF
   10 CONTINUE

  
      RETURN
      END


