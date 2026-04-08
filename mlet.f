










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
C*STARLET***************************************************************
C
C
C  UEBERSICHT ITERATIVE  NAVIER-STOKES LOESUNG IM U,V,W,P-
C             SYSTEM.
C             LAUFFAEHIG UNTER NOS/VE MIT FTN 5 (CYBER 990)
C                        UNTER UNIX (CRAY, CONVEX, SILICON GRAPHICS, KSR)
C
C  STICHWORTE QUADERFOERMIGE BERECHNUNGSGEBIETE
C             EINGEBETTETE GITTER
C             BELIEBIG GEFORMTE EINBAUTEN
C
C  BELEGUNG
C  D. KANAELE KANAL 1: UNFORMATIERTES LESEN DER 3-D STARTFELDER (SUBR.
C                      DIB)
C             KANAL 2: UNFORMATIERTES SCHREIBEN DER 3-D FELDER AM ENDE
C                      EINES LAUFES (SUBR. DOB)
C             KANAL 3: FORMATIERTES LESEN DER 3-D STARTFELDER (SUBR.
C                      DIC)
C             KANAL 4: FORMATIERTES SCHREIBEN DER 3-D FELDER AM ENDE
C                      EINES LAUFES (SUBR. DOC)
C             KANAL 5: INPUT
C             KANAL 6: OUTPUT
C
C             KANAL 8: CHECK-OUTPUT, WIRD MIT ITPRIN GESTEUERT
C             KANAL 9: CHECK-OUTPUT, WIRD MIT ITPRIN GESTEUERT
C
C             KANAL11: UNFORMATIERTES LESEN DER EINTRITTS-PROFILE
C                      U,V,W UND G  IN SUBR. DEIBI. (NUR F. LES)
C             KANAL12: UNFORMATIERTES SCHREIBEN DER EINTRITTS-PROFILE
C                      U,V,W UND G IN SUBR. DEOBI. (NUR F. LES)
C             KANAL13: FORMATIERTES LESEN DER EINTRITTS-PROFILE
C                      U,V,W UND G  IN SUBR. DEICI. (NUR F. LES)
C             KANAL14: FORMATIERTES SCHREIBEN DER EINTRITTS-PROFILE
C                      U,V,W UND G IN SUBR. DEOCI. (NUR F. LES)
C             KANAL21: LESEN DER PARTIKEL-POSITIONEN UND -ORIENTIERUNGEN
C             KANAL22: SCHREIBEN DER ""
C             KANAL23: ZEITSCHREIBE DER ""
C             KANAL24: FORMATIERTES LESEN U. SCHREIBEN F. PARTIKEL-
C                      PROGRAMM (SIEHE D.A. GEORG EDER)
C             KANAL31: INPUT MEHRGITTERVERFAHREN
C             KANAL36: OUTPUT MEHRGITTERVERFAHREN
C             KANAL26: FORMATIERTES LESEN U. SCHREIBEN F. PARTIKEL-
C                      PROGRAMM (SIEHE D.A. GEORG EDER)
C             KANAL27: FORMATIERTES SCHREIBEN EINZELNER GESCHWINDIG-
C                      KEITSPROFILE (FUER ZEITREIHEN)
C             KANAL44: MELDUNGEN ZUR KONTROLLE DES RAUSSCHREIBENS
C                      VON ZEITRECORDS
C             KANAL52: STEUERDATEN FUER RAUSSCHREIBEN VON ZEITRECORDS
C
C             KANAL57: LESEN DER CPU-ZEIT-SCHRANKE
C
C  DEFINE-DIREKTIVEN     : ADBA,   CRAY,   CYBER,  EULER,  LEAPF,
C                          XHOMOG, YHOMOG,  RECORD, HTMG 
C
C  VERS:  09.01.87 (HW)  : D00990 (VOELLIG NEUE REVISION)
C         18. 2.92 (MM)  : REFERENZVERSION FUER KOMPLETT NEUE VERSION
C         02.04.92 (MM)  : VOELLIG NEUE VERSION, WANDERKENNUNG GESCHIEHT
C                          B-FELD, SCHLEIFEN LAUFEN DURCH KOERPER DURCH,
C                          BELIEBIG GEFORMTER KOERPER DADURCH MOEGLICH
C        10. 6.92  (WA)  : MEHRGITTERVERFAHREN FUER DRUCKKORREKTUR 
C                          EINGEFUEHRT
C            2.93  (MM)  : MIT UMBAU FUER BELIEBIG VIELE GITTER BEGONNEN
C         01.12.03 (TB)  : SCALAR TRANSPORT IMPLEMENTED
C         29.01.03 (TB)  : _KSR_ REMOVED
C
C*STARLET***************************************************************
C
C                                       GESAMTSPEICHER FUER EIN 3D-FELD
      PARAMETER   ( IDIM3D =   3000000 )
C                                       GESAMTSPEICHER FUER EIN 2D-FELD
      PARAMETER   ( IDIM2D =    200000 )
C                                       GESAMTSPEICHER FUER EIN 1D-FELD
      PARAMETER   ( IDIM1D =    6000 )
C                                       LAENGSTE 1D-DIMENSIONIERUNG 
C                                       EINES GITTERS
      PARAMETER   ( IDIMF  =    6000 )
C                                       GESAMTSPEICHER FUER EIN 
C                                       AUSWERTE-FELD
      PARAMETER   ( IDIMA  =  3000000 )
C                                       DIM. FUER LINIENFELDER
      PARAMETER   ( IDIM1L = 10 )
      PARAMETER   ( IDIM2L = 10 )
C                                       ANZAHL DER RANDSCHICHTEN
      PARAMETER   ( NBND   =  2     )
C                                       ANZAHL DER BUFFERSCHICHTEN
      PARAMETER   ( NBUF   =  3     )
C
C                                 MULTIGRID-VERWALTUNG
C
      

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


      COMMON /COMGRID/
     &               NGRID,    NGRDOLD,   NGRDSET, NGRDDFD,
     &               NAUFP,    NAUFOLD,
     &                 KMX,  JMX,  IMX,
     &                 KMXA, JMXA, IMXA,
     &                IP3D, IP2D, IP1D, IPBB,IPB3, IPBU,
     &                NOF3D,NOF2D,NOF1D,NOFBB,NOFB3,NOFBU,
     &                IPA,   IP1L,  IP2L,
     &                NOFA,  NOF1L, NOF2L,  
     &                 IC1,  IC2,  JC1,  JC2,  KC1,  KC2

      INTEGER
     &         KMX(MAXGRIDS),     JMX(MAXGRIDS),       IMX(MAXGRIDS),
     &        KMXA(MAXGRIDS),    JMXA(MAXGRIDS),      IMXA(MAXGRIDS),
     &        IP3D(MAXGRIDS),    IP2D(MAXGRIDS),      IP1D(MAXGRIDS),
     &        IPBB(MAXGRIDS),    IPB3(MAXGRIDS),      IPBU(MAXGRIDS),
     &        NOF3D,   NOF2D,   NOF1D,    NOFBB,   NOFBU,
     &         IPA(MAXGRIDS),    IP1L(MAXGRIDS),      IP2L(MAXGRIDS),
     &         NOFA,   NOF1L,   NOF2L,    NAUFP,
     &         IC1(MAXGRIDS),     IC2(MAXGRIDS),
     &         JC1(MAXGRIDS),     JC2(MAXGRIDS),
     &         KC1(MAXGRIDS),     KC2(MAXGRIDS)
 


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

C
C                            DARF NUR NACH COMGRID STEHEN !!!!
C
      COMMON /COMGVP/
     &               IPCORR, IPCGES,   NDIVLEPS,  LDIVLEPS, DIVGMX

      INTEGER
     &         IPCORR(MAXGRIDS),   IPCGES(MAXGRIDS),   
     &       NDIVLEPS(MAXGRIDS), LDIVLEPS(MAXGRIDS)

      REAL   DIVGMX(MAXGRIDS)
 
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                                 PARAMETER-STATEMENTS F. PARTIKEL-PROG.
C
      PARAMETER   (NDGL =  3, NMAX =   20, NPRUZL = 50000)
      PARAMETER   (NRZUL = 2, NWTVGL = 1)
      PARAMETER   (EPS = 0.1)
C
C
C                                 PARAMETER-STATEMENTS F. KORRELATIONS-
C                                 FUNKTIONEN, -KOEFFIZIENTEN, LEISTUNGS-
C                                 DICHTESPEKTREN UND HAEUFIGKEITSVER-
C                                 TEILUNGEN
C                                 NUR  ILIMXP  DARF VERAENDERT WERDEN !!
C
      PARAMETER   (JJ1L =  1, JJ2L =  2, ILIMXP =  30, ISLIDI = 500)
      PARAMETER   (LINFB= 25)
C
C                                 PARAMETER-STATEMENTS F. KORRELATIONS-
C                                 FUNKTIONEN UND -KOEFFIZIENTEN
C
CSR22      PARAMETER   (KKXL =  IIP + 1 - 2*NBND + LINFB)
CSR22      PARAMETER   (KKYL =  JJP + 1 - 2*NBND + LINFB)
CSR22      PARAMETER   (KKZL =  KKP + 1 - 2*NBND + LINFB)
      PARAMETER   (KKXL =  IDIMF + 1 - 2*NBND + LINFB)
      PARAMETER   (KKYL =  IDIMF + 1 - 2*NBND + LINFB)
      PARAMETER   (KKZL =  IDIMF + 1 - 2*NBND + LINFB)
C
C                                 PARAMETER-STATEMENTS F. LEISTUNGS-
C                                 DICHTESPEKTREN
C
CSR22      PARAMETER   (KSXL = (IIP + 1 - 2*NBND)/2 + 1 + LINFB)
CSR22      PARAMETER   (KSYL = (JJP + 1 - 2*NBND)/2 + 1 + LINFB)
CSR22      PARAMETER   (KSZL = (KKP + 1 - 2*NBND)/2 + 1 + LINFB)
      PARAMETER   (KSXL = (IDIMF + 1 - 2*NBND)/2 + 1 + LINFB)
      PARAMETER   (KSYL = (IDIMF + 1 - 2*NBND)/2 + 1 + LINFB)
      PARAMETER   (KSZL = (IDIMF + 1 - 2*NBND)/2 + 1 + LINFB)
C
C                                 PARAMETER-STATEMENTS F. HAEUFIGKEITS-
C                                 VERTEILUNGEN DER INKLINATIONSWINKEL
C                                 DER WIRBELVEKTOREN
C
      PARAMETER   (KLASS=  72)
      PARAMETER   (KH0L = KLASS + LINFB)
C
C                                 MAXIMALE DIMENSION DER "LINIEN-FELDER"
C                                 IN K-RICHTUNG
C
          PARAMETER   (KKML =  KKXL )

C
C                                 PARAMETER-STATEMENTS FUER DEN
C                                 GITTERGENERATOR  G R D F M I
C
C     PARAMETER    (IDIMF=       IDIM1D         )
      PARAMETER    (NGRP1  = IDIMF + 1,  NGRP3  = IDIMF + 3)
C

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


      PARAMETER (NPPHYS_MAX=1)
      COMMON /COPHYSPAR/
     $                  MTURB,  TU_LEVEL,   RHO,    GMOL,   UGRID,
     $                  VREF,   EXPON,  UTAUX,
     $                  CIDUFR, UFRCON, UFRFREQ, DELTA,  XREF,
     $                  IPP,    JPP,    KPP,
     $                  XPER,   YPER,   ZPER, 
     $                  NXPER,  NYPER,  NZPER,
     $                  NPPHYS, PPHYS, XPPHYS

      INTEGER           IPP,  JPP,  KPP,   NPPHYS
     
      REAL              TU_LEVEL,  RHO,   GMOL,   UGRID,  VREF,  
     $                  EXPON, UFRCON, 
     $                  PPHYS(NPPHYS_MAX), XPPHYS(NPPHYS_MAX)

      CHARACTER (LEN=16)      CIDUFR


      COMMON /COPRINT/
     &                IPRGRID,    PRINTFORMAT,
     &                ISELPREC,
     &                LLU,  LLV,  LLW,  LLP,  LLG,  LLB,
     &                NPRINTPOS,
     &                ISTPR,   JSTPR,   KSTPR,
     &                PRINTEBE,   IPRINTEBE,
     &                IP1,  IP2,  IPS,
     &                JP1,  JP2,  JPS,
     &                KP1,  KP2,  KPS,
     &                CIP1

      INTEGER
     &        IPRGRID,    ISELPRINT,
     &        LLU,  LLV,  LLW,  LLP,  LLG,  LLB,
     &        NPRINTPOS,
     &        ISTPR(8),   JSTPR(8),   KSTPR(8),
     &        IPRINTEBE,
     &        IP1,  IP2,  IPS,
     &        JP1,  JP2,  JPS,
     &        KP1,  KP2,  KPS,
     &        CIP1

      CHARACTER (LEN=16) PRINTFORMAT
      CHARACTER (LEN=8)  PRINTEBE


      PARAMETER ( NBSP  =  8 )
      
      COMMON /COGRDDEF/
     &                 NBX,      NBY,      NBZ,
     &                 NTXB,     NTYB,     NTZB,
     &                 DLX,      DLY,      DLZ,
     &                 NLXB,     NLYB,     NLZB,
     &                 SX,       SY,       SZ,
     &                 NRXB,     NRYB,     NRZB,
     &                 DRX,      DRY,      DRZ

      INTEGER
     &       NBX(MAXGRIDS),       NBY(MAXGRIDS),       NBZ(MAXGRIDS),
     & NTXB(NBSP,MAXGRIDS), NTYB(NBSP,MAXGRIDS), NTZB(NBSP,MAXGRIDS),
     & NLXB(NBSP,MAXGRIDS), NLYB(NBSP,MAXGRIDS), NLZB(NBSP,MAXGRIDS),
     & NRXB(NBSP,MAXGRIDS), NRYB(NBSP,MAXGRIDS), NRZB(NBSP,MAXGRIDS)

      REAL
     &  DLX(NBSP,MAXGRIDS),  DLY(NBSP,MAXGRIDS),  DLZ(NBSP,MAXGRIDS),
     &   SX(NBSP,MAXGRIDS),   SY(NBSP,MAXGRIDS),   SZ(NBSP,MAXGRIDS),
     &  DRX(NBSP,MAXGRIDS),  DRY(NBSP,MAXGRIDS),  DRZ(NBSP,MAXGRIDS)

      
      COMMON /COGRDPRO/ LEVEL,LCHILD,XMIN,YMIN,ZMIN,XTOT,YTOT,ZTOT
      COMMON /COGRDPRO/ XHOMOG,YHOMOG,ZHOMOG,NXGRAE,NYGRAE,NZGRAE
      COMMON /COGRDPRO/ GRADPX,UBULKX,LTST,LVP,LSCAI,LPLEVEL,LPOISSONDIR
      COMMON /COGRDPRO/ LSLICE,NXSLICE,NYSLICE,NZSLICE,NVPGRIDS
      COMMON /COGRDPRO/ CONV1SANF,CONV1SEND,TRANSLES1,TRANSLES2
      COMMON /COGRDPRO/ GRADPXOLD

      INTEGER LEVEL(MAXGRIDS),NVPGRIDS(MAXGRIDS)
      INTEGER NXGRAE(MAXGRIDS),NYGRAE(MAXGRIDS),NZGRAE(MAXGRIDS)
      INTEGER NXSLICE(MAXGRIDS),NYSLICE(MAXGRIDS),NZSLICE(MAXGRIDS)

      REAL XTOT(MAXGRIDS),YTOT(MAXGRIDS),ZTOT(MAXGRIDS)
      REAL XMIN(MAXGRIDS),YMIN(MAXGRIDS),ZMIN(MAXGRIDS)
      REAL GRADPX(MAXGRIDS),UBULKX(MAXGRIDS)
      REAL CONV1SANF(MAXGRIDS),CONV1SEND(MAXGRIDS)
      REAL TRANSLES1(MAXGRIDS),TRANSLES2(MAXGRIDS),GRADPXOLD(MAXGRIDS)

      LOGICAL LCHILD(MAXGRIDS),LSLICE(MAXGRIDS),LPOISSONDIR(MAXGRIDS)
      LOGICAL XHOMOG(MAXGRIDS),YHOMOG(MAXGRIDS),ZHOMOG(MAXGRIDS)
      LOGICAL LTST(MAXGRIDS),LVP(MAXGRIDS),LSCAI(MAXGRIDS)
      LOGICAL LPLEVEL(MAXGRIDS)
      
      COMMON /COGRDCON/
     &                 IVPCHILD,
     &                 IPARENT,ISLPAR,
     &                 IPOSITION, JPOSITION, KPOSITION,
     &                 NOFSLCHILDS,IGRDOFSLCHILD,
     &                 ISLPOS, JSLPOS, KSLPOS,
     &                 IFRNBR, IBANBR, IRINBR, ILENBR,
     &                 IBONBR, ITONBR


      INTEGER
     &       IVPCHILD(MAXGRIDS),
     &       IPARENT(MAXGRIDS),ISLPAR(MAXGRIDS),
     & IPOSITION(MAXGRIDS), JPOSITION(MAXGRIDS), KPOSITION(MAXGRIDS),
     & NOFSLCHILDS(MAXGRIDS),IGRDOFSLCHILD(MAXGRIDS,MAXGRIDS),
     &    ISLPOS(MAXGRIDS), JSLPOS(MAXGRIDS), KSLPOS(MAXGRIDS),
     &       IFRNBR(MAXBOCONDS,MAXGRIDS), IBANBR(MAXBOCONDS,MAXGRIDS),
     &       IRINBR(MAXBOCONDS,MAXGRIDS), ILENBR(MAXBOCONDS,MAXGRIDS),
     &       IBONBR(MAXBOCONDS,MAXGRIDS), ITONBR(MAXBOCONDS,MAXGRIDS)

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

      COMMON /CONLES/  CONV2S,CMUE,CAPPA,ECONST
      SAVE   /CONLES/

      COMMON /CONPOT/  CWA,CWB,CPO1,CPO2,CPO3,CPO4,CPO5,CPO6,CPO7,CPO8,
     $                 CPO9,CPO10,CPO11,CPO12
      SAVE   /CONPOT/

      COMMON /CLINOU/

     $                 LDIB,   LDIC,   LDIEIB, LDIEIC,
     $                 LDOB,   LDOC,   LDOEIB, LDOEIC

      SAVE   /CLINOU/
      LOGICAL

     $                 LDIB,   LDIC,   LDIEIB, LDIEIC,
     $                 LDOB,   LDOC,   LDOEIB, LDOEIC
C
      CHARACTER (LEN=8)   CIDENT(10),   CIDEND(10)
C

      REAL        U (   IDIM3D  ), V (   IDIM3D  ), W (   IDIM3D  ),
     $            UO(   IDIM3D  ), VO(   IDIM3D  ), WO(   IDIM3D  ),
     $            P (   IDIM3D  ), G (   IDIM3D  ), B (   IDIM3D  ),
     $            DP(   IDIM3D  ),
     $            UP(   IDIM3D  ), VP(   IDIM3D  ), WP(   IDIM3D  ),
     $           WCU(  IDIM3D   ),WCV(   IDIM3D  ),WCW(   IDIM3D  ),
     $           H3D1( IDIM3D   ),H3D2(  IDIM3D  ),H3D3(  IDIM3D  ),
     $            RDDX(IDIM1D),     RDDY(IDIM1D),     RDDZ(IDIM1D),
     $             RDX(IDIM1D),      RDY(IDIM1D),      RDZ(IDIM1D),
     $             DDX(IDIM1D),      DDY(IDIM1D),      DDZ(IDIM1D),
     $              DX(IDIM1D),       DY(IDIM1D),       DZ(IDIM1D),
     $               X(IDIM1D),        Y(IDIM1D),        Z(IDIM1D),

     $            RIDENT(100),     RIDEND(100),     DIVG(IDIM2D) 
     
      INTEGER  numxloc, numvar, numzloc, req_xloc
      PARAMETER ( numxloc = 5 )
      PARAMETER ( numvar = 4 )
      PARAMETER ( numzloc = 3 )
      PARAMETER ( numzloc_stream = 3 )

      INTEGER MYID_loc(numxloc),x_loc(numxloc),z_loc(numzloc)

      INTEGER yloc_stream, zloc_stream(numzloc_stream)

      INTEGER INUNIT, i, MYID_temp, x_loc_temp(numxloc)
     
C
      REAL SDIV(IDIM3D),BP(IDIM3D),BU(IDIM3D),BV(IDIM3D),BW(IDIM3D)


      REAL    TAU11( 1 ), TAU12( 1 ), TAU13( 1 ),
     $        TAU21( 1 ), TAU22( 1 ), TAU23( 1 ),
     $        TAU31( 1 ), TAU32( 1 ), TAU33( 1 )



      REAL     COEFFX(1),     COEFFY(1),     COEFFZ(1),
     $         COEFDX(1),     COEFDY(1),     COEFDZ(1),
     $         LCOL  (1),     DIAG  (1),     RCOL(1),
     $         UZ    (1),     RSP (1),
     $         RSGS3 (1),     FAKTOR(1)



C                             FELDER FUER DIE BERECHNUNG DER
C                             ABLEITUNG FUER DAS KOMPAKTVERFAHREN
c         oder 3D-Felder fuer koeffizienten A-Matrix bei iterative Verfahren
c                          AS-FDUI;AT-FDVI;AB-FDWI
      REAL        FUI(1), FVI(1), FWI(1)





      REAL FUJ(IDIM3D),FVJ(IDIM3D),FWJ(IDIM3D)


      REAL FUK(1),FVK(1),FWK(1)
C#if defined _PREPROC_
C      REAL HPX(IDIM3D)
C#endif

C--------------------------------------------------------

      REAL GEOVP ( IDIM3D )
      REAL GSAB( IDIM1D ),GSAT( IDIM1D ),GSAW( IDIM1D ),
     $     GSAE( IDIM1D ),GSAS( IDIM1D ),GSAN( IDIM1D )
      REAL RES  ( IDIM3D ), SIPLB ( IDIM3D ), SIPUT( IDIM3D ),
     $     SIPLW( IDIM3D ), SIPUE ( IDIM3D ), SIPLS( IDIM3D ),
     $     SIPUN( IDIM3D ), SIPLPR( IDIM3D )


      REAL        YHILF1(NDGL),    YHILF2(NDGL),    YMAT(NMAX,NDGL),
     $            K1(NDGL),        K2(NDGL),        K3(NDGL),
     $            K4(NDGL),        YV(NDGL)
C
      REAL        DIV(   IDIM3D  )

C                                 ACHTUNG !!! TK- UND TE-FELD WERDEN
C                                 MOMENTAN NICHT BENOETIGT
C
      REAL        TK( 1, 1, 1),    TE( 1, 1, 1)
C
      INTEGER     FINT,CHILDREN, 
     $            IIDENT(100),     IIDEND(100)
      INTEGER     IWTEND(NMAX),    ISTART(NRZUL),   MSTART(NRZUL)
C
      LOGICAL     LPRLE
      LOGICAL     PRUFUL, LPART

      REAL  UI1(IDIM2D*2),  VI1(IDIM2D*2),  WI1(IDIM2D*2),
     $      UI2(IDIM2D*2),  VI2(IDIM2D*2),  WI2(IDIM2D*2),
     $      GI1(IDIM2D),    GI2(IDIM2D),
     $      UFR(IDIM2D*2), VFR(IDIM2D*2), WFR(IDIM2D*2),
     $      PFR(IDIM2D*2), GFR(IDIM2D*2),
     $      UTO(IDIM2D*2),  VTO(IDIM2D*2),  WTO(IDIM2D*2),
     $      PTO(IDIM2D*2),  GTO(IDIM2D*2),
     $      UBACK(IDIM2D*2),  VBACK(IDIM2D*2),  WBACK(IDIM2D*2),
     $      VRI(IDIM2D*2),  
     $      UBO(IDIM2D*2),  VBO(IDIM2D*2),  WBO(IDIM2D*2),
     $      PBO(IDIM2D*2),  GBO(IDIM2D*2),
     $      UBA(IDIM2D*2), VBA(IDIM2D*2), WBA(IDIM2D*2),
     $      PBA(IDIM2D*2),  GBA(IDIM2D*2)

C                                            FIELDS FOR SCALAR TRANSPORT
C                                            FIELDS FOR DIV ANALYSIS
      INTEGER 
     $ IDIVMAX(MAXGRIDS),JDIVMAX(MAXGRIDS),KDIVMAX(MAXGRIDS),
     $ GRDIVMAX(MAXGRIDS)
      REAL    XDIVMAX(MAXGRIDS),YDIVMAX(MAXGRIDS),ZDIVMAX(MAXGRIDS)
C
C                   FELDER FUER FEEDBACK-ROUTINEN (Jens Neumann)
C
C      REAL        VAROPT(2,IDIM2D*8)

C
C                                 HILFSFELDER FUER EFVISC
C
      REAL        DUDY(IDIM2D),  DUDZ(IDIM2D),
     $            DVDX(IDIM2D),  DVDZ(IDIM2D),
     $            DWDX(IDIM2D),  DWDY(IDIM2D),
     $            CONV1S(IDIM1D)

 
      real  cdelta(1)

      real
     $  hilf0(1),hilf1(1),hilf2(1),
     $  hilf3(1),hilf4(1),hilf5(1),
     $  hilf6(1),
     $  lux(1),luy(1),luz(1),
     $  lvy(1),lvz(1),lwz(1),
     $  mux(1),muy(1),muz(1),
     $  mvy(1),mvz(1),mwz(1),
     $  fa11(1),fa12(1),fa13(1),
     $  fa21(1),fa22(1),fa23(1),
     $  fa31(1),fa32(1),fa33(1)

      real
     $ dudx3(1),dudy3(1),dudz3(1),
     $ dvdx3(1),dvdy3(1),dvdz3(1),
     $ dwdx3(1),dwdy3(1),dwdz3(1),
     $ norms(1),
     $ sdudx3(1),sdudy3(1),sdudz3(1),
     $ sdvdy3(1),sdvdz3(1),sdwdz3(1),
     $ uc(1),vc(1),wc(1),
     $ uuc(1),uvc(1),uwc(1),
     $ vvc(1),vwc(1),wwc(1)



C                                 ARRAYDIMENSIONIERUNG FUER DEN
C                                 GITTERGENERATOR  G R D F M I
C
	  INTEGER IGRHF(IDIMF,1)
C
C
      REAL    GRX0(IDIMF),     GRX1(IDIMF)
      REAL GRHF(IDIMF,100)
C
C                                 FELDER FUER DEN DIREKTEN POISSONLOESER
C
      REAL
     $           FPSFAK(IDIM2D),FCOSMY(IDIM1D),FCOSNY(IDIM1D)
C
C                                 FELDER FUER DIE FOURIERTRANSFORMATION
C
      INTEGER IFFTX(19,MAXGRIDS),IFFTY(19,MAXGRIDS)
      INTEGER IPERMUX(IDIM1D)	     ,IPERMUY(IDIM1D)
      REAL    RFFTX(IDIM1D*2)     ,RFFTY(IDIM1D*2)  
C
C                                 DIE FOLGENDEN FELDER WERDEN ZUR
C                                 AUSWERTUNG FUER DIE LARGE-EDDY-
C                                 SIMULATION BENOETIGT
C
      INTEGER  ISELEA (2,752)      , ISELEP (2,752),
     $         ISAMPA (752,MAXGRIDS),ISAMPP (752,MAXGRIDS),
     $         ISLINP (ISLIDI)     , ISLINA (ISLIDI),
     $         ISLINI (ISLIDI)

      REAL     HILF   (   IDIM3D     ),    HILFL  (KKML,JJ2L,ILIMXP ),
     $         FVT     (0:IDIMF),          HVK    (KLASS,JJ2L)

      REAL       FELD1(IDIM3D),FELD2(IDIM3D),FELD3(IDIM3D)

      REAL       HILF3D1(IDIM3D),HILF3D2(IDIM3D),HILF3D3(IDIM3D)


C-------------------------------- FELDER FUER DECONVOLVE SCA -TEST by FLORIAN 24.11.2003
C---------------------------------------------------------------------------------------


C                                 FUER DIE KONTROLLE DES RECHENVORGANGES
C
      REAL     EPSU(MAXGRIDS),  EPSV(MAXGRIDS),  EPSW(MAXGRIDS)
      REAL     ESUMG(MAXGRIDS),  ESUMS(MAXGRIDS)
      REAL     WSSX(MAXGRIDS),  WSSY(MAXGRIDS),  WSSZ(MAXGRIDS)
      REAL     WNSX(MAXGRIDS),  WNSY(MAXGRIDS),  WNSZ(MAXGRIDS)
      REAL     WALLSSX(6,MAXGRIDS),
     $         WALLSSY(6,MAXGRIDS),
     $         WALLSSZ(6,MAXGRIDS)
      REAL     UBULK(MAXGRIDS)

C
C
C
C           DIMENSIONIERUNGEN FUER DIE STATISTISCHE AUSWERTUNG
C
C                  LEVEL: <U>, <V>, <W>,<T> AND <P>
C
C
      REAL     AU     (   IDIMA         ), SU     (   IDIMA         ),
*     REAL     AU     (        1        ), SU     (        1        ),
*    $         UFG    (   IDIM3D     ),
     $         UFG    (       1      ),
*    $         AURG   (   IDIMA         ), SURG   (   IDIMA         ),
     $         AURG   (        1        ), SURG   (        1        ),
*    $         AUSKG  (   IDIMA         ), SUSKG  (   IDIMA         ),
     $         AUSKG  (        1        ), SUSKG  (        1        ),
*    $         AUFLG  (   IDIMA         ), SUFLG  (   IDIMA         ),
     $         AUFLG  (        1        ), SUFLG  (        1        ),
*    $         ARXUU  (   IDIM2L        ), SRXUU  (     IDIM1L      ),
     $         ARXUU  (        1        ), SRXUU  (        1        ),
*    $         ARYUU  (   IDIM2L        ), SRYUU  (     IDIM1L      ),
     $         ARYUU  (        1        ), SRYUU  (        1        ),
*    $         ARZUU  (   IDIM2L        ), SRZUU  (     IDIM1L      ),
     $         ARZUU  (        1        ), SRZUU  (        1        ),
*    $         ASXUU  (     IDIM1L      ), SSXUU  (     IDIM1L      ),
     $         ASXUU  (        1        ), SSXUU  (        1        ),
*    $         ASYUU  (     IDIM1L      ), SSYUU  (     IDIM1L      ),
     $         ASYUU  (        1        ), SSYUU  (        1        ),
*    $         ASZUU  (     IDIM1L      ), SSZUU  (     IDIM1L      )
     $         ASZUU  (        1        ), SSZUU  (        1        )
      REAL     AV     (   IDIMA         ), SV     (   IDIMA         ),
*     REAL     AV     (        1        ), SV     (        1        ),
*    $         VFG    (  IDIM3D      ),
     $         VFG    (       1      ),
*    $         AVRG   (   IDIMA         ), SVRG   (   IDIMA         ),
     $         AVRG   (        1        ), SVRG   (        1        ),
*    $         AVSKG  (   IDIMA         ), SVSKG  (   IDIMA         ),
     $         AVSKG  (        1        ), SVSKG  (        1        ),
*    $         AVFLG  (   IDIMA         ), SVFLG  (   IDIMA         ),
     $         AVFLG  (        1        ), SVFLG  (        1        ),
*    $         ARXVV  (   IDIM2L        ), SRXVV  (     IDIM1L      ),
     $         ARXVV  (        1        ), SRXVV  (        1        ),
*    $         ARYVV  (   IDIM2L        ), SRYVV  (     IDIM1L      ),
     $         ARYVV  (        1        ), SRYVV  (        1        ),
*    $         ARZVV  (   IDIM2L        ), SRZVV  (     IDIM1L      ),
     $         ARZVV  (        1        ), SRZVV  (        1        ),
*    $         ASXVV  (     IDIM1L      ), SSXVV  (     IDIM1L      ),
     $         ASXVV  (        1        ), SSXVV  (        1        ),
*    $         ASYVV  (     IDIM1L      ), SSYVV  (     IDIM1L      ),
     $         ASYVV  (        1        ), SSYVV  (        1        ),
*    $         ASZVV  (     IDIM1L      ), SSZVV  (     IDIM1L      )
     $         ASZVV  (        1        ), SSZVV  (        1        )
      REAL     AW     (   IDIMA         ), SW     (   IDIMA         ),
*     REAL     AW     (        1        ), SW     (        1        ),
*    $         WFG    (  IDIM3D      ),
     $         WFG    (       1      ),
*    $         AWRG   (   IDIMA         ), SWRG   (   IDIMA         ),
     $         AWRG   (        1        ), SWRG   (        1        ),
*    $         AWSKG  (   IDIMA         ), SWSKG  (   IDIMA         ),
     $         AWSKG  (        1        ), SWSKG  (        1        ),
*    $         AWFLG  (   IDIMA         ), SWFLG  (   IDIMA         ),
     $         AWFLG  (        1        ), SWFLG  (        1        ),
*    $         ARXWW  (   IDIM2L        ), SRXWW  (     IDIM1L      ),
     $         ARXWW  (        1        ), SRXWW  (        1        ),
*    $         ARYWW  (   IDIM2L        ), SRYWW  (     IDIM1L      ),
     $         ARYWW  (        1        ), SRYWW  (        1        ),
*    $         ARZWW  (   IDIM2L        ), SRZWW  (     IDIM1L      ),
     $         ARZWW  (        1        ), SRZWW  (        1        ),
*    $         ASXWW  (     IDIM1L      ), SSXWW  (     IDIM1L      ),
     $         ASXWW  (        1        ), SSXWW  (        1        ),
*    $         ASYWW  (     IDIM1L      ), SSYWW  (     IDIM1L      ),
     $         ASYWW  (        1        ), SSYWW  (        1        ),
*    $         ASZWW  (     IDIM1L      ), SSZWW  (     IDIM1L      )
     $         ASZWW  (        1        ), SSZWW  (        1        )
      REAL     AP     (   IDIMA         ), SP     (   IDIMA         ),
*     REAL     AP     (        1        ), SP     (        1        ),
*    $         PFG    (  IDIM3D      ),
     $         PFG    (       1      ),
*    $         APRG   (   IDIMA         ), SPRG   (   IDIMA         ),
     $         APRG   (        1        ), SPRG   (        1        ),
*    $         APSKG  (   IDIMA         ), SPSKG  (   IDIMA         ),
     $         APSKG  (        1        ), SPSKG  (        1        ),
*    $         APFLG  (   IDIMA         ), SPFLG  (   IDIMA         )
     $         APFLG  (        1        ), SPFLG  (        1        )
*     REAL     AEFG   (   IDIMA         ), SEFG   (   IDIMA         ),
      REAL     AEFG   (        1        ), SEFG   (        1        ),
*    $         AEFS   (   IDIMA         ), SEFS   (   IDIMA         )
     $         AEFS   (        1        ), SEFS   (        1        )
*     REAL     ADFG   (   IDIMA         ), SDFG   (   IDIMA         ),
      REAL     ADFG   (        1        ), SDFG   (        1        ),
*    $         ADUDX2 (   IDIMA         ), SDUDX2 (   IDIMA         ),
     $         ADUDX2 (        1        ), SDUDX2 (        1        ),
*    $         ADUDY2 (   IDIMA         ), SDUDY2 (   IDIMA         ),
     $         ADUDY2 (        1        ), SDUDY2 (        1        ),
*    $         ADUDZ2 (   IDIMA         ), SDUDZ2 (   IDIMA         ),
     $         ADUDZ2 (        1        ), SDUDZ2 (        1        ),
*    $         ADVDX2 (   IDIMA         ), SDVDX2 (   IDIMA         ),
     $         ADVDX2 (        1        ), SDVDX2 (        1        ),
*    $         ADVDY2 (   IDIMA         ), SDVDY2 (   IDIMA         ),
     $         ADVDY2 (        1        ), SDVDY2 (        1        ),
*    $         ADVDZ2 (   IDIMA         ), SDVDZ2 (   IDIMA         ),
     $         ADVDZ2 (        1        ), SDVDZ2 (        1        ),
*    $         ADWDX2 (   IDIMA         ), SDWDX2 (   IDIMA         ),
     $         ADWDX2 (        1        ), SDWDX2 (        1        ),
*    $         ADWDY2 (   IDIMA         ), SDWDY2 (   IDIMA         ),
     $         ADWDY2 (        1        ), SDWDY2 (        1        ),
*    $         ADWDZ2 (   IDIMA         ), SDWDZ2 (   IDIMA         )
     $         ADWDZ2 (        1        ), SDWDZ2 (        1        )
*     REAL     AUFWFG (   IDIMA         ), SUFWFG (   IDIMA         ),
      REAL     AUFWFG (        1        ), SUFWFG (        1        ),
*    $         AUFWFS (   IDIMA         ), SUFWFS (   IDIMA         ),
     $         AUFWFS (        1        ), SUFWFS (        1        ),
*    $         AUFWFM (   IDIMA         ), SUFWFM (   IDIMA         )
     $         AUFWFM (        1        ), SUFWFM (        1        )
*     REAL     AVFWFG (   IDIMA         ), SVFWFG (   IDIMA         ),
      REAL     AVFWFG (        1        ), SVFWFG (        1        ),
*    $         AVFWFS (   IDIMA         ), SVFWFS (   IDIMA         ),
     $         AVFWFS (        1        ), SVFWFS (        1        ),
*    $         AVFWFM (   IDIMA         ), SVFWFM (   IDIMA         ),
     $         AVFWFM (        1        ), SVFWFM (        1        ),
*    $         AUFVFG (   IDIMA         ), SUFVFG (   IDIMA         ),
     $         AUFVFG (        1        ), SUFVFG (        1        ),
*    $         AUFVFS (   IDIMA         ), SUFVFS (   IDIMA         ),
     $         AUFVFS (        1        ), SUFVFS (        1        ),
*    $         AUFVFM (   IDIMA         ), SUFVFM (   IDIMA         )
     $         AUFVFM (        1        ), SUFVFM (        1        )
*     REAL     OX     (  IDIM3D      ),
      REAL     OX     (       1      ),
*    $         OXFG   (  IDIM3D      ),
     $         OXFG   (       1      ),
*    $         AOX    (   IDIMA         ), SOX    (   IDIMA         ),
     $         AOX    (        1        ), SOX    (        1        ),
*    $         AOXRG  (   IDIMA         ), SOXRG  (   IDIMA         )
     $         AOXRG  (        1        ), SOXRG  (        1        )
*     REAL     OY     (  IDIM3D      ),
      REAL     OY     (       1      ),
*    $         OYFG   (  IDIM3D      ),
     $         OYFG   (       1      ),
*    $         AOY    (   IDIMA         ), SOY    (   IDIMA         ),
     $         AOY    (        1        ), SOY    (        1        ),
*    $         AOYRG  (   IDIMA         ), SOYRG  (   IDIMA         )
     $         AOYRG  (        1        ), SOYRG  (        1        )
*     REAL     OZ     (  IDIM3D      ),
      REAL     OZ     (       1      ),
*    $         OZFG   (  IDIM3D      ),
     $         OZFG   (       1      ),
*    $         AOZ    (   IDIMA         ), SOZ    (   IDIMA         ),
     $         AOZ    (        1        ), SOZ    (        1        ),
*    $         AOZRG  (   IDIMA         ), SOZRG  (   IDIMA         )
     $         AOZRG  (        1        ), SOZRG  (        1        )
*     REAL                                 O2FG   (  IDIM3D      ),
      REAL                                 O2FG   (       1      ),
*    $         AO2    (   IDIMA         ), SO2    (   IDIMA         ),
     $         AO2    (        1        ), SO2    (        1        ),
*    $         AO2RG  (   IDIMA         ), SO2RG  (   IDIMA         )
     $         AO2RG  (        1        ), SO2RG  (        1        )
*     REAL                                 HEFG   (  IDIM3D      ),
      REAL                                 HEFG   (       1      ),
*    $         AHE    (   IDIMA         ), SHE    (   IDIMA         ),
     $         AHE    (        1        ), SHE    (        1        ),
*    $         AHERG  (   IDIMA         ), SHERG  (   IDIMA         )
     $         AHERG  (        1        ), SHERG  (        1        )
      REAL
     $         ATAU11 (        1        ), STAU11 (        1        ),
     $         ATAU12 (        1        ), STAU12 (        1        ),
     $         ATAU13 (        1        ), STAU13 (        1        ),
     $         ATAU22 (        1        ), STAU22 (        1        ),
     $         ATAU23 (        1        ), STAU23 (        1        ),
     $         ATAU33 (        1        ), STAU33 (        1        )
C
C                                  KREUZKORRELATIONSFUNKTIONEN
C                                  ---------------------------
C
*     REAL     ARXUV  (   IDIM2L        ), SRXUV  (     IDIM1L      ),
      REAL     ARXUV  (        1        ), SRXUV  (        1        ),
*    $         ARYUV  (   IDIM2L        ), SRYUV  (     IDIM1L      ),
     $         ARYUV  (        1        ), SRYUV  (        1        ),
*    $         ARXUW  (   IDIM2L        ), SRXUW  (     IDIM1L      ),
     $         ARXUW  (        1        ), SRXUW  (        1        ),
*    $         ARYUW  (   IDIM2L        ), SRYUW  (     IDIM1L      ),
     $         ARYUW  (        1        ), SRYUW  (        1        ),
*    $         ARYVW  (   IDIM2L        ), SRYVW  (     IDIM1L      )
     $         ARYVW  (        1        ), SRYVW  (        1        )
C
C                                  KORRELATIONSKOEFFIZIENTEN
C                                  -------------------------
C
*     REAL     ACXUW  (   IDIM2L        ), SCXUW  (     IDIM1L      ),
      REAL     ACXUW  (        1        ), SCXUW  (        1        ),
*    $         ACZUW  (   IDIM2L        ), SCZUW  (     IDIM1L      )
     $         ACZUW  (        1        ), SCZUW  (        1        )
C
C                                  OBERE DREIECKSMATRIX DES TENSORS
C                                  DER KORRELATIONSFUNKTIONEN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  --------------------------------
C
*     REAL     ARXOXX (   IDIM2L        ), SRXOXX (     IDIM1L      ),
      REAL     ARXOXX (        1        ), SRXOXX (        1        ),
*    $         ARXOXY (   IDIM2L        ), SRXOXY (     IDIM1L      ),
     $         ARXOXY (        1        ), SRXOXY (        1        ),
*    $         ARXOXZ (   IDIM2L        ), SRXOXZ (     IDIM1L      ),
     $         ARXOXZ (        1        ), SRXOXZ (        1        ),
*    $         ARXOYY (   IDIM2L        ), SRXOYY (     IDIM1L      ),
     $         ARXOYY (        1        ), SRXOYY (        1        ),
*    $         ARXOYZ (   IDIM2L        ), SRXOYZ (     IDIM1L      ),
     $         ARXOYZ (        1        ), SRXOYZ (        1        ),
*    $         ARXOZZ (   IDIM2L        ), SRXOZZ (     IDIM1L      )
     $         ARXOZZ (        1        ), SRXOZZ (        1        )
*     REAL     ARYOXX (   IDIM2L        ), SRYOXX (     IDIM1L      ),
      REAL     ARYOXX (        1        ), SRYOXX (        1        ),
*    $         ARYOXY (   IDIM2L        ), SRYOXY (     IDIM1L      ),
     $         ARYOXY (        1        ), SRYOXY (        1        ),
*    $         ARYOXZ (   IDIM2L        ), SRYOXZ (     IDIM1L      ),
     $         ARYOXZ (        1        ), SRYOXZ (        1        ),
*    $         ARYOYY (   IDIM2L        ), SRYOYY (     IDIM1L      ),
     $         ARYOYY (        1        ), SRYOYY (        1        ),
*    $         ARYOYZ (   IDIM2L        ), SRYOYZ (     IDIM1L      ),
     $         ARYOYZ (        1        ), SRYOYZ (        1        ),
*    $         ARYOZZ (   IDIM2L        ), SRYOZZ (     IDIM1L      )
     $         ARYOZZ (        1        ), SRYOZZ (        1        )
*     REAL     ARZOXX (   IDIM2L        ), SRZOXX (     IDIM1L      ),
      REAL     ARZOXX (        1        ), SRZOXX (        1        ),
*    $         ARZOXY (   IDIM2L        ), SRZOXY (     IDIM1L      ),
     $         ARZOXY (        1        ), SRZOXY (        1        ),
*    $         ARZOXZ (   IDIM2L        ), SRZOXZ (     IDIM1L      ),
     $         ARZOXZ (        1        ), SRZOXZ (        1        ),
*    $         ARZOYY (   IDIM2L        ), SRZOYY (     IDIM1L      ),
     $         ARZOYY (        1        ), SRZOYY (        1        ),
*    $         ARZOYZ (   IDIM2L        ), SRZOYZ (     IDIM1L      ),
     $         ARZOYZ (        1        ), SRZOYZ (        1        ),
*    $         ARZOZZ (   IDIM2L        ), SRZOZZ (     IDIM1L      )
     $         ARZOZZ (        1        ), SRZOZZ (        1        )
C
C                                  LEISTUNGSDICHTESPEKTREN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  ----------------------------
C
*     REAL     ASXOXX (     IDIM1L      ), SSXOXX (     IDIM1L      ),
      REAL     ASXOXX (        1        ), SSXOXX (        1        ),
*    $         ASXOYY (     IDIM1L      ), SSXOYY (     IDIM1L      ),
     $         ASXOYY (        1        ), SSXOYY (        1        ),
*    $         ASXOZZ (     IDIM1L      ), SSXOZZ (     IDIM1L      ),
     $         ASXOZZ (        1        ), SSXOZZ (        1        ),
*    $         ASYOXX (     IDIM1L      ), SSYOXX (     IDIM1L      ),
     $         ASYOXX (        1        ), SSYOXX (        1        ),
*    $         ASYOYY (     IDIM1L      ), SSYOYY (     IDIM1L      ),
     $         ASYOYY (        1        ), SSYOYY (        1        ),
*    $         ASYOZZ (     IDIM1L      ), SSYOZZ (     IDIM1L      )
     $         ASYOZZ (        1        ), SSYOZZ (        1        )
C
C                                  HAEUFIGKEITSVERTEILUNG DER
C                                  INKLINATIONSWINKEL DER VORTICITY-
C                                  VEKTOREN
C                                  ---------------------------------
C
*     REAL     AHOZOY (     IDIM2L      ), SHOZOY (     IDIM2L      ),
      REAL     AHOZOY (        1        ), SHOZOY (        1        ),
*    $         AHOZOX (     IDIM2L      ), SHOZOX (     IDIM2L      ),
     $         AHOZOX (        1        ), SHOZOX (        1        ),
*    $         AHOYOX (     IDIM2L      ), SHOYOX (     IDIM2L      )
     $         AHOYOX (        1        ), SHOYOX (        1        )
C
C                                  SCALAR T FIELD STATISTICS
C                                  ---------------------------------


c
c Neue Statistik
c


C
C           DIMENSIONIERUNGEN FUER DIE STATISTISCHE AUSWERTUNG
C
C                LEVEL:
C
C               <UU> G   <VV> G   <WW> G <PP> G
C
C
      REAL     AUUM   (   IDIMA         ), SUUM   (   IDIMA         ),
     &         AVVM   (   IDIMA         ), SVVM   (   IDIMA         ),
     &         AWWM   (   IDIMA         ), SWWM   (   IDIMA         ),
     &         APPM   (   IDIMA         ), SPPM   (   IDIMA         )
c Included stat11.h

C
C           DIMENSIONIERUNGEN FUER DIE STATISTISCHE AUSWERTUNG
C
C                LEVEL:
C
C               <UV> G   <UW> G   <VW> G
C
      REAL     AUVM   (   IDIMA         ), SUVM   (   IDIMA         ),
     &         AUWM   (   IDIMA         ), SUWM   (   IDIMA         ),
     &         AVWM   (   IDIMA         ), SVWM   (   IDIMA         )
c Included stat12.h

C
C           DIMENSIONIERUNGEN FUER DIE STATISTISCHE AUSWERTUNG
C
C                LEVEL:
C
C               <(dU/dX)^2> G   <(dU/dY)^2> G   <(dU/dZ)^2> G
C               <(dV/dX)^2> G   <(dV/dY)^2> G   <(dV/dZ)^2> G
C               <(dW/dX)^2> G   <(dW/dY)^2> G   <(dW/dZ)^2> G
C
      REAL     AUXUXM (   IDIMA         ), SUXUXM (   IDIMA         ),
     &         AUYUYM (   IDIMA         ), SUYUYM (   IDIMA         ),
     &         AUZUZM (   IDIMA         ), SUZUZM (   IDIMA         ),
     &         AVXVXM (   IDIMA         ), SVXVXM (   IDIMA         ),
     &         AVYVYM (   IDIMA         ), SVYVYM (   IDIMA         ),
     &         AVZVZM (   IDIMA         ), SVZVZM (   IDIMA         ),
     &         AWXWXM (   IDIMA         ), SWXWXM (   IDIMA         ),
     &         AWYWYM (   IDIMA         ), SWYWYM (   IDIMA         ),
     &         AWZWZM (   IDIMA         ), SWZWZM (   IDIMA         )
c Included stat13.h





C
C
C
      EQUIVALENCE ( HILF  (1), HILFL (1,1,1))
      
CTBC  DIV NOW USED FOR DIVMAX OUTPUT - NO EQUIVALENCE HILF <-> DIV!!! 
C#ifdef _POISSONDIR_
C      EQUIVALENCE ( HILF  (1), DIV   (1))
C#else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC VORSICHT !!!!
CC     EQUIVALENCE ( HILF  (1), DIV   (1), DP(1))
C      EQUIVALENCE ( HILF  (1), DIV   (1))
C#endif
C
C                                 EQUIVALENCE-STATEMENTS FUER DEN
C                                 GITTERGENERATOR  G R D F M I
C
      EQUIVALENCE ( HILF  (1), GRHF  (1,1)  )
C
      DATA   PRUFUL                            /.FALSE./
C
C
C
C                         KANPR: AUSGABEKANAL FUER
C                         KONTROLLAUSDRUECKE 

      DATA       KANPR               /  8  /

      DATA   RKOMXP                            / 3.0    /
      DATA   IDOBD1, IDOBD2, IDOBD3, IDOBD4    /999, 999, 999, 999/
      DATA   RDOBD1, RDOBD2, RDOBD3, RDOBD4    /9.99,9.99,9.99,9.99/
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

       OPEN (KANPR,FILE='fort.8',FORM='FORMATTED')
C                                 FALLS  KKML "VON HAND"
C                                 DIMENSIONIERT WERDEN MUSS, IST
C                                 EINE UEBERPRUEFUNG UNERLAESSLICH !
C
      IF(KKML  .LT. (MAX0 (KKXL, KKYL, KKZL, KSXL, KSYL, KSZL,
     $                     KH0L)))         THEN
         WRITE (6,*) ' KKML  = ',KKML
         WRITE (6,*) ' KKXL  = ',KKXL,' KKYL = ',KKYL,' KKZL = ',KKZL
         WRITE (6,*) ' KSXL  = ',KSXL,' KSYL = ',KSYL,' KSZL = ',KSZL
         WRITE (6,*) ' KH0L  = ',KH0L
         CALL ERRR (502,' MLET   ')
      END IF      
C
C                                 CP-ZEITUEBERNAHME, EINLESEN VON
C                                 STEUERPARAMETERN, TEST AUF FELD-
C                                 GRENZENUEBERSCHREITUNG UND ERMITTLUNG
C                                 WICHTIGER KONSTANTEN
      CPSEC = GETSEC(0)
C                                 KONSTANTEN
      CALL SETKON
C
C                                 INITIALISIEREN DER MULTIGRID-VERW.
C
       CALL ICOMGRID
CTBCM 290103: LAMDA,PRMOL,PRTURB,TREF,EXPONT,CIDTFR,TFRCON,DELTAT added
       CALL ICOPHYSPAR
       CALL ICOLEVEL
       CALL ICOMGVP
       CALL ICOGRDCON
CTBCM 290103: Comments about new ITYPBOCONDS (15,16) in Header
       CALL ICOBOUND
       CALL ICOGRDDEF
       CALL ICOGRDPRO
       CALL ICOBODY
       DO I=1,MAXGRIDS
          EPSU(I) = 0.0
          EPSV(I) = 0.0
          EPSW(I) = 0.0
          ESUMG(I)= 0.0
          ESUMS(I)= 0.0
       ENDDO
C
C                                 STARTPARAMETER EINLESEN
      CALL STRLES
C
C                                 KONSTANTEN DES TURBULENZMODELLS
      CALL SKONLE (GMOL,RHO,1)
C
      CALL LINOUT (DREAD,  DWRITE, IWRB  )
C
C
C
C                                 AUSGABE EINIGER WICHTIGER KENNGROESSEN
C                                 ZU BEGINN DES LAUFES
C
C     CALL LISTI6 (KMX(1), JMX(1), IMX(1), CIDENT, IIDENT, RIDENT, 6)
C
C                                 VORBELEGUNG DER FELDER FUER DIE
C                                 STATISTISCHE AUSWERTUNG
C                                 -------------------------------
C
      IF(MTURB .GT. 0) THEN
C
C                                 BELEGUNG DES "INHALTSVERZEICHNISSES"
C                                 DER AUSZUWERTENDEN STATIST. GROESSEN
C
         NAUFP = ISLIDI

         CALL SEL0  (ISELEA,ISELEP,ISAMPA,ISAMPP,MAXGRIDS)

         CALL SEL1  (ISELEA,ISELEP,ISAMPA,ISAMPP,MAXGRIDS)

         CALL SEL11 (ISELEP)
         CALL SEL12 (ISELEP)
         CALL SEL13 (ISELEP)
      
C
C                                 BELEGUNG DER AUFPUNKTE FUER DIE BIL-
C                                 DUNG VON KORRELATIONEN, LEISTUNGS-
C                                 DICHTESPEKTREN ETC.
C
         CALL SELAUF  (ISLINP, ISLINA, ISLIDI, ISLINI, ILIMXP, ILIMX,
     $                 IB1, IB2, JB1, JB2, KB, NBND, ILINT0,
     $                 ILINTX, ILINTY, ILINTZ)
         IIDENT (75) = 1
      ENDIF
C
C                                 WENN DREAD DANN EINLESEN VON
C                                 STARTPARAMETERN UND FELDWERTEN, TEST
C                                 AUF PLAUSIBILITAET UND FELDGRENZEN
C
      IF(DREAD) THEN
         IF((LDIB).AND.(.NOT.LDIC)) THEN

            OPEN (1,FILE='fort.1',FORM='UNFORMATTED')
            CALL DIBHEAD (CIDENT,IIDENT,RIDENT,
     $                    ISELEA,ISAMPA,ISLINA,ISLIDI)
         ELSEIF((LDIC).AND.(.NOT.LDIB)) THEN

            OPEN (3,FILE='fort.3',FORM='FORMATTED')

            CALL DICHEAD (CIDENT,IIDENT,RIDENT,
     $                    ISELEA,ISAMPA,ISLINA,ISLIDI)
         ELSE
            CALL ERRR (503,' MLET ')
         ENDIF
      ELSE
         DCONT = .FALSE.
         NGRDOLD = 0
      ENDIF
C


C                                 IDENTIFIKATION:
C
      CALL SETID8 (CIDENT,IIDENT,RIDENT,'AERO WENGLE',VERS,MTSTEP,
     $             NRRUN,MPCORR,KB,JB1,JB2,IB1,IB2,RE,DT,EPCORR,
     $             ZTOT(1),YTOT(1),XTOT(1),
     $             BETA,VREF,EXPON,GRADPX,RHO,GMOL,DREAD,ZCUB
     $		   )
C
      CALL STRSET (IDIM3D,IDIM2D,IDIM1D,NBUF,NBND,
     $               IDIMA,IDIM1L,IDIM2L,IDIMF,IIDENT,
     $             X,Y,Z,ISAMPA)
      CALL BODYCHECK

C050696C
C050696C                           GITTER FUER GEBIETSZERLEGUNG
C050696C
C050696      DO ILEVEL = MINLEVEL,MAXLEVEL
C050696        DO I = 1,NOFSLCED(ILEVEL)
C050696           IGRID = IGRDOFSLCED(I,ILEVEL)
C050696
C050696           CALL SETSLICE (IGRID,IDIM3D,IDIM2D,IDIM1D,
C050696     $                    NBUF,NBND,IDIMA,IDIM1L,IDIM2L,IDIMF,X,Y,Z)
C050696
C050696        ENDDO
C050696      ENDDO
C050696
C050696      CALL SETMPI
C050696
      IF(DREAD) THEN     	
         DO IGRID = 1,NGRDOLD 
             CALL DIGRID (CIDENT,IIDENT,RIDENT,
     $                   DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,
     $              U,V,W,P,G,B,BP,BU,BV,BW,	
     $              ITSTEP,ITTOT,TIMEPH,
     $              IDIMF,NGRP1, NGRP3, IGRHF, GRX0, GRX1, HILF,
     $              IDIM3D,IDIM2D,IDIM1D,NBUF,NBND,IDIMA,IDIM1L,IDIM2L,
     $              IGRID,CONV1S)



CTBC  DO NOT READ OLD STATISTICS WHILE EXTENDING GRID -> OLD STAT. LOST!        
            IF(MTURB .GT. 0) THEN
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   BEI ERGEBNISSEN AUS MLET
C                                      WIRD KEINE STATISTIK EINGELESEN
            IF (IIDENT(5) .GT. 8) THEN

C                                 EINLESEN DER ENSEMBLE-MITTELWERTE
C                                 DES VORANGEGANGENEN LAUFES 
            CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $                   NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID) 
            CALL MGDIMA  (KKA,JJA,IIA,IGRID)
            CALL MGPOINA (IAV,I2L,I1L,IGRID)
                IF(LDIB  ) THEN
            CALL DIBCA   (  1  ,'BINAER  ',IIDENT,
     $                    ILIMXA,IDIBD1,IDIBD2,IDIBD3,IDIBD4,RKOMXA,
     $                    RDIBD1,RDIBD2,RDIBD3,RDIBD4,IGRID,IDIM3D,  

     $ KK,JJ,II,KKA,JJA,IIA,KK,JJ,II,KKA,JJA,IIA,NBND,
C-------------------------- AB HIER, WIE cstaca.h:
     $ ISELEA,ISELEP,
     $ ISAMPA,ISAMPP,ISLINP,ISLINA,ISLIDI,
     $ RIDENT,HILF,B,U,UO,AU,SU,UFG,AURG,SURG,AUSKG,SUSKG,AUFLG,SUFLG,V,
     $ VO,AV,SV,VFG,AVRG,SVRG,AVSKG,SVSKG,AVFLG,SVFLG,W,WO,AW,SW,WFG,
     $ AWRG,SWRG,AWSKG,SWSKG,AWFLG,SWFLG,P,AP,SP,PFG,APRG,SPRG,
     $ APSKG,SPSKG,APFLG,SPFLG,
     $ G,AEFG,SEFG,AEFS,SEFS,ADFG,SDFG,
     $ ADUDX2,SDUDX2,ADUDY2,SDUDY2,ADUDZ2,SDUDZ2,ADVDX2,SDVDX2,ADVDY2,
     $ SDVDY2,ADVDZ2,SDVDZ2,ADWDX2,SDWDX2,ADWDY2,SDWDY2,ADWDZ2,SDWDZ2,
     $ AUFWFG,SUFWFG,AUFWFS,SUFWFS,AUFWFM,SUFWFM,AVFWFG,SVFWFG,AVFWFS,
     $ SVFWFS,AVFWFM,SVFWFM,AUFVFG,SUFVFG,AUFVFS,SUFVFS,AUFVFM,
     $ SUFVFM,OX,AOX,SOX,OXFG,AOXRG,SOXRG,OY,AOY,SOY,OYFG,
     $ AOYRG,SOYRG,OZ,AOZ,SOZ,OZFG,AOZRG,SOZRG,AO2,SO2,O2FG,AO2RG,
     $ SO2RG,AHE,SHE,HEFG,AHERG,SHERG
     $ ,BP,BU,BV,BW
     $ ,HILF3D1,HILF3D2,HILF3D3
     $ ,AUUM,SUUM,AVVM,SVVM,AWWM,SWWM,APPM,SPPM
     $ ,AUVM,SUVM,AUWM,SUWM,AVWM,SVWM
     $ ,AUXUXM,SUXUXM,AUYUYM,SUYUYM,AUZUZM,SUZUZM,AVXVXM,SVXVXM,AVYVYM
     $ ,SVYVYM,AVZVZM,SVZVZM,AWXWXM,SWXWXM,AWYWYM,SWYWYM,AWZWZM,SWZWZM
     $ )

            CALL DIBCA1  (  1  ,'BINAER  ',HILFL,KKML,ILIMX ,

     $ KK,JJ,II,KKA,JJA,IIA,KK,JJ,II,KKA,JJA,IIA,JJ1L,JJ2L,ILIMXP,
     $ ISLIDI,LINFB,KKXL,KKYL,KKZL,KSXL,KSYL,KSZL,KH0L,IDIMF,KLASS,FVT,
     $ HVK,ISELEP,ISELEA,ISLINP,ISLINA,KSTW,RIDENT,
     $   HILF( 1 ),   UFG( 1 ),  AURG( 1 ),   VFG( 1 ),  AVRG( 1 ),
     $    WFG( 1 ),  AWRG( 1 ),    OX( 1 ),   AOX( 1 ),  OXFG( 1 ),
     $  AOXRG( 1 ),    OY( 1 ),   AOY( 1 ),  OYFG( 1 ), AOYRG( 1 ),
     $     OZ( 1 ),   AOZ( 1 ),  OZFG( 1 ), AOZRG( 1 ), ARXUU(I2L),
     $  SRXUU(I1L), ARYUU(I2L), SRYUU(I1L), ARZUU(I2L), SRZUU(I1L),
     $  ASXUU(I1L), SSXUU(I1L), ASYUU(I1L), SSYUU(I1L), ASZUU(I1L),
     $  SSZUU(I1L), ARXVV(I2L), SRXVV(I1L), ARYVV(I2L), SRYVV(I1L),
     $  ARZVV(I2L), SRZVV(I1L), ASXVV(I1L), SSXVV(I1L), ASYVV(I1L),
     $  SSYVV(I1L), ASZVV(I1L), SSZVV(I1L), ARXWW(I2L), SRXWW(I1L),
     $  ARYWW(I2L), SRYWW(I1L), ARZWW(I2L), SRZWW(I1L), ASXWW(I1L),
     $  SSXWW(I1L), ASYWW(I1L), SSYWW(I1L), ASZWW(I1L), SSZWW(I1L),
     $  ARXUV(I2L), SRXUV(I1L), ARYUV(I2L), SRYUV(I1L), ARXUW(I2L),
     $  SRXUW(I1L), ARYUW(I2L), SRYUW(I1L), ARYVW(I2L), SRYVW(I1L),
     $  ACXUW(I2L), SCXUW(I1L), ACZUW(I2L), SCZUW(I1L),ARXOXX(I2L),
     $ SRXOXX(I1L),ARXOXY(I2L),SRXOXY(I1L),ARXOXZ(I2L),SRXOXZ(I1L),
     $ ARXOYY(I2L),SRXOYY(I1L),ARXOYZ(I2L),SRXOYZ(I1L),ARXOZZ(I2L),
     $ SRXOZZ(I1L),ARYOXX(I2L),SRYOXX(I1L),ARYOXY(I2L),SRYOXY(I1L),
     $ ARYOXZ(I2L),SRYOXZ(I1L),ARYOYY(I2L),SRYOYY(I1L),ARYOYZ(I2L),
     $ SRYOYZ(I1L),ARYOZZ(I2L),SRYOZZ(I1L),ARZOXX(I2L),SRZOXX(I1L),
     $ ARZOXY(I2L),SRZOXY(I1L),ARZOXZ(I2L),SRZOXZ(I1L),ARZOYY(I2L),
     $ SRZOYY(I1L),ARZOYZ(I2L),SRZOYZ(I1L),ARZOZZ(I2L),SRZOZZ(I1L),
     $ ASXOXX(I1L),SSXOXX(I1L),ASXOYY(I1L),SSXOYY(I1L),ASXOZZ(I1L),
     $ SSXOZZ(I1L),ASYOXX(I1L),SSYOXX(I1L),ASYOYY(I1L),SSYOYY(I1L),
     $ ASYOZZ(I1L),SSYOZZ(I1L),AHOZOY(I2L),SHOZOY(I2L),AHOZOX(I2L),
     $ SHOZOX(I2L),AHOYOX(I2L),SHOYOX(I2L))
 
            END IF
            IF(LDIC  ) THEN
            CALL DIBCA   (  3  ,'CODIERT ',IIDENT,
     $                    ILIMXA,IDIBD1,IDIBD2,IDIBD3,IDIBD4,RKOMXA,
     $                    RDIBD1,RDIBD2,RDIBD3,RDIBD4,IGRID,IDIM3D,

     $ KK,JJ,II,KKA,JJA,IIA,KK,JJ,II,KKA,JJA,IIA,NBND,
C-------------------------- AB HIER, WIE cstaca.h:
     $ ISELEA,ISELEP,
     $ ISAMPA,ISAMPP,ISLINP,ISLINA,ISLIDI,
     $ RIDENT,HILF,B,U,UO,AU,SU,UFG,AURG,SURG,AUSKG,SUSKG,AUFLG,SUFLG,V,
     $ VO,AV,SV,VFG,AVRG,SVRG,AVSKG,SVSKG,AVFLG,SVFLG,W,WO,AW,SW,WFG,
     $ AWRG,SWRG,AWSKG,SWSKG,AWFLG,SWFLG,P,AP,SP,PFG,APRG,SPRG,
     $ APSKG,SPSKG,APFLG,SPFLG,
     $ G,AEFG,SEFG,AEFS,SEFS,ADFG,SDFG,
     $ ADUDX2,SDUDX2,ADUDY2,SDUDY2,ADUDZ2,SDUDZ2,ADVDX2,SDVDX2,ADVDY2,
     $ SDVDY2,ADVDZ2,SDVDZ2,ADWDX2,SDWDX2,ADWDY2,SDWDY2,ADWDZ2,SDWDZ2,
     $ AUFWFG,SUFWFG,AUFWFS,SUFWFS,AUFWFM,SUFWFM,AVFWFG,SVFWFG,AVFWFS,
     $ SVFWFS,AVFWFM,SVFWFM,AUFVFG,SUFVFG,AUFVFS,SUFVFS,AUFVFM,
     $ SUFVFM,OX,AOX,SOX,OXFG,AOXRG,SOXRG,OY,AOY,SOY,OYFG,
     $ AOYRG,SOYRG,OZ,AOZ,SOZ,OZFG,AOZRG,SOZRG,AO2,SO2,O2FG,AO2RG,
     $ SO2RG,AHE,SHE,HEFG,AHERG,SHERG
     $ ,BP,BU,BV,BW
     $ ,HILF3D1,HILF3D2,HILF3D3
     $ ,AUUM,SUUM,AVVM,SVVM,AWWM,SWWM,APPM,SPPM
     $ ,AUVM,SUVM,AUWM,SUWM,AVWM,SVWM
     $ ,AUXUXM,SUXUXM,AUYUYM,SUYUYM,AUZUZM,SUZUZM,AVXVXM,SVXVXM,AVYVYM
     $ ,SVYVYM,AVZVZM,SVZVZM,AWXWXM,SWXWXM,AWYWYM,SWYWYM,AWZWZM,SWZWZM
     $ )

            CALL DIBCA1  (  3  ,'CODIERT ',HILFL,KKML,ILIMX ,

     $ KK,JJ,II,KKA,JJA,IIA,KK,JJ,II,KKA,JJA,IIA,JJ1L,JJ2L,ILIMXP,
     $ ISLIDI,LINFB,KKXL,KKYL,KKZL,KSXL,KSYL,KSZL,KH0L,IDIMF,KLASS,FVT,
     $ HVK,ISELEP,ISELEA,ISLINP,ISLINA,KSTW,RIDENT,
     $   HILF( 1 ),   UFG( 1 ),  AURG( 1 ),   VFG( 1 ),  AVRG( 1 ),
     $    WFG( 1 ),  AWRG( 1 ),    OX( 1 ),   AOX( 1 ),  OXFG( 1 ),
     $  AOXRG( 1 ),    OY( 1 ),   AOY( 1 ),  OYFG( 1 ), AOYRG( 1 ),
     $     OZ( 1 ),   AOZ( 1 ),  OZFG( 1 ), AOZRG( 1 ), ARXUU(I2L),
     $  SRXUU(I1L), ARYUU(I2L), SRYUU(I1L), ARZUU(I2L), SRZUU(I1L),
     $  ASXUU(I1L), SSXUU(I1L), ASYUU(I1L), SSYUU(I1L), ASZUU(I1L),
     $  SSZUU(I1L), ARXVV(I2L), SRXVV(I1L), ARYVV(I2L), SRYVV(I1L),
     $  ARZVV(I2L), SRZVV(I1L), ASXVV(I1L), SSXVV(I1L), ASYVV(I1L),
     $  SSYVV(I1L), ASZVV(I1L), SSZVV(I1L), ARXWW(I2L), SRXWW(I1L),
     $  ARYWW(I2L), SRYWW(I1L), ARZWW(I2L), SRZWW(I1L), ASXWW(I1L),
     $  SSXWW(I1L), ASYWW(I1L), SSYWW(I1L), ASZWW(I1L), SSZWW(I1L),
     $  ARXUV(I2L), SRXUV(I1L), ARYUV(I2L), SRYUV(I1L), ARXUW(I2L),
     $  SRXUW(I1L), ARYUW(I2L), SRYUW(I1L), ARYVW(I2L), SRYVW(I1L),
     $  ACXUW(I2L), SCXUW(I1L), ACZUW(I2L), SCZUW(I1L),ARXOXX(I2L),
     $ SRXOXX(I1L),ARXOXY(I2L),SRXOXY(I1L),ARXOXZ(I2L),SRXOXZ(I1L),
     $ ARXOYY(I2L),SRXOYY(I1L),ARXOYZ(I2L),SRXOYZ(I1L),ARXOZZ(I2L),
     $ SRXOZZ(I1L),ARYOXX(I2L),SRYOXX(I1L),ARYOXY(I2L),SRYOXY(I1L),
     $ ARYOXZ(I2L),SRYOXZ(I1L),ARYOYY(I2L),SRYOYY(I1L),ARYOYZ(I2L),
     $ SRYOYZ(I1L),ARYOZZ(I2L),SRYOZZ(I1L),ARZOXX(I2L),SRZOXX(I1L),
     $ ARZOXY(I2L),SRZOXY(I1L),ARZOXZ(I2L),SRZOXZ(I1L),ARZOYY(I2L),
     $ SRZOYY(I1L),ARZOYZ(I2L),SRZOYZ(I1L),ARZOZZ(I2L),SRZOZZ(I1L),
     $ ASXOXX(I1L),SSXOXX(I1L),ASXOYY(I1L),SSXOYY(I1L),ASXOZZ(I1L),
     $ SSXOZZ(I1L),ASYOXX(I1L),SSYOXX(I1L),ASYOYY(I1L),SSYOYY(I1L),
     $ ASYOZZ(I1L),SSYOZZ(I1L),AHOZOY(I2L),SHOZOY(I2L),AHOZOX(I2L),
     $ SHOZOX(I2L),AHOYOX(I2L),SHOYOX(I2L))
            END IF

            END IF        
            END IF 
 

         ENDDO
C                                 GESAMTANZAHL DER ITERATIONEN
C                                 PHYS. GESAMTZEIT 
         ITTOT  = IIDENT(30)
         TIMEPH = RIDENT(30)        
            
CTBA1 070203: NOW INITIALISE SCALAR FIELD IF IIDENT(50) = 0 
C             THEN SET IIDENT(50) TO 1 
C             -> T-FIELD IS SET AND WILL BE WRITTEN TO RESULT 
C
C***********************************************************************
C123     CALL HPRORD (LPART,ITBGN,IVGL)
C***********************************************************************
C
      END IF
CTBC  END DREAD CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C---- ------------------------------------ EINLESEN DER PARTIKEL 
C---- ------------------------------------ EINLESEN DER PARTIKEL FERTIG
C
C                            BELEGUNG VON WEITEREN GITTERN
C                            FALLS LGRIDNEW = .TRUE.
C


      DO IGRID=NGRDOLD+1,NGRDDFD
        CALL INIGRIDGEO
     $             (CIDENT,IIDENT,RIDENT,
     $              DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,
     $              U,V,W,P,G,B,BP,BU,BV,BW,
     $              ITSTEP,ITTOT,TIMEPH,
     $              IDIMF,NGRP1, NGRP3, IGRHF, GRX0, GRX1, GRHF,
     $              IDIM3D,IDIM2D,IDIM1D,NBUF,NBND,IDIMA,IDIM1L,IDIM2L,
     $              IGRID,CONV1S)

      ENDDO              

C     jk 28.7.2005
C     Gitternachbarn stezten
      WRITE(6,*)'VOR BLOCKBP!'
      CALL MGNBRSET (DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,HILF,
     $     IDIM1D,IDIM2D,NBND)
C     BP Feld erzeugen
      CALL BLOCKBP (
     $     DDX,DDY,DDZ,
     $     DX,DY,DZ,
     $     X,Y,Z,
     $     BP, HILF,
     $     IDIM3D,IDIM2D,IDIM1D)

      DO IGRID=NGRDOLD+1,NGRDDFD
        CALL INIGRIDVAR
     $             (CIDENT,IIDENT,RIDENT,
     $              DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,
     $              U,V,W,P,G,B,BP,BU,BV,BW,
     $              ITSTEP,ITTOT,TIMEPH,
     $              IDIMF,NGRP1, NGRP3, IGRHF, GRX0, GRX1, GRHF,
     $              IDIM3D,IDIM2D,IDIM1D,NBUF,NBND,IDIMA,IDIM1L,IDIM2L,
     $              IGRID,CONV1S)

      ENDDO              

C                           GITTER FUER GEBIETSZERLEGUNG
C
      DO ILEVEL = MINLEVEL,MAXLEVEL
        DO I = 1,NOFSLCED(ILEVEL)
           IGRID = IGRDOFSLCED(I,ILEVEL)
C         CALL SLICEGRD (U,V,W,P,G,B,HILF,IDIM3D,IGRID)
        ENDDO
      ENDDO
C
C                       JETZT SIND BERECHNINGSGITTER ZERSCHNITTEN
C                       GEBIETSZERLEGUNG ABGESCHLOSSEN
C
C
C                            GROBGITTER FUER MULTIGRID-DRUCKKORREKTUR
C

      DO ILEVEL = MINLEVEL,MAXLEVEL
      DO I = 1,NOFVPIT(ILEVEL)
         IGRID = IGRDOFVPIT( I , ILEVEL )
           
         CALL MGVPSET
     $             (CIDENT,IIDENT,RIDENT,
     $              DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,
     $              U,V,W,P,G,B,BU,BV,BW,BP,	
     $              ITSTEP,ITTOT,TIMEPH,
     $              IDIMF,NGRP1, NGRP3, IGRHF, GRX0, GRX1, GRHF,
     $              IDIM3D,IDIM2D,IDIM1D,NBUF,NBND,IDIMA,IDIM1L,IDIM2L,
     $              IDIMF,IGRID,CONV1S
     $             )     
         ENDDO
      ENDDO

C                       SETZEN DER INDIZES FUER DIE BLOWING/SUCTION-
C                       RANDBEDINGUNG
CTBC1                   SET INDICES FOR SCALAR BOUNDARY CONDITIONS
      CALL MGBLOINDSET (DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,HILF,
     $               IDIM1D,IDIMF,NBND,UBO,VBO,WBO)

C                       SETZEN DER PARENTGITTER, DIE NOCH NICHT
C                       UEBER STEUERFILE SPEZIFIZIERT WURDEN
C
      CALL MGPARSET (DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,HILF,
     $               IDIM1D,IDIMF,NBND)

C                       SETZEN DER NACHBARGITTER, DIE NOCH NICHT
C                       UEBER STEUERFILE SPEZIFIZIERT WURDEN
C
      CALL MGNBRSET (DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,HILF,
     $               IDIM1D,IDIM2D,NBND)
C
C                       INFORMATION UEBER VERNETZUNG DER GITTER
C
      CALL MGCONINF (DX,DY,DZ,X,Y,Z,IDIM1D,NBND)      
C
C                       EINIGE GITTERSPEZIFISCHE INITIALISIERUNGEN
C
      DO ILEVEL = MINLEVEL,MAXLEVEL
        DO I = 1,NOFLEVEL(ILEVEL)
           IGRID = IGRDOFLEVEL(I,ILEVEL)
C
C                                NUR, FALLS GITTER NICHT GEBIETSZERLEGT
C
         IF ( .NOT. LSLICE(IGRID) ) THEN
         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $           NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
         CALL MGBASBSCA (NFROSCA,NBACSCA,NRGTSCA,NLFTSCA,NBOTSCA,
     $        NTOPSCA,NCUBSCA,IGRID)
C         CALL WRITE3DDIAGY(KK,JJ,II,W(IP3),31,3)
C         CALL WRITE3DDIAGY(KK,JJ,II,U(IP3),41,4)
C
C                                 FESTSTELLEN IN WELCHEN KOORD.-RI.
C                                 DAS GITTER AEQUIDISTANT IST
C
         CALL AEQGRD (KK,JJ,II,KK,JJ,II,DX(IP1),DY(IP1),DZ(IP1),NBND,
     &                NXGRAE(IGRID),NYGRAE(IGRID),NZGRAE(IGRID),IGRID)
C
C                  BERECHNEN DER REZIPROKWERTE DER GEOMETRIEFAKTOREN
C
C
         CALL REZIPD(KK,JJ,II,KK,JJ,II,DX(IP1),DY(IP1),DZ(IP1),
     $               DDX(IP1),DDY(IP1),DDZ(IP1),
     $               RDX(IP1),RDY(IP1),RDZ(IP1),
     $               RDDX(IP1),RDDY(IP1),RDDZ(IP1))     

C
C                  BERECHNEN DER FAKTOREN FUER DIREKTEN POISSONLOESER
C
         IF (LPOISSONDIR(IGRID)) THEN
            CALL CALPSFAK 
     $                   (KK,JJ,II,NBND,
     $                     DX(IP1), DY(IP1), DZ(IP1),
     $                    DDX(IP1),DDY(IP1),DDZ(IP1),
     $                    FPSFAK(IP2),FCOSMY(IP1),FCOSNY(IP1),
     &                NXGRAE(IGRID),NYGRAE(IGRID),NZGRAE(IGRID),
     $                     NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,IGRID)
C
C                 SINUS- UND COSINUS-TABELLEN FUER FOURIERTRANSFORMATION
C                 IN X- UND Y-RICHTUNG
C

            CALL FFTSET 
     $     ( II-2*NBND , IFFTX(1,IGRID) , IPERMUX(IP1) , RFFTX(IP1*2))
            CALL FFTSET 
     $     ( JJ-2*NBND , IFFTY(1,IGRID) , IPERMUY(IP1) , RFFTY(IP1*2))


 

         ENDIF
C
C
C************************** BERECHNUNG DER KOEFFIZIENTEN DES
C*********************** KOMPAKTVERFAHRENS FUER JEDES GITTER
C


      CIP1 = 12*3*(IP1-1)+1



C********************************************************************
C   Berechnung der Ubertragunsfunktion des Filters fuer ADM-Scalar
C*********************************************************************
C

C************************************SET POISSON FAKTORS**************
      CALL GITEIG     (KK,JJ,II, DX(IP1),DY(IP1),DZ(IP1),
     $                 GSAW(IP1),GSAE(IP1),GSAN(IP1),
     $                 GSAS(IP1),GSAT(IP1),GSAB(IP1),
     $                 BP(IP3),
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,
     $                 GEOVP(IP3)
     $                 )


      CALL    SIPLU    (KK,JJ,II,
     $                 GSAW(IP1),GSAE(IP1),GSAN(IP1),
     $                 GSAS(IP1),GSAT(IP1),GSAB(IP1),
     $                 BP(IP3),
     $                 GEOVP(IP3),SIPLW(IP3),SIPLS(IP3),
     $                 SIPLB(IP3),SIPLPR(IP3),SIPUE(IP3),
     $                 SIPUN(IP3),SIPUT(IP3))
CTBC     ENDIF (.NOT. LSLICE)
         ENDIF
        ENDDO
      ENDDO
C
C
C
C

      ITSTEP = 0

      IF(.NOT. DCONT) THEN
         ITTOT  = 0
         TIMEPH = 0.0
      END IF
C
C
C
C                       EINIGE GITTERSPEZIFISCHE INITIALISIERUNGEN
C
C-101296      DO IGRID = 1,NGRDSET
      DO ILEVEL = MAXLEVEL,MINLEVEL,-1
        DO I = 1,NOFTST(ILEVEL)
           IGRID = IGRDOFTST(I,ILEVEL)
         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $           NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)


C
C
C                             BELEGEN DES FRONT-BUFFERS
C
C
         IF ( NFRO .EQ. 2 .OR. NFRO .EQ. 11) THEN

           IF (IVPCHILD(IGRID).EQ.0) THEN
C                                            REGULAERES GITTER

             IF(LDIEIB) THEN
C
C                                 EINLESEN DES G-PROFILS UND DER GESCHW.
C                                 PROFILE AM EINTRITTSRAND (UNFOR-
C                                 MATIERT AUF KANAL 11)
C
CTBC11 260203 READING OF T PROFILES NOT YET IMPLEMENTED CCCCCCCCCCCCCCCC
                CALL DEIBI   (KMX(1),JMX(1),IMX(1),KMX(1),JMX(1),IMX(1),
     $                        ZTOT(1),YTOT(1),XTOT(1),
     $                        CIDEND,IIDEND,RIDEND,
C#ifdef _TSCAL_
C     $            UFR(IBB),VFR(IBB),WFR(IBB),TFR(IBB),PFR(IBB),GFR(IBB),
C     $               UI1,VI1,WI1,TI1,UI2,VI2,WI2,TI2,GI1,GI2,
C#else
     $                     UFR(IBB),VFR(IBB),WFR(IBB),PFR(IBB),GFR(IBB),
     $                        UI1,VI1,WI1,UI2,VI2,WI2,GI1,GI2,
C#endif 
     $                        UGRID,ITSTEP,ITTOT,TIMEPH,DT)
C
             END IF
             IF(LDIEIC) THEN
C
C                                 EINLESEN DES G-PROFILS UND DER GESCHW.
C                                 PROFILE AM EINTRITTSRAND (FORMATIERT
C                                 AUF KANAL 13)
C
CTBC11 260203 READING OF T PROFILES NOT YET IMPLEMENTED CCCCCCCCCCCCCCCC
                CALL DEICI   (KMX(1),JMX(1),IMX(1),KMX(1),JMX(1),IMX(1),
     $                        ZTOT(1),YTOT(1),XTOT(1),
     $                        CIDEND,IIDEND,RIDEND,
C#ifdef _TSCAL_
C     $            UFR(IBB),VFR(IBB),WFR(IBB),TFR(IBB),PFR(IBB),GFR(IBB),
C     $               UI1,VI1,WI1,TI1,UI2,VI2,WI2,TI2,GI1,GI2,
C#else
     $                     UFR(IBB),VFR(IBB),WFR(IBB),PFR(IBB),GFR(IBB),
     $                        UI1,VI1,WI1,UI2,VI2,WI2,GI1,GI2,
C#endif 
     $                        UGRID,ITSTEP,ITTOT,TIMEPH,DT)
C  
        END IF
      
             IF((.NOT.LDIEIB).AND.(.NOT.LDIEIC)) THEN
               CALL SVEIPR (KK,JJ, 2,KK,JJ, 2,
     $                     UFR(IBB),VFR(IBB),WFR(IBB),PFR(IBB),
     $                     GFR(IBB),Y(IP1),Z(IP1),DZ(IP1),VCON,WCON,
     $                     YBANF(1,1,IGRID),YBEND(1,1,IGRID),
     $                     ZBANF(1,1,IGRID),ZBEND(1,1,IGRID),
     $                     FREQB(1,1,IGRID),ANIVEAU(1,1,IGRID),
     $                     AUB(1,1,IGRID),AVB(1,1,IGRID),AWB(1,1,IGRID),
     $                     TIMEPH,DT)
             ENDIF
           ENDIF
         ENDIF

C hier Randschicht noch ok
C       CALL WRITE3DDIAGY (KMX(1),JMX(1),IMX(1),U,99,3)
C       stop

         IF ( NFRO .EQ. 11 .OR. NFRO .EQ. 12) THEN
            CALL MGSVFLU (U,V,W,P,G,AU,AV,AW,
     $                    HILF,UFR,VFR,WFR,PFR,GFR,X,Y,Z,
     $                    IDIM3D,IDIM2D,IDIM1D,IDIMA,
     $                    NBND,IGRID,ITSTEP,1,NFRO
     $                   )
            
         ENDIF

C
C
C                             BELEGEN DES BOTTOM-BUFFERS
C
C
         IF ( NBOT .EQ. 2 .OR. NBOT .EQ. 11) THEN

            IF (IVPCHILD(IGRID).EQ.0) THEN
C                                            REGULAERES GITTER
               IF((.NOT.LDIEIB).AND.(.NOT.LDIEIC)) THEN
                  CALL SV3D (JJ,II, 2,JJ,II, 2,
     $                 WBO(IBB),UBO(IBB),VBO(IBB),PBO(IBB),
     $                 GBO(IBB),HILF(IBB),
     $                   Z(IP1),  X(IP1),  Y(IP1),
     $                  DZ(IP1), DX(IP1), DY(IP1),
     $                 DDZ(IP1),DDX(IP1),DDY(IP1),
     $                 VCON,WCON,
     $                 YBANF(1,5,IGRID),YBEND(1,5,IGRID),
     $                 XBANF(1,5,IGRID),XBEND(1,5,IGRID),
     $                 ZBANF(1,5,IGRID),ZBEND(1,5,IGRID),
     $                 FREQB(1,5,IGRID),ANIVEAU(1,5,IGRID),
     $                 AUB(1,5,IGRID),AVB(1,5,IGRID),AWB(1,5,IGRID),
     $                 TIMEPH,DT,0)
               ENDIF
            ENDIF
         ENDIF
         IF ( NBOT .EQ. 11 .OR. NBOT .EQ. 12 ) THEN
            CALL MGBOFLU (U,V,W,P,G,AU,AV,AW,
     $                    HILF,UBO,VBO,WBO,PBO,GBO,X,Y,Z,
     $                    IDIM3D,IDIM2D,IDIM1D,IDIMA,
     $                    NBND,IGRID,ITSTEP,1,NBOT)
            
         ENDIF
      ENDDO
      ENDDO
C                                  LEVELS FROM FINE TO COARSE, IN
C                                  ORDER TO PROPAGATE BOUNDARY-COND.
C                                  FROM FINE TO COARSE GRIDS
      DO ILEVEL = MAXLEVEL,MINLEVEL,-1
         DO I = 1,NOFVPIT(ILEVEL)
            IGRID = IGRDOFVPIT(I,ILEVEL)
            IF (IPARENT(IGRID) .NE. 0) THEN
            IPROCF = 0
            IPROCC = 0

               CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $              NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
C
               IF ( LCHILD(IGRID) ) THEN
C
C                                 FRONT BUFFER
C
                  IF (      NFRO .EQ. 2 
     $                 .OR. NFRO .EQ. 11 
     $                 .OR. NFRO .EQ. 12 ) THEN

C                                 BEREITSTELLUNG DER POINTER
C                                 DIMENSIONIERUNGEN DES FEINGITTERS

                     CALL MGDIMS  (KKC,JJC,IIC,IPARENT(IGRID))
                     CALL MGPOINT 
     &              (IP3C,IP2C,IP1C,IBBC,IB3C,IBUC,IPARENT(IGRID))
                
                
                     CALL MGBFTC (KK,JJ, 2 ,DX(IP1),DY(IP1),DZ(IP1),
     &                    DDX(IP1),DDY(IP1),DDZ(IP1),UFR(IBB),
     &                    KKC,JJC, 2 ,DX(IP1C),DY(IP1C),DZ(IP1C),
     &                    DDX(IP1C),DDY(IP1C),DDZ(IP1C),UFR(IBBC),
     &                    KPOSITION(IGRID),JPOSITION(IGRID),
     &                    IPOSITION(IGRID),'U',IPROCF,IPROCC)

                  END IF
C
C                                BOTTOM BUFFER
C
                  IF (       NBOT .EQ. 2 
     $                  .OR. NBOT .EQ. 11 
     $                  .OR. NBOT .EQ. 12 ) THEN     	
                     CALL MGDIMS  (KKC,JJC,IIC,IPARENT(IGRID))
                     CALL MGPOINT 
     &                 (IP3C,IP2C,IP1C,IBBC,IB3C,IBUC,IPARENT(IGRID))
                     CALL MGBFTC (JJ,II, 2 ,DZ(IP1),DX(IP1),DY(IP1),
     &                    DDZ(IP1),DDX(IP1),DDY(IP1),WBO(IBB),
     &                    JJC,IIC, 2 ,DZ(IP1C),DX(IP1C),DY(IP1C),
     &                    DDZ(IP1C),DDX(IP1C),DDY(IP1C),WBO(IBBC),
     &                    JPOSITION(IGRID),IPOSITION(IGRID),
     &                    KPOSITION(IGRID),'U',IPROCF,IPROCC)
                  END IF
               END IF
         ENDIF
         ENDDO
      ENDDO
C
C
C
C                                 SHIFTEN DES U-FELDES FUER DIE GALILEI-
C                                 TRANSFORMATION. DAS GITTER WIRD MIT
C                                 DER GESCHW.  UGRID  BEWEGT.
      DO ILEVEL = MINLEVEL,MAXLEVEL
         DO I = 1,NOFTST(ILEVEL)
            IGRID = IGRDOFTST(I,ILEVEL)

               CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $              NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
C
C                                 BELEGUNG DER AUSTRITTSRANDFELDER
C                                 UBA,VBA,WBA FUER DIE KONVEKTIVE
C                                 AUSTRITTSRANDBEDINGUNG
C
C
               CALL VSHIFT  (KK,JJ,II,KK,JJ,II,U(IP3),UGRID,'S')
C
CTBC6 041202 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL SETBACKOLD   (KK,JJ,II,KK,JJ,II,UBA(IBB), VBA(IBB), WBA(IBB),
     $                   U(IP3),V(IP3),W(IP3)
     $                  )
C
         ENDDO
      ENDDO

C
C                                 SETZEN DER RANDBEDINGUNGEN
C

      DO ILEVEL = MINLEVEL,MAXLEVEL
Cifdef _MPI_
C
C                                 BEI MESSAGE-PASSING, KOMMUNIKATION
C
         CALL CONNECTMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,UTO,VTO,WTO,PTO,GTO,
     $              VFR,WFR,PFR,GFR,
     $              'T',ILEVEL,0,0,0,
     $              NOFTST(ILEVEL),IGRDOFTST(1,ILEVEL))

C
C                               FILL THE BOUNDARY BUFFER FOR LOCAL GRIDS
C
         CALL BPARMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,DUDY,DUDZ,DVDX,
     $              UTO,VTO,WTO,PTO,GTO,UBA,VBA,WBA,PBA,GBA,
     $              VFR,WFR,PFR,GFR,
     $              'T',ILEVEL,0,0,0,
     $              NOFTST(ILEVEL),IGRDOFTST(1,ILEVEL)
     $             )      
Cendif
        DO I = 1,NOFTST(ILEVEL)
           IGRID = IGRDOFTST(I,ILEVEL)

         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $           NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

         CALL BOUNDMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $              'T',IGRID,0,0,0,
     $              TIMEPH,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $              PFR,GFR
     $             )     
        ENDDO
      ENDDO
CCCCCCCCCCCCCCC  	AM 18.3.1996 AUSGESCHALTET,
CCCCCCCCCCCCCCC         DA REYNOLDSZAHLWECHSEL BEI LAMINAREN
CCCCCCCCCCCCCCC         LAEUFEN SONST NICHT MOEGLICH!
CCCCCCCCCCCCCCC         IF(MTURB .EQ. 1) THEN

      DO ILEVEL = MINLEVEL,MAXLEVEL
Cifdef _MPI_
C
C                                 BEI MESSAGE-PASSING, KOMMUNIKATION
C
         CALL CONNECTMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,UTO,VTO,WTO,PTO,GTO,
     $              VFR,WFR,PFR,GFR,
     $              'T',ILEVEL,0,0,0,
     $              NOFTST(ILEVEL),IGRDOFTST(1,ILEVEL))

Cendif
        DO I = 1,NOFTST(ILEVEL)
           IGRID = IGRDOFTST(I,ILEVEL)

         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $           NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
           CALL MGPOINT (IP3,IP2,IP1,IBB,IB3,IBU,IGRID)
C
C
C                                 SETZEN DER RANDBEDINGUNGEN
C
         CALL BOUNDMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $              'T',IGRID,0,0,0,
     $              TIMEPH,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $              PFR,GFR
     $             )     
C                                 CALCULATION OF THE LES CONSTANT CONV1S
C
      CALL LESCONST (KK,JJ,II,KK,JJ,II,
     $                              X(IP1),Y(IP1),Z(IP1),
     $                              DX(IP1),DY(IP1),DZ(IP1),
     $                              DDX(IP1),DDY(IP1),DDZ(IP1),
     $                              CONV1S(IP1),CONV1SANF(IGRID),
     $                              CONV1SEND(IGRID),
     $                              TRANSLES1(IGRID),TRANSLES2(IGRID))

C                                 BELEGUNG DES G-FELDES UND SETZEN
C                                 DER RANDBEDINGUNGEN FUER G.
C
         CALL EFVISC (KK,JJ,II,KK,JJ,II,
     $                              X(IP1),Y(IP1),Z(IP1),
     $                              DX(IP1),DY(IP1),DZ(IP1),
     $                              DDX(IP1),DDY(IP1),DDZ(IP1),
     $                              RDX(IP1),RDY(IP1),RDZ(IP1),
     $                              RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                              U(IP3),V(IP3),W(IP3),G(IP3),BP(IP3),
     $                              DUDY(IP2),DUDZ(IP2),DVDX(IP2),
     $                              DVDZ(IP2),DWDX(IP2),DWDY(IP2),
     $                              GMOL,RHO,UGRID,
     $                              IC1(IGRID),IC2(IGRID),
     $                              JC1(IGRID),JC2(IGRID),
     $                              KC1(IGRID),KC2(IGRID),
     $                              NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,
     $                              CONV1S(IP1),
     $   HILF0(1),HILF1(1),HILF2(1),HILF3(1),HILF4(1),
     $   HILF5(1),HILF6(1),
     $   LUX(1),LUY(1),LUZ(1),LVY(1),LVZ(1),LWZ(1),MUX(1),
     $   MUY(1),MUZ(1),MVY(1),MVZ(1),MWZ(1),
     $   FA11(1),FA12(1),FA13(1),FA21(1),FA22(1),FA23(1),
     $   FA31(1),FA32(1),FA33(1),
     $   DUDX3(1),DUDY3(1),DUDZ3(1),DVDX3(1),DVDY3(1),
     $   DVDZ3(1),
     $   DWDX3(1),DWDY3(1),DWDZ3(1),NORMS(1),SDUDX3(1),
     $   SDUDY3(1),SDUDZ3(1),
     $   SDVDY3(1),SDVDZ3(1),SDWDZ3(1),
     $   UC(1),VC(1),WC(1),UUC(1),UVC(1),UWC(1),VVC(1),
     $   VWC(1),WWC(1),CDELTA(1))


         ENDDO
       ENDDO

      DO ILEVEL = MINLEVEL,MAXLEVEL
Cifdef _MPI_
C
C                                 BEI MESSAGE-PASSING, KOMMUNIKATION
C
         CALL CONNECTMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,UTO,VTO,WTO,PTO,GTO,
     $              VFR,WFR,PFR,GFR,
     $              'T',ILEVEL,0,0,0,
     $              NOFTST(ILEVEL),IGRDOFTST(1,ILEVEL))

Cendif
        DO I = 1,NOFTST(ILEVEL)
           IGRID = IGRDOFTST(I,ILEVEL)

           CALL BOUNDMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $              'T',IGRID,0,0,0,
     $              TIMEPH,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $              PFR,GFR
     $             )       
         ENDDO
       ENDDO

CCCCCCCCCCCCCCC  	AM 18.3.1996 AUSGESCHALTET      ENDIF
C
C
      DO ILEVEL = MINLEVEL,MAXLEVEL
Cifdef _MPI_
C
C                                 BEI MESSAGE-PASSING, KOMMUNIKATION
C
         CALL CONNECTMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,UTO,VTO,WTO,PTO,GTO,
     $              VFR,WFR,PFR,GFR,
     $              'T',ILEVEL,0,0,0,
     $              NOFTST(ILEVEL),IGRDOFTST(1,ILEVEL))

Cendif
        DO I = 1,NOFTST(ILEVEL)
           IGRID = IGRDOFTST(I,ILEVEL)

         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $           NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
C
         CALL BOUNDMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $              'T',IGRID,0,0,0,
     $              TIMEPH,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $              PFR,GFR
     $             )       

C                                 BELEGUNG DER UO, VO, WO -FELDER MIT
C                                 DEN WERTEN DER U, V, W -FELDER
c        write(0,*) 'before boundmg    myid = ', myid
c        stop
         CALL BOUNDMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $              'Z',IGRID,0,0,0,
     $              TIMEPH,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $              PFR,GFR
     $             )       
c        write(0,*) 'after first boundmg with Z  myid=', myid
c        stop
C
        CALL EXCHAN  (KK,JJ,II,KK,JJ,II,
     $                U(IP3),V(IP3),W(IP3),UO(IP3),VO(IP3),WO(IP3),2)
CTBA2 19.11.02 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        ENDDO
      ENDDO
C
C
      IF (MTURB .EQ. 1) THEN
C                                 VORBELEGUNG DER STATISTIK

      DO ILEVEL = MINLEVEL,MAXLEVEL
C----------------------------------------------------------------------
       IF (LSCAI(ILEVEL)) THEN 
          NGRIDPERLEVEL = NOFSCAI(ILEVEL)
        DO I=1,NGRIDPERLEVEL
          NOFTHISGRID(I) = IGRDOFSCAI(I,ILEVEL)
        ENDDO
       ELSE
          NGRIDPERLEVEL = NOFTST(ILEVEL)
        DO I=1,NGRIDPERLEVEL
          NOFTHISGRID(I) = IGRDOFTST(I,ILEVEL)
        ENDDO
       ENDIF
C----------------------------------------------------------------------
        DO I = 1,NGRIDPERLEVEL
           IGRID = NOFTHISGRID(I)
C
C                                NUR, FALLS GITTER NICHT GEBIETSZERLEGT
C
          IF ( .NOT. LSLICE(IGRID) ) THEN
             
             CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $            NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
             CALL MGDIMA  (KKA,JJA,IIA,IGRID)
             CALL MGPOINA (IAV,I2L,I1L,IGRID)
C
C                                 STATISTIK WIRD AUF REGULAEREN GITTERN 
C                                 GEMACHT
C
C
             CALL SETSTA  (KPP,JPP,IPP,IC1(IGRID),IC2(IGRID),JC1(IGRID),
     $                    JC2(IGRID),KC1(IGRID),KC2(IGRID),
     $                          NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,
     $                        XHOMOG(IGRID),YHOMOG(IGRID),ZHOMOG(IGRID),
     $                 GMOL,RHO,UGRID,DDX(IP1),DDY(IP1),DDZ(IP1),
     $                 DX(IP1),DY(IP1),DZ(IP1),X(IP1),Y(IP1),Z(IP1),
     $                 RDX(IP1),RDY(IP1),RDZ(IP1),
     $                 RDDX(IP1),RDDY(IP1),RDDZ(IP1),

     $ KK,JJ,II,KKA,JJA,IIA,KK,JJ,II,KKA,JJA,IIA,NBND,ISELEA,ISELEP,
     $ ISAMPA(1,IGRID),ISAMPP(1,IGRID),ISLINP,ISLINA,ISLIDI,RIDENT,
     $   HILF(IP3),     B(IP3),     U(IP3),    UO(IP3),    AU(IAV),
     $     SU(IAV),   UFG(IP3),  AURG(IAV),  SURG(IAV), AUSKG(IAV),
     $  SUSKG(IAV), AUFLG(IAV), SUFLG(IAV),     V(IP3),    VO(IP3),
     $     AV(IAV),    SV(IAV),   VFG(IP3),  AVRG(IAV),  SVRG(IAV),
     $  AVSKG(IAV), SVSKG(IAV), AVFLG(IAV), SVFLG(IAV),     W(IP3),
     $     WO(IP3),    AW(IAV),    SW(IAV),   WFG(IP3),  AWRG(IAV),
     $   SWRG(IAV), AWSKG(IAV), SWSKG(IAV), AWFLG(IAV), SWFLG(IAV),
     $      P(IP3),    AP(IAV),    SP(IAV),   PFG(IP3),  APRG(IAV),
     $   SPRG(IAV), APSKG(IAV), SPSKG(IAV), APFLG(IAV), SPFLG(IAV),
     $      G(IP3),  AEFG(IAV),  SEFG(IAV),  AEFS(IAV),  SEFS(IAV),      
     $   ADFG(IAV),  SDFG(IAV),        
     $ ADUDX2(IAV),SDUDX2(IAV),ADUDY2(IAV),
     $ SDUDY2(IAV),ADUDZ2(IAV),SDUDZ2(IAV),ADVDX2(IAV),SDVDX2(IAV),
     $ ADVDY2(IAV),SDVDY2(IAV),ADVDZ2(IAV),SDVDZ2(IAV),ADWDX2(IAV),
     $ SDWDX2(IAV),ADWDY2(IAV),SDWDY2(IAV),ADWDZ2(IAV),SDWDZ2(IAV),
     $ AUFWFG(IAV),SUFWFG(IAV),AUFWFS(IAV),SUFWFS(IAV),AUFWFM(IAV),
     $ SUFWFM(IAV),AVFWFG(IAV),SVFWFG(IAV),AVFWFS(IAV),SVFWFS(IAV),
     $ AVFWFM(IAV),SVFWFM(IAV),AUFVFG(IAV),SUFVFG(IAV),AUFVFS(IAV),
     $ SUFVFS(IAV),AUFVFM(IAV),SUFVFM(IAV),    OX(IP3),   AOX(IAV),
     $    SOX(IAV),  OXFG(IP3), AOXRG(IAV), SOXRG(IAV),    OY(IP3),
     $    AOY(IAV),   SOY(IAV),  OYFG(IP3), AOYRG(IAV), SOYRG(IAV),
     $     OZ(IP3),   AOZ(IAV),   SOZ(IAV),  OZFG(IP3), AOZRG(IAV),
     $  SOZRG(IAV),   AO2(IAV),   SO2(IAV),  O2FG(IP3), AO2RG(IAV),
     $  SO2RG(IAV),   AHE(IAV),   SHE(IAV),  HEFG(IP3), AHERG(IAV),
     $  SHERG(IAV)
     $   ,BP(IP3),     BU(IP3),     BV(IP3),     BW(IP3)
     $ ,HILF3D1(IP3),HILF3D2(IP3),HILF3D3(IP3)
     $ ,AUUM(IAV),SUUM(IAV),AVVM(IAV),SVVM(IAV),AWWM(IAV),SWWM(IAV)
     $ ,APPM(IAV),SPPM(IAV)
     $ ,AUVM(IAV),SUVM(IAV),AUWM(IAV),SUWM(IAV),AVWM(IAV),SVWM(IAV)
     $ ,AUXUXM(IAV),SUXUXM(IAV),AUYUYM(IAV),SUYUYM(IAV)
     $ ,AUZUZM(IAV),SUZUZM(IAV),AVXVXM(IAV),SVXVXM(IAV)
     $ ,AVYVYM(IAV),SVYVYM(IAV),AVZVZM(IAV),SVZVZM(IAV)
     $ ,AWXWXM(IAV),SWXWXM(IAV),AWYWYM(IAV),SWYWYM(IAV)
     $ ,AWZWZM(IAV),SWZWZM(IAV)

     $ )
 
            CALL SETST1  (HILFL,KKML,RKOMXP,RKOMXA,

     $ KK,JJ,II,KKA,JJA,IIA,KK,JJ,II,KKA,JJA,IIA,JJ1L,JJ2L,ILIMXP,
     $ ISLIDI,LINFB,KKXL,KKYL,KKZL,KSXL,KSYL,KSZL,KH0L,IDIMF,KLASS,FVT,
     $ HVK,ISELEP,ISELEA,ISLINP,ISLINA,KSTW,RIDENT,
     $   HILF( 1 ),   UFG( 1 ),  AURG( 1 ),   VFG( 1 ),  AVRG( 1 ),
     $    WFG( 1 ),  AWRG( 1 ),    OX( 1 ),   AOX( 1 ),  OXFG( 1 ),
     $  AOXRG( 1 ),    OY( 1 ),   AOY( 1 ),  OYFG( 1 ), AOYRG( 1 ),
     $     OZ( 1 ),   AOZ( 1 ),  OZFG( 1 ), AOZRG( 1 ), ARXUU(I2L),
     $  SRXUU(I1L), ARYUU(I2L), SRYUU(I1L), ARZUU(I2L), SRZUU(I1L),
     $  ASXUU(I1L), SSXUU(I1L), ASYUU(I1L), SSYUU(I1L), ASZUU(I1L),
     $  SSZUU(I1L), ARXVV(I2L), SRXVV(I1L), ARYVV(I2L), SRYVV(I1L),
     $  ARZVV(I2L), SRZVV(I1L), ASXVV(I1L), SSXVV(I1L), ASYVV(I1L),
     $  SSYVV(I1L), ASZVV(I1L), SSZVV(I1L), ARXWW(I2L), SRXWW(I1L),
     $  ARYWW(I2L), SRYWW(I1L), ARZWW(I2L), SRZWW(I1L), ASXWW(I1L),
     $  SSXWW(I1L), ASYWW(I1L), SSYWW(I1L), ASZWW(I1L), SSZWW(I1L),
     $  ARXUV(I2L), SRXUV(I1L), ARYUV(I2L), SRYUV(I1L), ARXUW(I2L),
     $  SRXUW(I1L), ARYUW(I2L), SRYUW(I1L), ARYVW(I2L), SRYVW(I1L),
     $  ACXUW(I2L), SCXUW(I1L), ACZUW(I2L), SCZUW(I1L),ARXOXX(I2L),
     $ SRXOXX(I1L),ARXOXY(I2L),SRXOXY(I1L),ARXOXZ(I2L),SRXOXZ(I1L),
     $ ARXOYY(I2L),SRXOYY(I1L),ARXOYZ(I2L),SRXOYZ(I1L),ARXOZZ(I2L),
     $ SRXOZZ(I1L),ARYOXX(I2L),SRYOXX(I1L),ARYOXY(I2L),SRYOXY(I1L),
     $ ARYOXZ(I2L),SRYOXZ(I1L),ARYOYY(I2L),SRYOYY(I1L),ARYOYZ(I2L),
     $ SRYOYZ(I1L),ARYOZZ(I2L),SRYOZZ(I1L),ARZOXX(I2L),SRZOXX(I1L),
     $ ARZOXY(I2L),SRZOXY(I1L),ARZOXZ(I2L),SRZOXZ(I1L),ARZOYY(I2L),
     $ SRZOYY(I1L),ARZOYZ(I2L),SRZOYZ(I1L),ARZOZZ(I2L),SRZOZZ(I1L),
     $ ASXOXX(I1L),SSXOXX(I1L),ASXOYY(I1L),SSXOYY(I1L),ASXOZZ(I1L),
     $ SSXOZZ(I1L),ASYOXX(I1L),SSYOXX(I1L),ASYOYY(I1L),SSYOYY(I1L),
     $ ASYOZZ(I1L),SSYOZZ(I1L),AHOZOY(I2L),SHOZOY(I2L),AHOZOX(I2L),
     $ SHOZOX(I2L),AHOYOX(I2L),SHOYOX(I2L))

C

         ENDIF

        ENDDO
       ENDDO

      ENDIF
C
C                               ENDE DER GITTERSPEZIFISCHEN VORBELEGUNG
C
C
C
      IF(DREAD) THEN
                    CLOSE(1)
                    CLOSE(3)
      ENDIF
C
C                                                 KONTROLLAUSDRUCK
C


      IF(ITPRIN.EQ.0) THEN
        IPG = IPRGRID
        IP3 = IP3D(IPG)
        IP1 = IP1D(IPG)
      CALL PRLE3M (1,1,4,'J',KMX(IPG),JMX(IPG),IMX(IPG),
     $             U(IP3),1,V(IP3),1,W(IP3),1,P(IP3),1,
     $             G(IP3),1,B(IP3),0,DIV(IP3),0,  1,  1,  4,  1,  1,  4,
     $             DX(IP1),DY(IP1),DZ(IP1),X(IP1),Y(IP1),Z(IP1),
     $             KANPR,ITSTEP)
      ENDIF
C
C
C                                 ZUERST WIRD EIN DIVERGENZFREIES FELD
C                                 ERZEUGT !!
C
                               OMBETA = OMG*RHO/(-2.0*DT)
                               WSOR = 1.0



        CALL MGPOISL1    (IDIM3D,IDIM2D,IDIM1D,NBND,
     $                    X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                    RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                    U,V,W,P,DP,G,B,DIV,RES,
     $                    UFR,VFR,WFR,
     $                    UTO,VTO,WTO,PTO,GTO,VRI,HILF,
     $                    WSOR,RHO,DIVG,IPRGRID,
     $                    MINLEVEL,MAXLEVEL,MAXLEVEL,TIMEPH,
     $                    UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                    PFR,GFR,
     $                    GSAW,GSAE,GSAN,GSAS,GSAT,GSAB,
     $                    BP,BU,BV,BW,SDIV
     $                   ,GEOVP,SIPLW,SIPLS,
     $                    SIPLB,SIPLPR,SIPUE,SIPUN,SIPUT
     $                   )
C
C
C               REDUZIERUNG DES DRUCKNIVEAUS FUER GEWUENSCHTE GITTER
C
      DO ILEVEL = MINLEVEL,MAXLEVEL
        DO I = 1,NOFTST(ILEVEL)
           IGRID = IGRDOFTST(I,ILEVEL)
        IF (LPLEVEL(IGRID)) THEN
           CALL MGDIMS  (KK,JJ,II,IGRID)
           CALL MGPOINT (IP3,IP2,IP1,IBB,IB3,IBU,IGRID)
           CALL MGBASB  (NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
           CALL PLEVEL  (KK,JJ,II,KK,JJ,II,P(IP3),BP(IP3),
     $             DDX(IP1), DDY(IP1), DDZ(IP1),NFRO,NRGT,NBOT)

        ENDIF
        ENDDO
      ENDDO

C                                 **************************************
C                                 ENDE DER INITIALISIERUNGSPHASE
C                                 **************************************
C

      DO ILEVEL = MINLEVEL,MAXLEVEL
C----------------------------------------------------------------------
       IF (LSCAI(ILEVEL)) THEN
          NGRIDPERLEVEL = NOFSCAI(ILEVEL)
        DO I=1,NGRIDPERLEVEL
          NOFTHISGRID(I) = IGRDOFSCAI(I,ILEVEL)
        ENDDO
       ELSE
          NGRIDPERLEVEL = NOFTST(ILEVEL)
        DO I=1,NGRIDPERLEVEL
          NOFTHISGRID(I) = IGRDOFTST(I,ILEVEL)
        ENDDO
       ENDIF
C----------------------------------------------------------------------

        DO I = 1,NGRIDPERLEVEL
           IGRID = NOFTHISGRID(I)

               CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $                     NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
C         CALL WRITE3DDIAGY(KK,JJ,II,W(IP3),33,3)
C         CALL WRITE3DDIAGY(KK,JJ,II,U(IP3),43,4)
               IF (LPOISSONDIR(IGRID)) THEN
                  CALL DIVCAL
     $                 (KK,JJ,II,KK,JJ,II,
     $                 RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                 U(IP3),V(IP3),W(IP3),BP(IP3),DIV(IP3),
     $                 1.0,-2,DIVGMX(IGRID),BU(IP3),BV(IP3),
     $                 BW(IP3),SDIV(IP3))
               ENDIF

                  
               IF(IPRINT_WSS .EQ. 1) THEN              	
                  CALL WSSINT(KK,JJ,II,NBND,
     $                 DX(IP1),DY(IP1),DZ(IP1),
     $                 DDX(IP1),DDY(IP1),DDZ(IP1),
     $                 U(IP3),V(IP3),W(IP3),G(IP3),
     $                 WSSX(IGRID),WSSY(IGRID),WSSZ(IGRID),
     $                 GMOL,RHO,UGRID,IC1(IGRID),IC2(IGRID),
     $                 JC1(IGRID),JC2(IGRID),KC1(IGRID),KC2(IGRID),
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP)
               ENDIF
               IF(IPRINT_WNS .EQ. 1) THEN
                  CALL WNSINT(KK,JJ,II,NBND,
     $                 DX(IP1),DY(IP1),DZ(IP1),
     $                 DDX(IP1),DDY(IP1),DDZ(IP1),
     $                 P(IP3),G(IP3),
     $                 WNSX(IGRID),WNSY(IGRID),WNSZ(IGRID),
     $                 GMOL,RHO,UGRID,IC1(IGRID),IC2(IGRID),
     $                 JC1(IGRID),JC2(IGRID),KC1(IGRID),KC2(IGRID),
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP)
               ENDIF
               IF(MTURB .EQ. 0) THEN
                 CALL ERNORM(KK,JJ,II,KK,JJ,II,NBND,
     $                       U(IP3),UO(IP3),EPSU(IGRID))
                 CALL ERNORM(KK,JJ,II,KK,JJ,II,NBND,
     $                       V(IP3),VO(IP3),EPSV(IGRID))
                 CALL ERNORM(KK,JJ,II,KK,JJ,II,NBND,
     $                       W(IP3),WO(IP3),EPSW(IGRID))
               ELSE

                 CALL ENERFG  (KK,JJ,II,KK,JJ,II,
     $           DDX(IP1),DDY(IP1),DDZ(IP1),DX(IP1), DY(IP1), DZ(IP1),
     $             X(IP1),  Y(IP1),  Z(IP1),U(IP3),V(IP3),W(IP3),
     $           HILF(IP3),ESUMG(IGRID),
     $           NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,1,GRADPX(IGRID))
             ENDIF
CTBA4 070503
         
        ENDDO
      ENDDO

       DO IGRID = 1,NGRID
          IF (LSLICE(IGRID)) THEN
          	            
             CALL ITSAMPLE (IGRID,IPCGES,DIVGMX,EPSU,EPSV,EPSW,
     $            ESUMG,ESUMS,WSSX,WSSY,WSSZ,WNSX,WNSY,WNSZ,
     $            UBULK,WALLSSX,WALLSSY,WALLSSZ
     $    ,IDIVMAX,JDIVMAX,KDIVMAX,XDIVMAX,YDIVMAX,ZDIVMAX,GRDIVMAX       
     $                     )                  
          ENDIF


      
             CALL ITINF (ITTOT,TIMEPH,IPCGES(IGRID),DIVGMX(IGRID),IGRID,
     $            CPSEC,MTURB,EPSU(IGRID),EPSV(IGRID),EPSW(IGRID),
     $            ESUMG(IGRID),ESUMS(IGRID),
     $            IPRINT_WSS,WSSX(IGRID),WSSY(IGRID),WSSZ(IGRID),
     $            IPRINT_WNS,WNSX(IGRID),WNSY(IGRID),WNSZ(IGRID)
     $            ,UBULK(IGRID),GRADPX(IGRID),
     $            WALLSSX(1,IGRID),WALLSSY(1,IGRID),WALLSSZ(1,IGRID)
     $            ,IDIVMAX(IGRID),JDIVMAX(IGRID),KDIVMAX(IGRID)
     $            ,XDIVMAX(IGRID),YDIVMAX(IGRID),ZDIVMAX(IGRID)
     $            ,GRDIVMAX(IGRID)      
     $                 )
          
      ENDDO

C#ifdef _MPI_
C                       EINLESEN FUER FREQUENZOPTIMIERUNG (Jens Neumann)
C      CALL READVAROPT(VAROPT,IDIM2D,NDIMVAROPT)
C        IF (MYID .EQ. 0) THEN
C      WRITE(6,*)'OPTF: READ: TIMEOPT,FREQOPT',TIMEALT,FREQOPT
C        ENDIF
C#else
C      CALL READVAROPT(VAROPT,IDIM2D,NDIMVAROPT)
C      WRITE(6,*)'OPTF: READ: TIMEOPT,FREQOPT',TIMEALT,FREQOPT
C#endif

    
C 
C                                 **************************************
C                                 BEGINN DER ZEITITERATIONEN 
C                                 **************************************
C

      
CCCCCCCCCCC ******** Read the PARAMETERS from printdata.dat ******** CCCCCCC

      INUNIT = 1000
      OPEN(UNIT=INUNIT,FILE='printdata.dat',FORM='FORMATTED')

      DO i = 1, numxloc
        READ(INUNIT,*)x_loc_temp(i)
      END DO

      DO i = 1, numzloc
        READ(INUNIT,*)z_loc(i)
      END DO

      req_xloc = 0 
      DO i = 1, numxloc
        READ(INUNIT,*)MYID_temp
        if (MYID_temp .NE. 999) then
          req_xloc = req_xloc + 1
          MYID_loc(req_xloc) = MYID_temp 
          x_loc(req_xloc) = x_loc_temp(i)
      end if
      END DO

      READ(INUNIT,*)yloc_stream

      DO I = 1, numzloc_stream
        READ(INUNIT,*)zloc_stream(I)
      END DO

      CLOSE(INUNIT)
      
      DO 2000 ISTEP=1,MTSTEP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  TESTWEISE VOR AUFRUF ALLER ZEITSCHRITT-
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  UND RANDBED.-ROUTINEN GESCHALTET



CCCCCCCCCCCCCCCCCCCCCCC RANDOME ZAHL PRO ZEITSCHRITT
      RANNUM= RANF()
CCCCCCCCCCCCCCCCCCCCCCC


      DO ILEVEL = MINLEVEL,MAXLEVEL
Cifdef _MPI_
C
C                                 BEI MESSAGE-PASSING, KOMMUNIKATION
C
Cendif

        DO I = 1,NOFTST(ILEVEL)
           IGRID = IGRDOFTST(I,ILEVEL)

         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $           NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
C
C
C
C                                 BELEGUNG DER UO, VO, WO -FELDER MIT
C                                 DEN WERTEN DER U, V, W -FELDER
c        write(0,*) 'before boundmg    myid = ', myid
c        stop
         CALL BOUNDMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $              'X',IGRID,0,0,0,
     $              TIMEPH,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $              PFR,GFR
     $             )      


            CALL BOUNDMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $              'Z',IGRID,0,0,0,TIMEPH,UBA,VBA,WBA,
     $              PBA,GBA,UBO,VBO,WBO,
     $              PFR,GFR
     $             )       
      
c        write(0,*) 'after first boundmg with Z  myid=', myid
c        stop
C

        ENDDO
      ENDDO

           CALL  TST3RK 
     $                   (IDIM3D,IDIM2D,IDIM1D,IDIMA,
     $                         X,     Y,     Z,
     $                        DX,    DY,    DZ,    DDX,    DDY,    DDZ,
     $                       RDX,   RDY,   RDZ,   RDDX,   RDDY,   RDDZ,
     $                    FPSFAK,FCOSMY,FCOSNY,
     $                     IFFTX,IPERMUX,RFFTX,IFFTY,IPERMUY,RFFTY,
     $                         U,     V,     W,     UO,     VO,     WO,
     $                         P,    DP,     G,      B,   BP,     BU,
     $                        BV,    BW,   SDIV,   HILF,    DIV,
     $                      DIVG,    AU,    AV,     AW,
     $                       WCU,    WCV,    WCW,
     $                       UFR,   VFR,   WFR,    PFR,    GFR,
     $                       UTO,VTO,WTO,PTO,GTO,
     $                       UBACK,    VBACK,    WBACK,
     $                       VRI, 
     $                       UI1,   VI1,   WI1,    UI2,    VI2,    WI2,
     $                       GI1,   GI2,
     $                      DUDY,  DUDZ,  DVDX,   DVDZ,   DWDX,   DWDY,
     $                      NBUF,ITSTEP, ITTOT,TIMEPH,
     $                      NBND, CIDEND, IIDEND, RIDEND,
     $                      UBA,VBA,WBA,PBA,GBA,
     $                      GEOVP,UBO,VBO,WBO,PBO,GBO,CONV1S,
     $                      RSGS3,FAKTOR,
     $                      FUI,FVI,FWI,
     $                      FUJ,FVJ,FWJ,
     $                      FUK,FVK,FWK,
     $                      COEFFX,COEFFY,COEFFZ,COEFDX,COEFDY,COEFDZ,
     $                      LCOL,DIAG,RCOL,UZ,RSP,
     $   HILF0,HILF1,HILF2,HILF3,HILF4,HILF5,HILF6,
     $   LUX,LUY,LUZ,LVY,LVZ,LWZ,MUX,MUY,MUZ,MVY,MVZ,MWZ,
     $   FA11,FA12,FA13,FA21,FA22,FA23,FA31,FA32,FA33,
     $   DUDX3,DUDY3,DUDZ3,DVDX3,DVDY3,DVDZ3,
     $   DWDX3,DWDY3,DWDZ3,NORMS,SDUDX3,SDUDY3,SDUDZ3,
     $   SDVDY3,SDVDZ3,SDWDZ3,
     $   UC,VC,WC,UUC,UVC,UWC,VVC,VWC,WWC,CDELTA,
     $   H3D1,H3D2,H3D3,
     $   FELD1,FELD2,FELD3,
     $   TAU11,TAU12,TAU13,TAU21,TAU22,TAU23,TAU31,TAU32,TAU33
     $ ,RES
     $ ,GSAW,GSAE,GSAN,GSAS,GSAT,GSAB
     $ ,SIPLW,SIPLS,SIPLB,SIPLPR,SIPUE,SIPUN,SIPUT
C#if defined _PREPROC_
C     $ ,HPX
C#endif
     $  )


C#ifdef _MPI_
C     WTIME2 = MPI_WTIME()  - WTIME1
C     MTIME2 = FLOAT(MCLOCK())/100.  - MTIME1
C     WRITE (6,*)
C    $'TST1G, WTIME:',WTIME2,'  MTIME:',MTIME2,' DIFF:',WTIME2-MTIME2
C      WTIME1 = MPI_WTIME()
C      MTIME1 = FLOAT(MCLOCK())/100.
Cendif
C
      IF(ITPRIN.GT.0) THEN
      IF(MOD(ITTOT,ITPRIN) .EQ. 0) THEN
        IPG = IPRGRID
        IP3 = IP3D(IPG)
        IP1 = IP1D(IPG)
      CALL PRLE3M (1,1,4,'J',KMX(IPG),JMX(IPG),IMX(IPG),
     $             U(IP3),1,V(IP3),1,W(IP3),1,P(IP3),1,
     $             G(IP3),1,B(IP3),0,DIV(IP3),0,  1,  1,  4,  1,  1,  4,
     $             DX(IP1),DY(IP1),DZ(IP1),X(IP1),Y(IP1),Z(IP1),
     $             KANPR,ITSTEP)
      ENDIF
      ENDIF


C***********************************************************************
C                                 AUFRUF DES PROGRAMMPAKETES FUER
C                                 DEN TRANSPORT MASSELOSER PARTIKEL
C
C123               IF (.NOT.LPART .OR. PRUFUL) GOTO 2210
C
C123               CALL HPRO (DCONT,ITSTEP,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
C123 $                        KSTW,U,UO,V,VO,W,WO,ITBGN,ITTOT,T0,
C123 $                        TIMEPH,DT,IVGL,NBND,II,JJ,KK,
C123 $                        IMXS,JMXS,KMXS,NOUT1S,NOUT2S,
C123 $                        NOUT3S,NCOUNT,ICALL,IWTEND,YV,YHILF1,
C123 $                        YHILF2,K1,K2,K3,K4,YMAT,NPA,
C123 $                        ISTART,MSTART,NWTVGL,RINDEF,NMAX,
C123 $                        NDGL,NREIHE,EPS,NRZUL,PRUFUL,NPRUZL,NPRU)
C***********************************************************************
C
C123  2210 CONTINUE
C
                   IF(LDOEIB) THEN
C
C                                 AUSGABE DER GESCHWINDIGKEITSKOMP.
C                                 UND DES G-FELDES DER EBENEN I=IWRB
C                                 UND I=IWRB+1 (UNFORMATIERTE AUS-
C                                 GABE AUF KANAL 12)
C
                   CALL DEOBI  (KK,JJ,II,KMXS,JMXS,IMXS,
     $                          XTOT(IGRID),YTOT(IGRID),ZTOT(IGRID),
     $                          CIDEND,IIDEND,RIDEND,CIDENT,IIDENT,
     $                          RIDENT,U,V,W,P,G,P1,UGRID,ITSTEP,
     $                          ITTOT,ITINT,MTSTEP,TIMEPH,DT,IWRB)
                   END IF
                   IF(LDOEIC) THEN
C
C                                 AUSGABE DER GESCHWINDIGKEITSKOMP.
C                                 UND DES G-FELDES DER EBENEN I=IWRB
C                                 UND I=IWRB+1 (FORMATIERTE AUSGABE
C                                 AUF KANAL 14)
C
                   CALL DEOCI  (KK,JJ,II,KMXS,JMXS,IMXS,
     $                          XTOT(IGRID),YTOT(IGRID),ZTOT(IGRID),
     $                          CIDEND,IIDEND,RIDEND,CIDENT,IIDENT,
     $                          RIDENT,U,V,W,P,G,P1,UGRID,ITSTEP,
     $                          ITTOT,ITINT,MTSTEP,TIMEPH,DT,IWRB)
                   END IF


c#ifdef  
c#ifdef _INTERPOL_

!!!!     HERE WE CAN EXTRACT A PART OF THE SOLUTION AND INTERPOLATE ON ANOTHER GRID

c      DO ILEVEL = MINLEVEL,MAXLEVEL
c        DO I = 1,NOFTST(ILEVEL)
c           IGRID = IGRDOFTST(I,ILEVEL)
c#ifdef _MPI_
c       IF (MYID .EQ. IDPROCOFGRD(IGRID))  THEN
c#endif
cC
cC                                NUR, FALLS GITTER NICHT GEBIETSZERLEGT
cC
c         IF ( .NOT. LSLICE(IGRID) ) THEN
c
c                     CALL MGDIMS  (KK,JJ,II,IGRID)
c                     CALL MGPOINT (IP3,IP2,IP1,IBB,IB3,IBU,IGRID)
c                     CALL MGBASB  (NFRO,NBAC,NRGT,NLFT,
c     $                             NBOT,NTOP,NCUB,IGRID)
c                     CALL MGDIMA  (KKA,JJA,IIA,IGRID)
c                     CALL MGPOINA (IAV,I2L,I1L,IGRID)
c
cC
c                   IF(ITSTEP .GT. 0 .AND. MOD(ITTOT,ITFLUC) .EQ. 0)
c     $                   THEN
cC                                 STATISTIK WIRD AUF REGULAEREN GITTERN
cC          
cC                      GEMACHT
c      write(*,*) 'taking sample'
cC     
cC     B muss wahrscheinlich durch BP ersetzt werden...
cC
cC


c      call solint (kk,jj,ii,u(ip3),v(ip3),w(ip3),p(ip3),BP(ip3),
c     $x(ip1),y(ip1),z(ip1),dx(ip1),dy(ip1),dz(ip1),timeph
c     $,myid,ddx(ip1),ddy(ip1),ddz(ip1),0,numprocs)
c#ifdef _MPI_
c      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
c      if (myid.eq.0) then
c 
c      call solint (kk,jj,ii,u(ip3),v(ip3),w(ip3),p(ip3),BP(ip3),
c     $x(ip1),y(ip1),z(ip1),dx(ip1),dy(ip1),dz(ip1),timeph
c     $,myid,ddx(ip1),ddy(ip1),ddz(ip1),1,numprocs)
c         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
c      else
c         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
c      end if
c
c#else
c      call solint (kk,jj,ii,u(ip3),v(ip3),w(ip3),p(ip3),BP(ip3),
c     $x(ip1),y(ip1),z(ip1),dx(ip1),dy(ip1),dz(ip1),timeph
c     $,myid,ddx(ip1),ddy(ip1),ddz(ip1),1,numprocs)
c#endif
c
c#ifdef _ANALYZ_
c#ifdef _MPI_ 
c      if (myid.eq.0) then
c#endif
c      call analyz
c#ifdef _MPI_ 
c      endif
c#endif
c#endif
c
c                      END IF
c                      END IF
c#ifdef _MPI_
c               ENDIF
c#endif
c                    ENDDO
c                    ENDDO
c#endif
c#endif

                   IF(MTURB .GT. 0) THEN
C
C                                 STATISTISCHE AUSWERTUNG ZUR SPAETEREN
C                                 BILDUNG VON ENSEMBLE-MITTELWERTEN
C
C                   DO IGRID = 1,NGRID
      DO ILEVEL = MINLEVEL,MAXLEVEL
C----------------------------------------------------------------------
       IF (LSCAI(ILEVEL)) THEN
          NGRIDPERLEVEL = NOFSCAI(ILEVEL)
        DO I=1,NGRIDPERLEVEL
          NOFTHISGRID(I) = IGRDOFSCAI(I,ILEVEL)
        ENDDO
       ELSE
          NGRIDPERLEVEL = NOFTST(ILEVEL)
        DO I=1,NGRIDPERLEVEL
          NOFTHISGRID(I) = IGRDOFTST(I,ILEVEL)
        ENDDO
       ENDIF
C----------------------------------------------------------------------

        DO I = 1,NGRIDPERLEVEL
           IGRID = NOFTHISGRID(I)
C
C                                NUR, FALLS GITTER NICHT GEBIETSZERLEGT
C
         IF ( .NOT. LSLICE(IGRID) ) THEN

                     CALL MGDIMS  (KK,JJ,II,IGRID)
                     CALL MGPOINT (IP3,IP2,IP1,IBB,IB3,IBU,IGRID)
                     CALL MGBASB  (NFRO,NBAC,NRGT,NLFT,
     $                             NBOT,NTOP,NCUB,IGRID)
                     CALL MGDIMA  (KKA,JJA,IIA,IGRID)
                     CALL MGPOINA (IAV,I2L,I1L,IGRID)

C
                      IF(ITSTEP .GT. 0 .AND. MOD(ITTOT,ITFLUC) .EQ. 0)
     $                   THEN
C                                 STATISTIK WIRD AUF REGULAEREN GITTERN 
C                                 GEMACHT
C
                         CALL SUMSTA  (KPP,JPP,IPP,IC1(IGRID),
     $                              IC2(IGRID),JC1(IGRID),JC2(IGRID),
     $                              KC1(IGRID),KC2(IGRID),NFRO,NBAC,
     $                              NRGT,NLFT,NBOT,NTOP,NCUB,
     $                        XHOMOG(IGRID),YHOMOG(IGRID),ZHOMOG(IGRID),
     $                 GMOL,RHO,UGRID,DDX(IP1),DDY(IP1),DDZ(IP1),
     $                 DX(IP1),DY(IP1),DZ(IP1),X(IP1),Y(IP1),Z(IP1),
     $                 RDX(IP1),RDY(IP1),RDZ(IP1),RDDX(IP1),
     $                 RDDY(IP1),RDDZ(IP1),

     $ KK,JJ,II,KKA,JJA,IIA,KK,JJ,II,KKA,JJA,IIA,NBND,ISELEA,ISELEP,
     $ ISAMPA(1,IGRID),ISAMPP(1,IGRID),ISLINP,ISLINA,ISLIDI,RIDENT,
     $   HILF(IP3),     B(IP3),     U(IP3),    UO(IP3),    AU(IAV),
     $     SU(IAV),   UFG(IP3),  AURG(IAV),  SURG(IAV), AUSKG(IAV),
     $  SUSKG(IAV), AUFLG(IAV), SUFLG(IAV),     V(IP3),    VO(IP3),
     $     AV(IAV),    SV(IAV),   VFG(IP3),  AVRG(IAV),  SVRG(IAV),
     $  AVSKG(IAV), SVSKG(IAV), AVFLG(IAV), SVFLG(IAV),     W(IP3),
     $     WO(IP3),    AW(IAV),    SW(IAV),   WFG(IP3),  AWRG(IAV),
     $   SWRG(IAV), AWSKG(IAV), SWSKG(IAV), AWFLG(IAV), SWFLG(IAV),
     $      P(IP3),    AP(IAV),    SP(IAV),   PFG(IP3),  APRG(IAV),
     $   SPRG(IAV), APSKG(IAV), SPSKG(IAV), APFLG(IAV), SPFLG(IAV),
     $      G(IP3),  AEFG(IAV),  SEFG(IAV),  AEFS(IAV),  SEFS(IAV),      
     $   ADFG(IAV),  SDFG(IAV),        
     $ ADUDX2(IAV),SDUDX2(IAV),ADUDY2(IAV),
     $ SDUDY2(IAV),ADUDZ2(IAV),SDUDZ2(IAV),ADVDX2(IAV),SDVDX2(IAV),
     $ ADVDY2(IAV),SDVDY2(IAV),ADVDZ2(IAV),SDVDZ2(IAV),ADWDX2(IAV),
     $ SDWDX2(IAV),ADWDY2(IAV),SDWDY2(IAV),ADWDZ2(IAV),SDWDZ2(IAV),
     $ AUFWFG(IAV),SUFWFG(IAV),AUFWFS(IAV),SUFWFS(IAV),AUFWFM(IAV),
     $ SUFWFM(IAV),AVFWFG(IAV),SVFWFG(IAV),AVFWFS(IAV),SVFWFS(IAV),
     $ AVFWFM(IAV),SVFWFM(IAV),AUFVFG(IAV),SUFVFG(IAV),AUFVFS(IAV),
     $ SUFVFS(IAV),AUFVFM(IAV),SUFVFM(IAV),    OX(IP3),   AOX(IAV),
     $    SOX(IAV),  OXFG(IP3), AOXRG(IAV), SOXRG(IAV),    OY(IP3),
     $    AOY(IAV),   SOY(IAV),  OYFG(IP3), AOYRG(IAV), SOYRG(IAV),
     $     OZ(IP3),   AOZ(IAV),   SOZ(IAV),  OZFG(IP3), AOZRG(IAV),
     $  SOZRG(IAV),   AO2(IAV),   SO2(IAV),  O2FG(IP3), AO2RG(IAV),
     $  SO2RG(IAV),   AHE(IAV),   SHE(IAV),  HEFG(IP3), AHERG(IAV),
     $  SHERG(IAV)
     $   ,BP(IP3),     BU(IP3),     BV(IP3),     BW(IP3)
     $ ,HILF3D1(IP3),HILF3D2(IP3),HILF3D3(IP3)
     $ ,AUUM(IAV),SUUM(IAV),AVVM(IAV),SVVM(IAV),AWWM(IAV),SWWM(IAV)
     $ ,APPM(IAV),SPPM(IAV)
     $ ,AUVM(IAV),SUVM(IAV),AUWM(IAV),SUWM(IAV),AVWM(IAV),SVWM(IAV)
     $ ,AUXUXM(IAV),SUXUXM(IAV),AUYUYM(IAV),SUYUYM(IAV)
     $ ,AUZUZM(IAV),SUZUZM(IAV),AVXVXM(IAV),SVXVXM(IAV)
     $ ,AVYVYM(IAV),SVYVYM(IAV),AVZVZM(IAV),SVZVZM(IAV)
     $ ,AWXWXM(IAV),SWXWXM(IAV),AWYWYM(IAV),SWYWYM(IAV)
     $ ,AWZWZM(IAV),SWZWZM(IAV)

     $ )
                         CALL SUMST1  (ILIMX,RKOMXP,NBND,
     $                                 X(IP1),Y(IP1),Z(IP1),
     $                                 DX(IP1),DY(IP1),DZ(IP1),
     $                                 DDX(IP1),DDY(IP1),DDZ(IP1),
     %                    NXGRAE(IGRID),NYGRAE(IGRID),NZGRAE(IGRID),
     $                        XHOMOG(IGRID),YHOMOG(IGRID),ZHOMOG(IGRID),

     $ KK,JJ,II,KKA,JJA,IIA,KK,JJ,II,KKA,JJA,IIA,JJ1L,JJ2L,ILIMXP,
     $ ISLIDI,LINFB,KKXL,KKYL,KKZL,KSXL,KSYL,KSZL,KH0L,IDIMF,KLASS,FVT,
     $ HVK,ISELEP,ISELEA,ISLINP,ISLINA,KSTW,RIDENT,
     $   HILF( 1 ),   UFG( 1 ),  AURG( 1 ),   VFG( 1 ),  AVRG( 1 ),
     $    WFG( 1 ),  AWRG( 1 ),    OX( 1 ),   AOX( 1 ),  OXFG( 1 ),
     $  AOXRG( 1 ),    OY( 1 ),   AOY( 1 ),  OYFG( 1 ), AOYRG( 1 ),
     $     OZ( 1 ),   AOZ( 1 ),  OZFG( 1 ), AOZRG( 1 ), ARXUU(I2L),
     $  SRXUU(I1L), ARYUU(I2L), SRYUU(I1L), ARZUU(I2L), SRZUU(I1L),
     $  ASXUU(I1L), SSXUU(I1L), ASYUU(I1L), SSYUU(I1L), ASZUU(I1L),
     $  SSZUU(I1L), ARXVV(I2L), SRXVV(I1L), ARYVV(I2L), SRYVV(I1L),
     $  ARZVV(I2L), SRZVV(I1L), ASXVV(I1L), SSXVV(I1L), ASYVV(I1L),
     $  SSYVV(I1L), ASZVV(I1L), SSZVV(I1L), ARXWW(I2L), SRXWW(I1L),
     $  ARYWW(I2L), SRYWW(I1L), ARZWW(I2L), SRZWW(I1L), ASXWW(I1L),
     $  SSXWW(I1L), ASYWW(I1L), SSYWW(I1L), ASZWW(I1L), SSZWW(I1L),
     $  ARXUV(I2L), SRXUV(I1L), ARYUV(I2L), SRYUV(I1L), ARXUW(I2L),
     $  SRXUW(I1L), ARYUW(I2L), SRYUW(I1L), ARYVW(I2L), SRYVW(I1L),
     $  ACXUW(I2L), SCXUW(I1L), ACZUW(I2L), SCZUW(I1L),ARXOXX(I2L),
     $ SRXOXX(I1L),ARXOXY(I2L),SRXOXY(I1L),ARXOXZ(I2L),SRXOXZ(I1L),
     $ ARXOYY(I2L),SRXOYY(I1L),ARXOYZ(I2L),SRXOYZ(I1L),ARXOZZ(I2L),
     $ SRXOZZ(I1L),ARYOXX(I2L),SRYOXX(I1L),ARYOXY(I2L),SRYOXY(I1L),
     $ ARYOXZ(I2L),SRYOXZ(I1L),ARYOYY(I2L),SRYOYY(I1L),ARYOYZ(I2L),
     $ SRYOYZ(I1L),ARYOZZ(I2L),SRYOZZ(I1L),ARZOXX(I2L),SRZOXX(I1L),
     $ ARZOXY(I2L),SRZOXY(I1L),ARZOXZ(I2L),SRZOXZ(I1L),ARZOYY(I2L),
     $ SRZOYY(I1L),ARZOYZ(I2L),SRZOYZ(I1L),ARZOZZ(I2L),SRZOZZ(I1L),
     $ ASXOXX(I1L),SSXOXX(I1L),ASXOYY(I1L),SSXOYY(I1L),ASXOZZ(I1L),
     $ SSXOZZ(I1L),ASYOXX(I1L),SSYOXX(I1L),ASYOYY(I1L),SSYOYY(I1L),
     $ ASYOZZ(I1L),SSYOZZ(I1L),AHOZOY(I2L),SHOZOY(I2L),AHOZOX(I2L),
     $ SHOZOX(I2L),AHOYOX(I2L),SHOYOX(I2L))

C
C        IF (LSLICE(IGRID)) 
C    $       CALL SLICEGRD (U,V,W,P,G,B,HILF,IDIM3D,IGRID)
                      END IF
                      END IF
                    ENDDO
                    ENDDO
                   END IF

CCCCC (NP) TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
c       CALL  WRITE3DMPI(KK,JJ,II,X(IP1),Y(IP1),Z(IP1),
c     $                              W(IP3),10,3,1,ITTOT)


CCCCC (NP) TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST 


C#ifdef _MPI_
C      WTIME2 = MPI_WTIME()  - WTIME1
C      MTIME2 = FLOAT(MCLOCK())/100.  - MTIME1
C      WRITE (6,*)
C     $'BIS ERRNORM, WTIME:',
C     $WTIME2,'  MTIME:',MTIME2,' DIFF:',WTIME2-MTIME2
C       WTIME1 = MPI_WTIME()
C       MTIME1 = FLOAT(MCLOCK())/100.
C#endif
C
C                                 AUSGABE DER ZWISCHENERGEBNISSE
C
         IF(ITSTEP            .EQ. 0       .OR.
     $      ITSTEP            .EQ. MTSTEP  .OR.
     $      ITTOT             .EQ. ITINT   .OR.
     $      MOD (ITTOT,IPINF) .EQ. 0           ) THEN

      DO ILEVEL = MINLEVEL,MAXLEVEL
C----------------------------------------------------------------------
       IF (LSCAI(ILEVEL)) THEN
          NGRIDPERLEVEL = NOFSCAI(ILEVEL)
        DO I=1,NGRIDPERLEVEL
          NOFTHISGRID(I) = IGRDOFSCAI(I,ILEVEL)
        ENDDO
       ELSE
          NGRIDPERLEVEL = NOFTST(ILEVEL)
        DO I=1,NGRIDPERLEVEL
          NOFTHISGRID(I) = IGRDOFTST(I,ILEVEL)
        ENDDO
       ENDIF
C----------------------------------------------------------------------
      
        DO I = 1,NGRIDPERLEVEL
           IGRID = NOFTHISGRID(I)

               CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $                     NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
               IF (LPOISSONDIR(IGRID)) THEN
                  CALL DIVCAL
     $                 (KK,JJ,II,KK,JJ,II,
     $                 RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                 U(IP3),V(IP3),W(IP3),BP(IP3),DIV(IP3),
     $                 1.0,-2,DIVGMX(IGRID),BU(IP3),BV(IP3),
     $                 BW(IP3),SDIV(IP3))
               ENDIF
               IF(IPRINT_WSS .EQ. 1) THEN
                CALL WSSINT(KK,JJ,II,NBND,
     $                    DX(IP1),DY(IP1),DZ(IP1),
     $                    DDX(IP1),DDY(IP1),DDZ(IP1),
     $                    U(IP3),V(IP3),W(IP3),G(IP3),
     $                    WSSX(IGRID),WSSY(IGRID),WSSZ(IGRID),
     $                    GMOL,RHO,UGRID,IC1(IGRID),IC2(IGRID),
     $                    JC1(IGRID),JC2(IGRID),KC1(IGRID),KC2(IGRID),
     $                     NFRO,NBAC,NRGT,NLFT,NBOT,NTOP)
               ENDIF
               IF(IPRINT_WNS .EQ. 1) THEN
                CALL WNSINT(KK,JJ,II,NBND,
     $                    DX(IP1),DY(IP1),DZ(IP1),
     $                    DDX(IP1),DDY(IP1),DDZ(IP1),
     $                    P(IP3),G(IP3),
     $                    WNSX(IGRID),WNSY(IGRID),WNSZ(IGRID),
     $                    GMOL,RHO,UGRID,IC1(IGRID),IC2(IGRID),
     $                    JC1(IGRID),JC2(IGRID),KC1(IGRID),KC2(IGRID),
     $                     NFRO,NBAC,NRGT,NLFT,NBOT,NTOP)
               ENDIF       



               IF(MTURB .EQ. 0) THEN
                 CALL ERNORM(KK,JJ,II,KK,JJ,II,NBND,
     $                       U(IP3),UO(IP3),EPSU(IGRID))
                 CALL ERNORM(KK,JJ,II,KK,JJ,II,NBND,
     $                       V(IP3),VO(IP3),EPSV(IGRID))
                 CALL ERNORM(KK,JJ,II,KK,JJ,II,NBND,
     $                       W(IP3),WO(IP3),EPSW(IGRID))
               ELSE   
                 CALL ENERFG  (KK,JJ,II,KK,JJ,II,
     $           DDX(IP1),DDY(IP1),DDZ(IP1),DX(IP1), DY(IP1), DZ(IP1),
     $             X(IP1),  Y(IP1),  Z(IP1),U(IP3),V(IP3),W(IP3),
     $           HILF(IP3),ESUMG(IGRID),
     $           NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,1,GRADPX(IGRID))
               ENDIF 
CTBA4 070503
     
            ENDDO
            ENDDO
        
C#ifdef _MPI_
C      WTIME2 = MPI_WTIME()  - WTIME1
C      MTIME2 = FLOAT(MCLOCK())/100.  - MTIME1
C      WRITE (6,*)
C     $'ENER, WTIME:',WTIME2,'  MTIME:',MTIME2,' DIFF:',WTIME2-MTIME2
C       WTIME1 = MPI_WTIME()
C       MTIME1 = FLOAT(MCLOCK())/100.
C#endif

       DO IGRID = 1,NGRID
C          write (0,*) 'igrid:',igrid
          IF (LSLICE(IGRID)) THEN               
             CALL ITSAMPLE (IGRID,IPCGES,DIVGMX,EPSU,EPSV,EPSW,
     $            ESUMG,ESUMS,WSSX,WSSY,WSSZ,WNSX,WNSY,WNSZ,
     $            UBULK,WALLSSX,WALLSSY,WALLSSZ
     $    ,IDIVMAX,JDIVMAX,KDIVMAX,XDIVMAX,YDIVMAX,ZDIVMAX,GRDIVMAX     
     $                     )                  
          ENDIF

            CALL ITINF (ITTOT,TIMEPH,IPCGES(IGRID),DIVGMX(IGRID),IGRID,
     $              CPSEC,MTURB,EPSU(IGRID),EPSV(IGRID),EPSW(IGRID),
     $                  ESUMG(IGRID),ESUMS(IGRID),
     $                 IPRINT_WSS,WSSX(IGRID),WSSY(IGRID),WSSZ(IGRID),
     $                 IPRINT_WNS,WNSX(IGRID),WNSY(IGRID),WNSZ(IGRID)
     $                 ,UBULK(IGRID),GRADPX(IGRID),
     $              WALLSSX(1,IGRID),WALLSSY(1,IGRID),WALLSSZ(1,IGRID)
     $             ,IDIVMAX(IGRID),JDIVMAX(IGRID),KDIVMAX(IGRID)
     $             ,XDIVMAX(IGRID),YDIVMAX(IGRID),ZDIVMAX(IGRID)
     $             ,GRDIVMAX(IGRID)      
     $                 )
C
            ENDDO
         ENDIF
         
      if (req_xloc .GT. 0) then

      CALL PRINTDATA(MYID,KK,JJ,II,NBND, U(IP3), V(IP3),
     $	W(IP3), P(IP3), ISTEP, NXSLICE, MTSTEP, req_xloc,
     $  numzloc, numvar, x_loc, z_loc, MYID_loc, ITTOT,
     $  yloc_stream, zloc_stream, numzloc_stream)

      end if
C
C                   AUSGABE DER ZWISCHENERGEBNISSE FERTIG
C
C#ifdef _MPI_
C      WTIME2 = MPI_WTIME()  - WTIME1
C      MTIME2 = FLOAT(MCLOCK())/100.  - MTIME1
C      WRITE (6,*)
C     $'ITINF, WTIME:',WTIME2,'  MTIME:',MTIME2,' DIFF:',WTIME2-MTIME2
C       WTIME1 = MPI_WTIME()
C       MTIME1 = FLOAT(MCLOCK())/100.
C#endif


         IIDENT(30) = ITTOT
         RIDENT(30) = TIMEPH
C
         IF(ITTOT  .EQ. ITINT)  GOTO 2400
C
C         OPEN(57,file='fort.57')
C         REWIND(57)
C         READ(57,*,END=1189) TSTOP
C        WRITE(6,*)TSTOP
         OPEN(57,file="fort.57")
         READ(57,*,iostat=IO57) TSTOP
         CLOSE(57)

         IF((GETSEC(0)-CPSEC .GE. TSTOP).AND.(IO57.EQ.0)) THEN
           WRITE(6,*) TSTOP, GETSEC(0)-CPSEC
           
 2301      FORMAT(1H1,10X,'TSTOP =',F10.2,'Processor time =',F10.2)

           GO TO 2400
         ENDIF
 1189    CONTINUE
C
C
C#ifdef _MPI_
C      WTIME2 = MPI_WTIME()  - WTIME1
C      MTIME2 = FLOAT(MCLOCK())/100.  - MTIME1
C      WRITE (6,*)
C     $'FORT57, WTIME:',WTIME2,'  MTIME:',MTIME2,' DIFF:',WTIME2-MTIME2
C#endif
C        
C
C                                 **************************************
C                                 VOLLSTAENDIGER ZEITSCHRITT ABGESCHLOSS
C                                 **************************************
C                                 HIER ZEITSCHRITTSCHLEIFE BEENDET. AB-
C                                 SPEICHERN BZW AUSDRUCK DER VARIABLEN
C                                 UND PROGRAMMSTOP
 2000 CONTINUE
 2400 CONTINUE



c#ifdef _MPI_
c                                 RAUSSCHREIBEN FUER FREQUENZOPTIMIERUNG
c      IF (MYID .EQ. 0) THEN
c         CALL WRIVAROPT(VAROPT,IDIM2D,NDIMVAROPT)
c      ENDIF
c#else
c      CALL WRIVAROPT(VAROPT,IDIM2D,NDIMVAROPT)
c#endif

      IF(ITPRIN.EQ.-8) THEN
        IPG = IPRGRID
        IP3 = IP3D(IPG)
        IP1 = IP1D(IPG)
      CALL PRLE3M (1,1,4,'J',KMX(IPG),JMX(IPG),IMX(IPG),
     $             U(IP3),1,V(IP3),1,W(IP3),1,P(IP3),1,
     $             G(IP3),1,B(IP3),0,DIV(IP3),0,  1,  1,  4,  1,  1,  4,
     $             DX(IP1),DY(IP1),DZ(IP1),X(IP1),Y(IP1),Z(IP1),
     $             KANPR,ITSTEP)
      ENDIF


C      DO IGRID = 1,NGRID
      DO ILEVEL = MINLEVEL,MAXLEVEL
C----------------------------------------------------------------------
       IF (LSCAI(ILEVEL)) THEN
          NGRIDPERLEVEL = NOFSCAI(ILEVEL)
        DO I=1,NGRIDPERLEVEL
          NOFTHISGRID(I) = IGRDOFSCAI(I,ILEVEL)
        ENDDO
       ELSE
          NGRIDPERLEVEL = NOFTST(ILEVEL)
        DO I=1,NGRIDPERLEVEL
          NOFTHISGRID(I) = IGRDOFTST(I,ILEVEL)
        ENDDO
       ENDIF
C----------------------------------------------------------------------

        DO I = 1,NGRIDPERLEVEL
           IGRID = NOFTHISGRID(I)

C
C                                NUR, FALLS GITTER NICHT GEBIETSZERLEGT
C
         IF ( .NOT. LSLICE(IGRID) ) THEN

C                                 AB HIER WIEDER ALLE GITTER 
C                                 ZUSAMMENGESETZT

                     CALL MGDIMS  (KK,JJ,II,IGRID)
                     CALL MGPOINT (IP3,IP2,IP1,IBB,IB3,IBU,IGRID)
                     CALL MGBASB  (NFRO,NBAC,NRGT,NLFT,
     $                             NBOT,NTOP,NCUB,IGRID)
       CALL MGDIMA  (KKA,JJA,IIA,IGRID)
       CALL MGPOINA (IAV,I2L,I1L,IGRID)
  



      IF(MTURB         .GT. 0) THEN
C
C                                 NEUE STATIST. MITTELWERTE
C
         CALL STMIMP  (DREAD,DCONT,NPRNEU,FPRNEU,
     $     DX(IP1),DY(IP1),DZ(IP1),DDX(IP1),DDY(IP1),DDZ(IP1)
     $     ,X(IP1),Y(IP1),Z(IP1),

     $ KK,JJ,II,KKA,JJA,IIA,KK,JJ,II,KKA,JJA,IIA,NBND,ISELEA,ISELEP,
     $ ISAMPA(1,IGRID),ISAMPP(1,IGRID),ISLINP,ISLINA,ISLIDI,RIDENT,
     $   HILF(IP3),     B(IP3),     U(IP3),    UO(IP3),    AU(IAV),
     $     SU(IAV),   UFG(IP3),  AURG(IAV),  SURG(IAV), AUSKG(IAV),
     $  SUSKG(IAV), AUFLG(IAV), SUFLG(IAV),     V(IP3),    VO(IP3),
     $     AV(IAV),    SV(IAV),   VFG(IP3),  AVRG(IAV),  SVRG(IAV),
     $  AVSKG(IAV), SVSKG(IAV), AVFLG(IAV), SVFLG(IAV),     W(IP3),
     $     WO(IP3),    AW(IAV),    SW(IAV),   WFG(IP3),  AWRG(IAV),
     $   SWRG(IAV), AWSKG(IAV), SWSKG(IAV), AWFLG(IAV), SWFLG(IAV),
     $      P(IP3),    AP(IAV),    SP(IAV),   PFG(IP3),  APRG(IAV),
     $   SPRG(IAV), APSKG(IAV), SPSKG(IAV), APFLG(IAV), SPFLG(IAV),
     $      G(IP3),  AEFG(IAV),  SEFG(IAV),  AEFS(IAV),  SEFS(IAV),      
     $   ADFG(IAV),  SDFG(IAV),        
     $ ADUDX2(IAV),SDUDX2(IAV),ADUDY2(IAV),
     $ SDUDY2(IAV),ADUDZ2(IAV),SDUDZ2(IAV),ADVDX2(IAV),SDVDX2(IAV),
     $ ADVDY2(IAV),SDVDY2(IAV),ADVDZ2(IAV),SDVDZ2(IAV),ADWDX2(IAV),
     $ SDWDX2(IAV),ADWDY2(IAV),SDWDY2(IAV),ADWDZ2(IAV),SDWDZ2(IAV),
     $ AUFWFG(IAV),SUFWFG(IAV),AUFWFS(IAV),SUFWFS(IAV),AUFWFM(IAV),
     $ SUFWFM(IAV),AVFWFG(IAV),SVFWFG(IAV),AVFWFS(IAV),SVFWFS(IAV),
     $ AVFWFM(IAV),SVFWFM(IAV),AUFVFG(IAV),SUFVFG(IAV),AUFVFS(IAV),
     $ SUFVFS(IAV),AUFVFM(IAV),SUFVFM(IAV),    OX(IP3),   AOX(IAV),
     $    SOX(IAV),  OXFG(IP3), AOXRG(IAV), SOXRG(IAV),    OY(IP3),
     $    AOY(IAV),   SOY(IAV),  OYFG(IP3), AOYRG(IAV), SOYRG(IAV),
     $     OZ(IP3),   AOZ(IAV),   SOZ(IAV),  OZFG(IP3), AOZRG(IAV),
     $  SOZRG(IAV),   AO2(IAV),   SO2(IAV),  O2FG(IP3), AO2RG(IAV),
     $  SO2RG(IAV),   AHE(IAV),   SHE(IAV),  HEFG(IP3), AHERG(IAV),
     $  SHERG(IAV)
     $   ,BP(IP3),     BU(IP3),     BV(IP3),     BW(IP3)
     $ ,HILF3D1(IP3),HILF3D2(IP3),HILF3D3(IP3)
     $ ,AUUM(IAV),SUUM(IAV),AVVM(IAV),SVVM(IAV),AWWM(IAV),SWWM(IAV)
     $ ,APPM(IAV),SPPM(IAV)
     $ ,AUVM(IAV),SUVM(IAV),AUWM(IAV),SUWM(IAV),AVWM(IAV),SVWM(IAV)
     $ ,AUXUXM(IAV),SUXUXM(IAV),AUYUYM(IAV),SUYUYM(IAV)
     $ ,AUZUZM(IAV),SUZUZM(IAV),AVXVXM(IAV),SVXVXM(IAV)
     $ ,AVYVYM(IAV),SVYVYM(IAV),AVZVZM(IAV),SVZVZM(IAV)
     $ ,AWXWXM(IAV),SWXWXM(IAV),AWYWYM(IAV),SWYWYM(IAV)
     $ ,AWZWZM(IAV),SWZWZM(IAV)

     $ )

         CALL STMIM1  (DREAD,DCONT,NPRNEU,FPRNEU,ILINT0,ILINTX,ILINTY,
     $                 ILINTZ,

     $ KK,JJ,II,KKA,JJA,IIA,KK,JJ,II,KKA,JJA,IIA,JJ1L,JJ2L,ILIMXP,
     $ ISLIDI,LINFB,KKXL,KKYL,KKZL,KSXL,KSYL,KSZL,KH0L,IDIMF,KLASS,FVT,
     $ HVK,ISELEP,ISELEA,ISLINP,ISLINA,KSTW,RIDENT,
     $   HILF( 1 ),   UFG( 1 ),  AURG( 1 ),   VFG( 1 ),  AVRG( 1 ),
     $    WFG( 1 ),  AWRG( 1 ),    OX( 1 ),   AOX( 1 ),  OXFG( 1 ),
     $  AOXRG( 1 ),    OY( 1 ),   AOY( 1 ),  OYFG( 1 ), AOYRG( 1 ),
     $     OZ( 1 ),   AOZ( 1 ),  OZFG( 1 ), AOZRG( 1 ), ARXUU(I2L),
     $  SRXUU(I1L), ARYUU(I2L), SRYUU(I1L), ARZUU(I2L), SRZUU(I1L),
     $  ASXUU(I1L), SSXUU(I1L), ASYUU(I1L), SSYUU(I1L), ASZUU(I1L),
     $  SSZUU(I1L), ARXVV(I2L), SRXVV(I1L), ARYVV(I2L), SRYVV(I1L),
     $  ARZVV(I2L), SRZVV(I1L), ASXVV(I1L), SSXVV(I1L), ASYVV(I1L),
     $  SSYVV(I1L), ASZVV(I1L), SSZVV(I1L), ARXWW(I2L), SRXWW(I1L),
     $  ARYWW(I2L), SRYWW(I1L), ARZWW(I2L), SRZWW(I1L), ASXWW(I1L),
     $  SSXWW(I1L), ASYWW(I1L), SSYWW(I1L), ASZWW(I1L), SSZWW(I1L),
     $  ARXUV(I2L), SRXUV(I1L), ARYUV(I2L), SRYUV(I1L), ARXUW(I2L),
     $  SRXUW(I1L), ARYUW(I2L), SRYUW(I1L), ARYVW(I2L), SRYVW(I1L),
     $  ACXUW(I2L), SCXUW(I1L), ACZUW(I2L), SCZUW(I1L),ARXOXX(I2L),
     $ SRXOXX(I1L),ARXOXY(I2L),SRXOXY(I1L),ARXOXZ(I2L),SRXOXZ(I1L),
     $ ARXOYY(I2L),SRXOYY(I1L),ARXOYZ(I2L),SRXOYZ(I1L),ARXOZZ(I2L),
     $ SRXOZZ(I1L),ARYOXX(I2L),SRYOXX(I1L),ARYOXY(I2L),SRYOXY(I1L),
     $ ARYOXZ(I2L),SRYOXZ(I1L),ARYOYY(I2L),SRYOYY(I1L),ARYOYZ(I2L),
     $ SRYOYZ(I1L),ARYOZZ(I2L),SRYOZZ(I1L),ARZOXX(I2L),SRZOXX(I1L),
     $ ARZOXY(I2L),SRZOXY(I1L),ARZOXZ(I2L),SRZOXZ(I1L),ARZOYY(I2L),
     $ SRZOYY(I1L),ARZOYZ(I2L),SRZOYZ(I1L),ARZOZZ(I2L),SRZOZZ(I1L),
     $ ASXOXX(I1L),SSXOXX(I1L),ASXOYY(I1L),SSXOYY(I1L),ASXOZZ(I1L),
     $ SSXOZZ(I1L),ASYOXX(I1L),SSYOXX(I1L),ASYOYY(I1L),SSYOYY(I1L),
     $ ASYOZZ(I1L),SSYOZZ(I1L),AHOZOY(I2L),SHOZOY(I2L),AHOZOX(I2L),
     $ SHOZOX(I2L),AHOYOX(I2L),SHOYOX(I2L))

C
         CALL SUMAAN  (

     $ KK,JJ,II,KKA,JJA,IIA,KK,JJ,II,KKA,JJA,IIA,NBND,ISELEA,ISELEP,
     $ ISAMPA(1,IGRID),ISAMPP(1,IGRID),ISLINP,ISLINA,ISLIDI,RIDENT,
     $   HILF(IP3),     B(IP3),     U(IP3),    UO(IP3),    AU(IAV),
     $     SU(IAV),   UFG(IP3),  AURG(IAV),  SURG(IAV), AUSKG(IAV),
     $  SUSKG(IAV), AUFLG(IAV), SUFLG(IAV),     V(IP3),    VO(IP3),
     $     AV(IAV),    SV(IAV),   VFG(IP3),  AVRG(IAV),  SVRG(IAV),
     $  AVSKG(IAV), SVSKG(IAV), AVFLG(IAV), SVFLG(IAV),     W(IP3),
     $     WO(IP3),    AW(IAV),    SW(IAV),   WFG(IP3),  AWRG(IAV),
     $   SWRG(IAV), AWSKG(IAV), SWSKG(IAV), AWFLG(IAV), SWFLG(IAV),
     $      P(IP3),    AP(IAV),    SP(IAV),   PFG(IP3),  APRG(IAV),
     $   SPRG(IAV), APSKG(IAV), SPSKG(IAV), APFLG(IAV), SPFLG(IAV),
     $      G(IP3),  AEFG(IAV),  SEFG(IAV),  AEFS(IAV),  SEFS(IAV),      
     $   ADFG(IAV),  SDFG(IAV),        
     $ ADUDX2(IAV),SDUDX2(IAV),ADUDY2(IAV),
     $ SDUDY2(IAV),ADUDZ2(IAV),SDUDZ2(IAV),ADVDX2(IAV),SDVDX2(IAV),
     $ ADVDY2(IAV),SDVDY2(IAV),ADVDZ2(IAV),SDVDZ2(IAV),ADWDX2(IAV),
     $ SDWDX2(IAV),ADWDY2(IAV),SDWDY2(IAV),ADWDZ2(IAV),SDWDZ2(IAV),
     $ AUFWFG(IAV),SUFWFG(IAV),AUFWFS(IAV),SUFWFS(IAV),AUFWFM(IAV),
     $ SUFWFM(IAV),AVFWFG(IAV),SVFWFG(IAV),AVFWFS(IAV),SVFWFS(IAV),
     $ AVFWFM(IAV),SVFWFM(IAV),AUFVFG(IAV),SUFVFG(IAV),AUFVFS(IAV),
     $ SUFVFS(IAV),AUFVFM(IAV),SUFVFM(IAV),    OX(IP3),   AOX(IAV),
     $    SOX(IAV),  OXFG(IP3), AOXRG(IAV), SOXRG(IAV),    OY(IP3),
     $    AOY(IAV),   SOY(IAV),  OYFG(IP3), AOYRG(IAV), SOYRG(IAV),
     $     OZ(IP3),   AOZ(IAV),   SOZ(IAV),  OZFG(IP3), AOZRG(IAV),
     $  SOZRG(IAV),   AO2(IAV),   SO2(IAV),  O2FG(IP3), AO2RG(IAV),
     $  SO2RG(IAV),   AHE(IAV),   SHE(IAV),  HEFG(IP3), AHERG(IAV),
     $  SHERG(IAV)
     $   ,BP(IP3),     BU(IP3),     BV(IP3),     BW(IP3)
     $ ,HILF3D1(IP3),HILF3D2(IP3),HILF3D3(IP3)
     $ ,AUUM(IAV),SUUM(IAV),AVVM(IAV),SVVM(IAV),AWWM(IAV),SWWM(IAV)
     $ ,APPM(IAV),SPPM(IAV)
     $ ,AUVM(IAV),SUVM(IAV),AUWM(IAV),SUWM(IAV),AVWM(IAV),SVWM(IAV)
     $ ,AUXUXM(IAV),SUXUXM(IAV),AUYUYM(IAV),SUYUYM(IAV)
     $ ,AUZUZM(IAV),SUZUZM(IAV),AVXVXM(IAV),SVXVXM(IAV)
     $ ,AVYVYM(IAV),SVYVYM(IAV),AVZVZM(IAV),SVZVZM(IAV)
     $ ,AWXWXM(IAV),SWXWXM(IAV),AWYWYM(IAV),SWYWYM(IAV)
     $ ,AWZWZM(IAV),SWZWZM(IAV)

     $ )

      END IF
C
C                                 STABILITAETSANALYSE
C
      CALL STABED (KK,JJ,II,KK,JJ,II,DDX(IP1),DDY(IP1),DDZ(IP1),
     $             U(IP3),V(IP3),W(IP3),G(IP3),BP(IP3),RHO,
     $             DT,XTOT(IGRID),YTOT(IGRID),ZTOT(IGRID),IGRID)
C
C                                 RUECKTRANSFORMATION DES GESCHW.-
C                                 FELDES (GALILEI-TRANSFORMATION)
C
      CALL VSHIFT (KK,JJ,II,KMXS,JMXS,IMXS,U,UGRID,'A')
            CALL BOUNDMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $              'T',IGRID,0,0,0,
     $              TIMEPH,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $              PFR,GFR
     $             )       
C



      CALL SETREF (KK,JJ,II,KKA,JJA,IIA,KK,JJ,II,KKA,JJA,IIA,
     $             KPP,JPP,IPP,DDX,DDY,DDZ,
     $             ZTOT(IGRID),YTOT(IGRID),XTOT(IGRID),XREF,
     $             RHO,GMOL,GRADPX,IIDENT,RIDENT,AU(IAV),
     $             ALREF,UREF,TRF,ISETRE,
     $             XHOMOG(IGRID),YHOMOG(IGRID),ZHOMOG(IGRID),IGRID)

       ENDIF


      ENDDO
      ENDDO

      DO IGRID = 1,NGRID
        
        IF (LSLICE(IGRID)) THEN
         
        DO I = 1, 752
          ISAMPP (I,IGRID) = ISAMPP (I,IGRDOFSLCHILD(1,IGRID))
          ISAMPA (I,IGRID) = ISAMPA (I,IGRDOFSLCHILD(1,IGRID))
        ENDDO
   
        ENDIF
      ENDDO

C                      SCHREIBEN DES KOPFES DES DATENFILES
C
      IF(DWRITE) THEN
C
         IF((LDOB).AND.(.NOT.LDOC)) THEN
            OPEN (2,FILE='fort.2',FORM='UNFORMATTED')
            CALL DOBHEAD (CIDENT,IIDENT,RIDENT,
     $                    ISELEP,ISAMPP,ISLINP,ISLIDI)
         ELSEIF((LDOC).AND.(.NOT.LDOB)) THEN
            OPEN (4,FILE='fort.4',FORM='FORMATTED')
            CALL DOCHEAD (CIDENT,IIDENT,RIDENT,
     $                    ISELEP,ISAMPP,ISLINP,ISLIDI)
         ELSE
            CALL ERRR (503,' MLET ')
         ENDIF
C---------------------------- BEI EINFACHEM MPI MUSS JEDER PROZESSOR SCHREIBEN
C---------------------------- BEI _MPI_PARTICLES_ NUR DER MASTER

         DO IGRID = 1,NGRID
            CALL DOGRID (CIDENT,IIDENT,RIDENT,
     $                   DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,
     $              U,V,W,P,G,B,HILF,
     $              IDIM3D,IDIM2D,IDIM1D,NBUF,NBND,IDIMA,IDIM1L,IDIM2L,
     $              IGRID)

            IF(MTURB .GT. 0) THEN

               CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $                   NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
               CALL MGDIMA  (KKA,JJA,IIA,IGRID)
               CALL MGPOINA (IAV,I2L,I1L,IGRID)

c hier Fehler bei Aufruf dobca... scheinbar BP Feld nicht ok
c      call WRITE3DDIAGY (KK,JJ,II,BP,20,1)
c hier BP Feld noch ok

               IF(LDOB  ) THEN
                  CALL DOBCA  (  2  ,'BINAER  ',IIDENT,
     $                ISTALM,ISTPLM,JSTALM,JSTPLM,
     $                ITTOT,ITFLUC,DT,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                ILIMX ,IDOBD1,IDOBD2,IDOBD3,IDOBD4,RKOMXP,RDOBD1,
     $                RDOBD2,RDOBD3,RDOBD4,IGRID,IDIM3D,
     $                XHOMOG(IGRID),YHOMOG(IGRID),ZHOMOG(IGRID),IDIMA,

     $ KK,JJ,II,KKA,JJA,IIA,KK,JJ,II,KKA,JJA,IIA,NBND,
C-------------------------- AB HIER, WIE cstaca.h:
     $ ISELEA,ISELEP,
     $ ISAMPA,ISAMPP,ISLINP,ISLINA,ISLIDI,
     $ RIDENT,HILF,B,U,UO,AU,SU,UFG,AURG,SURG,AUSKG,SUSKG,AUFLG,SUFLG,V,
     $ VO,AV,SV,VFG,AVRG,SVRG,AVSKG,SVSKG,AVFLG,SVFLG,W,WO,AW,SW,WFG,
     $ AWRG,SWRG,AWSKG,SWSKG,AWFLG,SWFLG,P,AP,SP,PFG,APRG,SPRG,
     $ APSKG,SPSKG,APFLG,SPFLG,
     $ G,AEFG,SEFG,AEFS,SEFS,ADFG,SDFG,
     $ ADUDX2,SDUDX2,ADUDY2,SDUDY2,ADUDZ2,SDUDZ2,ADVDX2,SDVDX2,ADVDY2,
     $ SDVDY2,ADVDZ2,SDVDZ2,ADWDX2,SDWDX2,ADWDY2,SDWDY2,ADWDZ2,SDWDZ2,
     $ AUFWFG,SUFWFG,AUFWFS,SUFWFS,AUFWFM,SUFWFM,AVFWFG,SVFWFG,AVFWFS,
     $ SVFWFS,AVFWFM,SVFWFM,AUFVFG,SUFVFG,AUFVFS,SUFVFS,AUFVFM,
     $ SUFVFM,OX,AOX,SOX,OXFG,AOXRG,SOXRG,OY,AOY,SOY,OYFG,
     $ AOYRG,SOYRG,OZ,AOZ,SOZ,OZFG,AOZRG,SOZRG,AO2,SO2,O2FG,AO2RG,
     $ SO2RG,AHE,SHE,HEFG,AHERG,SHERG
     $ ,BP,BU,BV,BW
     $ ,HILF3D1,HILF3D2,HILF3D3
     $ ,AUUM,SUUM,AVVM,SVVM,AWWM,SWWM,APPM,SPPM
     $ ,AUVM,SUVM,AUWM,SUWM,AVWM,SVWM
     $ ,AUXUXM,SUXUXM,AUYUYM,SUYUYM,AUZUZM,SUZUZM,AVXVXM,SVXVXM,AVYVYM
     $ ,SVYVYM,AVZVZM,SVZVZM,AWXWXM,SWXWXM,AWYWYM,SWYWYM,AWZWZM,SWZWZM
     $ )

                  CALL DOBCA1 (2,'BINAER  ',ISTALM,ISTPLM,JSTALM,JSTPLM,
     $                ITTOT,ITFLUC,DT,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                ILINT0,ILINTX,ILINTY,ILINTZ,

     $ KK,JJ,II,KKA,JJA,IIA,KK,JJ,II,KKA,JJA,IIA,JJ1L,JJ2L,ILIMXP,
     $ ISLIDI,LINFB,KKXL,KKYL,KKZL,KSXL,KSYL,KSZL,KH0L,IDIMF,KLASS,FVT,
     $ HVK,ISELEP,ISELEA,ISLINP,ISLINA,KSTW,RIDENT,
     $   HILF( 1 ),   UFG( 1 ),  AURG( 1 ),   VFG( 1 ),  AVRG( 1 ),
     $    WFG( 1 ),  AWRG( 1 ),    OX( 1 ),   AOX( 1 ),  OXFG( 1 ),
     $  AOXRG( 1 ),    OY( 1 ),   AOY( 1 ),  OYFG( 1 ), AOYRG( 1 ),
     $     OZ( 1 ),   AOZ( 1 ),  OZFG( 1 ), AOZRG( 1 ), ARXUU(I2L),
     $  SRXUU(I1L), ARYUU(I2L), SRYUU(I1L), ARZUU(I2L), SRZUU(I1L),
     $  ASXUU(I1L), SSXUU(I1L), ASYUU(I1L), SSYUU(I1L), ASZUU(I1L),
     $  SSZUU(I1L), ARXVV(I2L), SRXVV(I1L), ARYVV(I2L), SRYVV(I1L),
     $  ARZVV(I2L), SRZVV(I1L), ASXVV(I1L), SSXVV(I1L), ASYVV(I1L),
     $  SSYVV(I1L), ASZVV(I1L), SSZVV(I1L), ARXWW(I2L), SRXWW(I1L),
     $  ARYWW(I2L), SRYWW(I1L), ARZWW(I2L), SRZWW(I1L), ASXWW(I1L),
     $  SSXWW(I1L), ASYWW(I1L), SSYWW(I1L), ASZWW(I1L), SSZWW(I1L),
     $  ARXUV(I2L), SRXUV(I1L), ARYUV(I2L), SRYUV(I1L), ARXUW(I2L),
     $  SRXUW(I1L), ARYUW(I2L), SRYUW(I1L), ARYVW(I2L), SRYVW(I1L),
     $  ACXUW(I2L), SCXUW(I1L), ACZUW(I2L), SCZUW(I1L),ARXOXX(I2L),
     $ SRXOXX(I1L),ARXOXY(I2L),SRXOXY(I1L),ARXOXZ(I2L),SRXOXZ(I1L),
     $ ARXOYY(I2L),SRXOYY(I1L),ARXOYZ(I2L),SRXOYZ(I1L),ARXOZZ(I2L),
     $ SRXOZZ(I1L),ARYOXX(I2L),SRYOXX(I1L),ARYOXY(I2L),SRYOXY(I1L),
     $ ARYOXZ(I2L),SRYOXZ(I1L),ARYOYY(I2L),SRYOYY(I1L),ARYOYZ(I2L),
     $ SRYOYZ(I1L),ARYOZZ(I2L),SRYOZZ(I1L),ARZOXX(I2L),SRZOXX(I1L),
     $ ARZOXY(I2L),SRZOXY(I1L),ARZOXZ(I2L),SRZOXZ(I1L),ARZOYY(I2L),
     $ SRZOYY(I1L),ARZOYZ(I2L),SRZOYZ(I1L),ARZOZZ(I2L),SRZOZZ(I1L),
     $ ASXOXX(I1L),SSXOXX(I1L),ASXOYY(I1L),SSXOYY(I1L),ASXOZZ(I1L),
     $ SSXOZZ(I1L),ASYOXX(I1L),SSYOXX(I1L),ASYOYY(I1L),SSYOYY(I1L),
     $ ASYOZZ(I1L),SSYOZZ(I1L),AHOZOY(I2L),SHOZOY(I2L),AHOZOX(I2L),
     $ SHOZOX(I2L),AHOYOX(I2L),SHOYOX(I2L))

               END IF
               IF(LDOC  ) THEN
                  CALL DOBCA  (  4  ,'CODIERT ',IIDENT,
     $                ISTALM,ISTPLM,JSTALM,JSTPLM,
     $                ITTOT,ITFLUC,DT,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                ILIMX ,IDOBD1,IDOBD2,IDOBD3,IDOBD4,RKOMXP,RDOBD1,
     $                RDOBD2,RDOBD3,RDOBD4,IGRID,IDIM3D,
     $                XHOMOG(IGRID),YHOMOG(IGRID),ZHOMOG(IGRID),IDIMA,

     $ KK,JJ,II,KKA,JJA,IIA,KK,JJ,II,KKA,JJA,IIA,NBND,
C-------------------------- AB HIER, WIE cstaca.h:
     $ ISELEA,ISELEP,
     $ ISAMPA,ISAMPP,ISLINP,ISLINA,ISLIDI,
     $ RIDENT,HILF,B,U,UO,AU,SU,UFG,AURG,SURG,AUSKG,SUSKG,AUFLG,SUFLG,V,
     $ VO,AV,SV,VFG,AVRG,SVRG,AVSKG,SVSKG,AVFLG,SVFLG,W,WO,AW,SW,WFG,
     $ AWRG,SWRG,AWSKG,SWSKG,AWFLG,SWFLG,P,AP,SP,PFG,APRG,SPRG,
     $ APSKG,SPSKG,APFLG,SPFLG,
     $ G,AEFG,SEFG,AEFS,SEFS,ADFG,SDFG,
     $ ADUDX2,SDUDX2,ADUDY2,SDUDY2,ADUDZ2,SDUDZ2,ADVDX2,SDVDX2,ADVDY2,
     $ SDVDY2,ADVDZ2,SDVDZ2,ADWDX2,SDWDX2,ADWDY2,SDWDY2,ADWDZ2,SDWDZ2,
     $ AUFWFG,SUFWFG,AUFWFS,SUFWFS,AUFWFM,SUFWFM,AVFWFG,SVFWFG,AVFWFS,
     $ SVFWFS,AVFWFM,SVFWFM,AUFVFG,SUFVFG,AUFVFS,SUFVFS,AUFVFM,
     $ SUFVFM,OX,AOX,SOX,OXFG,AOXRG,SOXRG,OY,AOY,SOY,OYFG,
     $ AOYRG,SOYRG,OZ,AOZ,SOZ,OZFG,AOZRG,SOZRG,AO2,SO2,O2FG,AO2RG,
     $ SO2RG,AHE,SHE,HEFG,AHERG,SHERG
     $ ,BP,BU,BV,BW
     $ ,HILF3D1,HILF3D2,HILF3D3
     $ ,AUUM,SUUM,AVVM,SVVM,AWWM,SWWM,APPM,SPPM
     $ ,AUVM,SUVM,AUWM,SUWM,AVWM,SVWM
     $ ,AUXUXM,SUXUXM,AUYUYM,SUYUYM,AUZUZM,SUZUZM,AVXVXM,SVXVXM,AVYVYM
     $ ,SVYVYM,AVZVZM,SVZVZM,AWXWXM,SWXWXM,AWYWYM,SWYWYM,AWZWZM,SWZWZM
     $ )

                  CALL DOBCA1 (4,'CODIERT ',ISTALM,ISTPLM,JSTALM,JSTPLM,
     $                ITTOT,ITFLUC,DT,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                ILINT0,ILINTX,ILINTY,ILINTZ,

     $ KK,JJ,II,KKA,JJA,IIA,KK,JJ,II,KKA,JJA,IIA,JJ1L,JJ2L,ILIMXP,
     $ ISLIDI,LINFB,KKXL,KKYL,KKZL,KSXL,KSYL,KSZL,KH0L,IDIMF,KLASS,FVT,
     $ HVK,ISELEP,ISELEA,ISLINP,ISLINA,KSTW,RIDENT,
     $   HILF( 1 ),   UFG( 1 ),  AURG( 1 ),   VFG( 1 ),  AVRG( 1 ),
     $    WFG( 1 ),  AWRG( 1 ),    OX( 1 ),   AOX( 1 ),  OXFG( 1 ),
     $  AOXRG( 1 ),    OY( 1 ),   AOY( 1 ),  OYFG( 1 ), AOYRG( 1 ),
     $     OZ( 1 ),   AOZ( 1 ),  OZFG( 1 ), AOZRG( 1 ), ARXUU(I2L),
     $  SRXUU(I1L), ARYUU(I2L), SRYUU(I1L), ARZUU(I2L), SRZUU(I1L),
     $  ASXUU(I1L), SSXUU(I1L), ASYUU(I1L), SSYUU(I1L), ASZUU(I1L),
     $  SSZUU(I1L), ARXVV(I2L), SRXVV(I1L), ARYVV(I2L), SRYVV(I1L),
     $  ARZVV(I2L), SRZVV(I1L), ASXVV(I1L), SSXVV(I1L), ASYVV(I1L),
     $  SSYVV(I1L), ASZVV(I1L), SSZVV(I1L), ARXWW(I2L), SRXWW(I1L),
     $  ARYWW(I2L), SRYWW(I1L), ARZWW(I2L), SRZWW(I1L), ASXWW(I1L),
     $  SSXWW(I1L), ASYWW(I1L), SSYWW(I1L), ASZWW(I1L), SSZWW(I1L),
     $  ARXUV(I2L), SRXUV(I1L), ARYUV(I2L), SRYUV(I1L), ARXUW(I2L),
     $  SRXUW(I1L), ARYUW(I2L), SRYUW(I1L), ARYVW(I2L), SRYVW(I1L),
     $  ACXUW(I2L), SCXUW(I1L), ACZUW(I2L), SCZUW(I1L),ARXOXX(I2L),
     $ SRXOXX(I1L),ARXOXY(I2L),SRXOXY(I1L),ARXOXZ(I2L),SRXOXZ(I1L),
     $ ARXOYY(I2L),SRXOYY(I1L),ARXOYZ(I2L),SRXOYZ(I1L),ARXOZZ(I2L),
     $ SRXOZZ(I1L),ARYOXX(I2L),SRYOXX(I1L),ARYOXY(I2L),SRYOXY(I1L),
     $ ARYOXZ(I2L),SRYOXZ(I1L),ARYOYY(I2L),SRYOYY(I1L),ARYOYZ(I2L),
     $ SRYOYZ(I1L),ARYOZZ(I2L),SRYOZZ(I1L),ARZOXX(I2L),SRZOXX(I1L),
     $ ARZOXY(I2L),SRZOXY(I1L),ARZOXZ(I2L),SRZOXZ(I1L),ARZOYY(I2L),
     $ SRZOYY(I1L),ARZOYZ(I2L),SRZOYZ(I1L),ARZOZZ(I2L),SRZOZZ(I1L),
     $ ASXOXX(I1L),SSXOXX(I1L),ASXOYY(I1L),SSXOYY(I1L),ASXOZZ(I1L),
     $ SSXOZZ(I1L),ASYOXX(I1L),SSYOXX(I1L),ASYOYY(I1L),SSYOYY(I1L),
     $ ASYOZZ(I1L),SSYOZZ(I1L),AHOZOY(I2L),SHOZOY(I2L),AHOZOX(I2L),
     $ SHOZOX(I2L),AHOYOX(I2L),SHOYOX(I2L))

               END IF
            END IF
         ENDDO

         write (6,*) 'Closing channel 2'
         CLOSE (2)
         write (6,*) 'Closing channel 23'
         CLOSE(23)
         write (6,*) 'Closing channel 4'
         CLOSE (4)

      END IF
C


      IF (MTURB.EQ.1)THEN
C
C                                 ERMITTLUNG VON MAXIMA, MINIMA,
C                                 INTEGRALEN LAENGENMASSEN VON "LINIEN"-
C                                 FELDERN
C
         CALL LININF  (ILINT0, ILINTX, ILINTY, ILINTZ,

     $ KK,JJ,II,KKA,JJA,IIA,KK,JJ,II,KKA,JJA,IIA,JJ1L,JJ2L,ILIMXP,
     $ ISLIDI,LINFB,KKXL,KKYL,KKZL,KSXL,KSYL,KSZL,KH0L,IDIMF,KLASS,FVT,
     $ HVK,ISELEP,ISELEA,ISLINP,ISLINA,KSTW,RIDENT,
     $   HILF( 1 ),   UFG( 1 ),  AURG( 1 ),   VFG( 1 ),  AVRG( 1 ),
     $    WFG( 1 ),  AWRG( 1 ),    OX( 1 ),   AOX( 1 ),  OXFG( 1 ),
     $  AOXRG( 1 ),    OY( 1 ),   AOY( 1 ),  OYFG( 1 ), AOYRG( 1 ),
     $     OZ( 1 ),   AOZ( 1 ),  OZFG( 1 ), AOZRG( 1 ), ARXUU(I2L),
     $  SRXUU(I1L), ARYUU(I2L), SRYUU(I1L), ARZUU(I2L), SRZUU(I1L),
     $  ASXUU(I1L), SSXUU(I1L), ASYUU(I1L), SSYUU(I1L), ASZUU(I1L),
     $  SSZUU(I1L), ARXVV(I2L), SRXVV(I1L), ARYVV(I2L), SRYVV(I1L),
     $  ARZVV(I2L), SRZVV(I1L), ASXVV(I1L), SSXVV(I1L), ASYVV(I1L),
     $  SSYVV(I1L), ASZVV(I1L), SSZVV(I1L), ARXWW(I2L), SRXWW(I1L),
     $  ARYWW(I2L), SRYWW(I1L), ARZWW(I2L), SRZWW(I1L), ASXWW(I1L),
     $  SSXWW(I1L), ASYWW(I1L), SSYWW(I1L), ASZWW(I1L), SSZWW(I1L),
     $  ARXUV(I2L), SRXUV(I1L), ARYUV(I2L), SRYUV(I1L), ARXUW(I2L),
     $  SRXUW(I1L), ARYUW(I2L), SRYUW(I1L), ARYVW(I2L), SRYVW(I1L),
     $  ACXUW(I2L), SCXUW(I1L), ACZUW(I2L), SCZUW(I1L),ARXOXX(I2L),
     $ SRXOXX(I1L),ARXOXY(I2L),SRXOXY(I1L),ARXOXZ(I2L),SRXOXZ(I1L),
     $ ARXOYY(I2L),SRXOYY(I1L),ARXOYZ(I2L),SRXOYZ(I1L),ARXOZZ(I2L),
     $ SRXOZZ(I1L),ARYOXX(I2L),SRYOXX(I1L),ARYOXY(I2L),SRYOXY(I1L),
     $ ARYOXZ(I2L),SRYOXZ(I1L),ARYOYY(I2L),SRYOYY(I1L),ARYOYZ(I2L),
     $ SRYOYZ(I1L),ARYOZZ(I2L),SRYOZZ(I1L),ARZOXX(I2L),SRZOXX(I1L),
     $ ARZOXY(I2L),SRZOXY(I1L),ARZOXZ(I2L),SRZOXZ(I1L),ARZOYY(I2L),
     $ SRZOYY(I1L),ARZOYZ(I2L),SRZOYZ(I1L),ARZOZZ(I2L),SRZOZZ(I1L),
     $ ASXOXX(I1L),SSXOXX(I1L),ASXOYY(I1L),SSXOYY(I1L),ASXOZZ(I1L),
     $ SSXOZZ(I1L),ASYOXX(I1L),SSYOXX(I1L),ASYOYY(I1L),SSYOYY(I1L),
     $ ASYOZZ(I1L),SSYOZZ(I1L),AHOZOY(I2L),SHOZOY(I2L),AHOZOX(I2L),
     $ SHOZOX(I2L),AHOYOX(I2L),SHOYOX(I2L))
      END IF
C
C
C
C                                 AUSGABE DER STATISTISCHEN MITTELWERTE
C
      IF(ITPRIN.EQ.-9) THEN
         DO IGRID = 1,NGRID


               CALL MGDIMA  (KKA,JJA,IIA,IGRID)
               CALL MGPOINA (IAV,I2L,I1L,IGRID)
           CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $                   NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

      CALL PRLE3E (JPP,  1,JPP,'J',IP1,IPS,IP2,KP1,KPS,KP2,
     $             DX(IP1),DY(IP1),DZ(IP1),X(IP1),Y(IP1),Z(IP1),
     $             KANPR,NRRUN,ITFLUC,DT,TRF,
     $             XHOMOG(IGRID),YHOMOG(IGRID),ZHOMOG(IGRID),

     $ KK,JJ,II,KKA,JJA,IIA,KK,JJ,II,KKA,JJA,IIA,NBND,ISELEA,ISELEP,
     $ ISAMPA(1,IGRID),ISAMPP(1,IGRID),ISLINP,ISLINA,ISLIDI,RIDENT,
     $   HILF(IP3),     B(IP3),     U(IP3),    UO(IP3),    AU(IAV),
     $     SU(IAV),   UFG(IP3),  AURG(IAV),  SURG(IAV), AUSKG(IAV),
     $  SUSKG(IAV), AUFLG(IAV), SUFLG(IAV),     V(IP3),    VO(IP3),
     $     AV(IAV),    SV(IAV),   VFG(IP3),  AVRG(IAV),  SVRG(IAV),
     $  AVSKG(IAV), SVSKG(IAV), AVFLG(IAV), SVFLG(IAV),     W(IP3),
     $     WO(IP3),    AW(IAV),    SW(IAV),   WFG(IP3),  AWRG(IAV),
     $   SWRG(IAV), AWSKG(IAV), SWSKG(IAV), AWFLG(IAV), SWFLG(IAV),
     $      P(IP3),    AP(IAV),    SP(IAV),   PFG(IP3),  APRG(IAV),
     $   SPRG(IAV), APSKG(IAV), SPSKG(IAV), APFLG(IAV), SPFLG(IAV),
     $      G(IP3),  AEFG(IAV),  SEFG(IAV),  AEFS(IAV),  SEFS(IAV),      
     $   ADFG(IAV),  SDFG(IAV),        
     $ ADUDX2(IAV),SDUDX2(IAV),ADUDY2(IAV),
     $ SDUDY2(IAV),ADUDZ2(IAV),SDUDZ2(IAV),ADVDX2(IAV),SDVDX2(IAV),
     $ ADVDY2(IAV),SDVDY2(IAV),ADVDZ2(IAV),SDVDZ2(IAV),ADWDX2(IAV),
     $ SDWDX2(IAV),ADWDY2(IAV),SDWDY2(IAV),ADWDZ2(IAV),SDWDZ2(IAV),
     $ AUFWFG(IAV),SUFWFG(IAV),AUFWFS(IAV),SUFWFS(IAV),AUFWFM(IAV),
     $ SUFWFM(IAV),AVFWFG(IAV),SVFWFG(IAV),AVFWFS(IAV),SVFWFS(IAV),
     $ AVFWFM(IAV),SVFWFM(IAV),AUFVFG(IAV),SUFVFG(IAV),AUFVFS(IAV),
     $ SUFVFS(IAV),AUFVFM(IAV),SUFVFM(IAV),    OX(IP3),   AOX(IAV),
     $    SOX(IAV),  OXFG(IP3), AOXRG(IAV), SOXRG(IAV),    OY(IP3),
     $    AOY(IAV),   SOY(IAV),  OYFG(IP3), AOYRG(IAV), SOYRG(IAV),
     $     OZ(IP3),   AOZ(IAV),   SOZ(IAV),  OZFG(IP3), AOZRG(IAV),
     $  SOZRG(IAV),   AO2(IAV),   SO2(IAV),  O2FG(IP3), AO2RG(IAV),
     $  SO2RG(IAV),   AHE(IAV),   SHE(IAV),  HEFG(IP3), AHERG(IAV),
     $  SHERG(IAV)
     $   ,BP(IP3),     BU(IP3),     BV(IP3),     BW(IP3)
     $ ,HILF3D1(IP3),HILF3D2(IP3),HILF3D3(IP3)
     $ ,AUUM(IAV),SUUM(IAV),AVVM(IAV),SVVM(IAV),AWWM(IAV),SWWM(IAV)
     $ ,APPM(IAV),SPPM(IAV)
     $ ,AUVM(IAV),SUVM(IAV),AUWM(IAV),SUWM(IAV),AVWM(IAV),SVWM(IAV)
     $ ,AUXUXM(IAV),SUXUXM(IAV),AUYUYM(IAV),SUYUYM(IAV)
     $ ,AUZUZM(IAV),SUZUZM(IAV),AVXVXM(IAV),SVXVXM(IAV)
     $ ,AVYVYM(IAV),SVYVYM(IAV),AVZVZM(IAV),SVZVZM(IAV)
     $ ,AWXWXM(IAV),SWXWXM(IAV),AWYWYM(IAV),SWYWYM(IAV)
     $ ,AWZWZM(IAV),SWZWZM(IAV)

     $ )

C      CALL PRLE3L (ILINT0,ILINTX,ILINTY,ILINTZ, 1 ,MSLIN,
C     $             KANPR,NRRUN,LINPRN,
C         "cstmg1.h"


       ENDDO
C
      ENDIF
C***********************************************************************
C123  IF (LPART) THEN
C
C                                 ABSPEICHERN ALLER PARTIKELKOORDINATEN
C
C123     CALL WTCON (NRRUN,NPRU,IMXS,JMXS,KMXS,X,Y,Z,DX,DY,DZ,
C123 $               IB1,IB2,JB1,JB2,KB,DT,T0,TIMEPH,NWTVGL,
C123 $               IVGL,ITBGN,NREIHE,ITTOT,NCOUNT,ICALL,
C123 $               NOUT1S,NOUT2S,NOUT3S,IWTEND,YMAT,NPA,
C123 $               NMAX,NDGL,II,JJ,KK)
C123  END IF
C***********************************************************************
        WRITE (6,*) 'VOR WRITE_PARTICLES'

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      WRITE(6,*) ' MLET WURDE ORDNUNGSGEMAESS BEENDET.'
      STOP

      END
