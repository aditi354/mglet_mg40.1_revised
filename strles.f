










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
      SUBROUTINE STRLES
C*MGLET*****************************************************************
C        S T R L E S      STEUERPARAMETER EINLESEN
C*MGLET*****************************************************************
C
C VERS:  11.07.85 (HW)  : STRLES AUS STREAD2 ABGELEITET
C        23.08.85 (HW)  : ERWEITERUNG
C        15.11.85 (HW)  : ERWEITERUNG 'UGRID'
C        04.02.86 (HW)  : ERWEITERUNG 'DCONT'
C        24.04.89 (HW)  : FUER COMMON-BLOCK SAVE-STATEMENT
C        30.10.90 (HW)  : LOGICAL-GROESSEN FUER INPUT/OUTPUT WERDEN
C                         IN EIGENER NAMELIST "INPOUT" EINGELESEN
C        13. 4.92 (MM)  : ERWEITERUNG 'CIDUFR' ENTSCHEIDUNG UEBER
C                         ANSTROEMPROFIL
C        15. 3.93 (MM)  : VOELLIG NEUE VERSION
C        30.12.93 (MM)  : ERWEITERUNG DER GITTERDEFINITIONEN
C        06.11.01 (TB)  : INTRODUCTION OF SCALAR TEMPERATURE TRANSPORT
C                         PARAMETERS (PHYSPAR,BFRONT,BBOTTOM,BCUBE)
C        25.08.03 (FS)  : ERWEITERT FUER T-MISCHER
C        04.01.04 (FS)  : PARTIKEL TRANSPORT
C
C*STARLET***************************************************************
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


      COMMON /CLINOU/

     $                 LDIB,   LDIC,   LDIEIB, LDIEIC,
     $                 LDOB,   LDOC,   LDOEIB, LDOEIC

      SAVE   /CLINOU/
      LOGICAL

     $                 LDIB,   LDIC,   LDIEIB, LDIEIC,
     $                 LDOB,   LDOC,   LDOEIB, LDOEIC

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


      PARAMETER (NMBODY = 100)
      
      COMMON /COBODY/
     &               NBODY,  CTYP,
     &               IB1,   IB2,   JB1,   JB2,   KB1,   KB2,
     &               XB1,   XB2,   YB1,   YB2,   ZB1,   ZB2,
     &               NCOUN, XMIT,  HEIGHT, ALPHA, CDIR

      INTEGER
     &       NBODY,
     &       IB1(NMBODY),   IB2(NMBODY),
     &       JB1(NMBODY),   JB2(NMBODY),
     &       KB1(NMBODY),   KB2(NMBODY),
     &       NCOUN(NMBODY) 

      REAL
     &       XB1(NMBODY),   XB2(NMBODY),
     &       YB1(NMBODY),   YB2(NMBODY),
     &       ZB1(NMBODY),   ZB2(NMBODY),
     &      XMIT(NMBODY),HEIGHT(NMBODY),ALPHA(NMBODY)

      CHARACTER (LEN=16) CTYP(NMBODY),CDIR(NMBODY)
 

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

         PARAMETER   ( NPREC= 16500000 , NMREC= 1 )

      COMMON /CORECDEF/
     &                 NVREC,    ITREC,    MAXREC,
     &                 NGREC,
     &                 CIDREC,   IVEREC,   IVOREC,
     &                 KANREC,   KANGEO,
     &                 INXREC,   INYREC,   INZREC,
     &                 XUGREC,      XOGREC,
     &                 YUGREC,      YOGREC,
     &                 ZUGREC,      ZOGREC

       INTEGER
     +         NGREC(NMREC)
     +        ,KANREC(NMREC),KANGEO(NMREC)
     +        ,IVEREC(NMREC),IVOREC(NMREC)
     +        ,INXREC(NMREC),INYREC(NMREC),INZREC(NMREC)
     +        ,NTREC(NMREC)


       REAL
     +      XUGREC(NMREC),XOGREC(NMREC)
     +     ,YUGREC(NMREC),YOGREC(NMREC)
     +     ,ZUGREC(NMREC),ZOGREC(NMREC)
C
       CHARACTER (LEN=16) CIDREC(NMREC),CIDRE2
C



C

      NAMELIST /STPARM/
     $                  VERS,   NRRUN,
     $                  DREAD,  DWRITE, DCONT,  DCUB,   
     $                  LBODYINSERT,
     $                  LGRIDNEW,       LGRIDREMOVE,
     $                  NPRNEU, FPRNEU,  MTSTEP,
     $                  DT,     ITPRIN, IPINF,  ITINT,
     $                  IPRINT_WSS, IPRINT_WNS, ITFLUC, ITMIT,  
     $                  MPCORR, EPCORR, EPFAK,
     $                  MPCVOR, MPCNACH,  IVPINF,
     $                  OMG,    BETA,   LDIMLO, ISETRE,
     $                  LINPRN, MSLIN,  IWRB,
     $                  NBODY,  NGRID,  NAUFP,
     $                  NVREC,  ITREC,  MAXREC,  LREC 

      NAMELIST /INPOUT/
     $                 LDIB,   LDIC,   LDIEIB, LDIEIC,
     $                 LDOB,   LDOC,   LDOEIB, LDOEIC

      NAMELIST /PHYSPAR/
     $                  MTURB,  TU_LEVEL,   RHO,    GMOL,   UGRID,
     $                  VREF,   EXPON,     UTAUX,
     $                  CIDUFR, UFRCON, UFRFREQ, DELTA,  XREF,
     $                  IPP,    JPP,    KPP,
     $                  XPER,   YPER,   ZPER, 
     $                  NXPER,  NYPER,  NZPER,
     $                  NPPHYS
CTBA2 030203 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



      NAMELIST /PRINTPAR/
     &                IPRGRID,  PRINTFORMAT,
     &                PRINTEBE,   IPRINTEBE,
     &                IP1,  IP2,  IPS,
     &                JP1,  JP2,  JPS,
     &                KP1,  KP2,  KPS,
     &                ISELPREC,
     &                LLU,  LLV,  LLW,  LLP,  LLG,  LLB,	
     &                NPRINTPOS


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC ZU TESTZWECKEN EINLESEN
C                                    AUF KANAL 15
C        WRITE(6,*) 'EINLESEN AUF KANAL 15 !!!!!!!!!!!!!!!!!!!!!!!'
      CALL RMCOMMENT (5,15)

C     OPEN (15,FILE="mrgl",FORM="FORMATTED")
      
      READ (15,STPARM,END=10000)
10000 CONTINUE
      WRITE(6,STPARM)
C
      READ (15,INPOUT,END=10001)
10001 CONTINUE
      WRITE(6,INPOUT)
C
      READ (15,PHYSPAR,END=10002)
10002 CONTINUE
      WRITE(6,PHYSPAR)
C
      CALL READ_PPHYS(NPPHYS_MAX,NPPHYS,PPHYS,XPPHYS)
C
      DO I=1,NPPHYS
	WRITE(6,*)'PPHYS:',PPHYS(I),'XPPHYS:',XPPHYS(I)
      ENDDO



      READ (15,PRINTPAR,END=10003)
10003 CONTINUE
      WRITE(6,PRINTPAR)
C
      DO I=1,NPRINTPOS
                      CALL READ_PRINTPOS (ISTPR(I),JSTPR(I),KSTPR(I))
      ENDDO
C
      DO I=1,NBODY
             CALL READ_BODYDEF 
     &                  (CTYP(I),
     &                   IB1(I),IB2(I),JB1(I),JB2(I),KB1(I),KB2(I),
     &                   XB1(I),XB2(I),YB1(I),YB2(I),ZB1(I),ZB2(I),
     &                   NCOUN(I),XMIT(I),HEIGHT(I),ALPHA(I),CDIR(I))
      ENDDO
C
      DO I=1,NGRID
             CALL READ_GRDDEF
     &                (LEVEL(I),    LCHILD(I),  IPARENT(I),
     &                 IPOSITION(I), JPOSITION(I), KPOSITION(I),
     &                 LTST(I),  LVP(I), LSCAI(I), LPOISSONDIR(I),
     &                 LSLICE(I),NXSLICE(I),NYSLICE(I),NZSLICE(I),
     &                 NVPGRIDS(I),  GRADPX(I), UBULKX(I),
     &                 IMX(I),      JMX(I),      KMX(I),
     &                 XTOT(I),     YTOT(I),     ZTOT(I),
     &                 NBX(I),      NBY(I),      NBZ(I),
     &                 XMIN(I),    YMIN(I),    ZMIN(I),
     &                 NBOCD(1,I),
     &                 FRONT(1,I),    BACK(1,I), 
     &                 RIGHT(1,I),    LEFT(1,I), 
     &                 BOTTOM(1,I),   TOP(1,I),      CUBE(1,I),
     &                 IFRNBR(1,I), IBANBR(1,I), IRINBR(1,I), 
     &                 ILENBR(1,I), IBONBR(1,I), ITONBR(1,I),
     &                 IBPOS(1,1,I),  JBPOS(1,1,I), KBPOS(1,1,I),
     &                 IBANF(1,1,I),  JBANF(1,1,I), KBANF(1,1,I),
     &                 IBEND(1,1,I),  JBEND(1,1,I), KBEND(1,1,I),
     &                 XHOMOG(I),   YHOMOG(I),   ZHOMOG(I),
     &                 LPLEVEL(I),
     &                 NTXB(1,I),   NTYB(1,I),   NTZB(1,I),
     &                 DLX(1,I),    DLY(1,I),    DLZ(1,I),
     &                 NLXB(1,I),   NLYB(1,I),   NLZB(1,I),
     &                 SX(1,I),     SY(1,I),     SZ(1,I),
     &                 NRXB(1,I),   NRYB(1,I),   NRZB(1,I),
     &                 DRX(1,I),    DRY(1,I),    DRZ(1,I),
     &                 MAXBOCONDS,NBSP  , 7,
     &                 XBANF(1,1,I), XBEND(1,1,I), YBANF(1,1,I), 
     &                 YBEND(1,1,I), ZBANF(1,1,I), ZBEND(1,1,I),
     &                 FREQB(1,1,I), 
     &                 ANIVEAU(1,1,I), AUB(1,1,I), AVB(1,1,I), 
     &                 AWB(1,1,I),FLOWTYP(1,1,I),
     &                 CONV1SANF(I),CONV1SEND(I),
     &                 TRANSLES1(I),TRANSLES2(I),WAVENUMBER(1,1,I),
     &                 LOPT(1,I),NTOPT1(1,I),NTOPT2(1,I),XMPOS1(1,I),
     &                 YMPOS1(1,I),ZMPOS1(1,I),XMPOS2(1,I),
     &                 YMPOS2(1,I),ZMPOS2(1,I),PHASE(1,1,I)
CTBA1 030203 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     &                 )
            NGRDDFD=I
            NGRDSET=I
      ENDDO

      DO I=1,NVREC
            CALL READ_RECDEF
     &               (NGREC(I),    CIDREC(I),
     &                IVEREC(I),   IVOREC(I),
     &                KANREC(I),   KANGEO(I),
     &                INXREC(I),   INYREC(I),   INZREC(I),
     &                XUGREC(I),   XOGREC(I),
     &                YUGREC(I),   YOGREC(I),
     &                ZUGREC(I),   ZOGREC(I))
      ENDDO

      RETURN
      END


      SUBROUTINE READ_PRINTPOS (IST,JST,KST)


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

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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

      NAMELIST /PRINTPOS/
     &                ISTPR,  JSTPR,   KSTPR


      READ (15,PRINTPOS,END=10004)
10004 CONTINUE
      WRITE(6,PRINTPOS)

      IST = ISTPR
      JST = JSTPR
      KST = KSTPR

      RETURN
      END


      SUBROUTINE READ_BODYDEF
     &                (CTYPK,
     &                 IB1K,   IB2K,   JB1K,   JB2K,   KB1K,   KB2K,
     &                 XB1K,   XB2K,   YB1K,   YB2K,   ZB1K,   ZB2K,
     &                 NCOUNK, XMITK,  HEIGHTK, ALPHAK, CDIRK)


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

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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

      INTEGER IB1,   IB2,   JB1,   JB2,   KB1,   KB2,   NCOUN,
     &        IB1K,  IB2K,  JB1K,  JB2K,  KB1K,  KB2K,  NCOUNK

      REAL XB1, XB2, YB1, YB2, ZB1, ZB2, XMIT, HEIGHT, ALPHA,
     &     XB1K,XB2K,YB1K,YB2K,ZB1K,ZB2K,XMITK,HEIGHTK,ALPHAK

      CHARACTER (LEN=16) CTYP ,CDIR 
      CHARACTER (LEN=16) CTYPK,CDIRK

      NAMELIST /BODYDEF/
     &                  CTYP,
     &                  IB1,   IB2,   JB1,   JB2,   KB1,   KB2,
     &                  XB1,   XB2,   YB1,   YB2,   ZB1,   ZB2,
     &                  NCOUN, XMIT,  HEIGHT, ALPHA, CDIR


      READ (15,BODYDEF,END=10005)
10005 CONTINUE
      WRITE(6,BODYDEF)

      CTYPK = CTYP

      IB1K  = IB1
      IB2K  = IB2
      JB1K  = JB1
      JB2K  = JB2
      KB1K  = KB1
      KB2K  = KB2

      XB1K  = XB1
      XB2K  = XB2
      YB1K  = YB1
      YB2K  = YB2
      ZB1K  = ZB1
      ZB2K  = ZB2

      NCOUNK = NCOUN
      XMITK  = XMIT
      HEIGHTK = HEIGHT
      ALPHAK = ALPHA
      CDIRK  = CDIR

      RETURN
      END


      SUBROUTINE READ_GRDDEF
     &                (LEVELK,   LCHILDK,   IPARK,
     &                 IPOSK, JPOSK, KPOSK,
     &                 LTSTK, LVPK, LSCAIK, LPOISSONDIRK,  
     &                 LSLICEK,  NXSLICEK,  NYSLICEK,  NZSLICEK,
     &                 NVPK,  GRADPXK,  UBULKXK,
     &                 IMXK,     JMXK,     KMXK,
     &                 XTOTK,    YTOTK,    ZTOTK,
     &                 NBXK,     NBYK,     NBZK,
     &                 XMINK,   YMINK,   ZMINK,
     &                 NBOC,
     &                   FR,     BA,   RI,   LE,   BO,   TO,   CU,
     &                 IFRN,   IBAN, IRIN, ILEN, IBON, ITON,
     &                 IBPOS,  JBPOS,   KBPOS,
     &                 IBANF,  JBANF,   KBANF,
     &                 IBEND,  JBEND,   KBEND,
     &                 XHO,      YHO,      ZHO,   LPK,
     &                 NTXB,     NTYB,     NTZB,
     &                 DLX,      DLY,      DLZ,
     &                 NLXB,     NLYB,     NLZB,
     &                 SX,       SY,       SZ,
     &                 NRXB,     NRYB,     NRZB,
     &                 DRX,      DRY,      DRZ,
     &                 MB,     NBSP  , NDIR,
     &                 XBA2D, XBE2D, YBA2D, YBE2D, ZBA2D, ZBE2D, 
     &                 FRE2D, ANIV2D, AU2D, AV2D, AW2D,FLTYP,
     &                 CONV1SANFK,CONV1SENDK,TRANSLES1K,TRANSLES2K,
     &                 WAVENUM,LOP,NTOP1,NTOP2,XMPO1,YMPO1,ZMPO1,
     &                 XMPO2,YMPO2,ZMPO2,PHAS
CTBA1 030203 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     &                 )     


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

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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

      INTEGER
     &   NBX,          NBY,          NBZ,
     &   NTXB(NBSP),   NTYB(NBSP),   NTZB(NBSP),
     &   NLXB(NBSP),   NLYB(NBSP),   NLZB(NBSP),
     &   NRXB(NBSP),   NRYB(NBSP),   NRZB(NBSP)

      INTEGER
     &    NBOC(7),
     &   IFRN(MB),IBAN(MB),ILEN(MB),IRIN(MB),IBON(MB),ITON(MB),
     &  IBPOS(MB,7),      JBPOS(MB,7),       KBPOS(MB,7),
     &  IBANF(MB,7),      JBANF(MB,7),       KBANF(MB,7),
     &  IBEND(MB,7),      JBEND(MB,7),       KBEND(MB,7),
     &   LOP(MB),NTOP1(MB),NTOP2(MB)


      REAL
     &      XMIN,    YMIN,     ZMIN,
     &  DLX(NBSP),  DLY(NBSP),  DLZ(NBSP),
     &   SX(NBSP),   SY(NBSP),   SZ(NBSP),
     &  DRX(NBSP),  DRY(NBSP),  DRZ(NBSP),
     &  XBA2D(MB,NDIR),XBE2D(MB,NDIR),YBA2D(MB,NDIR),YBE2D(MB,NDIR),
     &  ZBA2D(MB,NDIR),ZBE2D(MB,NDIR),
     &  FRE2D(MB,NDIR),ANIV2D(MB,NDIR),
     &  AU2D(MB,NDIR),AV2D(MB,NDIR),AW2D(MB,NDIR),
     &  WAVENUM(MB,NDIR),
c     ,NTOP1(MB),NTOP2(MB),
     &  XMPO1(MB),YMPO1(MB),
     &  ZMPO1(MB),XMPO2(MB),YMPO2(MB),ZMPO2(MB),PHAS(MB,NDIR)

CTBA1 030203 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      CHARACTER (LEN=16)
     &    FR(MB), BA(MB), RI(MB),  LE(MB),  BO(MB), TO(MB), CU(MB),
     &    FLTYP(MB,NDIR)

      LOGICAL
     &      LCHILD,    LCHILDK,
     &      XHOMOG,   YHOMOG,   ZHOMOG,
     &      XHO,      YHO,      ZHO,
     &      LTSTK,    LVPK, LSCAIK, LPOISSONDIRK,  
     &      LSLICEK,
     &      LTST, LVP, LSCAI, LPOISSONDIR,  
     &      LSLICE,
     &      LPLEVEL,  LPK


      NAMELIST /GRDDEF/
     &                 LEVEL,    LCHILD,   IPARENT,
     &                 IPOSITION, JPOSITION, KPOSITION,
     &                 LTST, LVP, LSCAI, LPOISSONDIR,  
     &                 LSLICE,  NXSLICE,  NYSLICE,  NZSLICE,
     &                 NVPGRIDS,  GRADPX, UBULKX,
     &                 IMX,      JMX,      KMX,
     &                 XTOT,     YTOT,     ZTOT,
     &                 NBX,      NBY,      NBZ,
     &                 XMIN,    YMIN,    ZMIN,
     &                 NBFRONT, NBBACK,    NBRIGHT,    NBLEFT,
     &                 NBBOTTOM,    NBTOP, NBCUBE,
     &                 XHOMOG,   YHOMOG,   ZHOMOG,
     &                 LPLEVEL,
     &                 CONV1SANF,CONV1SEND,TRANSLES1,TRANSLES2

      DATA CONV1SANF,CONV1SEND,TRANSLES1,TRANSLES2 /0.1,0.1,-10000.,0.0/
      DATA FPRNEU /1.0/



      READ (15,GRDDEF,END=10006)
10006 CONTINUE
      WRITE(6,GRDDEF)

      LEVELK = LEVEL
      LCHILDK= LCHILD
      IPARK  = IPARENT
      IPOSK  = IPOSITION
      JPOSK  = JPOSITION
      KPOSK  = KPOSITION
      LTSTK  = LTST
      LVPK   = LVP
      LSCAIK = LSCAI
      LPOISSONDIRK=LPOISSONDIR
      LSLICEK= LSLICE
      NXSLICEK=NXSLICE
      NYSLICEK=NYSLICE
      NZSLICEK=NZSLICE
      NVPK   = NVPGRIDS
      GRADPXK= GRADPX
      UBULKXK= UBULKX
      IMXK   = IMX
      JMXK   = JMX
      KMXK   = KMX
      XTOTK = XTOT
      YTOTK = YTOT
      ZTOTK = ZTOT
      NBXK   = NBX
      NBYK   = NBY
      NBZK   = NBZ
      XMINK = XMIN
      YMINK = YMIN
      ZMINK = ZMIN
      NBOC(1)=NBFRONT
      NBOC(2)=NBBACK
      NBOC(3)=NBRIGHT
      NBOC(4)=NBLEFT
      NBOC(5)=NBBOTTOM
      NBOC(6)=NBTOP
      NBOC(7)=NBCUBE
      XHO  =  XHOMOG
      YHO  =  YHOMOG
      ZHO  =  ZHOMOG
      LPK  =  LPLEVEL
      CONV1SANFK=CONV1SANF
      CONV1SENDK=CONV1SEND
      TRANSLES1K=TRANSLES1
      TRANSLES2K=TRANSLES2

      IF (.NOT.LCHILD) THEN

        DO I=1,NBX
               CALL READ_GRDBLX
     &         (NTXB(I),DLX(I),NLXB(I),SX(I),NRXB(I),DRX(I))
C.....  Fuer Orrsommer stabilitaetsuntersuchung Lx=2*pi

        ENDDO

        DO I=1,NBY
               CALL READ_GRDBLY
     &         (NTYB(I),DLY(I),NLYB(I),SY(I),NRYB(I),DRY(I))
        ENDDO

        DO I=1,NBZ
               CALL READ_GRDBLZ
     &         (NTZB(I),DLZ(I),NLZB(I),SZ(I),NRZB(I),DRZ(I))
        ENDDO

      ENDIF

      DO I=1,NBFRONT
              CALL  READ_BFRDEF (   FR(I),    IFRN(I),
     &                           IBPOS(I,1), JBPOS(I,1), KBPOS(I,1),
     &                           IBANF(I,1), JBANF(I,1), KBANF(I,1),
     &                           IBEND(I,1), JBEND(I,1), KBEND(I,1),
     &                           XBA2D(I,1), XBE2D(I,1),
     &                           YBA2D(I,1), YBE2D(I,1), ZBA2D(I,1),
     &                           ZBE2D(I,1), FRE2D(I,1), ANIV2D(I,1),
     &                           AU2D(I,1), AV2D(I,1), AW2D(I,1),
     &                           FLTYP(I,1),WAVENUM(I,1)
CTBA1 030203 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     &				)
      ENDDO
      DO I=1,NBBACK
              CALL  READ_BBADEF (   BA(I),   IBAN(I),
     &                           IBPOS(I,2), JBPOS(I,2), KBPOS(I,2),
     &                           IBANF(I,2), JBANF(I,2), KBANF(I,2),
     &                           IBEND(I,2), JBEND(I,2), KBEND(I,2),
     &                           XBA2D(I,2), XBE2D(I,2),
     &                           YBA2D(I,2), YBE2D(I,2), ZBA2D(I,2),
     &                           ZBE2D(I,2), FRE2D(I,2), ANIV2D(I,2),
     &                           AU2D(I,2), AV2D(I,2), AW2D(I,2),
     &                           FLTYP(I,2),WAVENUM(I,2)
     &				)
      ENDDO
      DO I=1,NBRIGHT
              CALL  READ_BRIDEF (   RI(I),   IRIN(I),
     &                           IBPOS(I,3), JBPOS(I,3), KBPOS(I,3),
     &                           IBANF(I,3), JBANF(I,3), KBANF(I,3),
     &                           IBEND(I,3), JBEND(I,3), KBEND(I,3),
     &                           XBA2D(I,3), XBE2D(I,3),
     &                           YBA2D(I,3), YBE2D(I,3), ZBA2D(I,3),
     &                           ZBE2D(I,3), FRE2D(I,3), ANIV2D(I,3),
     &                           AU2D(I,3), AV2D(I,3), AW2D(I,3),
     &                           FLTYP(I,3),WAVENUM(I,3)
     &				)
      ENDDO
      DO I=1,NBLEFT
              CALL  READ_BLEDEF (   LE(I),   ILEN(I),
     &                           IBPOS(I,4), JBPOS(I,4), KBPOS(I,4),
     &                           IBANF(I,4), JBANF(I,4), KBANF(I,4),
     &                           IBEND(I,4), JBEND(I,4), KBEND(I,4),
     &                           XBA2D(I,4), XBE2D(I,4),
     &                           YBA2D(I,4), YBE2D(I,4), ZBA2D(I,4),
     &                           ZBE2D(I,4), FRE2D(I,4), ANIV2D(I,4),
     &                           AU2D(I,4), AV2D(I,4), AW2D(I,4),
     &                           FLTYP(I,4),WAVENUM(I,4)
     &				)
      ENDDO
      DO I=1,NBBOTTOM
              CALL  READ_BBODEF (   BO(I),   IBON(I),
     &                           IBPOS(I,5), JBPOS(I,5), KBPOS(I,5),
     &                           IBANF(I,5), JBANF(I,5), KBANF(I,5),
     &                           IBEND(I,5), JBEND(I,5), KBEND(I,5),
     &                           XBA2D(I,5), XBE2D(I,5), YBA2D(I,5),
     &                           YBE2D(I,5), ZBA2D(I,5), ZBE2D(I,5),
     &                           FRE2D(I,5), ANIV2D(I,5),
     &                           AU2D(I,5), AV2D(I,5), AW2D(I,5),
     &                           FLTYP(I,5), WAVENUM(I,5),
     &                           LOP(I),NTOP1(I),NTOP2(I),XMPO1(I),
     &                           YMPO1(I),ZMPO1(I),
     &                           XMPO2(I),YMPO2(I),ZMPO2(I),
     &                           PHAS(I,5)
CTBA1 030203 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     &				)
      ENDDO
      DO I=1,NBTOP
              CALL  READ_BTODEF (   TO(I),   ITON(I),
     &                           IBPOS(I,6), JBPOS(I,6), KBPOS(I,6),
     &                           IBANF(I,6), JBANF(I,6), KBANF(I,6),
     &                           IBEND(I,6), JBEND(I,6), KBEND(I,6),
     &                           XBA2D(I,6), XBE2D(I,6), YBA2D(I,6),
     &                           YBE2D(I,6), ZBA2D(I,6), ZBE2D(I,6),
     &                           FRE2D(I,6), ANIV2D(I,6),
     &                           AU2D(I,6), AV2D(I,6), AW2D(I,6),
     &                           FLTYP(I,6), WAVENUM(I,6),
     &                           LOP(I),NTOP1(I),NTOP2(I),XMPO1(I),
     &                           YMPO1(I),ZMPO1(I),
     &                           XMPO2(I),YMPO2(I),ZMPO2(I),
     &                           PHAS(I,6)
     &				)
      ENDDO
      DO I=1,NBCUBE
              CALL  READ_BCUDEF (   CU(I),
     &                           IBPOS(I,7), JBPOS(I,7), KBPOS(I,7),
     &                           IBANF(I,7), JBANF(I,7), KBANF(I,7),
     &                           IBEND(I,7), JBEND(I,7), KBEND(I,7)
CTBA1 030203 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     &				)
      ENDDO

      RETURN
      END


      SUBROUTINE READ_GRDBLX (NTXBK,DLXK,NLXBK,SXK,NRXBK,DRXK)


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

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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

               NAMELIST /GRDBLX/ NTXB,DLX,NLXB,SX,NRXB,DRX

               READ (15,GRDBLX,END=10007)
10007 CONTINUE
               WRITE(6,GRDBLX)

               NTXBK = NTXB
               DLXK  = DLX
               NLXBK = NLXB
               SXK   = SX
               NRXBK = NRXB
               DRXK  = DRX

               RETURN
      END

      SUBROUTINE READ_GRDBLY (NTYBK,DLYK,NLYBK,SYK,NRYBK,DRYK)


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

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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

               NAMELIST /GRDBLY/ NTYB,DLY,NLYB,SY,NRYB,DRY

               READ (15,GRDBLY,END=10007)
10007 CONTINUE
               WRITE(6,GRDBLY)

               NTYBK = NTYB
               DLYK  = DLY
               NLYBK = NLYB
               SYK   = SY
               NRYBK = NRYB
               DRYK  = DRY

               RETURN
      END

      SUBROUTINE READ_GRDBLZ (NTZBK,DLZK,NLZBK,SZK,NRZBK,DRZK)


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

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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

               NAMELIST /GRDBLZ/ NTZB,DLZ,NLZB,SZ,NRZB,DRZ

               READ (15,GRDBLZ,END=10007)
10007 CONTINUE
               WRITE(6,GRDBLZ)

               NTZBK = NTZB
               DLZK  = DLZ
               NLZBK = NLZB
               SZK   = SZ
               NRZBK = NRZB
               DRZK  = DRZ

               RETURN
      END
      SUBROUTINE READ_RECDEF
     &               (NGREH,   CIDREC,  
     &                IVEREC,   IVOREC,  IKANREC,  IKANGEO,
     &                NXREC,    NYREC,   NZREC,
     &                XUGREC,   XOGREC,
     &                YUGREC,   YOGREC,
     &                ZUGREC,   ZOGREC)


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

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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

      REAL   XUGREC, XOGREC, YUGREC, YOGREC, ZUGREC, ZOGREC,
     &       XUG,    XOG,    YUG,    YOG,    ZUG,    ZOG

      INTEGER NGREH,IVEREC, IVOREC,IKANREC,IKANGEO, NXREC, NYREC, NZREC,
     &        NGREC, IVEL,   IVOR,   KANREC, KANGEO, NX,    NY,    NZ   

       CHARACTER (LEN=16) CIDREC,CID


               NAMELIST /RECDEF/
     &                NGREC, CID,  IVEL,   IVOR,   KANREC,   KANGEO,
     &                NX,   NY,   NZ,  
     &                XUG,   YUG,  ZUG,   XOG,  YOG,   ZOG

               READ (15,RECDEF,END=10008)
10008 CONTINUE
               WRITE(6,RECDEF)


      NGREH  = NGREC
      CIDREC = CID
      IVEREC = IVEL
      IVOREC = IVOR
      IKANREC= KANREC
      IKANGEO= KANGEO
      NXREC  = NX
      NYREC  = NY
      NZREC  = NZ
      XUGREC = XUG
      XOGREC = XOG
      YUGREC = YUG
      YOGREC = YOG
      ZUGREC = ZUG
      ZOGREC = ZOG
      

      RETURN
      END


      SUBROUTINE  READ_BFRDEF (   FR,       IFRN,   
     & IB   ,JB   ,KB   ,IBA  ,JBA  ,KBA  ,IBE  ,JBE  ,KBE ,  
     & XA1D, XE1D, YA1D, YE1D, ZA1D, ZE1D,
     & FRE1D, ANIV1D, AU1D, AV1D, AW1D,FLTYP,WAVENUM
CTBA1 030203 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     & )     
      CHARACTER (LEN=16) FR,FRONT,FLOWTYP,FLTYP
	  DATA NEIGHBOUR /-1/
	  DATA IPOS, JPOS, KPOS / 3, 3, 3/
	  DATA IANF, IEND, JANF, JEND, KANF, KEND / 3, 3, 3, 3, 3, 3/

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

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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


      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

      NAMELIST /BFRONT/  FRONT,NEIGHBOUR,IPOS,JPOS,KPOS,
     &                   JANF,JEND,KANF,KEND,
     &                   XANF,XEND,YANF,YEND,ZANF,ZEND,FREQ,
     &                   ANIV,AU,AV,AW,FLOWTYP,WAVENUMBER
CTBA1 030203 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      XANF = -GREAT
      XEND =  GREAT
      YANF = -GREAT
      YEND =  GREAT
      ZANF = -GREAT
      ZEND =  GREAT

      READ (15,BFRONT,END=100)
  100 CONTINUE
      WRITE( 6,BFRONT)
        FR    =   FRONT
       IFRN   =   NEIGHBOUR
        IB    =   IPOS
        JB    =   JPOS
        KB    =   KPOS
        JBA   =   JANF
        JBE   =   JEND
        KBA   =   KANF
        KBE   =   KEND
        XA1D  =   XANF 
        XE1D  =   XEND 
        YA1D  =   YANF 
        YE1D  =   YEND 
        ZA1D  =   ZANF 
        ZE1D  =   ZEND 
        FRE1D =   FREQ
        ANIV1D =   ANIV
        AU1D   =   AU
        AV1D   =   AV
        AW1D   =   AW
        FLTYP  = FLOWTYP
        WAVENUM = WAVENUMBER
CTBA1 030203 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      RETURN
      END

      SUBROUTINE  READ_BBADEF (   BA,      IBAN,     
     & IB   ,JB   ,KB   ,IBA  ,JBA  ,KBA  ,IBE  ,JBE  ,KBE,
     & XA1D,XE1D,YA1D,YE1D,ZA1D,ZE1D,
     & FRE1D,ANIV1D,AU1D,AV1D,AW1D,FLTYP,WAVENUM
     & )
      CHARACTER (LEN=16) BA,BACK,FLOWTYP,FLTYP
	  DATA NEIGHBOUR /-1/
	  DATA IPOS, JPOS, KPOS / 3, 3, 3/
	  DATA IANF, IEND, JANF,JEND,KANF,KEND / 3, 3, 3, 3, 3, 3/

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

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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


      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

      NAMELIST /BBACK/  BACK,NEIGHBOUR,IPOS,JPOS,KPOS,
     &                   JANF,JEND,KANF,KEND,
     &                   XANF,XEND,YANF,YEND,ZANF,ZEND,FREQ,
     &                   ANIV,AU,AV,AW,FLOWTYP,WAVENUMBER

      XANF = -GREAT
      XEND =  GREAT
      YANF = -GREAT
      YEND =  GREAT
      ZANF = -GREAT
      ZEND =  GREAT

      READ (15,BBACK,END=100)
  100 CONTINUE
      WRITE( 6,BBACK)
        BA    =   BACK
       IBAN   =   NEIGHBOUR
        IB    =   IPOS
        JB    =   JPOS
        KB    =   KPOS
        JBA   =   JANF
        JBE   =   JEND
        KBA   =   KANF
        KBE   =   KEND
        XA1D  =   XANF 
        XE1D  =   XEND 
        YA1D  =   YANF 
        YE1D  =   YEND 
        ZA1D  =   ZANF 
        ZE1D  =   ZEND 
        FRE1D =   FREQ
        ANIV1D =   ANIV
        AU1D   =   AU
        AV1D   =   AV
        AW1D   =   AW
        FLTYP  = FLOWTYP
        WAVENUM = WAVENUMBER
      RETURN
      END

      SUBROUTINE  READ_BRIDEF (   RI,      IRIN,     
     & IB   ,JB   ,KB   ,IBA  ,JBA  ,KBA  ,IBE  ,JBE  ,KBE
     & ,XA1D,XE1D,YA1D,YE1D,ZA1D,ZE1D,
     & FRE1D,ANIV1D,AU1D,AV1D,AW1D,FLTYP,WAVENUM
     & )
      CHARACTER (LEN=16) RI,RIGHT,FLOWTYP,FLTYP
	  DATA NEIGHBOUR /-1/
	  DATA IPOS, JPOS, KPOS / 3, 3, 3/
	  DATA IANF, IEND, JANF,JEND,KANF,KEND / 3, 3, 3, 3, 3, 3/

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

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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


      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

      NAMELIST /BRIGHT/  RIGHT,NEIGHBOUR,IPOS,JPOS,KPOS,
     &                   IANF,IEND,KANF,KEND,
     &                   XANF,XEND,YANF,YEND,ZANF,ZEND,FREQ,
     &                   ANIV,AU,AV,AW,FLOWTYP,WAVENUMBER

      XANF = -GREAT
      XEND =  GREAT
      YANF = -GREAT
      YEND =  GREAT
      ZANF = -GREAT
      ZEND =  GREAT

      READ (15,BRIGHT,END=100)
  100 CONTINUE
      WRITE( 6,BRIGHT)
        RI    =   RIGHT
       IRIN   =   NEIGHBOUR
        IB    =   IPOS
        JB    =   JPOS
        KB    =   KPOS
        IBA   =   IANF
        IBE   =   IEND
        KBA   =   KANF
        KBE   =   KEND
        XA1D  =   XANF 
        XE1D  =   XEND 
        YA1D  =   YANF 
        YE1D  =   YEND 
        ZA1D  =   ZANF 
        ZE1D  =   ZEND 
        FRE1D =   FREQ
        ANIV1D =   ANIV
        AU1D   =   AU
        AV1D   =   AV
        AW1D   =   AW
        FLTYP  = FLOWTYP
        WAVENUM = WAVENUMBER


      RETURN
      END

      SUBROUTINE  READ_BLEDEF (   LE,      ILEN,     
     & IB   ,JB   ,KB   ,IBA  ,JBA  ,KBA  ,IBE  ,JBE  ,KBE,
     & XA1D,XE1D,YA1D,YE1D,ZA1D,ZE1D,
     & FRE1D,ANIV1D,AU1D,AV1D,AW1D,FLTYP,WAVENUM
     & )
      CHARACTER (LEN=16) LE,LEFT,FLOWTYP,FLTYP
	  DATA NEIGHBOUR /-1/
	  DATA IPOS, JPOS, KPOS / 3, 3, 3/
	  DATA IANF, IEND, JANF,JEND,KANF,KEND / 3, 3, 3, 3, 3, 3/

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

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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


      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

      NAMELIST /BLEFT/  LEFT,NEIGHBOUR,IPOS,JPOS,KPOS,
     &                   IANF,IEND,KANF,KEND,
     &                   XANF,XEND,YANF,YEND,ZANF,ZEND,FREQ,
     &                   ANIV,AU,AV,AW,FLOWTYP,WAVENUMBER

      XANF = -GREAT
      XEND =  GREAT
      YANF = -GREAT
      YEND =  GREAT
      ZANF = -GREAT
      ZEND =  GREAT

      READ (15,BLEFT,END=100)
  100 CONTINUE
      WRITE( 6,BLEFT)
        LE    =   LEFT
       ILEN   =   NEIGHBOUR
        IB    =   IPOS
        JB    =   JPOS
        KB    =   KPOS
        IBA   =   IANF
        IBE   =   IEND
        KBA   =   KANF
        KBE   =   KEND
        XA1D  =   XANF 
        XE1D  =   XEND 
        YA1D  =   YANF 
        YE1D  =   YEND 
        ZA1D  =   ZANF 
        ZE1D  =   ZEND 
        FRE1D =   FREQ
        ANIV1D =   ANIV
        AU1D   =   AU
        AV1D   =   AV
        AW1D   =   AW
        FLTYP  = FLOWTYP
        WAVENUM = WAVENUMBER


      RETURN
      END

      SUBROUTINE  READ_BBODEF (   BO,      IBON,     
     & IB   ,JB   ,KB   ,IBA  ,JBA  ,KBA  ,IBE  ,JBE  ,KBE ,
     & XA1D, XE1D, YA1D, YE1D, ZA1D, ZE1D,FRE1D, ANIV1D, AU1D, AV1D,
     & AW1D, FLTYP, WAVENUM,LOP,NTOP1,NTOP2,XMPO1,YMPO1,ZMPO1,
     & XMPO2,YMPO2,ZMPO2,PHAS
CTBA1 030203 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     & )     

	  DATA NEIGHBOUR /-1/
	  DATA IPOS, JPOS, KPOS / 3, 3, 3/
	  DATA IANF, IEND, JANF,JEND,KANF,KEND / 3, 3, 3, 3, 3, 3/

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

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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


      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

      CHARACTER (LEN=16) BO,BOTTOM,FLOWTYP,FLTYP
      NAMELIST /BBOTTOM/  BOTTOM,NEIGHBOUR,IPOS,JPOS,KPOS,
     &                   IANF,IEND,JANF,JEND,
     &                   XANF,XEND,YANF,YEND,ZANF,ZEND,FREQ,
     &                   ANIV,AU,AV,AW,FLOWTYP,WAVENUMBER,
     &                   LOPT,NTOPT1,NTOPT2,XMPOS1,YMPOS1,ZMPOS1,
     &                   XMPOS2,YMPOS2,ZMPOS2,PHASE
CTBA1 030203 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      XANF = -GREAT
      XEND =  GREAT
      YANF = -GREAT
      YEND =  GREAT
      ZANF = -GREAT
      ZEND =  GREAT

      READ (15,BBOTTOM,END=100)
  100 CONTINUE
      WRITE( 6,BBOTTOM)
        BO    =   BOTTOM
       IBON   =   NEIGHBOUR
        IB    =   IPOS
        JB    =   JPOS
        KB    =   KPOS
        IBA   =   IANF
        IBE   =   IEND
        JBA   =   JANF
        JBE   =   JEND
        XA1D  =   XANF 
        XE1D  =   XEND 
        YA1D  =   YANF 
        YE1D  =   YEND 
        ZA1D  =   ZANF 
        ZE1D  =   ZEND 
        FRE1D =   FREQ
        ANIV1D =   ANIV
        AU1D   =   AU
        AV1D   =   AV
        AW1D   =   AW
       FLTYP =   FLOWTYP
       WAVENUM = WAVENUMBER
       LOP    = LOPT
       NTOP1  = NTOPT1
       NTOP2  = NTOPT2
       XMPO1   = XMPOS1
       YMPO1   = YMPOS1
       ZMPO1   = ZMPOS1
       XMPO2   = XMPOS2
       YMPO2   = YMPOS2
       ZMPO2   = ZMPOS2
        PHAS = PHASE
CTBA2 030203 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      RETURN
      END

      SUBROUTINE  READ_BTODEF (   TO,      ITON,     
     & IB   ,JB   ,KB   ,IBA  ,JBA  ,KBA  ,IBE  ,JBE  ,KBE,
     & XA1D, XE1D, YA1D, YE1D, ZA1D, ZE1D,FRE1D, ANIV1D, AU1D, AV1D,
     & AW1D, FLTYP, WAVENUM,LOP,NTOP1,NTOP2,XMPO1,YMPO1,ZMPO1,
     & XMPO2,YMPO2,ZMPO2,PHAS
C-------------------------------------------------------------------
     & )
      CHARACTER (LEN=16) TO,TOP
	  DATA NEIGHBOUR /-1/
	  DATA IPOS, JPOS, KPOS / 3, 3, 3/
	  DATA IANF, IEND, JANF,JEND,KANF,KEND / 3, 3, 3, 3, 3, 3/

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

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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


      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

      NAMELIST /BTOP/  TOP,NEIGHBOUR,IPOS,JPOS,KPOS,
     &                 IANF,IEND,JANF,JEND,
     &                 XANF,XEND,YANF,YEND,ZANF,ZEND,FREQ,
     &                 ANIV,AU,AV,AW,FLOWTYP,WAVENUMBER,
     &                 LOPT,NTOPT1,NTOPT2,XMPOS1,YMPOS1,ZMPOS1,
     &                 XMPOS2,YMPOS2,ZMPOS2,PHASE
CTBA1 030203 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      XANF = -GREAT
      XEND =  GREAT
      YANF = -GREAT
      YEND =  GREAT
      ZANF = -GREAT
      ZEND =  GREAT

      READ (15,BTOP,END=100)
  100 CONTINUE
      WRITE( 6,BTOP)
        TO    =   TOP
       ITON   =   NEIGHBOUR
        IB    =   IPOS
        JB    =   JPOS
        KB    =   KPOS
        IBA   =   IANF
        IBE   =   IEND
        JBA   =   JANF
        JBE   =   JEND
        XA1D  =   XANF 
        XE1D  =   XEND 
        YA1D  =   YANF 
        YE1D  =   YEND 
        ZA1D  =   ZANF 
        ZE1D  =   ZEND 
        FRE1D =   FREQ
        ANIV1D =   ANIV
        AU1D   =   AU
        AV1D   =   AV
        AW1D   =   AW
       FLTYP =   FLOWTYP
       WAVENUM = WAVENUMBER
       LOP    = LOPT
       NTOP1  = NTOPT1
       NTOP2  = NTOPT2
       XMPO1   = XMPOS1
       YMPO1   = YMPOS1
       ZMPO1   = ZMPOS1
       XMPO2   = XMPOS2
       YMPO2   = YMPOS2
       ZMPO2   = ZMPOS2
        PHAS = PHASE
C----------------------------------------------------------------------

      RETURN
      END

      SUBROUTINE  READ_BCUDEF (   CU, 
     & IB   ,JB   ,KB   ,IBA  ,JBA  ,KBA  ,IBE  ,JBE  ,KBE
CTBA1 030203 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     & )     
	  DATA IPOS, JPOS, KPOS / 3, 3, 3/
	  DATA IANF, IEND, JANF,JEND,KANF,KEND / 3, 3, 3, 3, 3, 3/

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

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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

      CHARACTER (LEN=16) CU,CUBE
      NAMELIST /BCUBE/  CUBE,NEIGHBOUR,IPOS,JPOS,KPOS,
     &                   IANF,IEND,JANF,JEND,KANF,KEND
CTBA1 030203 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      READ (15,BCUBE,END=100)
  100 CONTINUE
      WRITE( 6,BCUBE)
        CU    =   CUBE 
        IB    =   IPOS
        JB    =   JPOS
        KB    =   KPOS
        IBA   =   IANF
        IBE   =   IEND
        JBA   =   JANF
        JBE   =   JEND
        KBA   =   KANF
        KBE   =   KEND
CTBA2 030203 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      RETURN
      END
      SUBROUTINE READ_PPHYS(NPPHYS_MAX,NPPHYS,PPHYS,XPPHYS)
      REAL PPHYS(NPPHYS_MAX),XPPHYS(NPPHYS_MAX)

      DO I=1,NPPHYS_MAX
         PPHYS(I)=0.0
         XPPHYS(I)=FLOAT(I)
      ENDDO

      IF (NPPHYS .GT. NPPHYS_MAX) CALL ERRR (501,'READ_PPHYS')
      IF (NPPHYS .LE. 0) RETURN

      DO I=1,NPPHYS
         READ(15,*,ERR=1) XPPHYS(I),PPHYS(I)
      ENDDO

    1 CONTINUE
      RETURN
      END
