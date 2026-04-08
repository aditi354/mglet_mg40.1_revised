










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
       SUBROUTINE TST3RK
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
     $     ,RES
     $     ,GSAW,GSAE,GSAN,GSAS,GSAT,GSAB
     $     ,SIPLW,SIPLS,SIPLB,SIPLPR,SIPUE,SIPUN,SIPUT
C#if defined _PREPROC_
C     $        ,HPX
C#endif
     $   )     
C*STARLET***************************************************************
C DEFINE STATEMENTS     : ADBA, DYNSGSM, EULER, FASTBUTBIGMEM, 
C                         FRED_BODY, INC, IRIS, LEAPF, 
C                         LPART_ORIENTATION, MPI_WRITE_TIMES, MPI, 
C                         TAU_NN, TSCAL
C 
C CALLED BY             : MLET 
C
C        22.05.02 (TB)  : SCALAR TRANSPORT IMPLEMENTED 
C                         RESTRICTED TO DEFINE STATEMENT _TSCAL_ 
C        02.06.02 (TB)  : FRED_BODY BASICS IMPLEMENTED
C        03.06.02 (TB)  : MPI IMPLEMENTED
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
C

      COMMON /CLINOU/

     $                 LDIB,   LDIC,   LDIEIB, LDIEIC,
     $                 LDOB,   LDOC,   LDOEIB, LDOEIC

      SAVE   /CLINOU/
      LOGICAL

     $                 LDIB,   LDIC,   LDIEIB, LDIEIC,
     $                 LDOB,   LDOC,   LDOEIB, LDOEIC

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


C
C
C
      REAL        RDDX(IDIM1D),     RDDY(IDIM1D),     RDDZ(IDIM1D),
     $             RDX(IDIM1D),      RDY(IDIM1D),      RDZ(IDIM1D),
     $             DDX(IDIM1D),      DDY(IDIM1D),      DDZ(IDIM1D),
     $              DX(IDIM1D),       DY(IDIM1D),       DZ(IDIM1D),
     $               X(IDIM1D),        Y(IDIM1D),        Z(IDIM1D)

      REAL      CONV1S(IDIM1D)

      REAL    TAU11( 1 ), TAU12( 1 ), TAU13( 1 ),
     $        TAU21( 1 ), TAU22( 1 ), TAU23( 1 ),
     $        TAU31( 1 ), TAU32( 1 ), TAU33( 1 )

c#ifdef _PASSIVE_PARTICLE_
c#include "copart2.h"
c#endif

      REAL        GEOVP(IDIM3D)

      REAL        U (   IDIM3D  ), V (   IDIM3D  ), W (   IDIM3D  ),
     $            UO(   IDIM3D  ), VO(   IDIM3D  ), WO(   IDIM3D  ),
     $            P (   IDIM3D  ), DP(   IDIM3D  ), G (   IDIM3D  ),
     $            B (   IDIM3D  ),
     $            H3D1(   IDIM3D ),H3D2(   IDIM3D  ), H3D3(   IDIM3D  ),
     $           WCU(  IDIM3D   ),WCV(   IDIM3D  ),WCW(   IDIM3D  )

      REAL      AU  (   IDIMA   ),AV (   IDIMA   ),AW (   IDIMA   )

      REAL        UI1(IDIM2D*2),  VI1(IDIM2D*2),  WI1(IDIM2D*2),
     $            UI2(IDIM2D*2),  VI2(IDIM2D*2),  WI2(IDIM2D*2),
     $            GI1(IDIM2D),    GI2(IDIM2D),
     $            UFR(IDIM2D*2), VFR(IDIM2D*2), WFR(IDIM2D*2),
     $            PFR(IDIM2D*2), GFR(IDIM2D*2),
     $            UBACK(IDIM2D*2),  VBACK(IDIM2D*2),  WBACK(IDIM2D*2),
     $            VRI(IDIM2D*2),  
     $            UBO(IDIM2D*2),VBO(IDIM2D*2),WBO(IDIM2D*2),
     $            PBO(IDIM2D*2),GBO(IDIM2D*2),
     $            UBA(IDIM2D*2), VBA(IDIM2D*2), WBA(IDIM2D*2),
     $            BP(   IDIM3D  ),BU(   IDIM3D  ), 
     $            BV(   IDIM3D  ), BW(   IDIM3D), SDIV(   IDIM3D  )
          
C                                 FELDER FUER DEN DIREKTEN POISSONLOESER
C
      REAL
     $           FPSFAK(IDIM2D),FCOSMY(IDIM1D),FCOSNY(IDIM1D)
C
c                                 MATRIX FUER ITERATIVE VERFAHREN
      REAL GSAB( IDIM1D ), GSAT( IDIM1D ), GSAW( IDIM1D ),
     $     GSAE( IDIM1D ), GSAS( IDIM1D ), GSAN( IDIM1D )
      REAL RES  ( IDIM3D ), SIPLB ( IDIM3D ), SIPUT( IDIM3D ),
     $     SIPLW( IDIM3D ), SIPUE ( IDIM3D ), SIPLS( IDIM3D ),
     $     SIPUN( IDIM3D ), SIPLPR( IDIM3D )
C                             FELDER FUER KOEFFIZENTEN DES KOMPAKTVER.

      REAL COEFFX(1),     COEFFY(1),     COEFFZ(1),
     &     COEFDX(1),     COEFDY(1),     COEFDZ(1),
     &     LCOL  (1),     DIAG  (1),     RCOL(1),
     &     FAKTOR(1),     UZ    (1),     RSP (1),
     &     RSGS3 (1)



      REAL  FUI(1),FVI(1),FWI(1)
      REAL FUJ(1),FVJ(1),FWJ(1)
      REAL FUK(1),FVK(1),FWK(1)
C#if defined _PREPROC_
C      REAL HPX(IDIM3D)
C#endif
C
C
C                                 FELDER FUER DIE FOURIERTRANSFORMATION
C
      INTEGER IFFTX(19,MAXGRIDS),IFFTY(19,MAXGRIDS)
      INTEGER IPERMUX(IDIM1D)	     ,IPERMUY(IDIM1D)
      REAL    RFFTX(IDIM1D*2)     ,RFFTY(IDIM1D*2)
C
C
      REAL     HILF( IDIM3D ), DIV( IDIM3D ),DIVG(IDIM2D)
C
C
C                                 FUER DIE KONTROLLE DES RECHENVORGANGES
C
      REAL     EPSU(MAXGRIDS),  EPSV(MAXGRIDS),  EPSW(MAXGRIDS)
C
C
C                                 HILFSFELDER FUER EFVISC
C
      REAL        DUDY(IDIM2D),  DUDZ(IDIM2D),
     $            DVDX(IDIM2D),  DVDZ(IDIM2D),
     $            DWDX(IDIM2D),  DWDY(IDIM2D)
C

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


CTBC  ATTENTION: CHANGED ACCORDING TO SE/FS
      CHARACTER (LEN=8) CIDEND(10)
C      CHARACTER (LEN=8)      CIDEND(10)

      INTEGER          IIDEND(100)
      REAL             RIDEND(100)

CCCCCCCCCC Gewichtungsfaktoren Runge Kutta
CC    Low-Storage: FELD1,.. und UO,VO,.. werden als
CC    Hilfsfelder verwendet

      REAL             M0, M01, M1, M2, TIMERK, DTHELP
      REAL       FELD1(IDIM3D),FELD2(IDIM3D),FELD3(IDIM3D)

C
      DATA       KANPR               /  8  /

      DTHELP = DT

C                          RUNGE KUTTA ZEITSCHRITT UEBER ALLE GITTER
C
C
C---------------------------- TRANPORT OF PARTICLE-ORIENTATIONS...
C---------------------------- END OF TRANPORT OF PARTICLE-ORIENTATIONS...
C
C                                 **************************************
C                                             ZEITSCHRITT
C                                 **************************************
C

CCCCCCCCCCCCCCCCCCC ende #ifndef DNS

      DO ILEVEL = MINLEVEL,MAXLEVEL
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

        DO I = 1,NOFTST(ILEVEL)
           IGRID = IGRDOFTST(I,ILEVEL)
         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $               NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

          CIP1 = 12*3*(IP1-1)+1

C
C                                 SETZEN DER RANDBEDINGUNGEN
C
            CALL BOUNDMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $              'T',IGRID,0,0,0,TIMEPH,UBA,VBA,WBA,
     $              PBA,GBA,UBO,VBO,WBO,
     $              PFR,GFR
     $             )     

C
C
C
CTBC               PRESET ALL FIELDS BY ZERO
                   CALL DPHI0(KK,JJ,II,KK,JJ,II,UO(IP3)) 
                   CALL DPHI0(KK,JJ,II,KK,JJ,II,VO(IP3))
                   CALL DPHI0(KK,JJ,II,KK,JJ,II,WO(IP3))
                   CALL DPHI0(KK,JJ,II,KK,JJ,II,FELD1(IP3))
                   CALL DPHI0(KK,JJ,II,KK,JJ,II,FELD2(IP3))
                   CALL DPHI0(KK,JJ,II,KK,JJ,II,FELD3(IP3))
                   IDUZ = 2

                   CALL WIFAK  (ITSTEP,ITMIT,WPHI,WKON,WDIF,WSOR,IDUZ)
                   CALL TSTLE4 (KK,JJ,II,KK,JJ,II,
     $                              X(IP1),Y(IP1),Z(IP1),
     $                              DX(IP1),DY(IP1),DZ(IP1),
     $                              DDX(IP1),DDY(IP1),DDZ(IP1),
     $                              RDX(IP1),RDY(IP1),RDZ(IP1),
     $                              RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                           U(IP3), V(IP3), W(IP3),
     $                          UO(IP3),VO(IP3),WO(IP3),
     $                          WCU(IP3),WCV(IP3),WCW(IP3),
     $                          P(IP3),G(IP3),B(IP3),
     $                          BP(IP3),BU(IP3),BV(IP3),BW(IP3),
     $                          DT,
     $                          GRADPX(IGRID),ZTOT,NBUF,ITSTEP,
     $                          WPHI,WKON,WDIF,WSOR,IDUZ,
     $                          IC1(IGRID),IC2(IGRID),
     $                          JC1(IGRID),JC2(IGRID),
     $                          KC1(IGRID),KC2(IGRID),
     $                          NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,
     $                       RSGS3(IP3),
     $                       FAKTOR(IP3),
     $                       FUI(IP3),FVI(IP3),FWI(IP3),
     $                       FUJ(IP3),FVJ(IP3),FWJ(IP3),
     $                       FUK(IP3),FVK(IP3),FWK(IP3),
     $                       COEFFX(CIP1),COEFFY(CIP1),
     $                       COEFFZ(CIP1),COEFDX(CIP1),
     $                       COEFDY(CIP1),COEFDZ(CIP1),
     $                       LCOL(IP1),DIAG(IP1),RCOL(IP1),
     $                       UZ(IP1),RSP(IP1),
     $                       H3D1(IP3),H3D2(IP3),H3D3(IP3)
C#if defined _PREPROC_
C     $        ,HPX(IP3)
C#endif
     $                       )


            M0 = 0.0
            M01 = 1.0
            M1 = (1.0/3.0)*DT
CTBC        ADD M01*U,V,W + M1*UO,VO,WO TO U,V,W
            CALL PHIADD (KK,JJ,II,U(IP3),UO(IP3),U(IP3),M0,M01,M1)
            CALL PHIADD (KK,JJ,II,V(IP3),VO(IP3),V(IP3),M0,M01,M1)
            CALL PHIADD (KK,JJ,II,W(IP3),WO(IP3),W(IP3),M0,M01,M1)
            TIMERK = TIMEPH + M1

      ENDDO
      ENDDO
C                             BELEGEN DES FRONT-BUFFERS
C

C                                  LEVELS FROM FINE TO COARSE, IN
C                                  ORDER TO PROPAGATE BOUNDARY-COND.
C                                  FROM FINE TO COARSE GRIDS
      DO ILEVEL = MAXLEVEL,MINLEVEL,-1
        DO I = 1,NOFTST(ILEVEL)
           IGRID = IGRDOFTST(I,ILEVEL)

         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $               NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)



C                      SETBACKOLD BELEGT DIE FELDER UBA,VBA,WBA MIT DEN
C                      GESCHW. DER I-EBENEN IM1,IM2,IM3, DIE FUER DIE
C                      KONVEKTIVE FORMULIERUNG DER AUSSTROEMRANDBEDINGUNG
C                      BENOETIGT WERDEN
CTBC                   HERE APPLIED TO U,V,W ONLY
      CALL SETBACKOLD   (KK,JJ,II,KK,JJ,II,UBA(IBB), VBA(IBB), WBA(IBB),
     $                  U(IP3),V(IP3),W(IP3)
     $                  )     

C
C
C
C
         IF ( NFRO .EQ. 2 .OR. NFRO .EQ. 11) THEN

           IF (IVPCHILD(IGRID).EQ.0) THEN


             IF(LDIEIB) THEN
C
C                                 EINLESEN DES G-PROFILS UND DER GESCHW.
C                                 PROFILE AM EINTRITTSRAND (UNFOR-
C                                 MATIERT AUF KANAL 11)
C
                CALL DEIBI   (KMX(1),JMX(1),IMX(1),KMX(1),JMX(1),IMX(1),
     $                        ZTOT(1),YTOT(1),XTOT(1),
     $                        CIDEND,IIDEND,RIDEND,
     $                     UFR(IBB),VFR(IBB),WFR(IBB),PFR(IBB),GFR(IBB),
     $                        UI1,VI1,WI1,UI2,VI2,WI2,GI1,GI2,
     $                        UGRID,ITSTEP,ITTOT,TIMERK,DT)
             END IF
             IF(LDIEIC) THEN
C
C                                 EINLESEN DES G-PROFILS UND DER GESCHW.
C                                 PROFILE AM EINTRITTSRAND (FORMATIERT
C                                 AUF KANAL 13)
C

                CALL DEICI   (KMX(1),JMX(1),IMX(1),KMX(1),JMX(1),IMX(1),
     $                        ZTOT(1),YTOT(1),XTOT(1),
     $                        CIDEND,IIDEND,RIDEND,
     $                     UFR(IBB),VFR(IBB),WFR(IBB),PFR(IBB),GFR(IBB),
     $                        UI1,VI1,WI1,UI2,VI2,WI2,GI1,GI2,
     $                        UGRID,ITSTEP,ITTOT,TIMERK,DT)
             END IF
             IF((.NOT.LDIEIB).AND.(.NOT.LDIEIC)) THEN
               CALL SVEIPR (KK,JJ, 2,KK,JJ, 2,
     $                     UFR(IBB),VFR(IBB),WFR(IBB),PFR(IBB),
     $                     GFR(IBB),Y(IP1),Z(IP1),DZ(IP1),VCON,WCON,
     $                     YBANF(1,1,IGRID),YBEND(1,1,IGRID),
     $                     ZBANF(1,1,IGRID),ZBEND(1,1,IGRID),
     $                     FREQB(1,1,IGRID),ANIVEAU(1,1,IGRID),
     $                     AUB(1,1,IGRID),AVB(1,1,IGRID),AWB(1,1,IGRID),
     $                     TIMERK,DT,UBO,VBO,WBO)
             ENDIF

           ENDIF
        ENDIF

         IF ( NFRO .EQ. 11 .OR. NFRO .EQ. 12) THEN
            CALL MGSVFLU (U,V,W,P,G,AU,AV,AW,
     $                    HILF,UFR,VFR,WFR,PFR,GFR,X,Y,Z,
     $                    IDIM3D,IDIM2D,IDIM1D,IDIMA,
     $                    NBND,IGRID,ITSTEP,1,NFRO
     $                   )     
            
         ENDIF

C
C                             BELEGEN DES BOTTOM-BUFFERS
C                             BEI FOLGENDEN AUFRUFEN IST I MIT K VERTAUSCHT!
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
     $                 TIMERK,DT,0)
             ENDIF
           ENDIF
         ENDIF
         IF ( NBOT .EQ. 11 .OR. NBOT .EQ. 12) THEN
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
     $         NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
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
                CALL MGPOINT (IP3C,IP2C,IP1C,IBBC,IB3C,
     $                        IBUC,IPARENT(IGRID))
                
                
                CALL MGBFTC (KK,JJ, 2 ,DX(IP1),DY(IP1),DZ(IP1),
     &               DDX(IP1),DDY(IP1),DDZ(IP1),UFR(IBB),
     &               KKC,JJC, 2 ,DX(IP1C),DY(IP1C),DZ(IP1C),
     &               DDX(IP1C),DDY(IP1C),DDZ(IP1C),UFR(IBBC),
     &               KPOSITION(IGRID),JPOSITION(IGRID),
     &               IPOSITION(IGRID),'U',IPROCF,IPROCC)

             END IF
C
C                                BOTTOM BUFFER
C
             IF (            NBOT .EQ. 2 
     $                  .OR. NBOT .EQ. 11 
     $                  .OR. NBOT .EQ. 12) THEN
          CALL MGDIMS  (KKC,JJC,IIC,IPARENT(IGRID))
          CALL MGPOINT (IP3C,IP2C,IP1C,IBBC,IB3C,IBUC,IPARENT(IGRID))
                CALL MGBFTC (JJ,II, 2 ,DZ(IP1),DX(IP1),DY(IP1),
     &               DDZ(IP1),DDX(IP1),DDY(IP1),WBO(IBB),
     &               JJC,IIC, 2 ,DZ(IP1C),DX(IP1C),DY(IP1C),
     &               DDZ(IP1C),DDX(IP1C),DDY(IP1C),WBO(IBBC),
     &               JPOSITION(IGRID),IPOSITION(IGRID),
     &               KPOSITION(IGRID),'U',IPROCF,IPROCC)
             ENDIF
          ENDIF
      ENDIF
      ENDDO
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  Ende erster Runge Kutta Zwischenschritt
C                                 IMPLIZITER TEIL UEBER ALLE GITTER
C                                 DRUCKKORREKTUR
C  

c*************TEST*****************+
      DO ILEVEL = MINLEVEL,MAXLEVEL
        DO I = 1,NOFTST(ILEVEL)
           IGRID = IGRDOFTST(I,ILEVEL)
         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $               NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

         CALL COP3DZERO(KK,JJ,II,U(IP3),V(IP3),W(IP3),
     $        U(IP3),V(IP3),W(IP3),P(IP3),BP(IP3)
     $        )

            CALL BOUNDMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $              'Y',IGRID,0,0,0,TIMEPH,UBA,VBA,WBA,
     $              PBA,GBA,UBO,VBO,WBO,
     $              PFR,GFR
     $             )



      ENDDO
      ENDDO
c*********************TEST**************************



      DO ILEVEL = MINLEVEL,MAXLEVEL
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
        DO I = 1,NOFTST(ILEVEL)
           IGRID = IGRDOFTST(I,ILEVEL)
         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $               NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

C                                 SETZEN DER RANDBEDINGUNGEN
C
            CALL BOUNDMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $              'T',IGRID,0,0,0,TIMERK,UBA,VBA,WBA,
     $              PBA,GBA,UBO,VBO,WBO,
     $              PFR,GFR
     $             )     
        ENDDO
      ENDDO 

                               DT = M1
                               WSOR = 1.0
         CALL MGPOISL1    (IDIM3D,IDIM2D,IDIM1D,NBND,
     $                    X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                    RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                    U,V,W,P,DP,G,B,DIV,RES,
     $                    UFR,VFR,WFR,
     $                    UTO,VTO,WTO,PTO,GTO,VRI,HILF,
     $                    WSOR,RHO,DIVG,IPRGRID,
     $                    MINLEVEL,MAXLEVEL,MAXLEVEL,TIMERK,
     $                    UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                    PFR,GFR,
     $                    GSAW,GSAE,GSAN,GSAS,GSAT,GSAB,
     $                    BP,BU,BV,BW,SDIV
     $                    ,GEOVP,SIPLW,SIPLS,
     $                    SIPLB,SIPLPR,SIPUE,SIPUN,SIPUT
     $                    )
                              DT = DTHELP                              
C                              
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  Beginn zweiter Runge Kutta Schritt
C

CCCCCCCCCCCCCCCCCCC ende #ifndef DNS

      DO ILEVEL = MINLEVEL,MAXLEVEL
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

        DO I = 1,NOFTST(ILEVEL)
           IGRID = IGRDOFTST(I,ILEVEL)
         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $               NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
          CIP1 = 12*3*(IP1-1)+1

C
C                                 SETZEN DER RANDBEDINGUNGEN
C
            CALL BOUNDMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $              'T',IGRID,0,0,0,TIMERK,UBA,VBA,WBA,
     $              PBA,GBA,UBO,VBO,WBO,
     $              PFR,GFR
     $             )     

C
C
C
CTBC               PRESET TEMPORARILY USED FELD-FIELDS BY ZERO AGAIN
                   CALL DPHI0(KK,JJ,II,KK,JJ,II,FELD1(IP3))
                   CALL DPHI0(KK,JJ,II,KK,JJ,II,FELD2(IP3))
                   CALL DPHI0(KK,JJ,II,KK,JJ,II,FELD3(IP3))
C                   
                   IDUZ = 2
                   CALL WIFAK  (ITSTEP,ITMIT,WPHI,WKON,WDIF,WSOR,IDUZ)
                   CALL TSTLE4 (KK,JJ,II,KK,JJ,II,
     $                              X(IP1),Y(IP1),Z(IP1),
     $                              DX(IP1),DY(IP1),DZ(IP1),
     $                              DDX(IP1),DDY(IP1),DDZ(IP1),
     $                              RDX(IP1),RDY(IP1),RDZ(IP1),
     $                              RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                           U(IP3), V(IP3), W(IP3),
     $                          FELD1(IP3),FELD2(IP3),FELD3(IP3),
     $                          WCU(IP3),WCV(IP3),WCW(IP3),
     $                          P(IP3),G(IP3),B(IP3),
     $                          BP(IP3),BU(IP3),BV(IP3),BW(IP3),
     $                          DT,
     $                          GRADPX(IGRID),ZTOT,NBUF,ITSTEP,
     $                          WPHI,WKON,WDIF,WSOR,IDUZ,
     $                          IC1(IGRID),IC2(IGRID),
     $                          JC1(IGRID),JC2(IGRID),
     $                          KC1(IGRID),KC2(IGRID),
     $                          NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,
     $                       RSGS3(IP3),
     $                       FAKTOR(IP3),
     $                       FUI(IP3),FVI(IP3),FWI(IP3),
     $                       FUJ(IP3),FVJ(IP3),FWJ(IP3),
     $                       FUK(IP3),FVK(IP3),FWK(IP3),
     $                       COEFFX(CIP1),COEFFY(CIP1),
     $                       COEFFZ(CIP1),COEFDX(CIP1),
     $                       COEFDY(CIP1),COEFDZ(CIP1),
     $                       LCOL(IP1),DIAG(IP1),RCOL(IP1),
     $                       UZ(IP1),RSP(IP1),
     $                       H3D1(IP3),H3D2(IP3),H3D3(IP3)
C#if defined _PREPROC_
C     $        ,HPX(IP3)
C#endif
     $                       )
 
            M0 = 0.0
            M01 = 1.0
            M1 = (-5.0/9.0)
            M2 = (15.0/16.0)*DT
CTBC  ADD M01*FELD1,2,3 + M1*UO,VO,WO TO FELD1,2,3
      CALL PHIADD (KK,JJ,II,FELD1(IP3),UO(IP3),FELD1(IP3),
     $      M0,M01,M1)
      CALL PHIADD (KK,JJ,II,FELD2(IP3),VO(IP3),FELD2(IP3),
     $      M0,M01,M1)       
      CALL PHIADD (KK,JJ,II,FELD3(IP3),WO(IP3),FELD3(IP3),
     $      M0,M01,M1)

CTBC       ADD M01*U,V,W + M2*FELD1,2,3 TO U,V,W
           CALL PHIADD (KK,JJ,II,U(IP3),FELD1(IP3),U(IP3),M0,M01,M2)
           CALL PHIADD (KK,JJ,II,V(IP3),FELD2(IP3),V(IP3),M0,M01,M2)
           CALL PHIADD (KK,JJ,II,W(IP3),FELD3(IP3),W(IP3),M0,M01,M2)
            TIMERK = TIMERK + 0.75*DT


      ENDDO
      ENDDO
C                             BELEGEN DES FRONT-BUFFERS
C

C                                  LEVELS FROM FINE TO COARSE, IN
C                                  ORDER TO PROPAGATE BOUNDARY-COND.
C                                  FROM FINE TO COARSE GRIDS
      DO ILEVEL = MAXLEVEL,MINLEVEL,-1
        DO I = 1,NOFTST(ILEVEL)
           IGRID = IGRDOFTST(I,ILEVEL)

         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $               NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)



C                      SETBACKOLD BELEGT DIE FELDER UBA,VBA,WBA MIT DEN
C                      GESCHW. DER I-EBENEN IM1,IM2,IM3, DIE FUER DIE
C                      KONVEKTIVE FORMULIERUNG DER AUSSTROEMRANDBEDINGUNG
C                      BENOETIGT WERDEN

      CALL SETBACKOLD   (KK,JJ,II,KK,JJ,II,UBA(IBB), VBA(IBB), WBA(IBB),
     $                  U(IP3),V(IP3),W(IP3)
     $                  )     

C
C
C
C
         IF ( NFRO .EQ. 2 .OR. NFRO .EQ. 11) THEN

           IF (IVPCHILD(IGRID).EQ.0) THEN


             IF(LDIEIB) THEN
C
C                                 EINLESEN DES G-PROFILS UND DER GESCHW.
C                                 PROFILE AM EINTRITTSRAND (UNFOR-
C                                 MATIERT AUF KANAL 11)
C
                CALL DEIBI   (KMX(1),JMX(1),IMX(1),KMX(1),JMX(1),IMX(1),
     $                        ZTOT(1),YTOT(1),XTOT(1),
     $                        CIDEND,IIDEND,RIDEND,
     $                     UFR(IBB),VFR(IBB),WFR(IBB),PFR(IBB),GFR(IBB),
     $                        UI1,VI1,WI1,UI2,VI2,WI2,GI1,GI2,
     $                        UGRID,ITSTEP,ITTOT,TIMERK,DT)
             END IF
             IF(LDIEIC) THEN
C
C                                 EINLESEN DES G-PROFILS UND DER GESCHW.
C                                 PROFILE AM EINTRITTSRAND (FORMATIERT
C                                 AUF KANAL 13)
C

                CALL DEICI   (KMX(1),JMX(1),IMX(1),KMX(1),JMX(1),IMX(1),
     $                        ZTOT(1),YTOT(1),XTOT(1),
     $                        CIDEND,IIDEND,RIDEND,
     $                     UFR(IBB),VFR(IBB),WFR(IBB),PFR(IBB),GFR(IBB),
     $                        UI1,VI1,WI1,UI2,VI2,WI2,GI1,GI2,
     $                        UGRID,ITSTEP,ITTOT,TIMERK,DT)
             END IF
             IF((.NOT.LDIEIB).AND.(.NOT.LDIEIC)) THEN
               CALL SVEIPR (KK,JJ, 2,KK,JJ, 2,
     $                     UFR(IBB),VFR(IBB),WFR(IBB),PFR(IBB),
     $                     GFR(IBB),Y(IP1),Z(IP1),DZ(IP1),VCON,WCON,
     $                     YBANF(1,1,IGRID),YBEND(1,1,IGRID),
     $                     ZBANF(1,1,IGRID),ZBEND(1,1,IGRID),
     $                     FREQB(1,1,IGRID),ANIVEAU(1,1,IGRID),
     $                     AUB(1,1,IGRID),AVB(1,1,IGRID),AWB(1,1,IGRID),
     $                     TIMERK,DT,UBO,VBO,WBO)
             ENDIF

           ENDIF
        ENDIF

         IF ( NFRO .EQ. 11 .OR. NFRO .EQ. 12) THEN
            CALL MGSVFLU (U,V,W,P,G,AU,AV,AW,
     $                    HILF,UFR,VFR,WFR,PFR,GFR,X,Y,Z,
     $                    IDIM3D,IDIM2D,IDIM1D,IDIMA,
     $                    NBND,IGRID,ITSTEP,1,NFRO
     $                   )     
            
         ENDIF

C
C                             BELEGEN DES BOTTOM-BUFFERS
C                             BEI FOLGENDEN AUFRUFEN IST I MIT K VERTAUSCHT!
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
     $                 TIMERK,DT,0)
             ENDIF
           ENDIF
         ENDIF

         IF ( NBOT .EQ. 11 .OR. NBOT .EQ. 12) THEN
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
     $         NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
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
                CALL MGPOINT (IP3C,IP2C,IP1C,IBBC,IB3C,
     $                        IBUC,IPARENT(IGRID))
                
                
                CALL MGBFTC (KK,JJ, 2 ,DX(IP1),DY(IP1),DZ(IP1),
     &               DDX(IP1),DDY(IP1),DDZ(IP1),UFR(IBB),
     &               KKC,JJC, 2 ,DX(IP1C),DY(IP1C),DZ(IP1C),
     &               DDX(IP1C),DDY(IP1C),DDZ(IP1C),UFR(IBBC),
     &               KPOSITION(IGRID),JPOSITION(IGRID),
     &               IPOSITION(IGRID),'U',IPROCF,IPROCC)

             ENDIF
C
C                                BOTTOM BUFFER
C
             IF (            NBOT .EQ. 2 
     $                  .OR. NBOT .EQ. 11 
     $                  .OR. NBOT .EQ. 12) THEN
          CALL MGDIMS  (KKC,JJC,IIC,IPARENT(IGRID))
          CALL MGPOINT (IP3C,IP2C,IP1C,IBBC,IB3C,IBUC,IPARENT(IGRID))
                CALL MGBFTC (JJ,II, 2 ,DZ(IP1),DX(IP1),DY(IP1),
     &               DDZ(IP1),DDX(IP1),DDY(IP1),WBO(IBB),
     &               JJC,IIC, 2 ,DZ(IP1C),DX(IP1C),DY(IP1C),
     &               DDZ(IP1C),DDX(IP1C),DDY(IP1C),WBO(IBBC),
     &               JPOSITION(IGRID),IPOSITION(IGRID),
     &               KPOSITION(IGRID),'U',IPROCF,IPROCC)
             ENDIF
          ENDIF
      ENDIF
      ENDDO
      ENDDO
CCCCCCCCCCCCCCCCC  Ende zweiter Zwischenschritt fuer Runge Kutta

c*************TEST*****************+
      DO ILEVEL = MINLEVEL,MAXLEVEL
        DO I = 1,NOFTST(ILEVEL)
           IGRID = IGRDOFTST(I,ILEVEL)
         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $               NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

         CALL COP3DZERO(KK,JJ,II,U(IP3),V(IP3),W(IP3),
     $        U(IP3),V(IP3),W(IP3),P(IP3),BP(IP3)
     $        )

            CALL BOUNDMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $              'Y',IGRID,0,0,0,TIMEPH,UBA,VBA,WBA,
     $              PBA,GBA,UBO,VBO,WBO,
     $              PFR,GFR
     $             )



      ENDDO
      ENDDO
c*********************TEST**************************
 

      DO ILEVEL = MINLEVEL,MAXLEVEL
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


        DO I = 1,NOFTST(ILEVEL)
           IGRID = IGRDOFTST(I,ILEVEL)
         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $               NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

C                                 SETZEN DER RANDBEDINGUNGEN
C
            CALL BOUNDMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $              'T',IGRID,0,0,0,TIMERK,UBA,VBA,WBA,
     $              PBA,GBA,UBO,VBO,WBO,
     $              PFR,GFR
     $             ) 


        ENDDO
      ENDDO

                               DT = DT*0.75
                               WSOR = 1.0


         CALL MGPOISL1    (IDIM3D,IDIM2D,IDIM1D,NBND,
     $                    X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                    RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                    U,V,W,P,DP,G,B,DIV,RES,
     $                    UFR,VFR,WFR,
     $                    UTO,VTO,WTO,PTO,GTO,VRI,HILF,
     $                    WSOR,RHO,DIVG,IPRGRID,
     $                    MINLEVEL,MAXLEVEL,MAXLEVEL,TIMERK,
     $                    UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                    PFR,GFR,
     $                    GSAW,GSAE,GSAN,GSAS,GSAT,GSAB,
     $                    BP,BU,BV,BW,SDIV
     $                    ,GEOVP,SIPLW,SIPLS,
     $                    SIPLB,SIPLPR,SIPUE,SIPUN,SIPUT
     $                    )
                     DT = DTHELP
C                     
     
C
C
C               REDUZIERUNG DES DRUCKNIVEAUS FUER GEWUENSCHTE GITTER    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  Beginn letzter Runge Kutta Schritt
C

CCCCCCCCCCCCCCCCCCC ende #ifndef DNS
C
      DO ILEVEL = MINLEVEL,MAXLEVEL
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

        DO I = 1,NOFTST(ILEVEL)
           IGRID = IGRDOFTST(I,ILEVEL)

         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $               NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
          CIP1 = 12*3*(IP1-1)+1

C
C                                 SETZEN DER RANDBEDINGUNGEN
C
            CALL BOUNDMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $              'T',IGRID,0,0,0,TIMERK,UBA,VBA,WBA,
     $              PBA,GBA,UBO,VBO,WBO,
     $              PFR,GFR
     $             )     

C
C
C
                CALL DPHI0(KK,JJ,II,KK,JJ,II,UO(IP3))
                CALL DPHI0(KK,JJ,II,KK,JJ,II,VO(IP3))
                CALL DPHI0(KK,JJ,II,KK,JJ,II,WO(IP3))              

                   IDUZ = 2
                   CALL WIFAK  (ITSTEP,ITMIT,WPHI,WKON,WDIF,WSOR,IDUZ)
                   CALL TSTLE4 (KK,JJ,II,KK,JJ,II,
     $                              X(IP1),Y(IP1),Z(IP1),
     $                              DX(IP1),DY(IP1),DZ(IP1),
     $                              DDX(IP1),DDY(IP1),DDZ(IP1),
     $                              RDX(IP1),RDY(IP1),RDZ(IP1),
     $                              RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                           U(IP3), V(IP3), W(IP3),
     $                          UO(IP3),VO(IP3),WO(IP3),
     $                          WCU(IP3),WCV(IP3),WCW(IP3),
     $                          P(IP3),G(IP3),B(IP3),
     $                          BP(IP3),BU(IP3),BV(IP3),BW(IP3),
     $                          DT,
     $                          GRADPX(IGRID),ZTOT,NBUF,ITSTEP,
     $                          WPHI,WKON,WDIF,WSOR,IDUZ,
     $                          IC1(IGRID),IC2(IGRID),
     $                          JC1(IGRID),JC2(IGRID),
     $                          KC1(IGRID),KC2(IGRID),
     $                          NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,
     $                       RSGS3(IP3),
     $                       FAKTOR(IP3),
     $                       FUI(IP3),FVI(IP3),FWI(IP3),
     $                       FUJ(IP3),FVJ(IP3),FWJ(IP3),
     $                       FUK(IP3),FVK(IP3),FWK(IP3),
     $                       COEFFX(CIP1),COEFFY(CIP1),
     $                       COEFFZ(CIP1),COEFDX(CIP1),
     $                       COEFDY(CIP1),COEFDZ(CIP1),
     $                       LCOL(IP1),DIAG(IP1),RCOL(IP1),
     $                       UZ(IP1),RSP(IP1),
     $                       H3D1(IP3),H3D2(IP3),H3D3(IP3)
C#if defined _PREPROC_
C     $        ,HPX(IP3)
C#endif
     $                       )
C
 
             M0 = 0.0
             M01 = 1.0
             M1 = (-153.0/128.0)
             M2 = (8.0/15.0)*DT

CTBC       ADD M01*UO,VO,WO + M1*FELD1,2,3 TO UO,VO,WO 
           CALL PHIADD (KK,JJ,II,UO(IP3),FELD1(IP3),UO(IP3),
     $      M0,M01,M1)
           CALL PHIADD (KK,JJ,II,VO(IP3),FELD2(IP3),VO(IP3),
     $      M0,M01,M1)       
           CALL PHIADD (KK,JJ,II,WO(IP3),FELD3(IP3),WO(IP3),
     $      M0,M01,M1)
CTBC         ADD M01*U,V,W + M2*UO,VO,WO TO U,V,W
             CALL PHIADD (KK,JJ,II,U(IP3),UO(IP3),U(IP3),M0,M01,M2)
             CALL PHIADD (KK,JJ,II,V(IP3),VO(IP3),V(IP3),M0,M01,M2)
             CALL PHIADD (KK,JJ,II,W(IP3),WO(IP3),W(IP3),M0,M01,M2)



      ENDDO
      ENDDO
C                             BELEGEN DES FRONT-BUFFERS
C
              TIMEPH     = TIMEPH + DT
C                                  LEVELS FROM FINE TO COARSE, IN
C                                  ORDER TO PROPAGATE BOUNDARY-COND.
C                                  FROM FINE TO COARSE GRIDS
      DO ILEVEL = MAXLEVEL,MINLEVEL,-1
        DO I = 1,NOFTST(ILEVEL)
           IGRID = IGRDOFTST(I,ILEVEL)

         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $               NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)



C                      SETBACKOLD BELEGT DIE FELDER UBA,VBA,WBA MIT DEN
C                      GESCHW. DER I-EBENEN IM1,IM2,IM3, DIE FUER DIE
C                      KONVEKTIVE FORMULIERUNG DER AUSSTROEMRANDBEDINGUNG
C                      BENOETIGT WERDEN

      CALL SETBACKOLD   (KK,JJ,II,KK,JJ,II,UBA(IBB), VBA(IBB), WBA(IBB),
     $                  U(IP3),V(IP3),W(IP3)
     $                  )     

C
C
C
C
         IF ( NFRO .EQ. 2 .OR. NFRO .EQ. 11) THEN

           IF (IVPCHILD(IGRID).EQ.0) THEN


             IF(LDIEIB) THEN
C
C                                 EINLESEN DES G-PROFILS UND DER GESCHW.
C                                 PROFILE AM EINTRITTSRAND (UNFOR-
C                                 MATIERT AUF KANAL 11)
C
                CALL DEIBI   (KMX(1),JMX(1),IMX(1),KMX(1),JMX(1),IMX(1),
     $                        ZTOT(1),YTOT(1),XTOT(1),
     $                        CIDEND,IIDEND,RIDEND,
     $                     UFR(IBB),VFR(IBB),WFR(IBB),PFR(IBB),GFR(IBB),
     $                        UI1,VI1,WI1,UI2,VI2,WI2,GI1,GI2,
     $                        UGRID,ITSTEP,ITTOT,TIMEPH,DT)
             END IF
             IF(LDIEIC) THEN
C
C                                 EINLESEN DES G-PROFILS UND DER GESCHW.
C                                 PROFILE AM EINTRITTSRAND (FORMATIERT
C                                 AUF KANAL 13)
C

                CALL DEICI   (KMX(1),JMX(1),IMX(1),KMX(1),JMX(1),IMX(1),
     $                        ZTOT(1),YTOT(1),XTOT(1),
     $                        CIDEND,IIDEND,RIDEND,
     $                     UFR(IBB),VFR(IBB),WFR(IBB),PFR(IBB),GFR(IBB),
     $                        UI1,VI1,WI1,UI2,VI2,WI2,GI1,GI2,
     $                        UGRID,ITSTEP,ITTOT,TIMEPH,DT)
             END IF
             IF((.NOT.LDIEIB).AND.(.NOT.LDIEIC)) THEN
               CALL SVEIPR (KK,JJ, 2,KK,JJ, 2,
     $                     UFR(IBB),VFR(IBB),WFR(IBB),PFR(IBB),
     $                     GFR(IBB),Y(IP1),Z(IP1),DZ(IP1),VCON,WCON,
     $                     YBANF(1,1,IGRID),YBEND(1,1,IGRID),
     $                     ZBANF(1,1,IGRID),ZBEND(1,1,IGRID),
     $                     FREQB(1,1,IGRID),ANIVEAU(1,1,IGRID),
     $                     AUB(1,1,IGRID),AVB(1,1,IGRID),AWB(1,1,IGRID),
     $                     TIMEPH,DT,UBO,VBO,WBO)
             ENDIF

           ENDIF
           ENDIF

         IF ( NFRO .EQ. 11 .OR. NFRO .EQ. 12) THEN
            CALL MGSVFLU (U,V,W,P,G,AU,AV,AW,
     $                    HILF,UFR,VFR,WFR,PFR,GFR,X,Y,Z,
     $                    IDIM3D,IDIM2D,IDIM1D,IDIMA,
     $                    NBND,IGRID,ITSTEP,1,NFRO
     $                   )     
            
         ENDIF
         
C
C                             BELEGEN DES BOTTOM-BUFFERS
C                             BEI FOLGENDEN AUFRUFEN IST I MIT K VERTAUSCHT!
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
         IF ( NBOT .EQ. 11 .OR. NBOT .EQ. 12) THEN
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
     $         NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
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
                CALL MGPOINT (IP3C,IP2C,IP1C,IBBC,IB3C,
     $                        IBUC,IPARENT(IGRID))
                
                
                CALL MGBFTC (KK,JJ, 2 ,DX(IP1),DY(IP1),DZ(IP1),
     &               DDX(IP1),DDY(IP1),DDZ(IP1),UFR(IBB),
     &               KKC,JJC, 2 ,DX(IP1C),DY(IP1C),DZ(IP1C),
     &               DDX(IP1C),DDY(IP1C),DDZ(IP1C),UFR(IBBC),
     &               KPOSITION(IGRID),JPOSITION(IGRID),
     &               IPOSITION(IGRID),'U',IPROCF,IPROCC)

             END IF
C
C                                BOTTOM BUFFER
C
             IF (            NBOT .EQ. 2 
     $                  .OR. NBOT .EQ. 11 
     $                  .OR. NBOT .EQ. 12) THEN
          CALL MGDIMS  (KKC,JJC,IIC,IPARENT(IGRID))
          CALL MGPOINT (IP3C,IP2C,IP1C,IBBC,IB3C,IBUC,IPARENT(IGRID))
                CALL MGBFTC (JJ,II, 2 ,DZ(IP1),DX(IP1),DY(IP1),
     &               DDZ(IP1),DDX(IP1),DDY(IP1),WBO(IBB),
     &               JJC,IIC, 2 ,DZ(IP1C),DX(IP1C),DY(IP1C),
     &               DDZ(IP1C),DDX(IP1C),DDY(IP1C),WBO(IBBC),
     &               JPOSITION(IGRID),IPOSITION(IGRID),
     &               KPOSITION(IGRID),'U',IPROCF,IPROCC)
             END IF
          END IF
      ENDIF
      ENDDO
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  Ende der Zeitintegration

c*************TEST*****************+
      DO ILEVEL = MINLEVEL,MAXLEVEL
        DO I = 1,NOFTST(ILEVEL)
           IGRID = IGRDOFTST(I,ILEVEL)
         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $               NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

         CALL COP3DZERO(KK,JJ,II,U(IP3),V(IP3),W(IP3),
     $        U(IP3),V(IP3),W(IP3),P(IP3),BP(IP3)
     $        )

            CALL BOUNDMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $              'Y',IGRID,0,0,0,TIMEPH,UBA,VBA,WBA,
     $              PBA,GBA,UBO,VBO,WBO,
     $              PFR,GFR
     $             )



      ENDDO
      ENDDO
c*********************TEST**************************


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  TESTWEISE HIER HOCHTICKERN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  DES ZEITSCHRITT-ZAEHLERS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  FUER KORREKTE ANSTEUERUNG
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  DER RANDBEDINGUNGEN

                   ITSTEP     = ITSTEP + 1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                           IF BODY IS IN THE COMPUTATIONAL DOMAIN
C                           THE BOUNDARY CONDITIONS FOR IT ARE NOT SET
C                           BEFORE THE FIRST PRESSURE CORRECTION.
C                           BCUBMG IS NOT INCLUDED IF BOUNDMG
C                           IS CALLED WITH ITYP 'P' IN CASE OF 
C                           DEFINE DIRECTIVE '_AVOIDLOADIMBALANCE_'
C                           SO BOUNDMG HAS TO BE CALLED WITH ITYP 'T'.
C                           NOCH EINMAL WERDEN RANDBEDINGUNGEN FUER 
C                           GITTER GESETZT
C
C                             BOUND und BPARMG ZWEIMAL
C      DO ICOUNT =1,2
C
C                               FILL THE BOUNDARY BUFFER FOR LOCAL GRIDS
C
      DO ILEVEL = MINLEVEL,MAXLEVEL
         CALL BPARMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,DUDY,DUDZ,DVDX,
     $              UTO,VTO,WTO,PTO,GTO,UBA,VBA,WBA,PBA,GBA,
     $              VFR,WFR,PFR,GFR,
     $              'T',ILEVEL,0,0,0,
     $              NOFTST(ILEVEL),IGRDOFTST(1,ILEVEL)
     $             )     
      ENDDO



      DO ILEVEL = MINLEVEL,MAXLEVEL
C
         CALL CONNECTMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,UTO,VTO,WTO,PTO,GTO,
     $              VFR,WFR,PFR,GFR,
     $              'T',ILEVEL,0,0,0,
     $              NOFTST(ILEVEL),IGRDOFTST(1,ILEVEL))

        DO I = 1,NOFTST(ILEVEL)
           IGRID = IGRDOFTST(I,ILEVEL)

            CALL MGPOINT (IP3,IP2,IP1,IBB,IB3,IBU,IGRID)

            CALL BOUNDMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $              'T',IGRID,0,0,0,TIMEPH,UBA,VBA,WBA,
     $              PBA,GBA,UBO,VBO,WBO,
     $              PFR,GFR
     $             )    

 
         ENDDO
         ENDDO
C                             BOUND und BPARMG ZWEIMAL
C      ENDDO
C
C                                 IMPLIZITER TEIL UEBER ALLE GITTER
C                                 DRUCKKORREKTUR
C
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
     $                    ,GEOVP,SIPLW,SIPLS,
     $                    SIPLB,SIPLPR,SIPUE,SIPUN,SIPUT
     $                    )
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
     $                   DDX(IP1), DDY(IP1), DDZ(IP1),NFRO,NRGT,NBOT)
        ENDIF
      ENDDO
      ENDDO
C
C
C                           NOCH EINMAL WERDEN RANDBED. FUER ALLE 
C                           GITTER GESETZT
C

      DO ILEVEL = MINLEVEL,MAXLEVEL

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

            CALL MGPOINT (IP3,IP2,IP1,IBB,IB3,IBU,IGRID)

            CALL BOUNDMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $              'T',IGRID,0,0,0,TIMEPH,UBA,VBA,WBA,
     $              PBA,GBA,UBO,VBO,WBO,
     $              PFR,GFR
     $             )     

C      CALL WRITE3DDIAGX (KMX,JMX,IMX,T(IP3),20,IMX-1)
C      CALL WRITE3DDIAGX (KMX,JMX,IMX,T(IP3),21,IMX-2)
C      CALL WRITE3DDIAGX (KMX,JMX,IMX,T(IP3),22,IMX-3)
C      CALL WRITE3DDIAGX (KMX,JMX,IMX,T(IP3),30,2)
C      CALL WRITE3DDIAGX (KMX,JMX,IMX,T(IP3),31,3)
C      CALL WRITE3DDIAGX (KMX,JMX,IMX,T(IP3),32,4)
C      STOP 'TST3RK'

         ENDDO
         ENDDO


CCC TESTWEISE vor Aufruf TST1G in MLET GESCHALTET  TIMEPH     = TIMEPH + DT
CCC TESTWEISE vor Aufruf von MGSVFLU                   ITSTEP     = ITSTEP + 1
C                   TIMEPH     = TIMEPH + DT
                   ITTOT      = ITTOT  + 1
C
C                                 **************************************
C                                 VOLLSTAENDIGER ZEITSCHRITT ABGESCHLOSS
C                                 **************************************
                                
      RETURN
      END
