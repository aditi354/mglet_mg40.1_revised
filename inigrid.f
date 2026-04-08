










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
         SUBROUTINE INIGRID
     $             (CIDENT,IIDENT,RIDENT,
     $              DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,
     $              U,V,W,P,G,B,BP,BU,BV,BW,
     $              ITSTEP,ITTOT,TIMEPH,
     $              IDIMF,NGRP1, NGRP3, IGRHF, GRX0, GRX1, HILF,
     $              IDIM3D,IDIM2D,IDIM1D,NBUF,NBND,IDIMA,IDIM1L,IDIM2L,
     $              IGRID,CONV1S)
C--MGLET----------------------------------------------------------------
C
C                   SETZEN UND VORBELEGEN EINES GITTERS
C
C
C        17. 6.93 (MM)  : DUMMYGITTER EINGEFUEHRT
C        24. 6.93 (MM)  :     ""       WIEDER AUSGEFUEHRT
C        22. 5.95 (MM.AO): MPI
C        11.09.1996 (A.O.): SUBGRIDSCALE CONSTANT DEPENDENT OF X
C        11.03.03 (TB): SCALAR FIELD T AND BOUNDARY FIELD BT IMPLEMENTED
C
C--MGLET----------------------------------------------------------------



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
      CHARACTER (LEN=8)    CIDENT(10)
C
      INTEGER IIDENT(100),
     $        IGRHF(IDIMF,NGRP1)
C
C
C
C
      REAL    RIDENT(100),
     $            X(IDIM1D),     Y(IDIM1D),     Z(IDIM1D),
     $           DX(IDIM1D),    DY(IDIM1D),    DZ(IDIM1D),
     $          DDX(IDIM1D),   DDY(IDIM1D),   DDZ(IDIM1D),
     $      U( IDIM3D ), V( IDIM3D ), W( IDIM3D ),
     $      P( IDIM3D ), G( IDIM3D ), B( IDIM3D ), HILF( IDIM3D ),
     $      BP( IDIM3D ),BU( IDIM3D ),BV( IDIM3D ),BW( IDIM3D )
C
C          SUBGRIDSCALE CONSTANT CONV1S DEPENDENT OF X
      REAL    CONV1S(IDIM1D)
C

C
C
C
C                                  EINSTROEMBUFFER EXISTIEREN NUR
C                                  EINMAL
C
         IPI =   1 
C
C

C
C                                 GITTERGENERIERUNG UND VORBELEGUNG VON
C                                 U,V,W,(T),P UND G.
C
      WRITE (6,*) 'GRID GENERATING:',IGRID

         CALL MGGRDGEN (DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,
     $                IDIMF,NGRP1, NGRP3, IGRHF, GRX0, GRX1, HILF,
     $                IDIM3D,IDIM2D,IDIM1D,NBND,
     $                IGRID,0)

C
C                        VERSENDEN DER GEOMETRIEINFORMATION AN DIE
C                        SUB-GITTER DER GEBIETSZERLEGUNG
C
      IF (LSLICE(IGRID)) THEN
        DO ISUB = 1,NOFSLCHILDS(IGRID)
          I = IGRDOFSLCHILD(ISUB,IGRID)
C
         CALL SLICEGEO (DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,B,HILF,
     $              IDIM3D,IDIM2D,IDIM1D,NBUF,NBND,IDIMA,IDIM1L,IDIM2L,
     $              I)


C                       SUCHEN DER PARENTGITTER
C
      CALL PARSET (DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,HILF,
     $               IDIM1D,IDIMF,NBND,I)

        ENDDO
      ENDIF

C
C
C                                 BELEGEN DES WANDERKENNUNGSFELDES
C                                 IF ONE GRID IS A CHILD OF
C                                 ANOTHER, THE BODY-FIELD WILL
C                                 BE PROLONGATED FROM THE
C                                 PARENT- (COARSE) GRID IN MGCTOF
C
C050696      IF ( .NOT. LCHILD(IGRID)) THEN


C050696       IF (LSLICE(IGRID)) THEN
C050696        IF (MYID .EQ. 0) THEN

C050696      WRITE (6,*) 'INIGRID, calling ibfield fort grid',igrid
C050696         CALL IBFIELD (KK,JJ,II,KK,JJ,II,NBND,
C050696     $              DDX(IP1),DDY(IP1),DDZ(IP1),
C050696     $              DX(IP1),DY(IP1),DZ(IP1),
C050696     $              X(IP1),Y(IP1),Z(IP1),
C050696     $              BIGBUF(IP3),IDIMF,NGRP3,HILF,ZCUB,NCUB,IGRID)

C050696         ENDIF
C050696       ELSE
C050696        IF (MYID .EQ. IDPROCOFGRD(IGRID)) THEN


C050696         CALL IBFIELD (KK,JJ,II,KK,JJ,II,NBND,
C050696     $              DDX(IP1),DDY(IP1),DDZ(IP1),
C050696     $              DX(IP1),DY(IP1),DZ(IP1),
C050696     $              X(IP1),Y(IP1),Z(IP1),
C050696     $              B(IP3),IDIMF,NGRP3,HILF,ZCUB,NCUB,IGRID)


C050696        ENDIF
C050696       ENDIF

C
C050696       IF (LSLICE(IGRID)) THEN
 
C                        VERSENDEN DES WANDERKENNUNGSFELDES
C
C050696       CALL SLICEFIELD (B,HILF,IDIM3D,IGRID)
C
C050696       ENDIF

C050696      ENDIF
C
C
      DO ISUB = 1,MAX(1,NOFSLCHILDS(IGRID))

        IF (LSLICE(IGRID)) THEN
         I = IGRDOFSLCHILD(ISUB,IGRID)
        ELSE
         I = IGRID
        ENDIF

             MYID   = 0
             IPROC  = 0
             IPPAR  = 0

C                                 BEREITSTELLUNG DER DIMENSIONIERUNG
C                                 BASIS-RANDBEDINGUNGEN
C
       CALL MGDIMS  (KK,JJ,II,I)
       CALL MGBASB  (NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,I)

C                                          POINTER
       CALL MGPOINT (IP3,IP2,IP1,IBB,IB3,IBU,I)
C
C
C


C
C                               VORBELEGUNG DER VARIABLEN
C
       IF (LCHILD(I)) THEN


C                               VARIABLE WERDEN AUS PARENT-GITTER
C                               INTERPOLIERT ("PROLONGIERT")
C                                          POINTER
C

       CALL MGDIMS  (KKC,JJC,IIC,IPARENT(I))
       CALL MGBASB  (NFROC,NBACC,NRGTC,NLFTC,NBOTC,NTOPC,NCUBC,
     $               IPARENT(I))
       IF(MYID .EQ. IPPAR) CALL MGPOINT (IP3C,IP2C,IP1C,IBBC,IB3C,
     $               IBUC,IPARENT(I))

              CALL MGCTOF (KK,JJ,II,DX(IP1),DY(IP1),DZ(IP1),
     &                     DDX(IP1),DDY(IP1),DDZ(IP1),U(IP3),
     &                     KKC,JJC,IIC,DX(IP1C),DY(IP1C),DZ(IP1C),
     &                     DDX(IP1C),DDY(IP1C),DDZ(IP1C),U(IP3C),
     &               IPOSITION(I),JPOSITION(I),KPOSITION(I),
     &                     'X',IPROC,IPPAR)

              CALL MGCTOF (KK,JJ,II,DX(IP1),DY(IP1),DZ(IP1),
     &                     DDX(IP1),DDY(IP1),DDZ(IP1),V(IP3),
     &                     KKC,JJC,IIC,DX(IP1C),DY(IP1C),DZ(IP1C),
     &                     DDX(IP1C),DDY(IP1C),DDZ(IP1C),V(IP3C),
     &               IPOSITION(I),JPOSITION(I),KPOSITION(I),
     &                     'Y',IPROC,IPPAR)

              CALL MGCTOF (KK,JJ,II,DX(IP1),DY(IP1),DZ(IP1),
     &                     DDX(IP1),DDY(IP1),DDZ(IP1),W(IP3),
     &                     KKC,JJC,IIC,DX(IP1C),DY(IP1C),DZ(IP1C),
     &                     DDX(IP1C),DDY(IP1C),DDZ(IP1C),W(IP3C),
     &               IPOSITION(I),JPOSITION(I),KPOSITION(I),
     &                     'Z',IPROC,IPPAR)
              CALL MGCTOF (KK,JJ,II,DX(IP1),DY(IP1),DZ(IP1),
     &                     DDX(IP1),DDY(IP1),DDZ(IP1),P(IP3),
     &                     KKC,JJC,IIC,DX(IP1C),DY(IP1C),DZ(IP1C),
     &                     DDX(IP1C),DDY(IP1C),DDZ(IP1C),P(IP3C),
     &               IPOSITION(I),JPOSITION(I),KPOSITION(I),
     &                     'P',IPROC,IPPAR)

              CALL MGCTOF (KK,JJ,II,DX(IP1),DY(IP1),DZ(IP1),
     &                     DDX(IP1),DDY(IP1),DDZ(IP1),G(IP3),
     &                     KKC,JJC,IIC,DX(IP1C),DY(IP1C),DZ(IP1C),
     &                     DDX(IP1C),DDY(IP1C),DDZ(IP1C),G(IP3C),
     &               IPOSITION(I),JPOSITION(I),KPOSITION(I),
     &                     'B',IPROC,IPPAR)

              CALL MGCTOF (KK,JJ,II,DX(IP1),DY(IP1),DZ(IP1),
     &                     DDX(IP1),DDY(IP1),DDZ(IP1),B(IP3),
     &                     KKC,JJC,IIC,DX(IP1C),DY(IP1C),DZ(IP1C),
     &                     DDX(IP1C),DDY(IP1C),DDZ(IP1C),B(IP3C),
     &               IPOSITION(I),JPOSITION(I),KPOSITION(I),
     &                     'B',IPROC,IPPAR)

       ELSE

              CALL SVLE1  (KK,JJ,II,KK,JJ,II,
     $                     DX(IP1),DY(IP1),DZ(IP1),
     &                     U(IP3),V(IP3),W(IP3),P(IP3),
     &                     X(IP1),Y(IP1),Z(IP1),0.0,0.0,
     &                     YBANF,YBEND,ZBANF,ZBEND,HILF,GRADPX(IGRID))

C
              CALL SETS   (KK,JJ,II,KK,JJ,II,P(IP3),0.0)
              CALL SETS   (KK,JJ,II,KK,JJ,II,G(IP3),GMOL)

         CALL IBFIELD (KK,JJ,II,KK,JJ,II,NBND,
     $              DDX(IP1),DDY(IP1),DDZ(IP1),
     $              DX(IP1),DY(IP1),DZ(IP1),
     $              X(IP1),Y(IP1),Z(IP1),
     $              B(IP3),BP(IP3),BU(IP3),BV(IP3),BW(IP3),
     $              IDIMF,NGRP3,HILF,ZCUB,NCUB,I
     $                 )     




       ENDIF
C
C
C                                 BELEGUNG DER GRENZEN DER BOUNDING-BOX
C
       CALL SBBC12 (KK,JJ,II,KK,JJ,II,BP(IP3),
     $              KC1(I),KC2(I),JC1(I),JC2(I),IC1(I),IC2(I),
     $                  NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,I)
C
C
C                                 BELEGUNG VON B MIT LAENGENMASS
C
       CALL DMIXLE (KK,JJ,II,KK,JJ,II,B(IP3),HILF(IP3),
     $              X(IP1),Y(IP1),Z(IP1),
     $              DX(IP1),DY(IP1),DZ(IP1),
     $              DDX(IP1),DDY(IP1),DDZ(IP1),
     $              IC1(I),IC2(I),JC1(I),JC2(I),KC1(I),KC2(I),
     $              NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,
     $              GRADPX(I),RHO,GMOL,CONV1S(IP1),
     $              CONV1SANF(I),CONV1SEND(I),TRANSLES1(I),TRANSLES2(I),
     $              U(IP3))
C
C
       ENDDO

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      RETURN
      END


