










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
         SUBROUTINE MGVPC4 (IDIM3D,IDIM2D,IDIM1D,NBND,
     $                      X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                      RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                    FPSFAK,FCOSMY,FCOSNY,IFFTX,RFFTX,IFFTY,RFFTY,    
     $                      U,V,W,P,DP,G,B,BP,BU,BV,BW,SDIV,
     $                      DIV,UFR,VFR,WFR,
     $                      UTO,VTO,WTO,PTO,GTO,VRI,HILF,
     $                      H2D1,H2D2,H2D3,H2D4,H2D5,H2D6,
     $                      OMBETA,WSOR, RHO,DIVG,IPRGRID,
     $                      MINCALLVL,ILEVEL,MAXCALLVL,TIMEPH,
     $                      UBA,VBA,WBA,PBA,GBA,GEOVP,UBO,VBO,WBO,
     $                      PFR,GFR,IALGO
     $                     )      

C**MGLET****************************************************************
C
C   M G V P C 1      DRUCKKORREKTURZYKLUS FUER MULTIGRID-VARIANTE
C                    HIER WERDEN GITTER DES LEVELS ILEVEL KORRIGIERT
C                    ERSTE GROBGITTERSTUFE  DES DRUCKKORREKTURZYKLUS
C
C     6.10.93 (MM):  ORIGINAL
C    27. 2.94 (MM):  MGVPIT EINGEFUEHRT
C     7.12.95 (MM):  ITERATIONS ARE PERFORMED ON DP, NEW PRESSURE
C                    IS CALCULATED BY PHIFLU AT THE END
C    18.12.95 (MM):  DELTA_P IS NOT CALCULATED. IT WILL BE CALCULATED
C                    ON THE FINE-GRID-LEVEL, IF NEEDED
C     29.01.03 (TB)  : _KSR_ REMOVED
C
C**MGLET****************************************************************
C

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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
C                            DARF NUR NACH COMGRID STEHEN !!!!
C
      COMMON /COMGVP/
     &               IPCORR, IPCGES,   NDIVLEPS,  LDIVLEPS, DIVGMX

      INTEGER
     &         IPCORR(MAXGRIDS),   IPCGES(MAXGRIDS),   
     &       NDIVLEPS(MAXGRIDS), LDIVLEPS(MAXGRIDS)

      REAL   DIVGMX(MAXGRIDS)
 

C
C
      REAL        RDDX(IDIM1D),     RDDY(IDIM1D),     RDDZ(IDIM1D),
     $             RDX(IDIM1D),      RDY(IDIM1D),      RDZ(IDIM1D),
     $             DDX(IDIM1D),      DDY(IDIM1D),      DDZ(IDIM1D),
     $              DX(IDIM1D),       DY(IDIM1D),       DZ(IDIM1D),
     $               X(IDIM1D),        Y(IDIM1D),        Z(IDIM1D)

      REAL        U (   IDIM3D  ), V (   IDIM3D  ), W (   IDIM3D  ),
     $            P (   IDIM3D  ), G (   IDIM3D  ), B (   IDIM3D  ),
     $            DP(   IDIM3D  ), BP(   IDIM3D  ), BU(   IDIM3D  ),
     $            BV(   IDIM3D  ), BW(   IDIM3D  ), SDIV(IDIM3D)

      REAL        UFR(IDIM2D*2),VRI(IDIM2D*2)
      REAL        UBO(IDIM2D*2),VBO(IDIM2D*2),WBO(IDIM2D*2)
      REAL        VFR(IDIM2D*2),WFR(IDIM2D*2)
      REAL        UBA(IDIM2D*2),VBA(IDIM2D*2),WBA(IDIM2D*2)

C
C                                 FELDER FUER DEN DIREKTEN POISSONLOESER
C
      REAL
     $           FPSFAK(IDIM2D),FCOSMY(IDIM1D),FCOSNY(IDIM1D)
C
C                                 FELDER FUER DIE FOURIERTRANSFORMATION
C
      INTEGER IFFTX(19,MAXGRIDS),IFFTY(19,MAXGRIDS)
      REAL    RFFTX(IDIM1D)     ,RFFTY(IDIM1D)
C
C
      REAL     HILF( IDIM3D ), DIV( IDIM3D ),DIVG(IDIM2D)
C
      LOGICAL HOMOG
C     HOMOG = .FALSE.
C
C
C
C
C                                 KONSISTENZ ???
C
      IF (ILEVEL .LT. MINCALLVL) CALL ERRR(501," MGVPC4 ")
C
C
C                                ABBRUCHSCHRANKE
C
      EPC4 = EPCORR / (EPFAK*8.0)
      MPC4 = MPCORR
C

      DO I = 1,NOFVPIT(ILEVEL)

           IGRID = IGRDOFVPIT( I , ILEVEL )

                      IPCORR(IGRID) = 0
                      LDIVLEPS(IGRID) = 0
                      NDIVLEPS(IGRID) = 0
      ENDDO
C
C---- --------------------------------- GESAMTER MULTIGRID-DRUCKKORREKTUR
C---- ---------------------------------ZYKLUS WIRD NUR DURCHLAUFEN, FALLS
C---- --------------------------------- MINDESTENS EIN ITERATIVES GITTER
CCCCCC      IF (NOFVPIT(ILEVEL) .GE. 1) THEN



      DO 2100 IPCOUNT = 1,MPC4

C---- ---------------------------------------- ZUERST DIREKTE LOESUNG
C---- ---------------------------------------- DER POISSONGLEICHUNG
C---- ---------------------------------------- FALLS MOEGLICH UND GEWUENSCHT
C
      IF (NOFPSDIR(ILEVEL) .GE. 1) THEN


       CALL MGPSDIR (IDIM3D,IDIM2D,IDIM1D,NBND,
     $               X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $               RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $               FPSFAK,FCOSMY,FCOSNY,
     $               IFFTX,IPERMUX,RFFTX,IFFTY,IPERMUY,RFFTY,  
     $               U,V,W,P,DP,G,B,DIV,UFR,VFR,WFR,
     $               UTO,VTO,WTO,PTO,GTO,VRI,HILF,
     $               H2D1,H2D2,H2D3,H2D4,H2D5,H2D6,
     $               WSOR, RHO,DIVG,
     $               ILEVEL,EPCORR,MPCVOR,TIMEPH,
     $               UBA,VBA,WBA,UBO,VBO,WBO,
     $               PFR,GFR,BU,BV,BW,SDIV
     $              )     

      ENDIF

C
C                                VORGLAETTUNGEN
C
      IF (NOFVPIT(ILEVEL) .GE. 1) THEN
C---- ---------------------------------------- ITERATIONS ARE DONE ON DP
         CALL MGVPIT (IDIM3D,IDIM2D,IDIM1D,NBND,
     $                      X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                      RDX,RDY,RDZ,RDDX,RDDY,RDDZ,   
     $                      U,V,W,P,DP,G,B,BP,BU,BV,BW,SDIV,
     $                      DIV,UFR,VFR,WFR,
     $                      UTO,VTO,WTO,PTO,GTO,VRI,HILF,
     $                      OMBETA,WSOR, RHO,DIVG,
     $                      ILEVEL,EPCORR,MPCVOR,TIMEPH,
     $                      UBA,VBA,WBA,PBA,GBA,GEOVP,UBO,VBO,WBO,
     $                      PFR,GFR,IALGO
     $                     )     
      ENDIF
C
C
C
C                                 BEGINN GROBGITTERZYKLUS
C

         IF ( ILEVEL .GT. MINCALLVL ) THEN

C
C                                RESTRIKTION DER GESCHWINDIGKEITSFELDER,
C                                 INITIALISIERUNG DES 
C                                 KORREKTURDRUCKES AUS GROBGITTER
C
C
C---- -------------------------------------------- ALL GRIDS OF LEVEL 
         DO I = 1,NOFVPIT(ILEVEL)
           IGRID = IGRDOFVPIT( I , ILEVEL )

           IF ( .NOT. LSLICE(IGRID)) THEN
C
C                                          NUR FUER GITTER, DIE
C                                          "PARENTGITTER" BESITZEN
C                                          "PARENT" == GROBGITTER
C
           IF (IPARENT(IGRID) .NE. 0 ) THEN

             IPROCF  = 1
             IPROCC  = 1
              CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $                   NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

              CALL MGDPB (KKC,JJC,IIC,IP3C,IP2C,IP1C,IBBC,IBUC,NFROC,
     $               NBACC,NRGTC,NLFTC,NBOTC,NTOPC,NCUBC,IPARENT(IGRID))

        IF (IPCOUNT.EQ.1) THEN
        CALL SETS   (KKC,JJC,IIC,KKC,JJC,IIC,SDIV(IP3C),0.0)
              CALL MGFTOC (KK,JJ,II,DX(IP1),DY(IP1),DZ(IP1),
     &                      DDX(IP1),DDY(IP1),DDZ(IP1),SDIV(IP3),
     &                      KKC,JJC,IIC,DX(IP1C),DY(IP1C),DZ(IP1C),
     &                      DDX(IP1C),DDY(IP1C),DDZ(IP1C),SDIV(IP3C),
     &               IPOSITION(IGRID),JPOSITION(IGRID),KPOSITION(IGRID),
     &                      'P',IPROCF,IPROCC,
     &                      NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,BP(IP3))
        ENDIF
               CALL MGFTOC (KK,JJ,II,DX(IP1),DY(IP1),DZ(IP1),
     &                      DDX(IP1),DDY(IP1),DDZ(IP1),U(IP3),
     &                      KKC,JJC,IIC,DX(IP1C),DY(IP1C),DZ(IP1C),
     &                      DDX(IP1C),DDY(IP1C),DDZ(IP1C),U(IP3C),
     &               IPOSITION(IGRID),JPOSITION(IGRID),KPOSITION(IGRID),
     &                      'U',IPROCF,IPROCC,
     &                      NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,BU(IP3))

               CALL MGFTOC (KK,JJ,II,DX(IP1),DY(IP1),DZ(IP1),
     &                      DDX(IP1),DDY(IP1),DDZ(IP1),V(IP3),
     &                      KKC,JJC,IIC,DX(IP1C),DY(IP1C),DZ(IP1C),
     &                      DDX(IP1C),DDY(IP1C),DDZ(IP1C),V(IP3C),
     &               IPOSITION(IGRID),JPOSITION(IGRID),KPOSITION(IGRID),
     &                      'V',IPROCF,IPROCC,
     &                      NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,BV(IP3))

               CALL MGFTOC (KK,JJ,II,DX(IP1),DY(IP1),DZ(IP1),
     &                      DDX(IP1),DDY(IP1),DDZ(IP1),W(IP3),
     &                      KKC,JJC,IIC,DX(IP1C),DY(IP1C),DZ(IP1C),
     &                      DDX(IP1C),DDY(IP1C),DDZ(IP1C),W(IP3C),
     &               IPOSITION(IGRID),JPOSITION(IGRID),KPOSITION(IGRID),
     &                      'W',IPROCF,IPROCC,
     &                      NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,BW(IP3))

C
C
CTEST              CALL SETS   (KK,JJ,II,KK,JJ,II,DP(IP3),0.0)
               CALL MGFTOC (KK,JJ,II,DX(IP1),DY(IP1),DZ(IP1),
     &                      DDX(IP1),DDY(IP1),DDZ(IP1),P(IP3),
     &                      KKC,JJC,IIC,DX(IP1C),DY(IP1C),DZ(IP1C),
     &                      DDX(IP1C),DDY(IP1C),DDZ(IP1C),P(IP3C),
     &               IPOSITION(IGRID),JPOSITION(IGRID),KPOSITION(IGRID),
     &                      'P',IPROCF,IPROCC,
     &                      NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,BP(IP3))

           ENDIF
           ENDIF

         ENDDO
C
C
C---- -------------------------------------------- ALL GRIDS OF LEVEL 
C                                                  
C                           DRUCK DES GROBGITTERS WIRD, UM DELTA_P BERECHNEN
C                           ZU KOENNEN HIER AUF HILF(FEINGITTER) ABGESETZT
C
C
C                      CONNECTION CONDITIONS FOR PRESSURE ON COARSE LEVEL
C
         CALL CONNECTMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,UTO,VTO,WTO,PTO,GTO,
     $              VFR,WFR,PFR,GFR,
     $              'T',ILEVEL,0,0,0,
     $              NOFVPIT(ILEVEL-1),IGRDOFVPIT(1,ILEVEL-1))

C---- -------------------------------------------- ALL GRIDS OF LEVEL 
         DO I = 1,NOFVPIT(ILEVEL)
           IGRID = IGRDOFVPIT( I , ILEVEL )

           IF ( .NOT. LSLICE(IGRID)) THEN
C
C                                          NUR FUER GITTER, DIE
C                                          "PARENTGITTER" BESITZEN
C                                          "PARENT" == GROBGITTER
C
           IF (IPARENT(IGRID) .NE. 0 ) THEN
             IPROCF  = 1
             IPROCC  = 1

              CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $                   NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

              CALL MGDPB (KKC,JJC,IIC,IP3C,IP2C,IP1C,IBBC,IBUC,NFROC,
     $               NBACC,NRGTC,NLFTC,NBOTC,NTOPC,NCUBC,IPARENT(IGRID))


              CALL BOUNDMG
     $                (IDIM3D,IDIM2D,IDIM1D,
     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $                 'T',IPARENT(IGRID),0,0,0,TIMEPH,UBA,VBA,WBA,
     $                 PBA,GBA,UBO,VBO,WBO,
     $                 PFR,GFR
     $                )


              CALL MGCTOF (KK,JJ,II,DX(IP1),DY(IP1),DZ(IP1),
     &                     DDX(IP1),DDY(IP1),DDZ(IP1),HILF(IP3),
     &                     KKC,JJC,IIC,DX(IP1C),DY(IP1C),DZ(IP1C),
     &                     DDX(IP1C),DDY(IP1C),DDZ(IP1C),P(IP3C),
     &               IPOSITION(IGRID),JPOSITION(IGRID),KPOSITION(IGRID),
     &                     'P',IPROCF,IPROCC)
           ENDIF
           ENDIF

         ENDDO
C---- -------------------------------------------- ALL GRIDS OF LEVEL 
C                                                  
C
C                                      GROBGITTERRELAXATION
C
              CALL MGVPC5 (IDIM3D,IDIM2D,IDIM1D,NBND,
     $                     X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                     RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                    FPSFAK,FCOSMY,FCOSNY,IFFTX,RFFTX,IFFTY,RFFTY,    
     $                      U,V,W,P,DP,G,B,BP,BU,BV,BW,SDIV,    
     $                     DIV,UFR,VFR,WFR,
     $                      UTO,VTO,WTO,PTO,GTO,VRI,HILF,
     $                      H2D1,H2D2,H2D3,H2D4,H2D5,H2D6,
     $                     OMBETA,WSOR, RHO,DIVG,IPRGRID,
     $                      MINCALLVL,ILEVEL-1,MAXCALLVL,TIMEPH,
     $                     UBA,VBA,WBA,PBA,GBA,GEOVP,UBO,VBO,WBO,
     $                      PFR,GFR,IALGO
     $                    )     

C
C---- -------------------------------------------- ALL GRIDS OF LEVEL 
         DO I = 1,NOFVPIT(ILEVEL)

           IGRID = IGRDOFVPIT( I , ILEVEL )
           IF ( .NOT. LSLICE(IGRID)) THEN
C                                          NUR FUER GITTER, DIE
C                                          "PARENTGITTER" BESITZEN
C                                          "PARENT" == GROBGITTER
C
           IF (IPARENT(IGRID) .NE. 0 ) THEN

           CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $                   NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
           CALL MGDPB (KKC,JJC,IIC,IP3C,IP2C,IP1C,IBBC,IBUC,NFROC,
     $               NBACC,NRGTC,NLFTC,NBOTC,NTOPC,NCUBC,IPARENT(IGRID))

C---- -------------------------------------------- PRESSURE CORRECTION OF
C---- -------------------------------------------- COARSE GRID IS PUT ON
C---- -------------------------------------------- HILF

C---- ---------- FIRST STEP: PROLONGATION OF THE NEW PRESSURE OF
C                            THE COARSE GRID TO FINE GRID

              CALL MGCTOF (KK,JJ,II,DX(IP1),DY(IP1),DZ(IP1),
     &                     DDX(IP1),DDY(IP1),DDZ(IP1),DP(IP3),
     &                     KKC,JJC,IIC,DX(IP1C),DY(IP1C),DZ(IP1C),
     &                     DDX(IP1C),DDY(IP1C),DDZ(IP1C),P(IP3C),
     &               IPOSITION(IGRID),JPOSITION(IGRID),KPOSITION(IGRID),
     &                     'P',IPROC,IPPAR)
          ENDIF
           ENDIF
           ENDDO
C---- -------------------------------------------- ALL GRIDS OF LEVEL 
C
C                      CONNECTION CONDITIONS FOR PRESSURE
C
         CALL CONNECTMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,UTO,VTO,WTO,PTO,GTO,
     $              VFR,WFR,PFR,GFR,
     $              'T',ILEVEL,0,0,0,
     $              NOFVPIT(ILEVEL),IGRDOFVPIT(1,ILEVEL))
         CALL BPARMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,H2D1,H2D2,H2D3,
     $              UTO,VTO,WTO,PTO,GTO,UBA,VBA,WBA,PBA,GBA,
     $              VFR,WFR,PFR,GFR,
     $              'T',ILEVEL,0,0,0,
     $              NOFVPIT(ILEVEL),IGRDOFVPIT(1,ILEVEL)
     $             )      

C
C---- ---------- SECOND STEP: BUILDING PRESSURE DIFFERENCE ON FINE GRID
C
C---- -------------------------------------------- ALL GRIDS OF LEVEL 
         DO I = 1,NOFVPIT(ILEVEL)

           IGRID = IGRDOFVPIT( I , ILEVEL )
           IF ( .NOT. LSLICE(IGRID)) THEN

C
           IF (IPARENT(IGRID) .NE. 0 ) THEN
           CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $                   NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
              CALL PHIADD(KK,JJ,II,
     &                       DP(IP3),HILF(IP3),DP(IP3),0.0,1.0,-1.0)

              CALL BOUNDMG
     $                (IDIM3D,IDIM2D,IDIM1D,
     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,DP,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $                 'T',IGRID,0,0,0,TIMEPH,UBA,VBA,WBA,
     $                 PBA,GBA,UBO,VBO,WBO,
     $                 PFR,GFR
     $                )     

C
C                              AUSFUEHREN DER DRUCKKORREKTUR
C                              ONLY DP IS CORRECTED
              CALL MGPCORR
     $             (KK,JJ,II,KK,JJ,II,
     $             DX(IP1),DY(IP1),DZ(IP1),
     $             DDX(IP1),DDY(IP1),DDZ(IP1),
     $             RDX(IP1),RDY(IP1),RDZ(IP1),
     $             RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $             U(IP3),V(IP3),W(IP3),P(IP3),B(IP3),DP(IP3),
     $             BP(IP3),BU(IP3),BV(IP3),BW(IP3),
     $             RHO,DT,WSOR,NBND,
     $             NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,0)

C              CALL MGPCORR
C     $                    (KK,JJ,II,KK,JJ,II,
C     $                     DX(IP1),DY(IP1),DZ(IP1),
C     $                    DDX(IP1),DDY(IP1),DDZ(IP1),
C     $                    RDX(IP1),RDY(IP1),RDZ(IP1),
C     $                   RDDX(IP1),RDDY(IP1),RDDZ(IP1),
C     $                       IPCORR,OMBETA,DT,
C     $              U(IP3),V(IP3),W(IP3),P(IP3),B(IP3),DP(IP3),
C     $              BP(IP3),BU(IP3),BV(IP3),BW(IP3),
C     $              RHO,DIVGMX(IGRID),DIVG,WSOR,NBND,
C     $                       NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,0)


C             CALL PRLE3M (1,1,4,'J',KK,JJ,II,
C    $                U(IP3),1,V(IP3),1,W(IP3),1,P(IP3),1,
C    $                G(IP3),1,B(IP3),0,DIV(IP3),0, 1, 1, 4, 1, 1, 4,
C    $                DX(IP1),DY(IP1),DZ(IP1),X(IP1),Y(IP1),Z(IP1),
C    $                18,ITSTEP)

           ENDIF
           ENDIF
           ENDDO
C---- -------------------------------------------- ALL GRIDS OF LEVEL 
C---- -------------------------------------------- ALL GRIDS OF LEVEL 
        DO I = 1,NOFVPIT(ILEVEL)
           IGRID = IGRDOFVPIT(I,ILEVEL)
           IF ( .NOT. LSLICE(IGRID)) THEN

              CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $                      NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

              CALL BOUNDMG
     $                (IDIM3D,IDIM2D,IDIM1D,
     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $                 'T',IGRID,0,0,0,TIMEPH,UBA,VBA,WBA,
     $                 PBA,GBA,UBO,VBO,WBO,
     $                 PFR,GFR
     $                )     

           ENDIF

         ENDDO
C---- -------------------------------------------- ALL GRIDS OF LEVEL 


C
C                      GROBGITTERZYKLUS ZUENDE
C

           ENDIF
C---- -------------------------------------------- ENDIF COARSER LEVEL EXISTS 
C
C                                                  UEBERPRUEFEN DER DIVERGENZ
C                                                  ZUM DEBUGGEN
      IF (IVPINF .GE. 2) THEN
C---- -------------------------------------------- ALL GRIDS OF LEVEL 
        DO I = 1,NOFVPIT(ILEVEL)
           IGRID = IGRDOFVPIT(I,ILEVEL)
           IF ( .NOT. LSLICE(IGRID)) THEN


              CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $                      NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
              CALL DIVCAL
     $                    (KK,JJ,II,KK,JJ,II,
     $                    RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                     U(IP3),V(IP3),W(IP3),BP(IP3),DIV(IP3),
     $                    1.0,-30,DIVGMX(IGRID),BU(IP3),BV(IP3),
     $                 BW(IP3),SDIV(IP3))

           ENDIF

         ENDDO
C---- -------------------------------------------- ALL GRIDS OF LEVEL 
           ENDIF
C---- -------------------------------------------- ENDE DER DEBUG-SCHLEIFE

C---- -------------------------------------------- IF ITERATIVE GRIDS
      IF (NOFVPIT(ILEVEL) .GE. 1) THEN

      CALL MGVPIT (IDIM3D,IDIM2D,IDIM1D,NBND,
     $                      X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                      RDX,RDY,RDZ,RDDX,RDDY,RDDZ,    
     $                      U,V,W,P,DP,G,B,BP,BU,BV,BW,SDIV,
     $                      DIV,UFR,VFR,WFR,
     $                      UTO,VTO,WTO,PTO,GTO,VRI,HILF,
     $                      OMBETA,WSOR, RHO,DIVG,
     $                      ILEVEL,EPCORR,MPCNACH,TIMEPH,
     $                      UBA,VBA,WBA,PBA,GBA,GEOVP,UBO,VBO,WBO,
     $                      PFR,GFR,IALGO
     $                     )     

       ENDIF
C---- -------------------------------------------- ENDIF ITERATIVE GRIDS

C
C                                   UEBERPRUEFEN DES ABBRUCHKRITERIUMS
C
C---- -------------------------------------------- ALL GRIDS OF LEVEL 
        DO I = 1,NOFVPIT(ILEVEL)
           IGRID = IGRDOFVPIT(I,ILEVEL)
           IF ( .NOT. LSLICE(IGRID)) THEN


              CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $                      NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
CC            CALL DIVCAL
C    $                    (KK,JJ,II,KK,JJ,II,
C    $                    RDDX(IP1),RDDY(IP1),RDDZ(IP1),
C    $                     U(IP3),V(IP3),W(IP3),B(IP3),DIV(IP3),
C    $                    1.0,-2,DIVGMX(IGRID))

           ENDIF

         ENDDO
C---- -------------------------------------------- ALL GRIDS OF LEVEL 

C
C                   FALLS MAXIMALE DIVERGENZ ALLER GITTER 
C                   DES AKTUELLEN LEVELS KLEINER
C                   ALS DIE SCHRANKE, WIRD NACHITERATION ABGEBROCHEN
C

         DIVMAX = 0.0

         DO I = 1,NOFVPIT(ILEVEL)

           IGRID = IGRDOFVPIT( I , ILEVEL )

           DIVMAX = MAX(DIVMAX,ABS(DIVGMX(IGRID)))

         ENDDO

C
C

         IF(DIVMAX .LT. EPC4) GOTO 2200



 2100    CONTINUE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  END OF PRESSURE-ITERATION

 2200    CONTINUE


CTEST         DO I = 1,NOFVPIT(ILEVEL)
CTEST
CTEST           IGRID = IGRDOFVPIT( I , ILEVEL )

CTEST#ifdef _MPI_
CTEST           IF ( MYID .EQ. IDPROCOFGRD(IGRID) ) THEN
CTEST#endif

CTEST           CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
CTEST     $                   NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

CTEST           CALL PHIADD(KK,JJ,II,P,DP,P,0.0,1.0,1.0)

CTEST#ifdef _MPI_
CTEST           ENDIF
CTEST#endif

CTEST         ENDDO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC BOUNDARY-CONDITIONS FOR DP
C
C
C                                 BEI MESSAGE-PASSING, KOMMUNIKATION
C                                 FUER DRUCKRANDBEDINGUNG IM
C                                 GROBGITTER WICHTIG
C
         CALL CONNECTMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,UTO,VTO,WTO,PTO,GTO,
     $              VFR,WFR,PFR,GFR,
     $              'T',ILEVEL,0,0,0,
     $              NOFVPIT(ILEVEL),IGRDOFVPIT(1,ILEVEL))

        DO I = 1,NOFVPIT(ILEVEL)
           IGRID = IGRDOFVPIT(I,ILEVEL)

           IF ( .NOT. LSLICE(IGRID)) THEN
C

              CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $                      NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

              CALL BOUNDMG
     $                (IDIM3D,IDIM2D,IDIM1D,
     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $                 'T',IGRID,0,0,0,TIMEPH,UBA,VBA,WBA,
     $                 PBA,GBA,UBO,VBO,WBO,
     $                 PFR,GFR
     $                )     
       ENDIF

         ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  

         RETURN
 1000    FORMAT (1X,A,I4,3X,A,2X,I4,3X,A,2X,I4,2X,A,2X,E10.4)
         END
