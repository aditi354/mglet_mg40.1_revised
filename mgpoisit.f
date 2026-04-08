










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
      SUBROUTINE MGPOISIT (IDIM3D,IDIM2D,IDIM1D,NBND,
     $                      X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                      RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                      U,V,W,P,DP,RHS,RES,RHO,
     $                      GSAW,GSAE,GSAN,GSAS,GSAT,GSAB,
     $                      G,B,UFR,VFR,WFR,HILF,
     $                      UTO,VTO,WTO,PTO,GTO,VRI,ILEVEL,
     $                      PFR,GFR,TIMEPH,  
     $                      EPIT,MPIT,
     $                      UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                      BP,
     $                      GSAP,SIPLW,SIPLS,
     $                      SIPLB,SIPLPR,SIPUE,SIPUN,SIPUT
     $                     ) 
C**MGLET****************************************************************
C
C   M G P O I S I T      MULTI-GRID-VELOCITY-PRESSURE-ITERATION
C                    DRUCKKORREKTURZYKLUS FUER MULTIGRID-VARIANTE
C                    HIER WERDEN GITTER DES LEVELS ILEVEL KORRIGIERT
C
C     4.10.93 (MM):  ORIGINAL
C    23.05.95 (MM.AO:) MPI EINGEFUEHRT
C    13.08.03 (TB):  EXTENDED FOR SCALAR TRANSPORT
C
C**MGLET****************************************************************


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

C
C                            DARF NUR NACH COMGRID STEHEN !!!!
C
      COMMON /COMGVP/
     &               IPCORR, IPCGES,   NDIVLEPS,  LDIVLEPS, DIVGMX

      INTEGER
     &         IPCORR(MAXGRIDS),   IPCGES(MAXGRIDS),   
     &       NDIVLEPS(MAXGRIDS), LDIVLEPS(MAXGRIDS)

      REAL   DIVGMX(MAXGRIDS)
 


      REAL        RDDX(IDIM1D),     RDDY(IDIM1D),     RDDZ(IDIM1D),
     $             RDX(IDIM1D),      RDY(IDIM1D),      RDZ(IDIM1D),
     $             DDX(IDIM1D),      DDY(IDIM1D),      DDZ(IDIM1D),
     $              DX(IDIM1D),       DY(IDIM1D),       DZ(IDIM1D),
     $               X(IDIM1D),        Y(IDIM1D),        Z(IDIM1D)

      REAL        U (   IDIM3D  ), V (   IDIM3D  ), W (   IDIM3D  ),
     $            P (   IDIM3D  ), DP(   IDIM3D  ), G (   IDIM3D  ), 
     $            B (   IDIM3D  )

      REAL        UFR(IDIM2D*2),VRI(IDIM2D*2)
      REAL        VFR(IDIM2D*2),WFR(IDIM2D*2)
      REAL   UBO(IDIM2D*2),VBO(IDIM2D*2),WBO(IDIM2D*2)
      REAL   UBA(IDIM2D*2),VBA(IDIM2D*2),WBA(IDIM2D*2)
      REAL   PBA(IDIM2D*2),GBA(IDIM2D*2)
      REAL   PFR(IDIM2D*2),GFR(IDIM2D*2)

      REAL   RHS (IDIM3D), HILF( IDIM3D )

      REAL   GSAW(IDIM1D), GSAE(IDIM1D), GSAN(IDIM1D),
     $       GSAS(IDIM1D), GSAT(IDIM1D), GSAB(IDIM1D)
      REAL   BP  ( IDIM3D )
      REAL RES   ( IDIM3D ), SIPLB ( IDIM3D ), SIPUT( IDIM3D ),
     $     SIPLW ( IDIM3D ), SIPUE ( IDIM3D ), SIPLS( IDIM3D ),
     $     SIPUN ( IDIM3D ), SIPLPR( IDIM3D ), GSAP ( IDIM3D )


C
      IF (MPIT.LE.0) RETURN
C
C                               ES WERDEN MPIT ITERATIONEN DURCHGEFUEHRT
C

      DO IPIT = 1,MPIT

      IF (NOFVPIT(ILEVEL) .GT. 1) THEN

         CALL CONNECTMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,DP,G,B,UFR,VRI,WBO,HILF,UTO,VTO,WTO,PTO,GTO,
     $              VFR,WFR,PFR,GFR,
     $              'S',ILEVEL,0,0,0,
     $              NOFVPIT(ILEVEL),IGRDOFVPIT(1,ILEVEL))

      ENDIF


C
C
         DO I = 1,NOFVPIT(ILEVEL)
           IGRID = IGRDOFVPIT( I , ILEVEL )
           DIVGMX(IGRID) = 0.0
         ENDDO
C
C               START OF MULTIGRIDCYCLE
C
C
         DO I = 1,NOFVPIT(ILEVEL)
         IGRID = IGRDOFVPIT( I , ILEVEL )


c              CALL BOUNDMG Ausgeschaltet wegen BTOOP1
          CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $         NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
          CALL SIPITER1 (KK,JJ,II,DP(IP3),RHS(IP3),RES(IP3),
     $         GSAW(IP1),GSAE(IP1),GSAN(IP1),
     $         GSAS(IP1),GSAT(IP1),GSAB(IP1),GSAP(IP3),
     $         BP(IP3),
     $         NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,
     $         SIPLW(IP3),SIPLS(IP3),SIPLB(IP3),
     $         SIPLPR(IP3)) 

       ENDDO

      IF (NOFVPIT(ILEVEL) .GT. 1) THEN

         CALL CONNECTMG
     $        (IDIM3D,IDIM2D,IDIM1D,
     $        X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $        U,V,W,RES,G,B,UFR,VRI,WBO,HILF,UTO,VTO,WTO,PTO,GTO,
     $        VFR,WFR,PFR,GFR,
     $        'S',ILEVEL,0,0,0,
     $        NOFVPIT(ILEVEL),IGRDOFVPIT(1,ILEVEL))

      ENDIF

         DO I = 1,NOFVPIT(ILEVEL)
         IGRID = IGRDOFVPIT( I , ILEVEL )


          CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $         NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

          CALL SIPITER2  (KK,JJ,II,DP(IP3),RES(IP3),
     $         NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,DIVGMX(IGRID),
     $         SIPUE(IP3),SIPUN(IP3),SIPUT(IP3)) 



          IPCORR(IGRID) = IPCORR(IGRID) + 1
          IPCGES(IGRID) = IPCGES(IGRID) + 1

          IF (IVPINF.EQ.2) WRITE (6,1000)
     $         'MGPOISIT, GRD',IGRID,'LEVEL',ILEVEL,
     $         'IPCORR',IPCORR(IGRID),'RESMAX',DIVGMX(IGRID)

          IF(DIVGMX(IGRID) .LT. EPIT) THEN
             LDIVLEPS(IGRID) = 1
             NDIVLEPS(IGRID) = NDIVLEPS(IGRID) + 1
          ENDIF



         ENDDO
      ENDDO
C                             ITERATIONEN SIND FERTIG
C

      IF (NOFVPIT(ILEVEL) .GT. 1) THEN

         CALL CONNECTMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,DP,G,B,UFR,VRI,WBO,HILF,UTO,VTO,WTO,PTO,GTO,
     $              VFR,WFR,PFR,GFR,
     $              'S',ILEVEL,0,0,0,
     $              NOFVPIT(ILEVEL),IGRDOFVPIT(1,ILEVEL))

      ENDIF

         RETURN
 1000    FORMAT (1X,A,I4,3X,A,I4,3X,A,2X,I4,2X,A,2X,E10.4)
         END
