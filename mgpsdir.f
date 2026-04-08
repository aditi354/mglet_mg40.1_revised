










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
         SUBROUTINE MGPSDIR (IDIM3D,IDIM2D,IDIM1D,NBND,
     $                      X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                      RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                      FPSFAK,FCOSMY,FCOSNY,
     $                      IFFTX,IPERMUX,RFFTX,IFFTY,IPERMUY,RFFTY,    
     $                      U,V,W,P,DP,G,B,DIV,UFR,VFR,WFR,
     $                      UTO,VTO,WTO,PTO,GTO,VRI,HILF,
     $                      H2D1,H2D2,H2D3,H2D4,H2D5,H2D6,
     $                      WSOR, RHO,DIVG,
     $                      ILEVEL,EPIT,MPIT,TIMEPH,
     $                      UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                      PFR,GFR,BU,BV,BW,SDIV
     $                      )     

C**MGLET****************************************************************
C
C   M G P S D I R    DIRECT-SOLUTION OF POIISON-EQUATION
C
C     4.10.93  (MM):  ORIGINAL
C     29.01.03 (TB): _KSR_ REMOVED
C
C**MGLET****************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

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
 
C
      REAL        RDDX(IDIM1D),     RDDY(IDIM1D),     RDDZ(IDIM1D),
     $             RDX(IDIM1D),      RDY(IDIM1D),      RDZ(IDIM1D),
     $             DDX(IDIM1D),      DDY(IDIM1D),      DDZ(IDIM1D),
     $              DX(IDIM1D),       DY(IDIM1D),       DZ(IDIM1D),
     $               X(IDIM1D),        Y(IDIM1D),        Z(IDIM1D)

      REAL        U (   IDIM3D  ), V (   IDIM3D  ), W (   IDIM3D  ),
     $            P (   IDIM3D  ), G (   IDIM3D  ), B (   IDIM3D  ),
     $            DP(   IDIM3D  ),BU (   IDIM3D  ),BV (   IDIM3D  ),
     $            BW(   IDIM3D  ),SDIV(   IDIM3D  )
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

      REAL        UFR(IDIM2D*2),VRI(IDIM2D*2)
      REAL        UBO(IDIM2D*2),VBO(IDIM2D*2),WBO(IDIM2D*2)
      REAL        VFR(IDIM2D*2),WFR(IDIM2D*2)
      REAL        PFR(IDIM2D*2),GFR(IDIM2D*2)
      REAL        UBA(IDIM2D*2),VBA(IDIM2D*2),WBA(IDIM2D*2)
      REAL      UTO(IDIM2D*2),  VTO(IDIM2D*2),  WTO(IDIM2D*2),
     $           PTO(IDIM2D*2),  GTO(IDIM2D*2)
C
      REAL     HILF( IDIM3D ),  DIV( IDIM3D ), DIVG( IDIM2D )
      REAL     H2D1( IDIM2D ), H2D2( IDIM2D ), H2D3( IDIM2D ),
     $         H2D4( IDIM2D ), H2D5( IDIM2D ), H2D6( IDIM2D )
C
C
C
C                                 VORBELEGUNG DER MAX. DIVERGENZ
         F = 1
         DO I = 1,NOFPSDIR(ILEVEL)
           IGRID = IGRDOFPSDIR( I , ILEVEL )
           DIVGMX(IGRID) = 0.0
         ENDDO
C
C---- -------------------------------------- CONNECTIVITY-BEDINGUNGEN
C
         CALL CONNECTMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,UTO,VTO,WTO,PTO,GTO,
     $              VFR,WFR,PFR,GFR, 
     $              'T',ILEVEL,0,0,0,
     $              NOFPSDIR(ILEVEL),IGRDOFPSDIR(1,ILEVEL))
         F = 2.0
C
C---- --------------------------------------------------------------------
C
         DO I = 1,NOFPSDIR(ILEVEL)
           IGRID = IGRDOFPSDIR( I , ILEVEL )

C
C---- --------------------------------------- POINTER UND DIMENSIONIERUNGEN
C
              CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $                      NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

C
C                              ZUERST EINMAL RANDBEDINGUNGEN SICHERSTELLEN
C
              CALL BOUNDMG
     $                (IDIM3D,IDIM2D,IDIM1D,
     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $                 'P',IGRID,0,0,0,
     $                 TIMEPH,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                 PFR,GFR
     $                )       
              F = 3.0
C
C                              DP-FELD AUF NULL SETZEN
C
              CALL SETS   (KK,JJ,II,KK,JJ,II,DP(IP3),0.0)
C
C---- ----------------------------------------- BERECHNUNG DER DIVERGENZ
C
C                           DER TRIDIAGONALLOESER VERLANGT EINEN
C                           VORFAKTOR 

              PREFAK = FPSFAK ( IP2 + 3 )/DT*RHO

C              WRITE (6,*) 'PREFAK IN MGPSDIR: ',PREFAK

              IREZIP = -1
              IF (IVPINF .EQ. 2) IREZIP = -3 
              CALL DIVCAL
     $                    (KK,JJ,II,KK,JJ,II,
     $                    RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                    U(IP3),V(IP3),W(IP3),B(IP3),DP(IP3),
     $                    PREFAK,IREZIP,DIVGMX(IGRID),
     $                    BU(IP3),BV(IP3),BW(IP3),SDIV(IP3))

C
C---- -------------------------------------- LOESUNG DER POISSONGLEICHUNG
C
            CALL CAP1252 (FPSFAK(IP2),FCOSMY(IP1),FCOSNY(IP1)
     $                    ,RFFTX(IP1*2),IPERMUX(IP1),IFFTX(1,IGRID)
     $                    ,RFFTY(IP1*2),IPERMUY(IP1),IFFTY(1,IGRID)
     $                    ,KK,JJ,II,NBND,DP(IP3),HILF(IP3)
     $                    ,H2D1(IP2),H2D2(IP2),H2D3(IP2)
     $                    ,H2D4(IP2),H2D5(IP2),H2D6(IP2)
     $                    ,NFRO,NBAC,NRGT,NLFT,NBOT,NTOP)
C
C---- -------------------------------------- SCHREIBEN D. KORREKTURDRUCKES
C
C
C---- -------------------------------------- SETZEN DER RANDBED. FUER DRUCK
C
              CALL BOUNDMG
     $                (IDIM3D,IDIM2D,IDIM1D,
     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,DP,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $                 'T',IGRID,0,0,0,
     $                 TIMEPH,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                 PFR,GFR
     $                )
C
C---- -------------------------------------- AUSFUEHREN DER DRUCKKORREKTUR
C
              CALL MGPCORR
     $                    (KK,JJ,II,KK,JJ,II,
     $                     DX(IP1),DY(IP1),DZ(IP1),
     $                    DDX(IP1),DDY(IP1),DDZ(IP1),
     $                    RDX(IP1),RDY(IP1),RDZ(IP1),
     $                   RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                     U(IP3),V(IP3),W(IP3),P(IP3),B(IP3),DP(IP3),
     $                     B(IP3),BU(IP3),BV(IP3),BW(IP3),
     $                       RHO,DT,WSOR,NBND,
     $                       NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,0)
C
C---- ------------------------------ SETZEN DER RANDBED. FUER GESCHWINDIGKEITEN
C
              CALL BOUNDMG
     $                (IDIM3D,IDIM2D,IDIM1D,
     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $                 'T',IGRID,0,0,0,
     $                 TIMEPH,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                 PFR,GFR
     $                )

C        DO F=IP1,(IP1+JJ)
C           WRITE(6,*)'IPERMUY: ',IPERMUY(F),F
C        ENDDO     
C
C
C---- ----------------------------------------
      ENDDO


C
C---- ----------------------------------------------------------------------
C
       DO I = 1,NOFPSDIR(ILEVEL)
         IGRID = IGRDOFPSDIR( I , ILEVEL )

              IPCORR(IGRID) = IPCORR(IGRID) + 1
              IPCGES(IGRID) = IPCGES(IGRID) + 1

C
C---- ----------------------------------------- KONTROLLE DER DIVERGENZ
C
      IF (IVPINF.EQ.2) THEN

              CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $                      NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

              PREFAK = 1.0

              CALL DIVCAL
     $                    (KK,JJ,II,KK,JJ,II,
     $                    RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                    U(IP3),V(IP3),W(IP3),B(IP3),DIV(IP3),
     $                    PREFAK,-3,DIVGMX(IGRID),
     $                    BU(IP3),BV(IP3),BW(IP3),SDIV(IP3))

      ENDIF
           IF (IVPINF.EQ.2) WRITE (6,1000)
     $        'MGPSDIR, GRD',IGRID,'LEVEL',ILEVEL,
     $        'IPCORR',IPCORR(IGRID),'DIV',DIVGMX(IGRID)


       ENDDO
 

         RETURN
 1000    FORMAT (1X,A,I4,3X,A,I4,3X,A,2X,I4,2X,A,2X,E10.4)
         END
