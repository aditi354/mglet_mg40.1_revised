










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
         SUBROUTINE MGVPIT (IDIM3D,IDIM2D,IDIM1D,NBND,
     $                      X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                      RDX,RDY,RDZ,RDDX,RDDY,RDDZ,     
     $                      U,V,W,P,DP,G,B,BP,BU,BV,BW,SDIV,
     $                      DIV,UFR,VFR,WFR,
     $                      UTO,VTO,WTO,PTO,GTO,VRI,HILF,
     $                      OMBETA,WSOR, RHO,DIVG,
     $                      ILEVEL,EPIT,MPIT,TIMEPH,
     $                      UBA,VBA,WBA,PBA,GBA,GEOVP,UBO,VBO,WBO,
     $                      PFR,GFR,IALGO
     $                     )          

C**MGLET****************************************************************
C
C   M G V P I T      MULTI-GRID-VELOCITY-PRESSURE-ITERATION
C                    DRUCKKORREKTURZYKLUS FUER MULTIGRID-VARIANTE
C                    HIER WERDEN GITTER DES LEVELS ILEVEL KORRIGIERT
C
C     4.10.93 (MM):  ORIGINAL
C    23.05.95 (MM.AO:) MPI EINGEFUEHRT
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
     $            DP(   IDIM3D  ), BP(   IDIM3D  ), BU(   IDIM3D  ),
     $            BV(   IDIM3D  ), BW(   IDIM3D  ),SDIV(  IDIM3D  )

      REAL        UFR(IDIM2D*2),VRI(IDIM2D*2)
      REAL        VFR(IDIM2D*2),WFR(IDIM2D*2)
      REAL   UBO(IDIM2D*2),VBO(IDIM2D*2),WBO(IDIM2D*2)
      REAL   UBA(IDIM2D*2),VBA(IDIM2D*2),WBA(IDIM2D*2)
      REAL   PBA(IDIM2D*2),GBA(IDIM2D*2)
      REAL   PFR(IDIM2D*2),GFR(IDIM2D*2)

C
      REAL     HILF( IDIM3D ), DIV( IDIM3D ),DIVG(IDIM2D)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IF (MPIT.LE.0) RETURN
C
C                             RANDBEDINGUNG WIRD  EINMAL AUFGEPRAEGT
C
C                             NOTWENDIG  IM RED-BLACK VERFAHREN
C
      IF (NOFVPIT(ILEVEL) .GT. 1) THEN

         DO IRB = 1,2
         DO JRB = 1,2
         DO KRB = 1,2

C
C                                 BEI MESSAGE-PASSING, KOMMUNIKATION
C
         CALL CONNECTMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,UTO,VTO,WTO,PTO,GTO,
     $              VFR,WFR,PFR,GFR,
     $              'P',ILEVEL,IRB,JRB,KRB,
     $              NOFVPIT(ILEVEL),IGRDOFVPIT(1,ILEVEL))

        DO I = 1,NOFVPIT(ILEVEL)
           IGRID = IGRDOFVPIT(I,ILEVEL)

              CALL BOUNDMG
     $                (IDIM3D,IDIM2D,IDIM1D,
     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $                 'P',IGRID,IRB,JRB,KRB,
     $                 TIMEPH,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                 PFR,GFR
     $                )     
              IF (IRB*JRB*KRB.EQ.1.AND.IALGO.NE.1) THEN

              CALL BOUNDMG
     $                (IDIM3D,IDIM2D,IDIM1D,
     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $                 'Z',IGRID,IRB,JRB,KRB,
     $                 TIMEPH,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                 PFR,GFR
     $                )     

              END IF

         ENDDO
         ENDDO
         ENDDO
         ENDDO
C
      ENDIF

      DO IPIT = 1,MPIT
C
C
C                                VARIATION OF OMG THROUGH OMBETA
C                                 VORBELEGUNG DER MAX. DIVERGENZ
CTEST      IF(MOD(IPIT,2) .EQ. 0) THEN
CTEST        OMBETA1=OMBETA*0.90
CTEST      ELSE
        OMBETA1=OMBETA
CTEST      ENDIF
C     WRITE(6,*) 'OMBETA=',OMBETA1
         DO I = 1,NOFVPIT(ILEVEL)
           IGRID = IGRDOFVPIT( I , ILEVEL )
           DIVGMX(IGRID) = 0.0
         ENDDO

C
C                                FALLS MEHR ALS EIN GITTER AUF DEM 
C                                LEVEL EXISTIEREN, WERDEN ITERATIONEN
C                                IM RED-BLACK-VERFAHREN DURCHLAUFEN
C
      IF (NOFVPIT(ILEVEL) .GT. 1) THEN

         DO IRB = 1,2
         DO JRB = 1,2
         DO KRB = 1,2
         CALL CONNECTMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,UTO,VTO,WTO,PTO,GTO,
     $              VFR,WFR,PFR,GFR,
     $              'P',ILEVEL,IRB,JRB,KRB,
     $              NOFVPIT(ILEVEL),IGRDOFVPIT(1,ILEVEL))

           DO I = 1,NOFVPIT(ILEVEL)
           IGRID = IGRDOFVPIT(I,ILEVEL)

              CALL BOUNDMG
     $                (IDIM3D,IDIM2D,IDIM1D,
     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $                 'P',IGRID,IRB,JRB,KRB,
     $                 TIMEPH,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                 PFR,GFR
     $                )     

              CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $                      NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
              CALL VPRBC
     $                    (KK,JJ,II,KK,JJ,II,
     $                     DX(IP1),DY(IP1),DZ(IP1),
     $                    DDX(IP1),DDY(IP1),DDZ(IP1),
     $                    RDX(IP1),RDY(IP1),RDZ(IP1),
     $                   RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                       IPCORR,OMBETA1,DT,
     $                     U(IP3),V(IP3),W(IP3),P(IP3),B(IP3),
     $                     BP(IP3),BU(IP3),BV(IP3),BW(IP3),SDIV(IP3),
     $                     HILF(IP3),
     $                       RHO,DIVGMX(IGRID),DIVG(IP2),WSOR,NBND,
     $                       NFRO,NBAC,NRGT,NLFT,IRB,JRB,KRB,GEOVP
     $                       ,IALGO)
              IF (IRB*JRB*KRB.eq.8.AND.IALGO.NE.1) THEN
              CALL BOUNDMG
     $                (IDIM3D,IDIM2D,IDIM1D,
     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $                 'Y',IGRID,IRB,JRB,KRB,
     $                 TIMEPH,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                 PFR,GFR
     $                )     
              END IF


           ENDDO
         ENDDO
         ENDDO
         ENDDO


C
C                                       NUR EIN GITTER AUF LEVEL
C
      ELSEIF (NOFVPIT(ILEVEL) .EQ. 1) THEN
         IGRID = IGRDOFVPIT( 1 , ILEVEL )
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  ZUM TEST FUER PARENT-RANDB.
CTEST         CALL CONNECTMG
CTEST     $             (IDIM3D,IDIM2D,IDIM1D,
CTEST     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
CTEST     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,UTO,VTO,WTO,PTO,GTO,
CTEST     $              VFR,WFR,PFR,GFR,
CTEST     $              'P',ILEVEL,0,0,0,
CTEST     $              NOFVPIT(ILEVEL),IGRDOFVPIT(1,ILEVEL))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  ZUM TEST FUER PARENT-RANDB.
  

              CALL BOUNDMG
     $                (IDIM3D,IDIM2D,IDIM1D,
     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $                 'P',IGRID,0,0,0,
     $                 TIMEPH,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                 PFR,GFR
     $                )     
              IF (IALGO.NE.1) THEN
              CALL BOUNDMG
     $                (IDIM3D,IDIM2D,IDIM1D,
     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $                 'Z',IGRID,IRB,JRB,KRB,
     $                 TIMEPH,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                 PFR,GFR
     $                )     
              END IF


             
              CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $                      NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
              CALL VPLESC
     $                    (KK,JJ,II,KK,JJ,II,
     $                     DX(IP1),DY(IP1),DZ(IP1),
     $                    DDX(IP1),DDY(IP1),DDZ(IP1),
     $                    RDX(IP1),RDY(IP1),RDZ(IP1),
     $                   RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                       IPCORR,OMBETA1,DT,
     $                     U(IP3),V(IP3),W(IP3),P(IP3),B(IP3),
     $                     BP(IP3),BU(IP3),BV(IP3),BW(IP3),SDIV(IP3),
     $                     HILF(IP3),
     $                       RHO,DIVGMX(IGRID),DIVG(IP2),WSOR,NBND,
     $                       NFRO,NBAC,NRGT,NLFT,GEOVP,
     $                       IALGO)
             IF (IALGO.NE.1) THEN
              CALL BOUNDMG
     $                (IDIM3D,IDIM2D,IDIM1D,
     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $                 'Y',IGRID,0,0,0,
     $                 TIMEPH,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                 PFR,GFR
     $                )     
             ENDIF


      ELSE

        CALL ERRR (501,'MGVPIT')

      ENDIF

         DO I = 1,NOFVPIT(ILEVEL)

           IGRID = IGRDOFVPIT( I , ILEVEL )


              IPCORR(IGRID) = IPCORR(IGRID) + 1
              IPCGES(IGRID) = IPCGES(IGRID) + 1

           IF (IVPINF.EQ.2) WRITE (6,1000)
     $        'MGVPIT, GRD',IGRID,'LEVEL',ILEVEL,
     $        'IPCORR',IPCORR(IGRID),'DIV',DIVGMX(IGRID)

              IF(DIVGMX(IGRID) .LT. EPIT) THEN
   
                LDIVLEPS(IGRID) = 1
                NDIVLEPS(IGRID) = NDIVLEPS(IGRID) + 1
   
              ENDIF



         ENDDO

      ENDDO
C                             ITERATIONEN SIND FERTIG
C
C                             RANDBEDINGUNG WIRD NOCH EINMAL AUFGEPRAEGT
C
C                             NOTWENDIG WIEDER IM RED-BLACK VERFAHREN
C
C     IF (NOFVPIT(ILEVEL) .GT. 1) THEN
         DO IRB = 1,2
         DO JRB = 1,2
         DO KRB = 1,2

C
C                                 BEI MESSAGE-PASSING, KOMMUNIKATION
C
         CALL CONNECTMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,UTO,VTO,WTO,PTO,GTO,
     $              VFR,WFR,PFR,GFR,
     $              'P',ILEVEL,IRB,JRB,KRB,
     $              NOFVPIT(ILEVEL),IGRDOFVPIT(1,ILEVEL))



        DO I = 1,NOFVPIT(ILEVEL)
           IGRID = IGRDOFVPIT(I,ILEVEL)


              CALL BOUNDMG
     $                (IDIM3D,IDIM2D,IDIM1D,
     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $                 'P',IGRID,IRB,JRB,KRB,
     $                 TIMEPH,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                 PFR,GFR
     $                )     

         ENDDO
         ENDDO
         ENDDO
         ENDDO
C      ENDIF

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


              CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $                      NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

              CALL BOUNDMG
     $                (IDIM3D,IDIM2D,IDIM1D,
     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $                 'T',IGRID,0,0,0,
     $                 TIMEPH,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                 PFR,GFR
     $                )     

         ENDDO
 

         RETURN
 1000    FORMAT (1X,A,I4,3X,A,I4,3X,A,2X,I4,2X,A,2X,E10.4)
         END
