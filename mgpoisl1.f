










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
       SUBROUTINE MGPOISL1 (IDIM3D,IDIM2D,IDIM1D,NBND,
     $                      X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                      RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                      U,V,W,P,DP,G,B,RHS,RES,UFR,VFR,WFR,
     $                      UTO,VTO,WTO,PTO,GTO,VRI,HILF,
     $                      WSOR,RHO,DIVG,IPRGRID,
     $                      MINCALLVL,ILEVEL,MAXCALLVL,TIMEPH,
     $                      UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                      PFR,GFR,
     $                      GSAW,GSAE,GSAN,GSAS,GSAT,GSAB
     $                      ,BP,BU,BV,BW,SRHS
     $                      ,GSAP,SIPLW,SIPLS,
     $                      SIPLB,SIPLPR,SIPUE,SIPUN,SIPUT
     $                     )
C**MGLET****************************************************************
C
C   M G P O I S L 1      MULTI-GRID-VELOCITY-PRESSURE-LEVEL-1
C                    DRUCKKORREKTURZYKLUS FUER MULTIGRID-VARIANTE
C                    HIER WERDEN GITTER DES LEVELS ILEVEL KORRIGIERT
C                    WURZEL DES DRUCKKORREKTURZYKLUS
C
C     7.08.02. (GT):  VON MGVPL1 ABGELEITET
C    24.02.03. (SE):  PARALLELE VERSION
C    13.08.03  (TB):  EXTENDED FOR SCALAR TRANSPORT
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

      INTEGER IDIM3D,IDIM2D,IDIM1D

      REAL        RDDX(IDIM1D),     RDDY(IDIM1D),     RDDZ(IDIM1D),
     $             RDX(IDIM1D),      RDY(IDIM1D),      RDZ(IDIM1D),
     $             DDX(IDIM1D),      DDY(IDIM1D),      DDZ(IDIM1D),
     $              DX(IDIM1D),       DY(IDIM1D),       DZ(IDIM1D),
     $               X(IDIM1D),        Y(IDIM1D),        Z(IDIM1D)

      REAL        U (   IDIM3D  ), V (   IDIM3D  ), W (   IDIM3D  ),
     $            P (   IDIM3D  ), G (   IDIM3D  ), B (   IDIM3D  ),
     $            DP(   IDIM3D  )


      REAL        UFR(IDIM2D*2),VRI(IDIM2D*2)
      REAL        UBO(IDIM2D*2),VBO(IDIM2D*2),WBO(IDIM2D*2)
      REAL        VFR(IDIM2D*2),WFR(IDIM2D*2)
      REAL        PFR(IDIM2D*2),GFR(IDIM2D*2)
      REAL        UBA(IDIM2D*2),VBA(IDIM2D*2),WBA(IDIM2D*2),
     $            UTO(IDIM2D*2),VTO(IDIM2D*2),WTO(IDIM2D*2),
     $            PTO(IDIM2D*2),GTO(IDIM2D*2)

      REAL     HILF( IDIM3D ), RHS( IDIM3D ), DIVG(IDIM2D)

      REAL MAXRES,MAXMAXRES,HELPRES


      REAL  BP  (IDIM3D),BU  (IDIM3D),BV  (IDIM3D),BW(IDIM3D),
     $      SRHS(IDIM3D)

      REAL  GSAW(IDIM1D),GSAE(IDIM1D),GSAN(IDIM1D),
     $      GSAS(IDIM1D),GSAT(IDIM1D),GSAB(IDIM1D)

      REAL RES   ( IDIM3D ), SIPLB ( IDIM3D ), SIPUT( IDIM3D ),
     $     SIPLW ( IDIM3D ), SIPUE ( IDIM3D ), SIPLS( IDIM3D ),
     $     SIPUN ( IDIM3D ), SIPLPR( IDIM3D ), GSAP ( IDIM3D )

      OMBETA = OMG*RHO/(-2.0*DT)

C                                 VORBELEGUNG DER MAX. DIVERGENZ
C                                 UND DER ANZAHL DER GESAMTEN
C                                 DRUCKKORREKTUREN/DT
      DO IGRID=1,MAXGRIDS
         DIVGMX(IGRID) = 0.0
         IPCGES(IGRID) = 0
         IPCORR(IGRID) = 0
      ENDDO

      DO I = 1,NOFVPIT(ILEVEL)
         IGRID = IGRDOFVPIT( I , ILEVEL )
         LDIVLEPS(IGRID) = 0
         NDIVLEPS(IGRID) = 0
      ENDDO

C
C---- --------------------------------- GESAMTER MULTIGRID-DRUCKKORREKTUR
C---- ---------------------------------ZYKLUS WIRD NUR DURCHLAUFEN, FALLS
C---- --------------------------------- MINDESTENS EIN ITERATIVES GITTER

      IF (NOFVPIT(ILEVEL) .GT. 1) THEN

         CALL CONNECTMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,UTO,VTO,WTO,PTO,GTO,
     $              VFR,WFR,PFR,GFR,
     $              'T',ILEVEL,0,0,0,
     $              NOFVPIT(ILEVEL),IGRDOFVPIT(1,ILEVEL))

      ENDIF

          DO I = 1,NOFVPIT(ILEVEL)
             IGRID = IGRDOFVPIT(I,ILEVEL)
             IF ( .NOT. LSLICE(IGRID)) THEN
          PREFAK = RHO/(WSOR*DT)

              CALL BOUNDMG
     $                (IDIM3D,IDIM2D,IDIM1D,
     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $                 'T',IGRID,0,0,0,TIMEPH,UBA,VBA,WBA,
     $                 PBA,GBA,UBO,VBO,WBO,
     $                 PFR,GFR
     $                ) 

          CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $         NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)


         CALL DIVCAL
     $        (KK,JJ,II,KK,JJ,II,
     $        RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $        U(IP3),V(IP3),W(IP3),BP(IP3),RHS(IP3),
     $        PREFAK,-3,DIVGMX(IGRID),
     $        U(IP3),V(IP3),W(IP3),RHS(IP3))    

         CALL DPHI0 (KK,JJ,II,KK,JJ,II,DP(IP3))


      ENDIF
      ENDDO

      IF (NOFVPIT(ILEVEL) .GT. 1) THEN

         CALL CONNECTMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,RHS,G,B,UFR,VRI,WBO,HILF,UTO,VTO,WTO,PTO,GTO,
     $              VFR,WFR,PFR,GFR,
     $              'S',ILEVEL,0,0,0,
     $              NOFVPIT(ILEVEL),IGRDOFVPIT(1,ILEVEL))

      ENDIF
C---- -------------------------------------------- ALL GRIDS OF LEVEL
      DO 2100 IPCOUNT = 1,MPCORR

C
C                                VORGLAETTUNGEN
C

         CALL MGPOISIT(IDIM3D,IDIM2D,IDIM1D,NBND,
     $        X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $        RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $        U,V,W,P,DP,RHS,RES,RHO,
     $        GSAW,GSAE,GSAN,GSAS,GSAT,GSAB,
     $        G,B,UFR,VFR,WFR,HILF,
     $        UTO,VTO,WTO,PTO,GTO,VRI,ILEVEL,
     $        PFR,GFR,TIMEPH,  
     $        EPCORR,MPCVOR,
     $        UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO
     $        ,BP
     $        ,GSAP,SIPLW,SIPLS,SIPLB,SIPLPR,SIPUE,SIPUN,SIPUT
     $                 )
C---- -------------------------------------------- Beginn Multigrid
C      Multigrid ist nicht implementiert.
C---- -------------------------------------------- End Multigrid

C                                NACHGLÄTTUNGEN

         CALL MGPOISIT(IDIM3D,IDIM2D,IDIM1D,NBND,
     $        X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $        RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $        U,V,W,P,DP,RHS,RES,RHO,
     $        GSAW,GSAE,GSAN,GSAS,GSAT,GSAB,
     $        G,B,UFR,VFR,WFR,HILF,
     $        UTO,VTO,WTO,PTO,GTO,VRI,ILEVEL,
     $        PFR,GFR,TIMEPH,  
     $        EPCORR,MPCNACH,
     $        UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO
     $        ,BP
     $        ,GSAP,SIPLW,SIPLS,SIPLB,SIPLPR,SIPUE,SIPUN,SIPUT
     $                 )
cc         DO I = 1,NOFVPIT(ILEVEL)
cc            IGRID = IGRDOFVPIT(I,ILEVEL)
cc            IF ( .NOT. LSLICE(IGRID)) THEN
cc#ifdef _MPI_
cc               IF (MYID .EQ. IDPROCOFGRD(IGRID))  THEN
cc#endif
C              CALL BOUNDMG
C     $                (IDIM3D,IDIM2D,IDIM1D,
C     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
C     $                 U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
C     $                 UTO,VTO,WTO,PTO,GTO,
C     $                 'T',IGRID,0,0,0,TIMEPH,UBA,VBA,WBA,
C     $                 PBA,GBA,UBO,VBO,WBO,
C     $                 PFR,GFR
C#ifdef _TSCAL_
C     $                ,T,TFR,TBA,TTO,TBO
C#endif
C     $                ) 
cc#ifdef _MPI_
cc              ENDIF
cc#endif
cc           ENDIF
cc         ENDDO

C                   FALLS MAXIMALE DIVERGENZ ALLER GITTER 
C                   DES AKTUELLEN LEVELS KLEINER
C                   ALS DIE SCHRANKE, WIRD NACHITERATION ABGEBROCHEN
C
         MAXRES = 0.0
         HELPRES = 0.0
         DO ILVL = ILEVEL,MINCALLVL,-1
            DO I = 1,NOFVPIT(ILVL)

               IGRID = IGRDOFVPIT( I , ILVL )

                MAXRES = MAX(MAXRES,ABS(DIVGMX(IGRID)))
                HELPRES = DIVGMX(IGRID)

               IF (IVPINF.EQ.1) THEN
                  WRITE (6,1000)'GRD',IGRID,'LEVEL',ILVL,
     &                 'MAXRES',IPCORR(IGRID),HELPRES,EPCORR
               ENDIF
            ENDDO
         ENDDO



c     ABBRUCHKRITERIUM

         IF(MAXRES .LT. EPCORR) GOTO 2200 

 2100 CONTINUE
C                             Druckkorrekturen zuende


 2200 CONTINUE

      DO ILVL = ILEVEL,MINCALLVL,-1
         DO I = 1,NOFTST(ILVL)

            IGRID = IGRDOFTST( I , ILVL )
            IF ( .NOT. LSLICE(IGRID)) THEN

              CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $             NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

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

              CALL BOUNDMG
     $                (IDIM3D,IDIM2D,IDIM1D,
     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $                 'P',IGRID,0,0,0,TIMEPH,UBA,VBA,WBA,
     $                 PBA,GBA,UBO,VBO,WBO,
     $                 PFR,GFR
     $                ) 


            ENDIF
         ENDDO
      ENDDO

      IF (NOFVPIT(ILEVEL) .GT. 1) THEN

         CALL CONNECTMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,UTO,VTO,WTO,PTO,GTO,
     $              VFR,WFR,PFR,GFR,
     $              'T',ILEVEL,0,0,0,
     $              NOFVPIT(ILEVEL),IGRDOFVPIT(1,ILEVEL))

      ENDIF

      DO ILVL = ILEVEL,MINCALLVL,-1
         DO I = 1,NOFTST(ILVL)

            IGRID = IGRDOFTST( I , ILVL )
            IF ( .NOT. LSLICE(IGRID)) THEN

                  CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)


                  CALL DIVCAL
     $                 (KK,JJ,II,KK,JJ,II,
     $                 RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                 U(IP3),V(IP3),W(IP3),BP(IP3),RHS(IP3),
     $                 1.0,-3,DIVGMX(IGRID),
     $                 BU(IP3),BV(IP3),BW(IP3),SRHS(IP3))


           ENDIF
           ENDDO
         ENDDO

      RETURN
 1000 FORMAT (A,2X,I4,3X,A,2X,I4,2X,A,2X,I3,2X,E10.4,2X,E10.4)
      END
