










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
       SUBROUTINE MGPOISC1 (IDIM3D,IDIM2D,IDIM1D,NBND,
     $                      X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                      RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                      FPSFAK,FCOSMY,FCOSNY,
     $                      U,V,W,P,DP,G,B,RES,RHS,UFR,VFR,WFR,
     $                      UTO,VTO,WTO,PTO,GTO,VRI,HILF,
     $                      H2D1,H2D2,H2D3,H2D4,H2D5,H2D6,
     $                      WSOR, RHO,DIVG,IPRGRID,
     $                      MINCALLVL,ILEVEL,MAXCALLVL,TIMEPH,
     $                      UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                      PFR,GFR,
     $                      GSAW,GSAE,GSAN,GSAS,GSAT,GSAB
     $                      ,GSAP,SIPLW,SIPLS,
     $                      SIPLB,SIPLPR,SIPUE,SIPUN,SIPUT
     $                     )
c      IMPLICIT NONE
C**MGLET****************************************************************
C
C   M G V P C 1      MULTI-GRID-VELOCITY-PRESSURE-LEVEL-2
C                    DRUCKKORREKTURZYKLUS FUER MULTIGRID-VARIANTE
C                    HIER WERDEN GITTER DES LEVELS ILEVEL KORRIGIERT
C                    WURZEL DES DRUCKKORREKTURZYKLUS
C
C     4.10.93 (MM):  ORIGINAL
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
     $            DP(   IDIM3D  )


      REAL        UFR(IDIM2D*2),VRI(IDIM2D*2)
      REAL        UBO(IDIM2D*2),VBO(IDIM2D*2),WBO(IDIM2D*2)
      REAL        VFR(IDIM2D*2),WFR(IDIM2D*2)
      REAL        PFR(IDIM2D*2),GFR(IDIM2D*2)
      REAL        UBA(IDIM2D*2),VBA(IDIM2D*2),WBA(IDIM2D*2)

C
C                                 FELDER FUER DEN DIREKTEN POISSONLOESER
C
      REAL
     $           FPSFAK(IDIM2D),FCOSMY(IDIM1D),FCOSNY(IDIM1D)
      REAL     H2D1( IDIM2D ), H2D2( IDIM2D ), H2D3( IDIM2D ),
     $         H2D4( IDIM2D ), H2D5( IDIM2D ), H2D6( IDIM2D )
C
C                                 FELDER FUER DIE FOURIERTRANSFORMATION
C
c      INTEGER IFFTX(19,MAXGRIDS),IFFTY(19,MAXGRIDS)
c      INTEGER IPERMUX(IDIM1D)	     ,IPERMUY(IDIM1D)
c      REAL    RFFTX(IDIM1D*2)     ,RFFTY(IDIM1D*2)
C
C
      REAL     HILF( IDIM3D ), DIV( IDIM3D ),DIVG(IDIM2D)
C
      INTEGER I,J,K
      REAL PREFAK 

      REAL GSAW(IDIM3D),GSAE(IDIM3D),GSAN(IDIM3D),
     $     GSAS(IDIM3D),GSAT(IDIM3D),GSAB(IDIM3D)
      REAL RHS( IDIM3D )
      REAL RES(IDIM3D )
      REAL SIPLB ( IDIM3D ), SIPUT( IDIM3D )
      REAL SIPLW ( IDIM3D ), SIPUE( IDIM3D )
      REAL SIPLS ( IDIM3D ), SIPUN( IDIM3D )
      REAL SIPLPR( IDIM3D ), GSAP ( IDIM3D )

C
C
c         OMBETA = OMG*RHO/(-2.0*DT)
c         write(6,*)'SIPNJK C1', SIPNJK(2),SIPNJK(1),KK,JJ,II
C
C

      wsor=2.0
         PREFAK = RHO/(WSOR*DT)
         EPC1 = EPCORR / (EPFAK*1.0)
      MPC1 = MPCORR
      WRITE(6,*)'MPC1',MPC1,EPC1
      DO I = 1,NOFVPIT(ILEVEL)
         IGRID = IGRDOFVPIT( I , ILEVEL )
         IPCORR(IGRID) = 0
         LDIVLEPS(IGRID) = 0
         NDIVLEPS(IGRID) = 0
      ENDDO

          CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $         NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

              CALL BOUNDMG
     $                (IDIM3D,IDIM2D,IDIM1D,
     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,DP,G,B,UFR,VFR,WFR,VRI,HILF,
     $                 UTO,VTO,WTO,PTO,GTO,
     $                 'T',IGRID,0,0,0,TIMEPH,UBA,VBA,WBA,
     $                 PBA,GBA,UBO,VBO,WBO,
     $                 PFR,GFR
     $                ) 
          

c         write(6,*)'C1- KK,JJ,II',KK,JJ,II
c          CALL DIVCAL
c     $         (KK,JJ,II,KK,JJ,II,
c     $         RDDX(IP1),RDDY(IP1),RDDZ(IP1),
c     $         U(IP3),V(IP3),W(IP3),B(IP3),RHS(IP3),
c     $         PREFAK,-3,DIVGMX(IGRID))
c          write(6,*)'MAX. DIVERGENZ vor korrektur C1',DIVGMX(IGRID)
c-------------------------------------------------------
C
C---- --------------------------------- GESAMTER MULTIGRID-DRUCKKORREKTUR
C---- ---------------------------------ZYKLUS WIRD NUR DURCHLAUFEN, FALLS
C---- --------------------------------- MINDESTENS EIN ITERATIVES GITTER
CCCCCC      IF (NOFVPIT(ILEVEL) .GE. 1) THEN

c      DO 2100 IPCOUNT = 1,MPCORR
      DO 2100 IPCOUNT = 1,MPC1
C
C                                   UEBERPRUEFEN DES ABBRUCHKRITERIUMS
C

C---- -------------------------------------------- ALL GRIDS OF LEVEL 

C
C                                VORGLAETTUNGEN
C


         CALL MGPOISIT      (IDIM3D,IDIM2D,IDIM1D,NBND,
     $                      X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                      RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                      U,V,W,P,DP,RHS,RES, RHO,
     $                      GSAW,GSAE,GSAN,GSAS,GSAT,GSAB,
     $                      G,B,UFR,VFR,WFR,HILF,
     $                      UTO,VTO,WTO,PTO,GTO,VRI,ILEVEL,
     $                      PFR,GFR,TIMEPH,  
     $                      EPCORR,MPCVOR,
     $                      UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO
     $                      ,GSAP,SIPLW,SIPLS,
     $                      SIPLB,SIPLPR,SIPUE,SIPUN,SIPUT
     $                      )

         CALL MGPOISIT      (IDIM3D,IDIM2D,IDIM1D,NBND,
     $                      X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                      RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                      U,V,W,P,DP,RHS,RES, RHO,
     $                      GSAW,GSAE,GSAN,GSAS,GSAT,GSAB,
     $                      G,B,UFR,VFR,WFR,HILF,
     $                      UTO,VTO,WTO,PTO,GTO,VRI,ILEVEL,
     $                      PFR,GFR,TIMEPH,  
     $                      EPCORR,MPCNACH,
     $                      UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO
     $                      ,GSAP,SIPLW,SIPLS,
     $                      SIPLB,SIPLPR,SIPUE,SIPUN,SIPUT
     $                      )
C

C                                   UEBERPRUEFEN DES ABBRUCHKRITERIUMS
C
C                        MAX. RESIDUUM LIEGT AUF FELD DIVGMX
C
C---- -------------------------------------------- ALL GRIDS OF LEVEL 

        DO I = 1,NOFVPIT(ILEVEL)
           IGRID = IGRDOFVPIT(I,ILEVEL)

           IF ( .NOT. LSLICE(IGRID)) THEN

c                  write(6,*)'RESMAX ',IGRID,DIVGMX(IGRID)

           ENDIF

         ENDDO

         RESMAX = 0.0

         DO I = 1,NOFVPIT(ILEVEL)

           IGRID = IGRDOFVPIT( I , ILEVEL )

           RESMAX = MAX(RESMAX,ABS(DIVGMX(IGRID)))

         ENDDO
         IF(RESMAX .LT. EPC1) GOTO 2200
 2100    CONTINUE
C                             Druckkorrekturen zuende


 2200    CONTINUE

C              CALL BOUNDMG
C     $                (IDIM3D,IDIM2D,IDIM1D,
C     $                 X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
C     $                 U,V,W,DP,G,B,UFR,VFR,WFR,VRI,HILF,
C     $                 UTO,VTO,WTO,PTO,GTO,
C     $                 'T',IGRID,0,0,0,TIMEPH,UBA,VBA,WBA,
C     $                 PBA,GBA,UBO,VBO,WBO,
C     $                 PFR,GFR
C#ifdef _TSCAL_
C     $                ,T,TFR,TBA,TTO,TBO
C#endif
C     $                )
 
         RETURN
 1000    FORMAT (1X,A,I4,3X,A,2X,I4,3X,A,2X,I4,2X,A,2X,E10.4)
         END
