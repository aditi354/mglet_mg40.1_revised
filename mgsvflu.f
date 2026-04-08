










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
      SUBROUTINE MGSVFLU (     
     $                    U,V,W,P,G,AU,AV,AW,
     $                    HILF,UFR,VFR,WFR,PFR,GFR,X,Y,Z,
     $                    IDIM3D,IDIM2D,IDIM1D,IDIMA,
     $                    NBND,IGRID,ITSTEP,NB,NFRO
     $                   )
C--MGLET----------------------------------------------------------
C
C     ADDS FLUCTUATIONS TO FRONT-BUFFERS FOR 
C     "STRUKTUR-PERIODISCHE" BOUNDARY CONDITION
C
C        UFR(X0)=UMEAN + (U(X) - <U>(X))
C
C     10.12.96 (MM.):    ORIGINAL
C     28.04.03 (TB) :    CHANGED FOR SCALAR TRANSPORT BOUNDARY TREATMENT
C
C          TSWITCH = 0 :  FLUCTUATION BC APPLIED TO U,V,W,(G) ONLY
C          TSWITCH = 1 :  FLUCTUATION BC APPLIED TO T ONLY
C          TSWITCH = 0 :  FLUCTUATION BC APPLIED TO ALL VARS U,V,W,(G),T
C
C-------10--------20--------30--------40--------50--------60--------7072

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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

C-------10--------20--------30--------40--------50--------60--------7072
      REAL U(IDIM3D),V(IDIM3D),W(IDIM3D),
     $     P(IDIM3D),G(IDIM3D),HILF(IDIM2D,6)
      REAL AU(IDIMA),AV(IDIMA),AW(IDIMA)
      REAL UFR(IDIM2D*2), VFR(IDIM2D*2), WFR(IDIM2D*2),
     $     PFR(IDIM2D*2), GFR(IDIM2D*2)
      REAL X(IDIM1D),Y(IDIM1D),Z(IDIM1D)

C
C-------10--------20--------30--------40--------50--------60--------7072
C                  IS IDIM3D SUFFICIENT?

      IF (IDIM3D .LT. IDIM2D*6) CALL ERRR (501,"MGSVFLU")

C
C---------------------------------------------------------------------72
C                      IS NFRO RIGHT?
      IF (NFRO .NE. 11  .AND. NFRO .NE. 12) CALL ERRR (502,"MGSVFLU")
C 
C-------10--------20--------30--------40--------50--------60--------7072
C
           CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $               NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
           CALL MGDIMA  (KKA,JJA,IIA,IGRID)
           CALL MGPOINA (IAV,I2L,I1L,IGRID)

           XPOS = XBANF(NB,1,IGRID)
C
C-------10--------20--------30--------40--------50--------60--------7072
C                              BEIM ERSTEN ZEITSCHRITT SONDERBEHANDLUNG
           IF (ITSTEP .EQ. 0) THEN

                IP=0
            CALL COPYREGION(KK,JJ,II,U(IP3),
     $           IP,1,KK,1,JJ,2,2,
     $           KK,JJ,2,UFR(IBB),IP,
     $           0,0,0,
     $           HILF(1,5),MYID)

            CALL COPYREGION(KK,JJ,II,V(IP3),
     $           IP,1,KK,1,JJ,2,2,
     $           KK,JJ,2,VFR(IBB),IP,
     $           0,0,0,
     $           HILF(1,5),MYID)

            CALL COPYREGION(KK,JJ,II,W(IP3),
     $           IP,1,KK,1,JJ,2,2,
     $           KK,JJ,2,WFR(IBB),IP,
     $           0,0,0,
     $           HILF(1,5),MYID)


              RETURN
           ENDIF
C
C-------10--------20--------30--------40--------50--------60--------7072
C                         SUCHE NACH DER POSITION

      ILEVEL = LEVEL(IGRID)
        DO I = 1,NOFTST(ILEVEL)
           IGSEARCH = IGRDOFTST(I,ILEVEL)

           CALL MGDIMS  (KK2,JJ2,II2,IGSEARCH)          
           CALL MGPOINT (IP32,IP22,IP12,IBB2,IB32,IBU2,IGSEARCH)
           CALL MGDIMA  (KKA2,JJA2,IIA2,IGSEARCH)
           CALL MGPOINA (IAV2,I2L2,I1L2,IGSEARCH)


           CALL FINDPOSITION(II2,X(IP12),XPOS,INDEX,2)

           IF (INDEX .GT. 0 .AND. INDEX .LT. II) THEN

C-------10--------20--------30--------40--------50--------60--------7072

              CALL MGOVERLAP (JJ,Y(IP1),JJ2,Y(IP12),NBND, 0.0 ,
     $                  JA1,JE1,JA2,JE2)     
              CALL MGOVERLAP (KK,Z(IP1),KK2,Z(IP12),NBND, 0.0 ,
     $                  KA1,KE1,KA2,KE2)   

          IF  (( JE1 .EQ. JJ ) .AND. ( JA1 .EQ. 1 )
     $   .AND. ( KE1 .EQ. KK ) .AND. ( KA1 .EQ. 1 ))
     $    THEN
C                                             POSITION AND GRID FOUND!!!
C---------------------------------------------------------------------72

                IP=0
                IP2=0

C
C---------------------------------------------------------------------72
C                   IF BOUND.-COND. IS "RECYCLE" THEN ANALYSE AVERAGE 
C                   FIELD FIRST
C
             IF (NFRO .EQ. 12) THEN
                CALL COPYREGION(KKA2,JJA2,IIA2,AU(IAV2),
     $               IP2,1,KKA2,1,JJA2,INDEX,INDEX,
     $               KKA2,JJA2,1,HILF(1,1),IP,
     $               0,0,1-INDEX,
     $               HILF(1,5),MYID)
             
                CALL BL_DIAGNOSE (KKA2,JJA2,1,NBND,HILF(1,1),Z(IP1),
     $               ZBANF(NB,1,IGRID),
     $               DELTA_ACT,UTAU_ACT,UMAX_ACT)

	write (6,*) 'delta...',DELTA_ACT,UTAU_ACT,UMAX_ACT

C
             ENDIF
C
C                   ANALYSIS OF AVERAGE FIELD FINISHED
C---------------------------------------------------------------------72
C                                   GETTING INFLOW PROFILE INTO HILF
C                                   AND PUTTING IT (WITH APPROPRIATE
C                                   SCALING)ONTO UFR, VFR, WFR, GFR, TFR
C
C---------------------------------------------------------------------72
C                                     HERE: U
             CALL COPYREGION(KK2,JJ2,II2,U(IP32),
     $            IP2,1,KK2,1,JJ2,INDEX,INDEX,
     $            KK,JJ,1,HILF(1,1),IP,
     $            0,0,1-INDEX,
     $            HILF(1,2),MYID)


                IF (NFRO .EQ. 11) THEN
                   CALL SVFLU(KK,JJ, 2,UFR(IBB),HILF(1,1),HILF(1,2),
     $                  HILF(1,3),JJA2,
     $                  NBND,Z(IP1),ZBANF(NB,1,IGRID),DELTA)
                ELSEIF (NFRO .EQ. 12) THEN
                   CALL SVREC(KK,JJ, 2,2,UFR(IBB),HILF(1,1),HILF(1,2),
     $                  HILF(1,3),HILF(1,4),Z(IP1),ZBANF(NB,1,IGRID),
     $                  DELTA,DELTA_ACT,UTAUX,UTAU_ACT,UFRCON,UMAX_ACT,
     $                  GMOL)
                ENDIF

             
C     
C---------------------------------------------------------------------72
C     HERE: V

             CALL COPYREGION(KK2,JJ2,II2,V(IP32),
     $            IP2,1,KK2,1,JJ2,INDEX,INDEX,
     $            KK,JJ,1,HILF(1,1),IP,
     $            0,0,1-INDEX,
     $            HILF(1,2),MYID)
             
                IF (NFRO .EQ. 11) THEN
                   CALL SVFLU(KK,JJ, 2,VFR(IBB),HILF(1,1),HILF(1,2),
     $                  HILF(1,3),JJA2,
     $                  NBND,Z(IP1),ZBANF(NB,1,IGRID),DELTA)
                ELSEIF (NFRO .EQ. 12) THEN
                   CALL SVREC(KK,JJ, 2,2,VFR(IBB),HILF(1,1),HILF(1,2),
     $                  HILF(1,3),HILF(1,4),Z(IP1),ZBANF(NB,1,IGRID),
     $                  DELTA,DELTA_ACT,UTAUX,UTAU_ACT,UFRCON,UMAX_ACT,
     $                  GMOL)
                ENDIF
C
C---------------------------------------------------------------------72
C                                     HERE: W

             CALL COPYREGION(KK2,JJ2,II2,W(IP32),
     $            IP2,1,KK2,1,JJ2,INDEX,INDEX,
     $            KK,JJ,1,HILF(1,1),IP,
     $            0,0,1-INDEX,
     $            HILF(1,2),MYID)

                IF (NFRO .EQ. 11) THEN
                   CALL SVFLU(KK,JJ, 2,WFR(IBB),HILF(1,1),HILF(1,2),
     $                  HILF(1,3),JJA2,
     $                  NBND,Z(IP1),ZBANF(NB,1,IGRID),DELTA)
                ELSEIF (NFRO .EQ. 12) THEN
                   CALL SVREC(KK,JJ, 2,2,WFR(IBB),HILF(1,1),HILF(1,2),
     $                  HILF(1,3),HILF(1,4),Z(IP1),ZBANF(NB,1,IGRID),
     $                  DELTA,DELTA_ACT,UTAUX,UTAU_ACT,UFRCON,UMAX_ACT,
     $                  GMOL)
                ENDIF
C
C


          GOTO 1000
       ENDIF
      ENDIF
      ENDDO
     
 1000 CONTINUE
C-------10--------20--------30--------40--------50--------60--------7072
      RETURN
      END
