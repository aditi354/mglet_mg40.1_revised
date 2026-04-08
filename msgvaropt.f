










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
      SUBROUTINE MSGVAROPT (U,V,W,P,G,
     $                    HILF,VAROPT,X,Y,Z,
     $                    IDIM3D,IDIM2D,IDIM1D,IDIMA,
     $                    NBND,IGRID,ITOPT,NB,ITTOT,NTOP1,NTOP2,TIMEPH)
C--MGLET----------------------------------------------------------
C
C     ADDS XRT = XR(TIMEPH) TO BUFFER VAROPT
C     WORKS ONLY IF YHOMOG=.T.
C
C     OUTPUT:
C        VAROPT(ITOPT)
C
C     VERSION 30.07.2000 (J.N.)
C
C-------10--------20--------30--------40--------50--------60--------7072

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )

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


      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

C-------10--------20--------30--------40--------50--------60--------7072
      REAL U(IDIM3D),V(IDIM3D),W(IDIM3D),
     $     P(IDIM3D),G(IDIM3D)
      REAL X(IDIM1D),Y(IDIM1D),Z(IDIM1D)
      REAL VAROPT(2,IDIM2D*8)
      REAL UWALL(IDIM1D,IDIM1D),HILF(IDIM1D)
      REAL XRT,XRTF,XRTL,XRTI
C
C-------10--------20--------30--------40--------50--------60--------7072
C                  IS IDIM3D SUFFICIENT?

      IF (NTOP2-NTOP1 .GT. IDIM2D*2) CALL ERRR (501,"MSGVAROPT")

C
C---------------------------------------------------------------------72
C
           CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $               NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
           CALL MGDIMA  (KKA,JJA,IIA,IGRID)
           CALL MGPOINA (IAV,I2L,I1L,IGRID)

           XPOS1 = XMPOS1(NB,IGRID)
           XPOS2 = XMPOS2(NB,IGRID)
           YPOS1 = YMPOS1(NB,IGRID)
           YPOS2 = YMPOS2(NB,IGRID)
           ZPOS  = ZMPOS1(NB,IGRID)

C
C ANALYSE ON PLANE Z=ZPOS=ZMPOS1 , ZMPOS2 EMPTY
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

           CALL FINDPOSITION(II2,X(IP12),XPOS1,INDX1,2)
           CALL FINDPOSITION(II2,X(IP12),XPOS2,INDX2,2)
           CALL FINDPOSITION(JJ2,Y(IP12),YPOS1,INDY1,2)
           CALL FINDPOSITION(JJ2,Y(IP12),YPOS2,INDY2,2)
           CALL FINDPOSITION(KK2,Z(IP12),ZPOS ,INDZ,2)

           CALL MGOVERLAP (II,X(IP1),II2,X(IP12),NBND, 0.0 ,
     $                  IA1,IE1,IA2,IE2)
           CALL MGOVERLAP (JJ,Y(IP1),JJ2,Y(IP12),NBND, 0.0 ,
     $                  JA1,JE1,JA2,JE2)

          IF  (( IE1 .EQ. II ) .AND. ( IA1 .EQ. 1 )
     $   .AND. ( JE1 .EQ. JJ ) .AND. ( JA1 .EQ. 1 ))
     $    THEN

          INDX1=3
          INDX2=II-2
          INDY1=3
          INDY2=JJ-2

C                                             POSITION AND GRID FOUND!!!
C---------------------------------------------------------------------72

                IP=0
                IP2=0

C---------------------------------------------------------------------72
C---------------------------------------------------------------------72
C                                     GETTING VARIABLE (OPTIMISE) 
C                                     AND PUTTING IT TO FIELD VAROPT
C                                     HERE: U
C---------------------------------------------------------------------72
C                                     1. STEP: AVERAGE SPANWISE


C---------------------------------------------------------------------72
C                                     3. STEP: PROCESS DATA 
C                                              ON MASTER PROC
C                                       3A: AVERAGE SPANWISE BLOCKS
C                                       3B: PUT U(I,K=3) INTO HILF
       IF ( MYID .EQ. 0 ) THEN
         IPRCOUNT = 0
         IGSEARCH = 1
         IPRLAST  = NYSLICE(IGSEARCH)*(NXSLICE(IGSEARCH)-1) + 1
         DO IPR = 1 , IPRLAST , NYSLICE(IGSEARCH)
           DO IPRO = IPR , IPR+NYSLICE(IGSEARCH)-1
             DO IV = 1+NBND , II-NBND
               UWALL(IPR,IV) = 0.5*(UWALL(IPR,IV) + UWALL(IPRO,IV))
             ENDDO
           ENDDO
           DO IW = IPRCOUNT*(II-2*NBND)+1,(IPRCOUNT+1)*(II-2*NBND)
             HILF(IW)=UWALL(IPR,IW-IPRCOUNT*(II-2*NBND)+2)
C            write(6,*)'HILF::',IW,IPR,IPRCOUNT,HILF(IW)
           ENDDO
           IPRCOUNT = IPRCOUNT + 1
         ENDDO

C---------------------------------------------------------------------72
C                                     4. STEP: CREATE VAROPT (3 METHODS)
C                                              COUNTING STARTS AT XPOS1
C                                              COUNTING ENDS AT XPOS2
C                                       4A: FIRST UPWARDS ZEROCROSSING (XRTF)
C                                       4B: LAST UPWARDS ZEROCROSSING (XRTL)
C                                       4C: INTEGRAL ZEROCROSSING (XRTI)
         XRTF=XPOS2
         XRTL=XPOS1
         XRTI=XPOS1

         DO IW=2,IMX(IGSEARCH)-2*NBND
           IF ( ABS(HILF(IW)-HILF(IW-1)) .LT. SMALL ) GOTO 10
           IF ( (X(IW+3).LE.XPOS1) .OR. (X(IW+2).GE.XPOS2) ) GOTO 10
C                                      CASE: 2 TIMES U<0
           IF ( (HILF(IW-1).LT.0.0)  .AND. (HILF(IW).LT.0.0) ) THEN
              IF (X(IW+2).LE.XPOS1) THEN
                XRTI = XRTI + (X(IW+3)-XPOS1)
              ELSE IF (X(IW+3).GE.XPOS2) THEN
                XRTI = XRTI + (XPOS2-X(IW+2))
              ELSE
                XRTI = XRTI + (X(IW+3)-X(IW+2))
              ENDIF
           ENDIF
C                                      CASE: ZEROCROSSING UPWARDS
           IF ( (HILF(IW-1).LT.0.0) .AND. (HILF(IW).GE.0.0) ) THEN
             XRT = X(IW+2) -
     $             HILF(IW-1)* 
     $             ( X(IW+3)-X(IW+2) ) / ( HILF(IW)-HILF(IW-1) )
             IF ( (XRT .LE. XPOS1) .OR. (XRT .GE. XPOS2) ) GOTO 10
             IF ( XRT .LT. XRTF ) XRTF = XRT
             IF ( XRT .GT. XRTL ) XRTL = XRT
             IF ( X(IW+2).GT.XPOS1 ) THEN
               XRTI = XRTI + (XRT-X(IW+2))
             ELSE
               XRTI = XRTI + (XRT-XPOS1)
             ENDIF
           ENDIF
C                                      CASE: ZEROCROSSING DOWNWARDS
           IF ( (HILF(IW-1).GE.0.0) .AND. (HILF(IW).LT.0.0) ) THEN
             XRT = X(IW+2) -
     $             HILF(IW-1)* 
     $             ( X(IW+3)-X(IW+2) ) / ( HILF(IW)-HILF(IW-1) )
             IF ( (XRT .LE. XPOS1) .OR. (XRT .GE. XPOS2) ) GOTO 10
             IF ( X(IW+3).LT.XPOS2 ) THEN
               XRTI = XRTI + (X(IW+3)-XRT)
             ELSE
               XRTI = XRTI + (XPOS2-XRT)
             ENDIF
           ENDIF

   10      CONTINUE
         ENDDO

         VAROPT(1,ITOPT) = TIMEPH
         VAROPT(2,ITOPT) = XRTI
         write(6,*)'XRT ( F , L , I )',TIMEPH,XRTF,XRTL,XRTI

       ENDIF
C---------------------------------------------------------------------72


      ENDIF
      ENDDO

      RETURN
      END
