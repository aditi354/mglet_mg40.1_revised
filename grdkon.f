










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
      SUBROUTINE GRDKON (IGRID,IFAIL,IIDENT)
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                    PRUEFT, OB DIE IN COGRDDEF LIEGENDEN
C                    GITTERDEFINITIONEN MIT DEN IN COGRDPRO
C                    LIEGENDEN DES VORANGEGANGENEN LAUFES 
C                    IDENTISCH SIND
C                    HIER WIRD DAS MIT IGRID SPEZIFIZIERTE GITTER 
C                    GEPRUEFT
C
C     21. 5.93 (MM) : ORIGINAL
C     07.07.03 (TB) : EXTENSION IN I,J IMPLEMENTED
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
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
 

      
      COMMON /COGRDOLD/
     &                 KMXOLD,  JMXOLD,  IMXOLD,
     &                 LEVELOLD,    LCHILDOLD, IPARENTOLD,
     &                 IPOSOLD, JPOSOLD, KPOSOLD,
     &                 XMINOLD,    YMINOLD,    ZMINOLD,
     &                 XTOTOLD,     YTOTOLD,     ZTOTOLD,   
     &                 NBOCOLD,
     &                 FRONTOLD,    BACKOLD,     RIGHTOLD,     LEFTOLD,
     &                 BOTTOMOLD,   TOPOLD,      CUBEOLD,
     &                 IFRNBROLD, IBANBROLD, IRINBROLD, ILENBROLD,
     &                 IBONBROLD, ITONBROLD,
     &                 IBPOSOLD,   JBPOSOLD,    KBPOSOLD,
     &                 IBANFOLD,  IBENDOLD,
     &                 JBANFOLD,  JBENDOLD,
     &                 KBANFOLD,  KBENDOLD,
     &                 XHOMOGOLD,   YHOMOGOLD,   ZHOMOGOLD,
     &                 LPLEVELOLD

      INTEGER
     &       KMXOLD(MAXGRIDS), JMXOLD(MAXGRIDS), IMXOLD(MAXGRIDS),
     &       LEVELOLD(MAXGRIDS), 
     &       IPARENTOLD(MAXGRIDS),
     & IPOSOLD(MAXGRIDS), JPOSOLD(MAXGRIDS), KPOSOLD(MAXGRIDS),

     &  NBOCOLD(7,MAXGRIDS),
     &  IFRNBROLD(MAXBOCONDS,MAXGRIDS), IBANBROLD(MAXBOCONDS,MAXGRIDS),
     &  IRINBROLD(MAXBOCONDS,MAXGRIDS), ILENBROLD(MAXBOCONDS,MAXGRIDS),
     &  IBONBROLD(MAXBOCONDS,MAXGRIDS), ITONBROLD(MAXBOCONDS,MAXGRIDS),
     &  IBPOSOLD(MAXBOCONDS,7,MAXGRIDS),JBPOSOLD(MAXBOCONDS,7,MAXGRIDS),
     &  KBPOSOLD(MAXBOCONDS,7,MAXGRIDS),
     &  IBANFOLD(MAXBOCONDS,7,MAXGRIDS),IBENDOLD(MAXBOCONDS,7,MAXGRIDS),
     &  JBANFOLD(MAXBOCONDS,7,MAXGRIDS),JBENDOLD(MAXBOCONDS,7,MAXGRIDS),
     &  KBANFOLD(MAXBOCONDS,7,MAXGRIDS),KBENDOLD(MAXBOCONDS,7,MAXGRIDS)


      REAL
     &       XTOTOLD(MAXGRIDS), YTOTOLD(MAXGRIDS), ZTOTOLD(MAXGRIDS),
     &       XMINOLD(MAXGRIDS), YMINOLD(MAXGRIDS), ZMINOLD(MAXGRIDS)

      CHARACTER (LEN=16)
     &   FRONTOLD(MAXBOCONDS,MAXGRIDS),   BACKOLD(MAXBOCONDS,MAXGRIDS),
     &   RIGHTOLD(MAXBOCONDS,MAXGRIDS),   LEFTOLD(MAXBOCONDS,MAXGRIDS),
     &  BOTTOMOLD(MAXBOCONDS,MAXGRIDS),    TOPOLD(MAXBOCONDS,MAXGRIDS),
     &    CUBEOLD(MAXBOCONDS,MAXGRIDS)

      LOGICAL
     &      LCHILDOLD(MAXGRIDS),
     &  XHOMOGOLD(MAXGRIDS), YHOMOGOLD(MAXGRIDS), ZHOMOGOLD(MAXGRIDS),
     &      LPLEVELOLD(MAXGRIDS)
      
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


      INTEGER IIDENT(100)
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      IFAIL = 0

C
C                           UM BUCHSTABEN ZU SPAREN !
      I = IGRID 

      IF (KMX(I).NE.KMXOLD(I)) THEN
         IFAIL = IFAIL + 1
         WRITE (6,6002) 'GRD',IGRID,'    KMX   FAILS:'
      ENDIF
      IF (JMX(I).NE.JMXOLD(I)) THEN
         IFAIL = IFAIL + 1
         WRITE (6,6002) 'GRD',IGRID,'    JMX   FAILS:'
      ENDIF
      IF (IMX(I).NE.IMXOLD(I)) THEN
         IFAIL = IFAIL + 1
         WRITE (6,6002) 'GRD',IGRID,'    IMX   FAILS:'
      ENDIF

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC FALLS ERGEBNIS VON MLET 
C                                         EINGELESEN WIRD, KOENNEN
C                                         FOLGENDE UEBERPRUEFUNGEN
C                                         NICHT GEMACHT WERDEN

      IF (IIDENT(5).GT.8) THEN

        IF (LEVEL(I).NE.LEVELOLD(I)) THEN
           IFAIL = IFAIL + 1
           WRITE (6,6002) 'GRD',IGRID,'    LEVEL   FAILS:'
        ENDIF
        IF (LCHILD(I) .AND.(.NOT.LCHILDOLD(I))) THEN
           IFAIL = IFAIL + 1
           WRITE (6,6002) 'GRD',IGRID,'    LCHILD   FAILS:'
        ENDIF
        IF (LCHILD(I)) THEN
          IF (IPARENT(I).NE.IPARENTOLD(I)) THEN
             IFAIL = IFAIL + 1
             WRITE (6,6002) 'GRD',IGRID,'    IPARENT   FAILS:'
          ENDIF
          IF (IPOSITION(I).NE.IPOSOLD(I)) THEN
             IFAIL = IFAIL + 1
             WRITE (6,6002) 'GRD',IGRID,'    IPOSITION   FAILS:'
          ENDIF
          IF (JPOSITION(I).NE.JPOSOLD(I)) THEN
             IFAIL = IFAIL + 1
             WRITE (6,6002) 'GRD',IGRID,'    JPOSITION   FAILS:'
          ENDIF
          IF (KPOSITION(I).NE.KPOSOLD(I)) THEN
             IFAIL = IFAIL + 1
             WRITE (6,6002) 'GRD',IGRID,'    KPOSITION   FAILS:'
          ENDIF
        ENDIF
C        IF (XTOT(I).NE.XTOTOLD(I)) THEN
C           IFAIL = IFAIL + 1
C           WRITE (6,6002) 'GRD',IGRID,'    XTOT   FAILS:'
C        ENDIF
C        IF (YTOT(I).NE.YTOTOLD(I)) THEN
C           IFAIL = IFAIL + 1
C           WRITE (6,6002) 'GRD',IGRID,'    YTOT   FAILS:'
C        ENDIF
C        IF (ZTOT(I).NE.ZTOTOLD(I)) THEN
C           IFAIL = IFAIL + 1
C           WRITE (6,6002) 'GRD',IGRID,'    ZTOT   FAILS:'
C        ENDIF

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC FALLS ERGEBNIS VON MLET
C                                             KEIN FORTSETZUNGSLAUF
      ELSE
        DCONT = .FALSE.
      ENDIF
C
C                                         RANDBEDINGUNGEN KOENNEN BEI
C                                         DCONT=.FALSE. GEAENDERT
C                                         WERDEN
C
      IF (DCONT) THEN

C                                         AB VERSION 34 STEHEN MEHR
C                                         ALS EINE RANDBED. ZUR VERFUEGUNG

         IF (IIDENT(5).GE.34) THEN

          DO IDIR = 1,7
               WRITE (6,*) 'GRD',IGRID,'NBOCONDS:'
     &         ,NBOCOLD(IDIR,I),'FRAGE',NBOCD(IDIR,I)
            IF (NBOCOLD(IDIR,I) .NE. NBOCD(IDIR,I)) THEN
               IFAIL = IFAIL + 1
                  WRITE (6,6002) 'GRD',IGRID,'    NBOCONDS  FAIL:'
     &            ,NBOCOLD(IDIR,I),NBOCD(IDIR,I)
            ENDIF
          ENDDO

         ELSE

          DO IDIR = 1,7
             NBOCOLD(IDIR,I) = 1
          ENDDO

         ENDIF

         DO IB = 1,NBOCD( 1 , I )
            IF (FRONT(IB,I).NE.FRONTOLD(IB,I)) THEN
               IFAIL = IFAIL + 1
               WRITE (6,6002) 'GRD',IGRID,'    FRONT   FAILS:'
            ENDIF
         ENDDO

         DO IB = 1,NBOCD( 2 , I )
            IF (BACK(IB,I).NE.BACKOLD(IB,I)) THEN
               IFAIL = IFAIL + 1
               WRITE (6,6002) 'GRD',IGRID,'    BACK   FAILS:'
            ENDIF
         ENDDO

         DO IB = 1,NBOCD( 3 , I )
            IF (RIGHT(IB,I).NE.RIGHTOLD(IB,I)) THEN
               IFAIL = IFAIL + 1
               WRITE (6,6002) 'GRD',IGRID,'    RIGHT   FAILS:'
            ENDIF
         ENDDO

         DO IB = 1,NBOCD( 4 , I )
            IF (LEFT(IB,I).NE.LEFTOLD(IB,I)) THEN
               IFAIL = IFAIL + 1
               WRITE (6,6002) 'GRD',IGRID,'    LEFT   FAILS:'
            ENDIF
         ENDDO

         DO IB = 1,NBOCD( 5 , I )
            IF (BOTTOM(IB,I).NE.BOTTOMOLD(IB,I)) THEN
               IFAIL = IFAIL + 1
               WRITE (6,6002) 'GRD',IGRID,'    BOTTOM   FAILS:'
            ENDIF
         ENDDO

         DO IB = 1,NBOCD( 6 , I )
            IF (TOP(IB,I).NE.TOPOLD(IB,I)) THEN
               IFAIL = IFAIL + 1
               WRITE (6,6002) 'GRD',IGRID,'    TOP   FAILS:'
            ENDIF
         ENDDO

         DO IB = 1,NBOCD( 7 , I )
            IF (CUBE(IB,I).NE.CUBEOLD(IB,I)) THEN
               IFAIL = IFAIL + 1
               WRITE (6,6002) 'GRD',IGRID,'    CUBE   FAILS:'
            ENDIF
         ENDDO

      ENDIF
C
C                                            FOLGENDE EINSTELLUNGEN MUESSEN
C                                            IMMER GLEICH BLEIBEN
C                                            AUSSER BEI ERG. AUS MLET
      IF (IIDENT(5) .GT. 8) THEN
C
        IF (XHOMOG(I).AND.(.NOT.XHOMOGOLD(I))) THEN
           IFAIL = IFAIL + 1
           WRITE (6,6002) 'GRD',IGRID,'    XHOMOG   FAILS:'
        ENDIF
        IF (YHOMOG(I).AND.(.NOT.YHOMOGOLD(I))) THEN
           IFAIL = IFAIL + 1
           WRITE (6,6002) 'GRD',IGRID,'    YHOMOG   FAILS:'
        ENDIF
        IF (ZHOMOG(I).AND.(.NOT.ZHOMOGOLD(I))) THEN
           IFAIL = IFAIL + 1
           WRITE (6,6002) 'GRD',IGRID,'    ZHOMOG   FAILS:'
        ENDIF
C       IF (LPLEVEL(I).AND.(.NOT.LPLEVELOLD(I))) THEN
C          IFAIL = IFAIL + 1
C          WRITE (6,6002) 'GRD',IGRID,'    LPLEVEL   FAILS:'
C       ENDIF

      ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
 6000 FORMAT(1X,A3,I3,A10,I3,A10,I3,A10,I3)
 6001 FORMAT(1X,A3,I3,A10,F8.3,A10,F8.3,A10,F8.3)
 6002 FORMAT(1X,A3,I3,A30,I8)
 6003 FORMAT(1X,A3,I3,A10,A10,A10,A10)

      RETURN
      END
