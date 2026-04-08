










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
      SUBROUTINE BLOCKBP (
     $             DDX,DDY,DDZ,
     $             DX,DY,DZ,
     $             X,Y,Z,
     $             BP, HILF,
     $     IDIM3D,IDIM2D,IDIM1D)

C*MGLET*****************************************************************
C   B L O C K B P    BELEGT DAS KOERPERERKENNUNGSFELD BEI FREDBODY
C                    (BP-FELD)
C*MGLET*****************************************************************
C
C  KREUZINGER UND MANHART TURBULENZ GmbH
C  JK      11. 8.2005
C
C  PARAMETER
C             X              - KOORDINATEN IN X-RICHTUNG
C             Y              - KOORDINATEN IN Y-RICHTUNG
C             Z              - KOORDINATEN IN Z-RICHTUNG
C             DDX            - X-KANTENLAENGE DER KONTROLLVOLUMINA
C             DDY            - Y-KANTENLAENGE DER KONTROLLVOLUMINA
C             DDZ            - Z-KANTENLAENGE DER KONTROLLVOLUMINA
C             DX             - X-ABSTAND DER GITTERPUNKTE
C             DY             - Y-ABSTAND DER GITTERPUNKTE
C             BP             - BLOCKINGFELD FUER DEN DRUCK
C             HILF           - HILFSFELD          
C             IDIM3D,IDIM2D,IDIM1D - DIMENSIONIERUNGEN
C
C***********************************************************************


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
 
      INTEGER     STARTGRID

      REAL        X(IDIM1D),         Y(IDIM1D),         Z(IDIM1D),
     $           DX(IDIM1D),        DY(IDIM1D),        DZ(IDIM1D),
     $          DDX(IDIM1D),       DDY(IDIM1D),       DDZ(IDIM1D),
     $          BP( IDIM3D ),   HILF( IDIM3D ),
     $          NSINGLE, NALL
C

C SCHLEIFE UEBER ALLE NEUEN GITTER, BZW. DEREN SLICES
C fs  BEI NEWBLOCK=TRUE schleife auch über altes gitter
C      IF(CTYP(1)(1:8).EQ.'USEBLOCK') STARTGRID = NGRDOLD
c      IF(CTYP(1)(1:8).EQ.'NEWBLOCK') STARTGRID = 0 
C fs
 
c      DO IGRID=STARTGRID+1,NGRDDFD
      DO IGRID=NGRDOLD+1,NGRDDFD
      DO ISUB = 1,MAX(1,NOFSLCHILDS(IGRID))
         IF (LSLICE(IGRID)) THEN
            IGRIDNUM = IGRDOFSLCHILD(ISUB,IGRID)
         ELSE
            IGRIDNUM = IGRID
         ENDIF
         MYID   = 0
         IPROC  = 0
         IPPAR  = 0
         IF(MYID .EQ. IPROC) THEN
            write(*,*) 'BLOCKBP: FIND CELLS INTERSECTED BY TRIANGLES'
            write(*,*) 'BLOCKBP: IGRIDNUM',IGRIDNUM
            CALL MGPOINT (IP3,IP2,IP1,IBB,IB3,IBU,IGRIDNUM)
            CALL MGDIMS  (KK,JJ,II,IGRIDNUM)
            CALL MGBASB(NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRIDNUM)
C  hier werden die einzelnen Gitter bearbeitet, ANFANG               
C           Feld umspeichern : KJI -> IJK
            DO K=1,KK
            DO J=1,JJ
            DO I=1,II
               HILF(IP3-1 + I + II*(J-1) + II*JJ*(K-1))
     $              = BP(IP3-1 + K + KK*(J-1) + KK*JJ*(I-1))
            ENDDO
            ENDDO
            ENDDO
C     call bp fuer einen Prozessor (BP(IP3), Koordinaten(IP1), KK,JJ,II)
            CALL BLOCKGRID(II,JJ,KK,0,X(IP1),Y(IP1),Z(IP1),
     $           DX(IP1),DY(IP1),DZ(IP1),
     $           DDX(IP1),DDY(IP1),DDZ(IP1),
     $           IGRIDNUM,NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,
     $           HILF(IP3),-1)

C           Feld umspeichern : IJK -> KJI
            DO K=1,KK
            DO J=1,JJ
            DO I=1,II
               BP(IP3-1 + K + KK*(J-1) + KK*JJ*(I-1))
     $              = HILF(IP3-1 + I + II*(J-1) + II*JJ*(K-1))
            ENDDO
            ENDDO
            ENDDO

            rmeanbp = 0.0
            do I = 1,II*JJ*KK
               rmeanbp = rmeanbp + abs(BP(IP3-1+I))
            enddo
            rmeanbp = rmeanbp/(II*JJ*KK)
            write(*,*) 'BLOCKBP rmeanbp A',rmeanbp,IGRIDNUM
C  einzelne Gitter, ENDE
         ENDIF
      ENDDO
      ENDDO              

      IF (NGRDOLD.LT.NGRDDFD) THEN
C      IF(STARTGRID.LT.NGRDDFD) THEN

      DO ILOOP = 1,10
C tausche BP Randschichten aus, damit Fuellalgorithmus ueber 
C Gittergrenzen wirken kann 
         write(*,*) 'BLOCKBP: ADDITIONAL BLOCKING'
         write(*,*) 'BLOCKBP: LOOP',ILOOP
         DO ILEVEL = MINLEVEL,MAXLEVEL

            CALL CONFROMGBP 
     $           (IDIM3D,IDIM2D,IDIM1D,
     $           X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $           BP,
     $           NOFTST(ILEVEL),IGRDOFTST(1,ILEVEL))
            CALL CONBACMGBP 
     $           (IDIM3D,IDIM2D,IDIM1D,
     $           X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $           BP,
     $           NOFTST(ILEVEL),IGRDOFTST(1,ILEVEL))
            CALL CONRGTMGBP 
     $           (IDIM3D,IDIM2D,IDIM1D,
     $           X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $           BP,
     $           NOFTST(ILEVEL),IGRDOFTST(1,ILEVEL))
            CALL CONLFTMGBP 
     $           (IDIM3D,IDIM2D,IDIM1D,
     $           X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $           BP,
     $           NOFTST(ILEVEL),IGRDOFTST(1,ILEVEL))
         END DO



C Schleife ueber alle feinen Gitter, 
C ob Fuellalgorithmus nochmal aktiv war
         NSINGLE = 0.0
         NALL = 0.0
         DO IGRID=NGRDOLD+1,NGRDDFD
C         DO IGRID=STARTGRID+1,NGRDDFD
         DO ISUB = 1,MAX(1,NOFSLCHILDS(IGRID))
            IF (LSLICE(IGRID)) THEN
               IGRIDNUM = IGRDOFSLCHILD(ISUB,IGRID)
            ELSE
               IGRIDNUM = IGRID
            ENDIF
            MYID   = 0
            IPROC  = 0
            IPPAR  = 0
            IF(MYID .EQ. IPROC) THEN
C     hier werden die einzelnen Gitter bearbeitet, ANFANG               
            CALL MGPOINT (IP3,IP2,IP1,IBB,IB3,IBU,IGRIDNUM)
            CALL MGDIMS  (KK,JJ,II,IGRIDNUM)
            CALL MGBASB(NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRIDNUM)

            rmeanbp = 0.0
            do I = 1,II*JJ*KK
               rmeanbp = rmeanbp + abs(BP(IP3-1+I))
            enddo
            rmeanbp = rmeanbp/(II*JJ*KK)
            write(*,*) 'BLOCKBP rmeanbp B',rmeanbp,IGRIDNUM

            DO K = 3,KK-2
               DO J = 3,JJ-2
                  I = 2
                  IF ((BP(IP3-1+K+KK*(J-1)+KK*JJ*(I-1)) .EQ. 1)
     $                 .AND. (BP(IP3-1+K+KK*(J-1)+KK*JJ*(I)) .EQ. 0))
     $                 THEN
                     NSINGLE = NSINGLE + 1.0
C                     write(*,*) 'BBB',IGRIDNUM, I,J,K
                  END IF

                  I = II-2
                  IF ((BP(IP3-1+K+KK*(J-1)+KK*JJ*(I-1)) .EQ. 0)
     $                 .AND. (BP(IP3-1+K+KK*(J-1)+KK*JJ*(I)) .EQ. 1))
     $                 THEN
                     NSINGLE = NSINGLE + 1.0
C                     write(*,*) 'BBB',IGRIDNUM, I,J,K
                  END IF
               ENDDO
            ENDDO
            DO K = 3,KK-2
               DO I = 3,II-2
                  J = 2
                  IF ((BP(IP3-1+K+KK*(J-1)+KK*JJ*(I-1)) .EQ. 1)
     $                 .AND. (BP(IP3-1+K+KK*(J)+KK*JJ*(I-1)) .EQ. 0))
     $                 THEN
                     NSINGLE = NSINGLE + 1.0
C                     write(*,*) 'BBB',IGRIDNUM, I,J,K
                  END IF
                  J = JJ-2
                  IF ((BP(IP3-1+K+KK*(J-1)+KK*JJ*(I-1)) .EQ. 0)
     $                 .AND. (BP(IP3-1+K+KK*(J)+KK*JJ*(I-1)) .EQ. 1))
     $                 THEN
                     NSINGLE = NSINGLE + 1.0
C                     write(*,*) 'BBB',IGRIDNUM, I,J,K
                  END IF
               ENDDO
            ENDDO

C  auskommentiert, da CONNECT in z-Richtung nicht implementiert
C            DO J = 3,JJ-2
C               DO I = 3,II-2
C                  K = 2
C                  IF ((BP(IP3-1+K+KK*(J-1)+KK*JJ*(I-1)) .EQ. 1)
C     $                 .AND.(BP(IP3-1+K+1+KK*(J-1)+KK*JJ*(I-1)) .EQ. 0))
C     $                 then
C                     NSINGLE = NSINGLE + 1.0
C                     end if
C                  K = KK-2
C                  IF ((BP(IP3-1+K+KK*(J-1)+KK*JJ*(I-1)) .EQ. 0)
C     $                 .AND.(BP(IP3-1+K+1+KK*(J-1)+KK*JJ*(I-1)) .EQ. 1))
C     $               then
C                     NSINGLE = NSINGLE + 1.0
C                     end if
C               ENDDO
C            ENDDO
               
C     einzelne Gitter, ENDE
            ENDIF
         ENDDO
         ENDDO              
          NALL = NSINGLE
          IF (NALL .LT. 0.1) EXIT


C Schleife über alle neuen feinen Gitter fuer Fuellalgorithmus
         DO IGRID=NGRDOLD+1,NGRDDFD
C         DO IGRID=STARTGRID+1,NGRDDFD
         DO ISUB = 1,MAX(1,NOFSLCHILDS(IGRID))
            IF (LSLICE(IGRID)) THEN
               IGRIDNUM = IGRDOFSLCHILD(ISUB,IGRID)
            ELSE
               IGRIDNUM = IGRID
            ENDIF
            MYID   = 0
            IPROC  = 0
            IPPAR  = 0
            IF(MYID .EQ. IPROC) THEN
C     hier werden die einzelnen Gitter bearbeitet, ANFANG               
            CALL MGPOINT (IP3,IP2,IP1,IBB,IB3,IBU,IGRIDNUM)
            CALL MGDIMS  (KK,JJ,II,IGRIDNUM)
            CALL MGBASB(NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRIDNUM)
C           Feld umspeichern : KJI -> IJK
            DO K=1,KK
            DO J=1,JJ
            DO I=1,II
               HILF(IP3-1 + I + II*(J-1) + II*JJ*(K-1))
     $              = BP(IP3-1 + K + KK*(J-1) + KK*JJ*(I-1))
            ENDDO
            ENDDO
            ENDDO
C     call bp fuer einen Prozessor (BP(IP3), Koordinaten(IP1), KK,JJ,II)
            CALL BLOCKGRID(II,JJ,KK,0,X(IP1),Y(IP1),Z(IP1),
     $           DX(IP1),DY(IP1),DZ(IP1),
     $           DDX(IP1),DDY(IP1),DDZ(IP1),
     $           IGRIDNUM,NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,
     $           HILF(IP3),0)

C           Feld umspeichern : IJK -> KJI
            DO K=1,KK
            DO J=1,JJ
            DO I=1,II
               BP(IP3-1 + K + KK*(J-1) + KK*JJ*(I-1))
     $              = HILF(IP3-1 + I + II*(J-1) + II*JJ*(K-1))
            ENDDO
            ENDDO
            ENDDO

               
C     einzelne Gitter, ENDE
            ENDIF
         ENDDO
         ENDDO              



      END DO
      WRITE(*,*) 'BLOCKBP: ADDITIONAL BLOCKING FINISHED'
      WRITE(*,*) 'BLOCKBP: Iterationen ueber Gitter', ILOOP
      IF (ILOOP.EQ.10+1) CALL ERRR (583,'BLOCKBP')


C Schleife über alle neuen feinen Gitter, um Werte -999 auf Null zu setzten
C (s. Subroutine blockingadd)

      DO IGRID=NGRDOLD+1,NGRDDFD
C       DO IGRID=STARTGRID+1,NGRDDFD
         DO ISUB = 1,MAX(1,NOFSLCHILDS(IGRID))
            IF (LSLICE(IGRID)) THEN
               IGRIDNUM = IGRDOFSLCHILD(ISUB,IGRID)
            ELSE
               IGRIDNUM = IGRID
            ENDIF
            MYID   = 0
            IPROC  = 0
            IPPAR  = 0
            IF(MYID .EQ. IPROC) THEN
C     hier werden die einzelnen Gitter bearbeitet, ANFANG               
            CALL MGPOINT (IP3,IP2,IP1,IBB,IB3,IBU,IGRIDNUM)
            CALL MGDIMS  (KK,JJ,II,IGRIDNUM)

            DO K=1,KK
            DO J=1,JJ
            DO I=1,II
               IF (BP(IP3-1 + K + KK*(J-1) + KK*JJ*(I-1)).LE.-990.0)
     $              BP(IP3-1 + K + KK*(J-1) + KK*JJ*(I-1)) = 0.0
            ENDDO
            ENDDO
            ENDDO
               
C     einzelne Gitter, ENDE
            ENDIF
         ENDDO
      ENDDO              

      END IF
      
      END
