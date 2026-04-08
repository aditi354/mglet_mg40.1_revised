










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
       SUBROUTINE DMIXLE (KK,JJ,II,KMX,JMX,IMX,
     $                  B,HILF,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                  IC1,IC2,JC1,JC2,KC1,KC2,
     $                  NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,
     $                  GRADPX,RHO,GMOL,CONV1S,
     $                  CONV1SANF,CONV1SEND,TRANSLES1,TRANSLES2,
     $                  U)
C*STARLET***************************************************************
C        D M I X L E      BERECHNUNG DES MISCUNGSWEGES
C                         UND ABLEGEN AUF B-FELD
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        X, Y, Z        - KOORDINATEN DER ZELLDEFINITIONSPUNKTE
C        DX,DY,DZ       - ABSTAND DER GITTERPUNKTE
C        DDX,DDY,DDZ    - KANTENLAENGEN DER KONTROLLVOLUMINA
C        B(KK,JJ,II)    + WANDABSTANDSFELD ENTHAELT IM BERECHNUNGSGEBIET
C                         DEN MISCHUNGSWEG
C        IC1,JC1,KC1    - LINKE  RAENDER DER BOUNDING BOX
C        IC2,JC2,KC2    - RECHTE RAENDER DER BOUNDING BOX
C
C UPROG                 : ERRR
C
C VERS:  13.03.92 (MM)  : ORIGINAL
C         1.04.92 (MM)  : VEKTORISIERBARE VERSION DER SCHLEIFE 600
C         4. 4.92 (MM)  : SCHLEIFEN ZUR WANDKORREKTUR AM KOERPER
C                         LAUFEN NUR INNERHALB DER
C                         BOUNDING BOX (IC1,IC2,JC1,JC2,KC1,KC2)
C        26. 3.93 (MM)  : RANDBEDINGUNGEN DEHEN UEBER KOPF EIN
C        24.10.93 (MM)  : VAN-DRIEST''SCHE WANDDAEMPFUNG EINGEFUEHRT
C                         GILT NUR FUER KANALSTROEMUNG MIT ZTOT=2.0
C
C DEFINE-DIREKTIVEN     : DDMIX2, DDMIX3, DDMIX4,  VAN_DRIEST
C
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

      COMMON /CONLES/  CONV2S,CMUE,CAPPA,ECONST
      SAVE   /CONLES/
C
C
      REAL        X(II),         Y(JJ),         Z(KK),
     $           DX(II),        DY(JJ),        DZ(KK),
     $          DDX(II),       DDY(JJ),       DDZ(KK),
     $      B(KK,JJ,II),   HILF(KK,JJ),       U(KK,JJ,II)
C
C       SUBGRIDSCALE CONSTANT
      REAL   CONV1S(II)
C
C
C

      KM1    = KMX-1
      JM1    = JMX-1
      IM1    = IMX-1
      KM2    = KMX-2
      JM2    = JMX-2
      IM2    = IMX-2
C
      CAP05  = 0.5 * CAPPA
      SMAL2  = 2.0 * SMALL
      THIRD  = 1.0 / 3.0
C
C       GEOMETRY CONSTANT TO CALCULATE THE 
C         TRANSITION INTERPOLATION OF LES

      TRANSCON = 1./(SMALL+TRANSLES2-TRANSLES1)

      DO 100 I  = 1,IMX
         IF (X(I) .LE. TRANSLES1) THEN
           CONV1S(I) = CONV1SANF
        
         ELSEIF (X(I).LT.TRANSLES2 .AND. X(I).GT.TRANSLES1) THEN
           CONV1S(I) = ((TRANSLES2-X(I))*CONV1SANF+
     &               (X(I)-TRANSLES1)*CONV1SEND)*TRANSCON
         
         ELSEIF (X(I) .GE. TRANSLES2) THEN
           CONV1S(I) = CONV1SEND

         ENDIF
 
         DDXI   = DDX(I)
         RDDXPL = 1.0/DDXI
C
         DO 200 J  = 1,JMX
            DDYJ   = DDY(J)
            RDDYPL = 1.0/DDYJ

C
C
C
            DO 300 K  = 1,KMX
               DDZK   = DDZ(K)
               RDDZPL = 1.0/DDZK
C
C
C
  300   CONTINUE
  200   CONTINUE
  100   CONTINUE
C
C                                 AN FESTEN WAENDEN (NUR IN DER WANDNAEC
C                                 STEN ZELLE !!!) WIRD DER MISCHUNGSWEG
C                                 DEM PRANDTL''SCHEN MISCHUNGSWEG ANGEPAS
C
      IF(NCUB.EQ.5) THEN
C
C
C
         DO  I=IC1,IC2

            DO  J=JC1,JC2
            DO  K=KC1,KC2

              CAPDZP =    GREAT * (SIGN(0.5,B(K,J,I)*B(K+1,J,I))+0.5)
     $           - CAP05 * DDZ(K) * (SIGN(0.5,B(K,J,I)*B(K+1,J,I))-0.5)

              CAPDZN =    GREAT * (SIGN(0.5,B(K,J,I)*B(K-1,J,I))+0.5)
     $           - CAP05 * DDZ(K) * (SIGN(0.5,B(K,J,I)*B(K-1,J,I))-0.5)

              CAPDYP =    GREAT * (SIGN(0.5,B(K,J,I)*B(K,J+1,I))+0.5)
     $           - CAP05 * DDY(J) * (SIGN(0.5,B(K,J,I)*B(K,J+1,I))-0.5)

              CAPDYN =    GREAT * (SIGN(0.5,B(K,J,I)*B(K,J-1,I))+0.5)
     $           - CAP05 * DDY(J) * (SIGN(0.5,B(K,J,I)*B(K,J-1,I))-0.5)

              CAPDXP =    GREAT * (SIGN(0.5,B(K,J,I)*B(K,J,I+1))+0.5)
     $           - CAP05 * DDX(I) * (SIGN(0.5,B(K,J,I)*B(K,J,I+1))-0.5)

              CAPDXN =    GREAT * (SIGN(0.5,B(K,J,I)*B(K,J,I-1))+0.5)
     $           - CAP05 * DDX(I) * (SIGN(0.5,B(K,J,I)*B(K,J,I-1))-0.5)

              HILF(K,J) = 
     $        MIN(CAPDZP,CAPDZN,CAPDYP,CAPDYN,CAPDXP,CAPDXN)/CONV1S(I)

            ENDDO
            ENDDO

            DO  J=JC1,JC2
            DO  K=KC1,KC2

               B(K,J,I) = MIN(B(K,J,I),HILF(K,J))

            ENDDO
            ENDDO

         ENDDO
C
      ENDIF


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  TOP NOSLIP

      IF(NTOP.EQ.5) THEN
         DO  I=2,IM1
         DO  J=2,JM1
         B(KMX-2,J,I) = MIN(B(KMX-2,J,I),CAP05*DDZ(KMX-2)/CONV1S(I))
         ENDDO
         ENDDO

      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  BOTTOM NOSLIP

      IF(NBOT.EQ.5) THEN
         DO  I=2,IM1
         DO  J=2,JM1
              B(3,J,I) = MIN(B(3,J,I),CAP05*DDZ(3)/CONV1S(I))
         ENDDO
         ENDDO

      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  RIGHT NOSLIP

      IF(NRGT.EQ.5) THEN

         DO  I=2,IM1
         DO  K=2,KM1
              B(K,3,I) = MIN(B(K,3,I),CAP05*DDY(3)/CONV1S(I))
         ENDDO
         ENDDO

      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  LEFT NOSLIP

      IF(NLFT.EQ.5) THEN

         DO  I=2,IM1
         DO  K=2,KM1
              B(K,JM2,I) = MIN(B(K,JM2,I),CAP05*DDY(JM2)/CONV1S(I))
         ENDDO
         ENDDO

      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  FRONT NOSLIP

      IF(NFRO.EQ.5) THEN

         DO  J=2,JM1
         DO  K=2,KM1
              B(K,J,3) = MIN(B(K,J,3),CAP05*DDX(3)/CONV1S(3))
         ENDDO
         ENDDO

      ENDIF


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  BACK NOSLIP

      IF(NBAC.EQ.5) THEN

         DO  J=2,JM1
         DO  K=2,KM1
              B(K,J,IM2) = MIN(B(K,J,IM2),CAP05*DDX(IM2)/CONV1S(3))
         ENDDO
         ENDDO

      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                                   BELEGUNG DER FESTEN WAENDE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  BOTTOM SLIP

      IF(NBOT.EQ.6) THEN

         DO  I=1,IMX
         DO  J=1,JMX
                    B(1,J,I) = -2.0
                    B(2,J,I) = -2.0
         ENDDO
         ENDDO

      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  TOP SLIP

      IF(NTOP.EQ.6) THEN

         DO  I=1,IMX
         DO  J=1,JMX
                    B(KMX,J,I) = -2.0
                    B(KM1,J,I) = -2.0
         ENDDO
         ENDDO

      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  LEFT SLIP

      IF(NLFT.EQ.6) THEN

         DO  I=1,IMX
         DO  K=1,KMX
                    B(K,JMX,I) = -2.0
                    B(K,JM1,I) = -2.0
         ENDDO
         ENDDO

      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  RIGHT SLIP

      IF(NRGT.EQ.6) THEN

         DO  I=1,IMX
         DO  K=1,KMX
                    B(K,  1,I) = -2.0
                    B(K,  2,I) = -2.0
         ENDDO
         ENDDO

      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  BOTTOM NOSLIP

      IF(NBOT.EQ.5) THEN

         DO  I=1,IMX
         DO  J=1,JMX
                    B(1,J,I) = -1.0
                    B(2,J,I) = -1.0
         ENDDO
         ENDDO

      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  TOP NOSLIP

      IF(NTOP.EQ.5) THEN

         DO  I=1,IMX
         DO  J=1,JMX
                    B(KMX,J,I) = -1.0
                    B(KM1,J,I) = -1.0
         ENDDO
         ENDDO

      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  LEFT NOSLIP

      IF(NLFT.EQ.5) THEN

         DO  I=1,IMX
         DO  K=1,KMX
                    B(K,JMX,I) = -1.0
                    B(K,JM1,I) = -1.0
         ENDDO
         ENDDO

      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  RIGHT NOSLIP

      IF(NRGT.EQ.5) THEN

         DO  I=1,IMX
         DO  K=1,KMX
                    B(K,  1,I) = -1.0
                    B(K,  2,I) = -1.0
         ENDDO
         ENDDO

      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  FRONT NOSLIP

      IF(NFRO.EQ.5) THEN

         DO  J=1,JMX
         DO  K=1,KMX
              B(K,J,1) = -1.0
              B(K,J,2) = -1.0
         ENDDO
         ENDDO

      ENDIF


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  BACK NOSLIP

      IF(NBAC.EQ.5) THEN

         DO  J=1,JMX
         DO  K=1,KMX
              B(K,J,IMX) = -1.0
              B(K,J,IM1) = -1.0
         ENDDO
         ENDDO

      ENDIF


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      RETURN
      END
