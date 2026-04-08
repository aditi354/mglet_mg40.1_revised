










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
      SUBROUTINE DEOBI   (KK,JJ,II,KMX,JMX,IMX,XTOT,YTOT,ZTOT,
     $                    CIDEND,IIDEND,RIDEND,CIDENT,IIDENT,RIDENT,
     $                    U,V,W,P,G,P1,UGRID,ITSTEP,ITTOT,ITINT,MTSTEP,
     $                    TIMEPH,DT,ISW)
C*STARLET***************************************************************
C        D E O B I        AUSGABE DER GESCHWINDIGKEITSKOMPONENTEN U,V,W
C                         DER EBENEN I=ISW UND I=ISW+1 (BINAER,KANAL 12)
C*STARLET***************************************************************
C
C PARAM: KK,JJ,II       - ARRAYGRENZEN
C
C VERS:  14.02.86 (HW)  : ORIGINAL
C        08.01.87 (HW)  : COPFI IN COPFIB UMBENANNT (WEGEN BINAERER
C                         KOPIE VON TAPE11)
C        30.12.88 (HW)  : DTEPS  WIRD MASCHINENABHAENGIG BELEGT
C
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
      CHARACTER (LEN=8)     CIDENT(10),   CIDEND(10)
      INTEGER    IIDENT(100),  IIDEND(100)
      REAL       RIDENT(100),  RIDEND(100)
C
      REAL       U(KK,JJ,II),    V(KK,JJ,II),    W(KK,JJ,II),
     $           P(KK,JJ,II),    G(KK,JJ,II)
C
      REAL       P1(KK,JJ, 7)
C
      ISWP    = ISW + 1
      IF(ISW   .LT.       1) CALL ERRR (501,' DEOBI    ')
      IF(ISWP  .GT. (IMX-2)) CALL ERRR (502,' DEOBI    ')
C
      IF(ITTOT .GT.       0) GOTO 2100
C
C                                  VORBELEGUNG VON "TCYCLE" MIT
C                                  MTSTEP*DT.
C
C
         RIDENT(31) = FLOAT(IIDENT(1)) * RIDENT(2)
C
C                                  EINLESEN WICHTIGER KENNDATEN
C
         REWIND 12
         WRITE(12) CIDENT
         WRITE(12) IIDENT
         WRITE(12) RIDENT
C
         WRITE(12) KK,JJ,II,KMX,JMX,IMX
C
 2100 IF((ITSTEP .NE. 0) .OR. (ITTOT .EQ. 0)) GOTO 2110
         DTEPS      = 10.0 * SMAONE
         CALL COPFIB  (KK,JJ,II,KMX,JMX,IMX,CIDEND,IIDEND,RIDEND,P1,
     $                 ITSTEP,ITTOT,ITINT,MTSTEP,DT)
C
C                                  DER ZEITSCHRITT DES VOHERGEHENDEN
C                                  LAUFES MUSS MIT DEM AKTUELLEN
C                                  ZEITSCHRITT UEBEREINSTIMMEN !
C
         IF(ABS(RIDEND(2) - DT) .GT. DTEPS)
     $      CALL ERRR (505,' DEOBI    ')
         RETURN
C
 2110 IF(ABS(UGRID) .LT. SMALL) GOTO 2200
C
C                                  HIER: GALILEI-RUECKTRANSFORMATION
C
         JM2    = JMX - 2
         KM2    = KMX - 2
C
         DO 100 I=ISW,ISWP
            DO 105 J=3,JM2
               DO 110 K=3,KM2
  110             U(K,J,I) = U(K,J,I) + UGRID
  105       CONTINUE
  100    CONTINUE
C
 2200 WRITE(12) (((U(K,J,I),K=1,KMX),J=1,JMX),I=ISW,ISWP),
     $          (((V(K,J,I),K=1,KMX),J=1,JMX),I=ISW,ISWP),
     $          (((W(K,J,I),K=1,KMX),J=1,JMX),I=ISW,ISWP),
     $           ((G(K,J,ISWP),K=1,KMX),J=1,JMX)
C
      IF(ABS(UGRID) .LT. SMALL) GOTO 2300
C
C                                  GALILEI-TRANSFORMATION
C
         DO 200 I=ISW,ISWP
            DO 205 J=3,JM2
               DO 210 K=3,KM2
  210             U(K,J,I) = U(K,J,I) - UGRID
  205       CONTINUE
  200    CONTINUE
C
 2300 RETURN
      END
