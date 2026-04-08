










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
      SUBROUTINE AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,KMXA,JMXA,
     $                    IMXA,PHI,SPHI,XHOMOG,YHOMOG,ZHOMOG)
C*STARLET***************************************************************
C        A R E A M        IN AREAM WERDEN LINIEN- ODER FLAECHEN-
C                         MITTELWERTE EINER VARIABLEN  PHI  GEBILDET
C                         FALLS  XHOMOG  UND/ODER  YHOMOG  DEFINIERT
C                         WURDE.
C                         DIE MITTELWERTE WERDEN JEWEILS IN DEN
C                         ERSTEN INDEX DER VARIABLE PHI ABGELEGT.
C                         DABEI MUSS ANGENOMMEN WERDEN, DASS:
C                             - PHI EIN TEMPORAERES HILFSFELD IST
C                             - MINDESTENS EINE RANDSCHICHT VORLIEGT.
C                         GEMITTELT WIRD IMMER VON 3 BIS IMX(JMX)-2.
C                         DIE ENTSTEHENDEN MITTELWERTE WERDEN IM
C                         SPHI-FELD  AUFSUMMIERT, SO DASS AM ENDE
C                         EINES LAUFES DIE VORHANDENEN <PHI> - WERTE
C                         VERBESSERT WERDEN KOENNEN.
C                   
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KKA,JJA,IIA    - ARRAYDIMENSIONEN DER AUSWERTEFELDER
C                         BEISPIEL: FALLS J-RICHTUNG NICHT HOMOGEN IST,
C                         GILT JJA = JJ. ANDERNFALLS IST JJA = 1.
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        KMXA,JMXA,IMXA - GRENZEN DER AUSWERTEFELDER (MIT BOUND)
C        PHI(KK,JJ,II)  - ALLGEMEINE VARIABLE, Z.B. 'U', 'V', 'W' ...
C        SPHI(KKA,JJA,IIA) + FELD ZUM SUMMIEREN DER MITTELWERTE
C
C UPROG                 : KEINE
C
C        23.07.86 (HW)  : ORIGINAL
C        15. 5.93 (MM)  : VOELLIG NEUE VERSION
C
C*STARLET***************************************************************
C
      INTEGER  KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,KMXA,JMXA,IMXA
      REAL     PHI(KK,JJ,II),       SPHI(KKA,JJA,IIA)
      LOGICAL  XHOMOG,YHOMOG,ZHOMOG
C
C                             KONTROLLE
C
      IF (XHOMOG) THEN
            IF( IIA.NE.1) CALL ERRR (501,' AREAM')
            IF(IMXA.NE.1) CALL ERRR (502,' AREAM')
      ELSE
            IF( IIA.NE. II) CALL ERRR (503,' AREAM')
            IF(IMXA.NE.IMX) CALL ERRR (504,' AREAM')
      ENDIF
      IF (YHOMOG) THEN
            IF( JJA.NE.1) CALL ERRR (505,' AREAM')
            IF(JMXA.NE.1) CALL ERRR (506,' AREAM')
      ELSE
            IF( JJA.NE. JJ) CALL ERRR (507,' AREAM')
            IF(JMXA.NE.JMX) CALL ERRR (508,' AREAM')
      ENDIF
      IF (ZHOMOG) CALL ERRR (509,' AREAM')
C
      XFAK = 1.
      YFAK = 1.
C
      ISTART = 2
      ISTOP  = IMX-1
      JSTART = 2
      JSTOP  = JMX-1
C
C                                 X-RICHTUNG IST HOMOGEN, DAHER BIL-
C                                 DUNG EINES LINIENMITTELWERTES.
C                                 (FALLS Y-RI. AUCH HOMOGEN, DANN
C                                 FLAECHENMITTELWERT)
C
      IF (XHOMOG) THEN
C
C                                   VORBELEGUNG
C
            DO J=2,JMX-1
            DO K=2,KMX-1

            PHI(K,J,1) = PHI(K,J,3)

            ENDDO
            ENDDO
C
C                                   SUMMATION
C
         DO I=4,IMX-2
            DO J=2,JMX-1
            DO K=2,KMX-1

            PHI(K,J,1) = PHI(K,J,1) + PHI(K,J,I)

            ENDDO
            ENDDO
         ENDDO
C
C                                 ANZAHL DER STICHPROBEN
C
         XFAK = 1./FLOAT(IMX-4)
         ISTART = 1
         ISTOP  = 1
      ENDIF

C
C
C                                 Y-RICHTUNG IST HOMOGEN, DAHER BIL-
C                                 DUNG EINES LINIENMITTELWERTES.
C                                 (FALLS X-RI. AUCH HOMOGEN, DANN
C                                 FLAECHENMITTELWERT)
C
      IF (YHOMOG) THEN
C
C                                   VORBELEGUNG
C
            DO I=ISTART,ISTOP
            DO K=2,KMX-1

            PHI(K,1,I) = PHI(K,3,I)

            ENDDO
            ENDDO
C
C                                   SUMMATION
C
         DO J=4,JMX-2
            DO I=ISTART,ISTOP
            DO K=1,KMX-1

            PHI(K,1,I) = PHI(K,1,I) + PHI(K,J,I)

            ENDDO
            ENDDO
         ENDDO
C
C                                 ANZAHL DER STICHPROBEN
C
         YFAK = 1./FLOAT(JMX-4)
         JSTART = 1
         JSTOP  = 1

      ENDIF
C
C                                BELEGUNG DES SPHI-FELDES
C
C                                KONTROLLE
C
      IF(ISTOP.GT.IIA) CALL ERRR (510,' AREAM')
      IF(JSTOP.GT.JJA) CALL ERRR (511,' AREAM')
C
      DO I = ISTART,ISTOP
      DO J = JSTART,JSTOP
      DO K = 2,KMX-1

         SPHI(K,J,I) = SPHI(K,J,I) + PHI(K,J,I)*XFAK*YFAK

      ENDDO
      ENDDO
      ENDDO
C
C
C
      RETURN
      END
