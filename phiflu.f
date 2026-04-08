










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
      SUBROUTINE PHIFLU  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,UGRID,PHI,
     $                    PHIF,APHI,IVAR,XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C*STARLET***************************************************************
C        P H I F L U      PHIFLU BERECHNET DIE FLUKTUATIONEN EINER
C                         GROESSE  PHI  GEGENUEBER  <PHI>.
C                         PHIFLU = PHI - <PHI>
C
C                         A C H T U N G: FALLS FLUKTUATIONEN DES DRUCKES
C                         BESTIMMT WERDEN SOLLEN, MUSS DAS DRUCKNIVEAU
C                         VOR AUFRUF DIESES UNTERPROGRAMMES NEU GESETZT
C                         WERDEN (SUBR. PLEVEL) !!
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KKA,JJA,IIA    - ARRAYDIMENSIONEN DER AUSWERTEFELDER
C                         BEISPIEL: FALLS J-RICHTUNG NICHT HOMOGEN IST,
C                         GILT JJA = JJ. ANDERNFALLS IST JJA = 1.
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        KMXA,JMXA,IMXA - GRENZEN DER AUSWERTEFELDER (MIT BOUND)
C        UGRID          - GESCHWINDIGKEIT DES KOORDINATENSYSTEMS
C                         IN X-RICHTUNG (GALILEI-TRANSFORMATION)
C        PHI(KK,JJ,II)  - ALLGEMEINE VARIABLE: 'U', 'V', 'W' ODER 'P'
C                         ZULAESSIG
C        PHIF(KK,JJ,II) + ENTHAELT DIE FLUKTUATIONEN
C        APHI(KKA,JJA,IIA) - ENTHAELT DEN STATIST. MITTELWERT <PHI>
C        IVAR           - CHARACTER-VARIABLE DER FORM 'U', 'V', 'W', 'P'
C
C UPROG                 : ERRR
C
C DEFINE-DIREKTIVEN     : FRFIX, RIFIX, XHOMOG, YHOMOG
C
C        23.07.86 (HW)  : ORIGINAL
C        01.07.88 (HW)  : FUER VARIABLEN DIE IN ZELLMITTE DEFINIERT SIND
C                         (IVAR = '  QP' -> QUASI PRESSURE) WIRD DIE I=2
C                         EBENE FUER FRFIX AUSGEWERTET. ANALOG F. RIFIX
C        21.08.88 (HW)  : MEHRFACH IM MASCHENGITTER VERSCHOBENE
C                         VARIABLEN KOENNEN BERUECKSICHTIGT WERDEN
C        12.02.03 (TB)  : SCALAR FIELD STATISTICS ADDED
C
C*STARLET***************************************************************
C
      CHARACTER (LEN=4)  IVAR
C
C
      REAL     PHI(KK,JJ,II),       PHIF(KK,JJ,II),
     $         APHI(KKA,JJA,IIA)
C
      LOGICAL  XHOMOG,YHOMOG,ZHOMOG
C
C                             KONTROLLE
C
      IF (XHOMOG) THEN
            IF( IIA.NE.1) CALL ERRR (501,' PHIFLU')
            IF(IMXA.NE.1) CALL ERRR (502,' PHIFLU')
      ELSE
            IF( IIA.NE. II) CALL ERRR (503,' PHIFLU')
            IF(IMXA.NE.IMX) CALL ERRR (504,' PHIFLU')
      ENDIF
      IF (YHOMOG) THEN
            IF( JJA.NE.1) CALL ERRR (505,' PHIFLU')
            IF(JMXA.NE.1) CALL ERRR (506,' PHIFLU')
      ELSE
            IF( JJA.NE. JJ) CALL ERRR (507,' PHIFLU')
            IF(JMXA.NE.JMX) CALL ERRR (508,' PHIFLU')
      ENDIF
      IF (ZHOMOG) CALL ERRR (509,' PHIFLU')
C
      IF(IVAR .EQ. '   U' .OR. IVAR .EQ. '   V' .OR. IVAR .EQ. '   W'
     $   .OR. IVAR .EQ. '   P' .OR. IVAR .EQ. '  QP' .OR.
     $   IVAR .EQ. '  WV' .OR. IVAR .EQ. '  UW' .OR. IVAR .EQ. '  VU'
     $   .OR. IVAR .EQ. ' UVW'
     $     ) GOTO 2000
         CALL ERRR (501,' PHIFLU   ')
C
 2000 KSTAG  = 0
      JSTAG  = 0
      ISTAG  = 0
      IFRFIX = 0
      JRIFIX = 0
C
      IF(IVAR .EQ. '   P') GOTO 2010
      IF (NFRO .EQ. 2) THEN
C                            FRONT-FIXED RANDBED.
         IFRFIX = 1
      ENDIF
      IF (NRGT .EQ. 2) THEN
C                            RIGHT-FIXED RANDBED.
         JRIFIX = 1
      ENDIF
C
 2010 UGALIL = 0.0
C
      IF(IVAR .EQ. '   U') THEN
         ISTAG  = 1
         UGALIL = UGRID
      ENDIF
      IF(IVAR .EQ. '   V') JSTAG  = 1
      IF(IVAR .EQ. '   W') KSTAG  = 1
      IF(IVAR .EQ. '  WV') THEN
         KSTAG  = 1
         JSTAG  = 1
      END IF
      IF(IVAR .EQ. '  UW') THEN
         ISTAG  = 1
         KSTAG  = 1
      END IF
      IF(IVAR .EQ. '  VU') THEN
         JSTAG  = 1
         ISTAG  = 1
      END IF
      IF(IVAR .EQ. ' UVW') THEN
         ISTAG  = 1
         JSTAG  = 1
         KSTAG  = 1
      END IF
C
C                                 START- UND STOPINDIZES
C
C OLD:      ISTART = 3-ISTAG-IFRFIX
      ISTART = 2
C OLD:      ISTOP  = IMX - 2
      ISTOP  = IMX - 1
      IF (XHOMOG) THEN
C
C                                 X-RICHTUNG IST HOMOGEN
C
      ISTO   = 1
      ENDIF
C
C OLD:      JSTART = 3-JSTAG-JRIFIX
C OLD:      JSTOP  = JMX - 2
      JSTART = 2
      JSTOP  = JMX - 1
      IF (YHOMOG) THEN
C
C                                 Y-RICHTUNG IST HOMOGEN
C
      JSTO   = 1
      ENDIF
C
      DO 10 I = ISTART,ISTOP
      IF (.NOT.XHOMOG) THEN
         ISTO    = I
      ENDIF
         DO 20 J = JSTART,JSTOP
      IF (.NOT.YHOMOG) THEN
            JSTO    = J
      ENDIF
C
C                                 **************************************
C                                 BEGINN DER K - SCHLEIFE
C                                 **************************************
C
C
            KSTART = 3 - KSTAG
            KSTOP  = KMX - 2
C
            DO 30 K = KSTART,KSTOP
               PHIF(K,J,I) = PHI(K,J,I) + UGALIL - APHI(K,JSTO,ISTO)
C
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
      RETURN
      END
