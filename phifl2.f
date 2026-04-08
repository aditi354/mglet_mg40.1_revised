










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
      SUBROUTINE PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,PHIF,PHIR,
     $                    IVAR)
C*STARLET***************************************************************
C        P H I F L 2      PHIFL2 BERECHNET DIE QUADRATE DER FLUK-
C                         TUATIONEN:
C                         PHIR = PHIF **2
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KKA,JJA,IIA    - ARRAYDIMENSIONEN DER AUSWERTEFELDER
C                         BEISPIEL: FALLS J-RICHTUNG NICHT HOMOGEN IST,
C                         GILT JJA = JJ. ANDERNFALLS IST JJA = 1.
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        KMXA,JMXA,IMXA - GRENZEN DER AUSWERTEFELDER (MIT BOUND)
C        PHIF(KK,JJ,II) - ENTHAELT DIE FLUKTUATIONEN DER GROESSE PHI
C        PHIR(KK,JJ,II) + ENTHAELT DIE QUADRATE DER FLUKTUATIONEN DER
C                         GROESSE PHI
C        IVAR           - CHARACTER-VARIABLE DER FORM 'U', 'V', 'W', 'P'
C
C UPROG                 : ERRR
C
C DEFINE-DIREKTIVEN     : FRFIX, RIFIX, XHOMOG, YHOMOG
C
C        12.09.86 (HW)  : ORIGINAL
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
      REAL     PHIF(KK,JJ,II),      PHIR(KK,JJ,II)
C
      IF(IVAR .EQ. '   U' .OR. IVAR .EQ. '   V' .OR. IVAR .EQ. '   W'
     $   .OR. IVAR .EQ. '   P' .OR. IVAR .EQ. '  QP' .OR.
     $   IVAR .EQ. '  WV' .OR. IVAR .EQ. '  UW' .OR. IVAR .EQ. '  VU'
     $   .OR. IVAR .EQ. ' UVW'
     $     ) GOTO 2000     
         CALL ERRR (501,' PHIFL2   ')
C
 2000 KSTAG  = 0
      JSTAG  = 0
      ISTAG  = 0
      IFRFIX = 0
      JRIFIX = 0
C
      IF(IVAR .EQ. '   P') GOTO 2010
C
C
 2010 IF(IVAR .EQ. '   U') ISTAG  = 1
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
      ISTART = 3-ISTAG-IFRFIX
      ISTOP  = IMX - 2
C
C OLD:      JSTART = 3-JSTAG-JRIFIX
C OLD:      JSTOP  = JMX - 2
      JSTART = 2
      JSTOP  = JMX - 1
C
      DO 10 I = ISTART,ISTOP
         DO 20 J = JSTART,JSTOP
C
C                                 BEGINN DER K - SCHLEIFE
C                                 -----------------------
C
            IKST   = MIN0 (MAX0(I      , 3),IMX-2)
            IKSTPS = MIN0 (MAX0(I+ISTAG, 3),IMX-2)
            JKST   = MIN0 (MAX0(J      , 3),JMX-2)
            JKSTPS = MIN0 (MAX0(J+JSTAG, 3),JMX-2)
C
            KSTART = 3 - KSTAG
            KSTOP  = KMX - 2
C
            DO 30 K = KSTART,KSTOP
   30          PHIR(K,J,I) = PHIF(K,J,I)**2
   20    CONTINUE
   10 CONTINUE
C
      RETURN
      END
