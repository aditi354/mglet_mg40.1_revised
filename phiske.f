










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
      SUBROUTINE PHISKE  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,PHIF,APHIR,
     $                    PHISK,IVAR,XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C*STARLET***************************************************************
C        P H I S K E      PHISKE BERECHNET DIE SKEWNESS (SCHIEFE)
C                         DER FLUKTUATIONEN DER GROESSE PHI.
C                         PHISKE = (PHIF / <PHIRMS>) **3
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KKA,JJA,IIA    - ARRAYDIMENSIONEN DER AUSWERTEFELDER
C                         BEISPIEL: FALLS J-RICHTUNG NICHT HOMOGEN IST,
C                         GILT JJA = JJ. ANDERNFALLS IST JJA = 1.
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        KMXA,JMXA,IMXA - GRENZEN DER AUSWERTEFELDER (MIT BOUND)
C        PHIF(KK,JJ,II) - ENTHAELT DIE FLUKTUATIONEN DER GROESSE PHI
C        APHIR(KKA,JJA,IIA) - ENTHAELT DIE STATISTISCHEN MITTELWERTE
C                             <PHIRMS> DER ROOT-MEAN-SQUARE-WERTE
C        PHISK(KK,JJ,II)+ ENTHAELT DIE SKEWNESS-WERTE DER GROESSE PHIF
C        IVAR           - CHARACTER-VARIABLE DER FORM 'U', 'V', 'W', 'P'
C                                                     'T'
C
C UPROG                 : ERRR
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        24.07.86 (HW)  : ORIGINAL
C        12.02.03 (TB)  : SCALAR FIELD STATISTICS ADDED
C
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
      CHARACTER (LEN=4)  IVAR
C
C
      REAL     PHIF(KK,JJ,II),      APHIR(KKA,JJA,IIA),
     $         PHISK(KK,JJ,II)
C
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
      IF(IVAR .EQ. '   U' .OR. IVAR .EQ. '   V' .OR. IVAR .EQ. '   W'
     $   .OR. IVAR .EQ. '   P'
     $     ) GOTO 2000     
         CALL ERRR (501,' PHISKE   ')
C
 2000 KSTAG  = 0
      JSTAG  = 0
      ISTAG  = 0
      IFRFIX = 0
      JRIFIX = 0
      GROSS  = SQRT(GREAT)
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
 2010 IF(IVAR .EQ. '   U') ISTAG  = 1
      IF(IVAR .EQ. '   V') JSTAG  = 1
      IF(IVAR .EQ. '   W') KSTAG  = 1
C
C                                 START- UND STOPINDIZES
C
      ISTART = 3-ISTAG-IFRFIX
      ISTOP  = IMX - 2
      IF (XHOMOG) THEN
C
C                                 X-RICHTUNG IST HOMOGEN
C
      ISTO   = 1
      ENDIF
C
      JSTART = 3-JSTAG-JRIFIX
      JSTOP  = JMX - 2
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
C                                 BEGINN DER K - SCHLEIFE
C                                 -----------------------
C
C
            KSTART = 3 - KSTAG
            KSTOP  = KMX - 2
C
            DO 30 K = KSTART,KSTOP
               ARMS         = AMAX1 (SMALL,APHIR(K,JSTO,ISTO))
               PHISK(K,J,I) = (AMIN1(GROSS,ABS(PHIF(K,J,I) / ARMS)))**3
     $                      * SIGN(1.0,PHIF(K,J,I))
C
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
      RETURN
      END
