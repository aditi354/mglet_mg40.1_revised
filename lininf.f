










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
      SUBROUTINE LININF  (ILINT0, ILINTX, ILINTY, ILINTZ,

     $ KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,KMXA,JMXA,IMXA,JJ1L,JJ2L,ILIMXP,
     $ ISLIDI,LINFB,KKXL,KKYL,KKZL,KSXL,KSYL,KSZL,KH0L,IDIMF,KLASS,FVT,
     $ HVK,ISELEP,ISELEA,ISLINP,ISLINA,KSTW,RIDENT,HILF,UFG,AURG,VFG,
     $ AVRG,WFG,AWRG,OX,AOX,OXFG,AOXRG,OY,AOY,OYFG,AOYRG,OZ,AOZ,OZFG,
     $ AOZRG,ARXUU,SRXUU,ARYUU,SRYUU,ARZUU,SRZUU,ASXUU,SSXUU,ASYUU,SSYUU
     $ ,ASZUU,SSZUU,ARXVV,SRXVV,ARYVV,SRYVV,ARZVV,SRZVV,ASXVV,SSXVV,
     $ ASYVV,SSYVV,ASZVV,SSZVV,ARXWW,SRXWW,ARYWW,SRYWW,ARZWW,SRZWW,
     $ ASXWW,SSXWW,ASYWW,SSYWW,ASZWW,SSZWW, ARXUV,SRXUV,ARYUV,SRYUV,
     $ ARXUW,SRXUW,ARYUW,SRYUW,ARYVW,SRYVW, ACXUW,SCXUW,ACZUW,SCZUW,
     $ ARXOXX,SRXOXX,ARXOXY,SRXOXY,ARXOXZ,SRXOXZ,ARXOYY,SRXOYY,ARXOYZ,
     $ SRXOYZ,ARXOZZ,SRXOZZ,ARYOXX,SRYOXX,ARYOXY,SRYOXY,ARYOXZ,SRYOXZ,
     $ ARYOYY,SRYOYY,ARYOYZ,SRYOYZ,ARYOZZ,SRYOZZ,ARZOXX,SRZOXX,ARZOXY,
     $ SRZOXY,ARZOXZ,SRZOXZ,ARZOYY,SRZOYY,ARZOYZ,SRZOYZ,ARZOZZ,SRZOZZ,
     $ ASXOXX,SSXOXX,ASXOYY,SSXOYY,ASXOZZ,SSXOZZ,ASYOXX,SSYOXX,ASYOYY,
     $ SSYOYY,ASYOZZ,SSYOZZ,AHOZOY,SHOZOY,AHOZOX,SHOZOX,AHOYOX,SHOYOX)
C*STARLET***************************************************************
C        L I N I N F      IN LININF WERDEN ZUR INFORMATION EINIGE
C                         WICHTIGE GROESSEN, WIE Z.B. MAXIMA, MINIMA
C                         INTEGRALE LAENGENMASSE USW. FUER "LINIEN"-
C                         FELDER ERMITTELT.
C  A C H T U N G:         DIESE ROUTINE MUSS ERWEITERT WERDEN, FALLS
C                         DIE STATIST. AUSWERTUNG ERWEITERT WIRD !
C*STARLET***************************************************************
C
C PARAM: ILINT0         - ANZAHL DER AUFPUNKTE MIT RICHTUNGSINFORMATION
C                         '0' + 1 (AUFPKT. IN DEN I-LINES 2 ... ILINT0)
C        ILINTX         - ANZAHL DER AUFPUNKTE MIT RICHTUNGSINFORMATION
C                         'X' + 1 (AUFPKT. IN DEN I-LINES 2 ... ILINTX)
C        ILINTY         - ANZAHL DER AUFPUNKTE MIT RICHTUNGSINFORMATION
C                         'Y' + 1 (AUFPKT. IN DEN I-LINES 2 ... ILINTY)
C        ILINTZ         - ANZAHL DER AUFPUNKTE MIT RICHTUNGSINFORMATION
C                         'Z' + 1 (AUFPKT. IN DEN I-LINES 2 ... ILINTZ)
C        KK ,JJ ,II                 ARRAYDIMENSIONEN
C        KKA,JJA,IIA                ARRAYDIMENSIONEN DER AUSWERTEFELDER
C        KMX,JMX,IMX                GRENZEN DES BERECHNUNGSGEBIETES
C                                   MIT BOUND.
C        KMXA,JMXA,IMXA             GRENZEN DER AUSWERTEFELDER
C        JJ1L                       ARRAYDIMENSION IN J-RI. FUER
C                                   "LINIEN"-FELDER (= 1)
C        JJ2L                       ARRAYDIMENSION IN J-RI. FUER
C                                   "LINIEN"-FELDER (= 2)
C        ILIMXP                     ARRAYDIMENSION IN I-RI. FUER
C                                   "LINIEN"-FELDER
C        ISLIDI                     ARRAYDIMENSION FUER
C                                   DIE ISLIN.  -FELDER
C        LINFB                      ANZAHL DER INFORMATIONSBLOECKE
C                                   (-FELDER) JEDER "LINIE"
C        K12L                       ARRAYDIMENSION IN K-RICHTUNG FUER
C                                   "LINIEN"-FELDER:
C                                   1 = K : K_ORRELATIONSFUNKTIONEN BZW.
C                                           -KOEFFIZIENTEN
C                                   1 = S : S_PEKTREN
C                                   1 = H : H_AEUFIGKEITSVERTEILUNGEN
C                                   2 = 0 : KEINE AUSGEZEICHNETE KOORD.-
C                                           RICHTUNG
C                                   2 = X : X-RICHTUNG IST AUSGEZEICHNET
C                                   2 = Y : Y-RICHTUNG IST AUSGEZEICHNET
C                                   2 = Z : Z-RICHTUNG IST AUSGEZEICHNET
C        IDIMF                      ZUR DIMENSIONIERUNG DES FVT-FELDES
C        KLASS                      ANZAHL DER KLASSEN FUER DIE RELATIVE
C                                   HAEUFIGKEITSVERTEILUNG
C        FVT (0:IDIMF)              PUFFERFELD FUER DIE JEWEILS EXTRA-
C                                   HIERTE "ZEIT"REIHE
C        HVK (KLASS,                HILFSFELD FUER DIE HAEUFIGKEITS-
C             JJ2L  )               VERTEILUNG
C        ISELEP(2,752)              STEUERFELD F. DIE STAT. AUSWERTUNG
C                                   (ENTHAELT DIE STEUERDATEN DES MOMEN-
C                                   TANEN LAUFES)
C        ISELEA(2,752)              STEUERFELD F. DIE STAT. AUSWERTUNG
C                                   (ENTHAELT DIE STEUERDATEN DES VORAN-
C                                   GEGANGENEN LAUFES)
C        ISLINP(ISLIDI)             ENTHAELT DIE AUFPUNKTE IN KODIERTER
C                                   FORM FUER DEN MOMENTANEN LAUF
C        ISLINA(ISLIDI)             ENTHAELT DIE AUFPUNKTE IN KODIERTER
C                                   FORM FUER DEN VORANGEGANGENEN LAUF
C        KSTW(JJ,II)                STARTINDIZES DER W-KOMP.
C        RIDENT (100)               IDENT-FELD FUER REAL-KONSTANTEN
C        HILF   (KK ,JJ ,II )       ALLGEMEIN VERWENDBARES HILFSFELD
C        PHI                        PHI STEHT IM FOLGENDEN FUER EINE
C                                   BELIEBIGE VARIABLE (Z.B. U, V, W,
C                                   OMEGA_X --> OX, OY, OZ)
C        PHI    (KK ,JJ ,II )       MOMENTANWERT DER VARIABLEN PHI
C        PHIO   (KK ,JJ ,II )       MOMENTANWERT DER VARIABLEN PHI
C                                   ZUM ALTEN ZEITSCHRITT
C        APHI   (KKA,JJA,IIA)       ENSEMBLE-MITTELWERT DER VARIABLEN
C                                   PHI
C        SPHI   (KKA,JJA,IIA)       SUMMATIONSFELD FUER DIE MOMENTAN-
C                                   WERTE VON PHI ZUR VERBESSERUNG
C                                   DER ALTEN ENSEMBLE-MITTELWERTE
C        PHIFG  (KK ,JJ, II )       FLUKTUATIONEN (ANTEIL DER GROBSTRUK-
C                                   TUR) VON PHI.  PHIFG = PHI-<PHI>
C        APHIRG (KKA,JJA,IIA)       ENSEMBLE-MITTELWERT DER ROOT-MEAN-
C                                   SQUARE - WERTE DER FLUKTUATIONEN
C                                   VON PHI
C        SPHIRG (KKA,JJA,IIA)       SUMMATIONSFELD FUER DIE MOMENTANEN
C                                   RMS-WERTE.
C        A12345 (K?NL,JJAL,         ENSEMBLE-MITTELWERT EINES "LINIEN"-
C                ILIMXP)            FELDES
C        S12345 (K?NL,JJSL,         SUMMATIONSFELD FUER DIE MOMENTAN-
C                ILIMXP)            WERTE EINES "LINIEN"-FLEDES ZUR
C                                   VERBESSERUNG DER ALTEN ENSEMBLE-
C                                   MITTELWERTE
C                                   1 = R : KORRELATIONSFUNKTIONEN
C                                   1 = C : KORRELATIONSKOEFFIZIENTEN
C                                   1 = S : S_PEKTREN
C                                   1 = H : H_AEUFIGKEITSVERTEILUNGEN
C                                   2 = X : X-RICHTUNG IST AUSGEZEICHNET
C                                   2 = Y : Y-RICHTUNG IST AUSGEZEICHNET
C                                   2 = Z : Z-RICHTUNG IST AUSGEZEICHNET
C                                   A U S N A H M E N  BILDEN DIE FELDER
C                                   FUER DIE HAEUFIGKEITSVERTEILUNG.
C                                   "2" ENTFAELLT HIER
C                                   3 + (4) : 1. KOMPONENTE
C                                   (4) + 5 : 2. KOMPONENTE
C                                   A C H T U N G: DA IN FORTRAN NUR 6
C                                   BUCHSTABEN ZULAESSIG SIND, WIRD
C                                   BEI DEN OMEGA-KOMPONENTEN DAS 2. "O"
C                                   GESTRICHEN. (Z. B.: A12OXOY -->
C                                                       A12OXY     )
C
C UPROG                 : MAMILI, INTLIN
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        20.12.88 (HW)  : ORIGINAL
C
C*STARLET***************************************************************
C

      REAL     FVT    (0:IDIMF),           HVK    (KLASS,JJ2L)
      INTEGER  ISELEA (2,752),             ISELEP (2,752),
     $         ISLINP (ISLIDI),            ISLINA (ISLIDI),
     $         KSTW   (JJ  ,II  )
      REAL     RIDENT (100)
      REAL     HILF   (KK  ,JJ  ,II     )
C
      REAL     UFG    (KK  ,JJ  ,II     ), AURG   (KKA ,JJA ,IIA    ),
     $         VFG    (KK  ,JJ  ,II     ), AVRG   (KKA ,JJA ,IIA    ),
     $         WFG    (KK  ,JJ  ,II     ), AWRG   (KKA ,JJA ,IIA    ),
     $         OX     (KK  ,JJ  ,II     ), AOX    (KKA ,JJA ,IIA    ),
     $         OXFG   (KK  ,JJ  ,II     ), AOXRG  (KKA ,JJA ,IIA    ),
     $         OY     (KK  ,JJ  ,II     ), AOY    (KKA ,JJA ,IIA    ),
     $         OYFG   (KK  ,JJ  ,II     ), AOYRG  (KKA ,JJA ,IIA    ),
     $         OZ     (KK  ,JJ  ,II     ), AOZ    (KKA ,JJA ,IIA    ),
     $         OZFG   (KK  ,JJ  ,II     ), AOZRG  (KKA ,JJA ,IIA    )
      REAL     ARXUU  (KKXL,JJ2L,ILIMXP ), SRXUU  (KKXL,JJ1L,ILIMXP ),
     $         ARYUU  (KKYL,JJ2L,ILIMXP ), SRYUU  (KKYL,JJ1L,ILIMXP ),
     $         ARZUU  (KKZL,JJ2L,ILIMXP ), SRZUU  (KKZL,JJ1L,ILIMXP ),
     $         ASXUU  (KSXL,JJ1L,ILIMXP ), SSXUU  (KSXL,JJ1L,ILIMXP ),
     $         ASYUU  (KSYL,JJ1L,ILIMXP ), SSYUU  (KSYL,JJ1L,ILIMXP ),
     $         ASZUU  (KSZL,JJ1L,ILIMXP ), SSZUU  (KSZL,JJ1L,ILIMXP )
      REAL     ARXVV  (KKXL,JJ2L,ILIMXP ), SRXVV  (KKXL,JJ1L,ILIMXP ),
     $         ARYVV  (KKYL,JJ2L,ILIMXP ), SRYVV  (KKYL,JJ1L,ILIMXP ),
     $         ARZVV  (KKZL,JJ2L,ILIMXP ), SRZVV  (KKZL,JJ1L,ILIMXP ),
     $         ASXVV  (KSXL,JJ1L,ILIMXP ), SSXVV  (KSXL,JJ1L,ILIMXP ),
     $         ASYVV  (KSYL,JJ1L,ILIMXP ), SSYVV  (KSYL,JJ1L,ILIMXP ),
     $         ASZVV  (KSZL,JJ1L,ILIMXP ), SSZVV  (KSZL,JJ1L,ILIMXP )
      REAL     ARXWW  (KKXL,JJ2L,ILIMXP ), SRXWW  (KKXL,JJ1L,ILIMXP ),
     $         ARYWW  (KKYL,JJ2L,ILIMXP ), SRYWW  (KKYL,JJ1L,ILIMXP ),
     $         ARZWW  (KKZL,JJ2L,ILIMXP ), SRZWW  (KKZL,JJ1L,ILIMXP ),
     $         ASXWW  (KSXL,JJ1L,ILIMXP ), SSXWW  (KSXL,JJ1L,ILIMXP ),
     $         ASYWW  (KSYL,JJ1L,ILIMXP ), SSYWW  (KSYL,JJ1L,ILIMXP ),
     $         ASZWW  (KSZL,JJ1L,ILIMXP ), SSZWW  (KSZL,JJ1L,ILIMXP )
C
C                                  KREUZKORRELATIONSFUNKTIONEN
C                                  ---------------------------
C
      REAL     ARXUV  (KKXL,JJ2L,ILIMXP ), SRXUV  (KKXL,JJ1L,ILIMXP ),
     $         ARYUV  (KKYL,JJ2L,ILIMXP ), SRYUV  (KKYL,JJ1L,ILIMXP ),
     $         ARXUW  (KKXL,JJ2L,ILIMXP ), SRXUW  (KKXL,JJ1L,ILIMXP ),
     $         ARYUW  (KKYL,JJ2L,ILIMXP ), SRYUW  (KKYL,JJ1L,ILIMXP ),
     $         ARYVW  (KKYL,JJ2L,ILIMXP ), SRYVW  (KKYL,JJ1L,ILIMXP )
C
C                                  KORRELATIONSKOEFFIZIENTEN
C                                  -------------------------
C
      REAL     ACXUW  (KKXL,JJ2L,ILIMXP ), SCXUW  (KKXL,JJ1L,ILIMXP ),
     $         ACZUW  (KKZL,JJ2L,ILIMXP ), SCZUW  (KKZL,JJ1L,ILIMXP )
C
C                                  OBERE DREIECKSMATRIX DES TENSORS
C                                  DER KORRELATIONSFUNKTIONEN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  --------------------------------
C
      REAL     ARXOXX (KKXL,JJ2L,ILIMXP ), SRXOXX (KKXL,JJ1L,ILIMXP ),
     $         ARXOXY (KKXL,JJ2L,ILIMXP ), SRXOXY (KKXL,JJ1L,ILIMXP ),
     $         ARXOXZ (KKXL,JJ2L,ILIMXP ), SRXOXZ (KKXL,JJ1L,ILIMXP ),
     $         ARXOYY (KKXL,JJ2L,ILIMXP ), SRXOYY (KKXL,JJ1L,ILIMXP ),
     $         ARXOYZ (KKXL,JJ2L,ILIMXP ), SRXOYZ (KKXL,JJ1L,ILIMXP ),
     $         ARXOZZ (KKXL,JJ2L,ILIMXP ), SRXOZZ (KKXL,JJ1L,ILIMXP )
      REAL     ARYOXX (KKYL,JJ2L,ILIMXP ), SRYOXX (KKYL,JJ1L,ILIMXP ),
     $         ARYOXY (KKYL,JJ2L,ILIMXP ), SRYOXY (KKYL,JJ1L,ILIMXP ),
     $         ARYOXZ (KKYL,JJ2L,ILIMXP ), SRYOXZ (KKYL,JJ1L,ILIMXP ),
     $         ARYOYY (KKYL,JJ2L,ILIMXP ), SRYOYY (KKYL,JJ1L,ILIMXP ),
     $         ARYOYZ (KKYL,JJ2L,ILIMXP ), SRYOYZ (KKYL,JJ1L,ILIMXP ),
     $         ARYOZZ (KKYL,JJ2L,ILIMXP ), SRYOZZ (KKYL,JJ1L,ILIMXP )
      REAL     ARZOXX (KKZL,JJ2L,ILIMXP ), SRZOXX (KKZL,JJ1L,ILIMXP ),
     $         ARZOXY (KKZL,JJ2L,ILIMXP ), SRZOXY (KKZL,JJ1L,ILIMXP ),
     $         ARZOXZ (KKZL,JJ2L,ILIMXP ), SRZOXZ (KKZL,JJ1L,ILIMXP ),
     $         ARZOYY (KKZL,JJ2L,ILIMXP ), SRZOYY (KKZL,JJ1L,ILIMXP ),
     $         ARZOYZ (KKZL,JJ2L,ILIMXP ), SRZOYZ (KKZL,JJ1L,ILIMXP ),
     $         ARZOZZ (KKZL,JJ2L,ILIMXP ), SRZOZZ (KKZL,JJ1L,ILIMXP )
C
C                                  LEISTUNGSDICHTESPEKTREN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  ----------------------------
C
      REAL     ASXOXX (KSXL,JJ1L,ILIMXP ), SSXOXX (KSXL,JJ1L,ILIMXP ),
     $         ASXOYY (KSXL,JJ1L,ILIMXP ), SSXOYY (KSXL,JJ1L,ILIMXP ),
     $         ASXOZZ (KSXL,JJ1L,ILIMXP ), SSXOZZ (KSXL,JJ1L,ILIMXP ),
     $         ASYOXX (KSYL,JJ1L,ILIMXP ), SSYOXX (KSYL,JJ1L,ILIMXP ),
     $         ASYOYY (KSYL,JJ1L,ILIMXP ), SSYOYY (KSYL,JJ1L,ILIMXP ),
     $         ASYOZZ (KSYL,JJ1L,ILIMXP ), SSYOZZ (KSYL,JJ1L,ILIMXP )
C
C                                  HAEUFIGKEITSVERTEILUNG DER
C                                  INKLINATIONSWINKEL DER VORTICITY-
C                                  VEKTOREN
C                                  ---------------------------------
C
      REAL     AHOZOY (KH0L,JJ2L,ILIMXP ), SHOZOY (KH0L,JJ2L,ILIMXP ),
     $         AHOZOX (KH0L,JJ2L,ILIMXP ), SHOZOX (KH0L,JJ2L,ILIMXP ),
     $         AHOYOX (KH0L,JJ2L,ILIMXP ), SHOYOX (KH0L,JJ2L,ILIMXP )
C
C                                 STATIST. GROESSEN DER U-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 14) .EQ. 1 .OR. ISELEP(1, 14) .EQ. 2) THEN
         CALL MAMILI  (ARXUU  , 2 ,KKXL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTX )
         CALL INTLIN  (ARXUU  ,KKXL,JJ2L,ILIMXP,  2  ,ILINTX )
      ENDIF
      IF(ISELEP(1, 16) .EQ. 1 .OR. ISELEP(1, 16) .EQ. 2) THEN
         CALL MAMILI  (ARYUU  , 2 ,KKYL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTY )
         CALL INTLIN  (ARYUU  ,KKYL,JJ2L,ILIMXP,  2  ,ILINTY )
      ENDIF
      IF(ISELEP(1, 18) .EQ. 1 .OR. ISELEP(1, 18) .EQ. 2) THEN
         CALL MAMILI  (ARZUU  , 2 ,KKZL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTZ )
         CALL INTLIN  (ARZUU  ,KKZL,JJ2L,ILIMXP,  2  ,ILINTZ )
      ENDIF
      IF(ISELEP(1, 20) .EQ. 1 .OR. ISELEP(1, 20) .EQ. 2) THEN
         CALL MAMILI  (ASXUU  , 1 ,KSXL,JJ1L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTX )
      ENDIF
      IF(ISELEP(1, 22) .EQ. 1 .OR. ISELEP(1, 22) .EQ. 2) THEN
         CALL MAMILI  (ASYUU  , 1 ,KSYL,JJ1L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTY )
      ENDIF
      IF(ISELEP(1, 24) .EQ. 1 .OR. ISELEP(1, 24) .EQ. 2) THEN
         CALL MAMILI  (ASZUU  , 1 ,KSZL,JJ1L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTZ )
      ENDIF
C
C                                 STATIST. GROESSEN DER V-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 39) .EQ. 1 .OR. ISELEP(1, 39) .EQ. 2) THEN
         CALL MAMILI  (ARXVV  , 2 ,KKXL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTX )
         CALL INTLIN  (ARXVV  ,KKXL,JJ2L,ILIMXP,  2  ,ILINTX )
      ENDIF
      IF(ISELEP(1, 41) .EQ. 1 .OR. ISELEP(1, 41) .EQ. 2) THEN
         CALL MAMILI  (ARYVV  , 2 ,KKYL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTY )
         CALL INTLIN  (ARYVV  ,KKYL,JJ2L,ILIMXP,  2  ,ILINTY )
      ENDIF
      IF(ISELEP(1, 43) .EQ. 1 .OR. ISELEP(1, 43) .EQ. 2) THEN
         CALL MAMILI  (ARZVV  , 2 ,KKZL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTZ )
         CALL INTLIN  (ARZVV  ,KKZL,JJ2L,ILIMXP,  2  ,ILINTZ )
      ENDIF
      IF(ISELEP(1, 45) .EQ. 1 .OR. ISELEP(1, 45) .EQ. 2) THEN
         CALL MAMILI  (ASXVV  , 1 ,KSXL,JJ1L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTX )
      ENDIF
      IF(ISELEP(1, 47) .EQ. 1 .OR. ISELEP(1, 47) .EQ. 2) THEN
         CALL MAMILI  (ASYVV  , 1 ,KSYL,JJ1L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTY )
      ENDIF
      IF(ISELEP(1, 49) .EQ. 1 .OR. ISELEP(1, 49) .EQ. 2) THEN
         CALL MAMILI  (ASZVV  , 1 ,KSZL,JJ1L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTZ )
      ENDIF
C
C                                 STATIST. GROESSEN DER W-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 64) .EQ. 1 .OR. ISELEP(1, 64) .EQ. 2) THEN
         CALL MAMILI  (ARXWW  , 2 ,KKXL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTX )
         CALL INTLIN  (ARXWW  ,KKXL,JJ2L,ILIMXP,  2  ,ILINTX )
      ENDIF
      IF(ISELEP(1, 66) .EQ. 1 .OR. ISELEP(1, 66) .EQ. 2) THEN
         CALL MAMILI  (ARYWW  , 2 ,KKYL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTY )
         CALL INTLIN  (ARYWW  ,KKYL,JJ2L,ILIMXP,  2  ,ILINTY )
      ENDIF
      IF(ISELEP(1, 68) .EQ. 1 .OR. ISELEP(1, 68) .EQ. 2) THEN
         CALL MAMILI  (ARZWW  , 2 ,KKZL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTZ )
         CALL INTLIN  (ARZWW  ,KKZL,JJ2L,ILIMXP,  2  ,ILINTZ )
      ENDIF
      IF(ISELEP(1, 70) .EQ. 1 .OR. ISELEP(1, 70) .EQ. 2) THEN
         CALL MAMILI  (ASXWW  , 1 ,KSXL,JJ1L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTX )
      ENDIF
      IF(ISELEP(1, 72) .EQ. 1 .OR. ISELEP(1, 72) .EQ. 2) THEN
         CALL MAMILI  (ASYWW  , 1 ,KSYL,JJ1L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTY )
      ENDIF
      IF(ISELEP(1, 74) .EQ. 1 .OR. ISELEP(1, 74) .EQ. 2) THEN
         CALL MAMILI  (ASZWW  , 1 ,KSZL,JJ1L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTZ )
      ENDIF
C
C                                  HIER: KREUZKORRELATIONSFUNKTIONEN
C                                  ---------------------------------
C
      IF(ISELEP(1,210) .EQ. 1 .OR. ISELEP(1,210) .EQ. 2) THEN
         CALL MAMILI  (ARXUV  , 2 ,KKXL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTX )
         CALL INTLIN  (ARXUV  ,KKXL,JJ2L,ILIMXP,  2  ,ILINTX )
      ENDIF
      IF(ISELEP(1,212) .EQ. 1 .OR. ISELEP(1,212) .EQ. 2) THEN
         CALL MAMILI  (ARYUV  , 2 ,KKYL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTY )
         CALL INTLIN  (ARYUV  ,KKYL,JJ2L,ILIMXP,  2  ,ILINTY )
      ENDIF
      IF(ISELEP(1,216) .EQ. 1 .OR. ISELEP(1,216) .EQ. 2) THEN
         CALL MAMILI  (ARXUW  , 2 ,KKXL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTX )
         CALL INTLIN  (ARXUW  ,KKXL,JJ2L,ILIMXP,  2  ,ILINTX )
      ENDIF
      IF(ISELEP(1,218) .EQ. 1 .OR. ISELEP(1,218) .EQ. 2) THEN
         CALL MAMILI  (ARYUW  , 2 ,KKYL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTY )
         CALL INTLIN  (ARYUW  ,KKYL,JJ2L,ILIMXP,  2  ,ILINTY )
      ENDIF
      IF(ISELEP(1,224) .EQ. 1 .OR. ISELEP(1,224) .EQ. 2) THEN
         CALL MAMILI  (ARYVW  , 2 ,KKYL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTY )
         CALL INTLIN  (ARYVW  ,KKYL,JJ2L,ILIMXP,  2  ,ILINTY )
      ENDIF
C
C                                  HIER: KORRELATIONSKOEFFIZIENTEN
C                                  -------------------------------
C
      IF(ISELEP(1,234) .EQ. 1 .OR. ISELEP(1,234) .EQ. 2) THEN
         CALL MAMILI  (ACXUW  , 2 ,KKXL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTX )
      ENDIF
      IF(ISELEP(1,238) .EQ. 1 .OR. ISELEP(1,238) .EQ. 2) THEN
         CALL MAMILI  (ACZUW  , 2 ,KKZL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTZ )
      ENDIF
C
C                                  OBERE DREIECKSMATRIX DES TENSORS
C                                  DER KORRELATIONSFUNKTIONEN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  --------------------------------
C
      IF(ISELEP(1,246) .EQ. 1 .OR. ISELEP(1,246) .EQ. 2) THEN
         CALL MAMILI  (ARXOXX , 2 ,KKXL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTX )
         CALL INTLIN  (ARXOXX ,KKXL,JJ2L,ILIMXP,  2  ,ILINTX )
      ENDIF
      IF(ISELEP(1,248) .EQ. 1 .OR. ISELEP(1,248) .EQ. 2) THEN
         CALL MAMILI  (ARXOXY , 2 ,KKXL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTX )
         CALL INTLIN  (ARXOXY ,KKXL,JJ2L,ILIMXP,  2  ,ILINTX )
      ENDIF
      IF(ISELEP(1,250) .EQ. 1 .OR. ISELEP(1,250) .EQ. 2) THEN
         CALL MAMILI  (ARXOXZ , 2 ,KKXL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTX )
         CALL INTLIN  (ARXOXZ ,KKXL,JJ2L,ILIMXP,  2  ,ILINTX )
      ENDIF
      IF(ISELEP(1,252) .EQ. 1 .OR. ISELEP(1,252) .EQ. 2) THEN
         CALL MAMILI  (ARXOYY , 2 ,KKXL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTX )
         CALL INTLIN  (ARXOYY ,KKXL,JJ2L,ILIMXP,  2  ,ILINTX )
      ENDIF
      IF(ISELEP(1,254) .EQ. 1 .OR. ISELEP(1,254) .EQ. 2) THEN
         CALL MAMILI  (ARXOYZ , 2 ,KKXL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTX )
         CALL INTLIN  (ARXOYZ ,KKXL,JJ2L,ILIMXP,  2  ,ILINTX )
      ENDIF
      IF(ISELEP(1,256) .EQ. 1 .OR. ISELEP(1,256) .EQ. 2) THEN
         CALL MAMILI  (ARXOZZ , 2 ,KKXL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTX )
         CALL INTLIN  (ARXOZZ ,KKXL,JJ2L,ILIMXP,  2  ,ILINTX )
      ENDIF
      IF(ISELEP(1,258) .EQ. 1 .OR. ISELEP(1,258) .EQ. 2) THEN
         CALL MAMILI  (ARYOXX , 2 ,KKYL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTY )
         CALL INTLIN  (ARYOXX ,KKYL,JJ2L,ILIMXP,  2  ,ILINTY )
      ENDIF
      IF(ISELEP(1,260) .EQ. 1 .OR. ISELEP(1,260) .EQ. 2) THEN
         CALL MAMILI  (ARYOXY , 2 ,KKYL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTY )
         CALL INTLIN  (ARYOXY ,KKYL,JJ2L,ILIMXP,  2  ,ILINTY )
      ENDIF
      IF(ISELEP(1,262) .EQ. 1 .OR. ISELEP(1,262) .EQ. 2) THEN
         CALL MAMILI  (ARYOXZ , 2 ,KKYL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTY )
         CALL INTLIN  (ARYOXZ ,KKYL,JJ2L,ILIMXP,  2  ,ILINTY )
      ENDIF
      IF(ISELEP(1,264) .EQ. 1 .OR. ISELEP(1,264) .EQ. 2) THEN
         CALL MAMILI  (ARYOYY , 2 ,KKYL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTY )
         CALL INTLIN  (ARYOYY ,KKYL,JJ2L,ILIMXP,  2  ,ILINTY )
      ENDIF
      IF(ISELEP(1,266) .EQ. 1 .OR. ISELEP(1,266) .EQ. 2) THEN
         CALL MAMILI  (ARYOYZ , 2 ,KKYL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTY )
         CALL INTLIN  (ARYOYZ ,KKYL,JJ2L,ILIMXP,  2  ,ILINTY )
      ENDIF
      IF(ISELEP(1,268) .EQ. 1 .OR. ISELEP(1,268) .EQ. 2) THEN
         CALL MAMILI  (ARYOZZ , 2 ,KKYL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTY )
         CALL INTLIN  (ARYOZZ ,KKYL,JJ2L,ILIMXP,  2  ,ILINTY )
      ENDIF
      IF(ISELEP(1,270) .EQ. 1 .OR. ISELEP(1,270) .EQ. 2) THEN
         CALL MAMILI  (ARZOXX , 2 ,KKZL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTZ )
         CALL INTLIN  (ARZOXX ,KKZL,JJ2L,ILIMXP,  2  ,ILINTZ )
      ENDIF
      IF(ISELEP(1,272) .EQ. 1 .OR. ISELEP(1,272) .EQ. 2) THEN
         CALL MAMILI  (ARZOXY , 2 ,KKZL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTZ )
         CALL INTLIN  (ARZOXY ,KKZL,JJ2L,ILIMXP,  2  ,ILINTZ )
      ENDIF
      IF(ISELEP(1,274) .EQ. 1 .OR. ISELEP(1,274) .EQ. 2) THEN
         CALL MAMILI  (ARZOXZ , 2 ,KKZL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTZ )
         CALL INTLIN  (ARZOXZ ,KKZL,JJ2L,ILIMXP,  2  ,ILINTZ )
      ENDIF
      IF(ISELEP(1,276) .EQ. 1 .OR. ISELEP(1,276) .EQ. 2) THEN
         CALL MAMILI  (ARZOYY , 2 ,KKZL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTZ )
         CALL INTLIN  (ARZOYY ,KKZL,JJ2L,ILIMXP,  2  ,ILINTZ )
      ENDIF
      IF(ISELEP(1,278) .EQ. 1 .OR. ISELEP(1,278) .EQ. 2) THEN
         CALL MAMILI  (ARZOYZ , 2 ,KKZL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTZ )
         CALL INTLIN  (ARZOYZ ,KKZL,JJ2L,ILIMXP,  2  ,ILINTZ )
      ENDIF
      IF(ISELEP(1,280) .EQ. 1 .OR. ISELEP(1,280) .EQ. 2) THEN
         CALL MAMILI  (ARZOZZ , 2 ,KKZL,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTZ )
         CALL INTLIN  (ARZOZZ ,KKZL,JJ2L,ILIMXP,  2  ,ILINTZ )
      ENDIF
C
C                                  LEISTUNGSDICHTESPEKTREN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  ----------------------------
C
      IF(ISELEP(1,282) .EQ. 1 .OR. ISELEP(1,282) .EQ. 2) THEN
         CALL MAMILI  (ASXOXX , 1 ,KSXL,JJ1L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTX )
      ENDIF
      IF(ISELEP(1,284) .EQ. 1 .OR. ISELEP(1,284) .EQ. 2) THEN
         CALL MAMILI  (ASXOYY , 1 ,KSXL,JJ1L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTX )
      ENDIF
      IF(ISELEP(1,286) .EQ. 1 .OR. ISELEP(1,286) .EQ. 2) THEN
         CALL MAMILI  (ASXOZZ , 1 ,KSXL,JJ1L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTX )
      ENDIF
      IF(ISELEP(1,288) .EQ. 1 .OR. ISELEP(1,288) .EQ. 2) THEN
         CALL MAMILI  (ASYOXX , 1 ,KSYL,JJ1L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTY )
      ENDIF
      IF(ISELEP(1,290) .EQ. 1 .OR. ISELEP(1,290) .EQ. 2) THEN
         CALL MAMILI  (ASYOYY , 1 ,KSYL,JJ1L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTY )
      ENDIF
      IF(ISELEP(1,292) .EQ. 1 .OR. ISELEP(1,292) .EQ. 2) THEN
         CALL MAMILI  (ASYOZZ , 1 ,KSYL,JJ1L,ILIMXP,
     $                 1   ,KEND,  1  ,  1  ,  2  ,ILINTY )
      ENDIF
C
C                                  HAEUFIGKEITSVERTEILUNG DER
C                                  INKLINATIONSWINKEL DER VORTICITY-
C                                  VEKTOREN
C                                  ---------------------------------
C
      IF(ISELEP(1,306) .EQ. 1 .OR. ISELEP(1,306) .EQ. 2) THEN
         CALL MAMILI  (AHOZOY , 1 ,KH0L,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  2  ,  2  ,ILINT0 )
      ENDIF
      IF(ISELEP(1,308) .EQ. 1 .OR. ISELEP(1,308) .EQ. 2) THEN
         CALL MAMILI  (AHOZOX , 1 ,KH0L,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  2  ,  2  ,ILINT0 )
      ENDIF
      IF(ISELEP(1,310) .EQ. 1 .OR. ISELEP(1,310) .EQ. 2) THEN
         CALL MAMILI  (AHOYOX , 1 ,KH0L,JJ2L,ILIMXP,
     $                 1   ,KEND,  1  ,  2  ,  2  ,ILINT0 )
      ENDIF
C
      RETURN
      END
