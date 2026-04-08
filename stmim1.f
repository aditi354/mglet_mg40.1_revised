










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
      SUBROUTINE STMIM1  (DREAD,DCONT,NPRNEU,FPRNEU,
     $                    ILINT0,ILINTX,ILINTY,ILINTZ,

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
C        S T M I M 1      IN STMIM1 WERDEN DIE ALTEN ENSEMBLE-
C                         MITTELWERTE VON KORRELATIONSFUNKTIONEN,
C                         -KOEFFIZIENTEN, LEISTUNGSDICHTESPEKTREN UND
C                         HAEUFIGKEITSVERTEILUNGEN (ALLGEMEIN: ALLE
C                         "LINIEN"-FELDER) VERBESSERT. DIE GEWICHTUNG
C                         DER ALTEN (WAEHREND DER VORANGEGANGENEN LAEUFE
C                         ERZEUGT) UND DER NEUEN ENSEMBLE-MITTELWERTE
C                         (WAEHREND DES MOMENTANEN LAUFES ERZEUGT)
C                         ERGIBT SICH AUS DEM VERHAELTNIS DER ANZAHL
C                         DER STICHPROBEN, DIE ZU DEN JEWEILIGEN
C                         MITTELWERTEN FUEHRTEN.
C
C                         STMIM1 DIENT DAZU, DIE FELDER AUSZUWAEHLEN,
C                         DEREN ENSEMBLE-MITTELWERT VERBESSERT WERDEN
C                         SOLL. DIE EIGENTLICHE ARITHMET. OPERATION
C                         FINDET IN SUBR. STMNLI STATT.
C  A C H T U N G:         DIESE ROUTINE MUSS ERWEITERT WERDEN, FALLS
C                         DIE STATIST. AUSWERTUNG ERWEITERT WIRD !
C*STARLET***************************************************************
C
C PARAM: DREAD          - DREAD = .T. : DIE WAEHREND DES VORANGEGANGENEN
C                         LAUFES ERZEUGTEN ENSEMBLE-MITTELWERTE WERDEN
C                         EINGELESEN.
C                         DREAD = .F. : ES FINDET EINE VORBELEGUNG MIT
C                         "BESTMOEGLICHEN" ENSEMBLE-MITTELWERTEN STATT.
C        DCONT          - DCONT = .T. .AND. DREAD = .T. : DIE GEWICH-
C                         TUNG DER "ALTEN" ENSEMBLE-MITTELWERTE ERFOLGT
C                         ENTSPRECHEND DER ANZAHL DER "ALTEN" STICH-
C                         PROBEN.
C
C                         DCONT = .F. .AND. DREAD = .T. : DIE "ALTEN"
C                         ENSEMBLE-MITTELWERTE WERDEN MIT NPRNEU
C                         STICHPROBEN GEWICHTET.
C        NPRNEU         - SIEHE BESCHREIBUNG  "DCONT"
C        ILINT0         - ANZAHL DER AUFPUNKTE MIT RICHTUNGSINFORMATION
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
C UPROG                 : STMNLI
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        14.12.88 (HW)  : ORIGINAL AUS STMIMP ABGELEITET
C
C*STARLET***************************************************************
C
      LOGICAL  DREAD,  DCONT
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
         NPROLD = ISELEA(2, 14)
         IF(ISELEP(2, 14) .GT. 0) THEN
            NPRRUS = ISELEP(2, 14)
            CALL STMNLI  (ARXUU  ,JJ2L,SRXUU  ,JJ1L,KKXL,ILIMXP,
     $                    2  ,ILINTX, 1 , 1 , 1 ,KKXL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2, 14) = NPROLD
         ISELEP(2, 14) = ISELEP(2, 14) + NPROLD
      ENDIF
      IF(ISELEP(1, 16) .EQ. 1 .OR. ISELEP(1, 16) .EQ. 2) THEN
         NPROLD = ISELEA(2, 16)
         IF(ISELEP(2, 16) .GT. 0) THEN
            NPRRUS = ISELEP(2, 16)
            CALL STMNLI  (ARYUU  ,JJ2L,SRYUU  ,JJ1L,KKYL,ILIMXP,
     $                    2  ,ILINTY, 1 , 1 , 1 ,KKYL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2, 16) = NPROLD
         ISELEP(2, 16) = ISELEP(2, 16) + NPROLD
      ENDIF
      IF(ISELEP(1, 18) .EQ. 1 .OR. ISELEP(1, 18) .EQ. 2) THEN
         NPROLD = ISELEA(2, 18)
         IF(ISELEP(2, 18) .GT. 0) THEN
            NPRRUS = ISELEP(2, 18)
            CALL STMNLI  (ARZUU  ,JJ2L,SRZUU  ,JJ1L,KKZL,ILIMXP,
     $                    2  ,ILINTZ, 1 , 1 , 1 ,KKZL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2, 18) = NPROLD
         ISELEP(2, 18) = ISELEP(2, 18) + NPROLD
      ENDIF
      IF(ISELEP(1, 20) .EQ. 1 .OR. ISELEP(1, 20) .EQ. 2) THEN
         NPROLD = ISELEA(2, 20)
         IF(ISELEP(2, 20) .GT. 0) THEN
            NPRRUS = ISELEP(2, 20)
            CALL STMNLI  (ASXUU  ,JJ1L,SSXUU  ,JJ1L,KSXL,ILIMXP,
     $                    2  ,ILINTX, 1 , 1 , 1 ,KSXL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2, 20) = NPROLD
         ISELEP(2, 20) = ISELEP(2, 20) + NPROLD
      ENDIF
      IF(ISELEP(1, 22) .EQ. 1 .OR. ISELEP(1, 22) .EQ. 2) THEN
         NPROLD = ISELEA(2, 22)
         IF(ISELEP(2, 22) .GT. 0) THEN
            NPRRUS = ISELEP(2, 22)
            CALL STMNLI  (ASYUU  ,JJ1L,SSYUU  ,JJ1L,KSYL,ILIMXP,
     $                    2  ,ILINTY, 1 , 1 , 1 ,KSYL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2, 22) = NPROLD
         ISELEP(2, 22) = ISELEP(2, 22) + NPROLD
      ENDIF
      IF(ISELEP(1, 24) .EQ. 1 .OR. ISELEP(1, 24) .EQ. 2) THEN
         NPROLD = ISELEA(2, 24)
         IF(ISELEP(2, 24) .GT. 0) THEN
            NPRRUS = ISELEP(2, 24)
            CALL STMNLI  (ASZUU  ,JJ1L,SSZUU  ,JJ1L,KSZL,ILIMXP,
     $                    2  ,ILINTZ, 1 , 1 , 1 ,KSZL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2, 24) = NPROLD
         ISELEP(2, 24) = ISELEP(2, 24) + NPROLD
      ENDIF
C
C                                 STATIST. GROESSEN DER V-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 39) .EQ. 1 .OR. ISELEP(1, 39) .EQ. 2) THEN
         NPROLD = ISELEA(2, 39)
         IF(ISELEP(2, 39) .GT. 0) THEN
            NPRRUS = ISELEP(2, 39)
            CALL STMNLI  (ARXVV  ,JJ2L,SRXVV  ,JJ1L,KKXL,ILIMXP,
     $                    2  ,ILINTX, 1 , 1 , 1 ,KKXL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2, 39) = NPROLD
         ISELEP(2, 39) = ISELEP(2, 39) + NPROLD
      ENDIF
      IF(ISELEP(1, 41) .EQ. 1 .OR. ISELEP(1, 41) .EQ. 2) THEN
         NPROLD = ISELEA(2, 41)
         IF(ISELEP(2, 41) .GT. 0) THEN
            NPRRUS = ISELEP(2, 41)
            CALL STMNLI  (ARYVV  ,JJ2L,SRYVV  ,JJ1L,KKYL,ILIMXP,
     $                    2  ,ILINTY, 1 , 1 , 1 ,KKYL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2, 41) = NPROLD
         ISELEP(2, 41) = ISELEP(2, 41) + NPROLD
      ENDIF
      IF(ISELEP(1, 43) .EQ. 1 .OR. ISELEP(1, 43) .EQ. 2) THEN
         NPROLD = ISELEA(2, 43)
         IF(ISELEP(2, 43) .GT. 0) THEN
            NPRRUS = ISELEP(2, 43)
            CALL STMNLI  (ARZVV  ,JJ2L,SRZVV  ,JJ1L,KKZL,ILIMXP,
     $                    2  ,ILINTZ, 1 , 1 , 1 ,KKZL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2, 43) = NPROLD
         ISELEP(2, 43) = ISELEP(2, 43) + NPROLD
      ENDIF
      IF(ISELEP(1, 45) .EQ. 1 .OR. ISELEP(1, 45) .EQ. 2) THEN
         NPROLD = ISELEA(2, 45)
         IF(ISELEP(2, 45) .GT. 0) THEN
            NPRRUS = ISELEP(2, 45)
            CALL STMNLI  (ASXVV  ,JJ1L,SSXVV  ,JJ1L,KSXL,ILIMXP,
     $                    2  ,ILINTX, 1 , 1 , 1 ,KSXL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2, 45) = NPROLD
         ISELEP(2, 45) = ISELEP(2, 45) + NPROLD
      ENDIF
      IF(ISELEP(1, 47) .EQ. 1 .OR. ISELEP(1, 47) .EQ. 2) THEN
         NPROLD = ISELEA(2, 47)
         IF(ISELEP(2, 47) .GT. 0) THEN
            NPRRUS = ISELEP(2, 47)
            CALL STMNLI  (ASYVV  ,JJ1L,SSYVV  ,JJ1L,KSYL,ILIMXP,
     $                    2  ,ILINTY, 1 , 1 , 1 ,KSYL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2, 47) = NPROLD
         ISELEP(2, 47) = ISELEP(2, 47) + NPROLD
      ENDIF
      IF(ISELEP(1, 49) .EQ. 1 .OR. ISELEP(1, 49) .EQ. 2) THEN
         NPROLD = ISELEA(2, 49)
         IF(ISELEP(2, 49) .GT. 0) THEN
            NPRRUS = ISELEP(2, 49)
            CALL STMNLI  (ASZVV  ,JJ1L,SSZVV  ,JJ1L,KSZL,ILIMXP,
     $                    2  ,ILINTZ, 1 , 1 , 1 ,KSZL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2, 49) = NPROLD
         ISELEP(2, 49) = ISELEP(2, 49) + NPROLD
      ENDIF
C
C                                 STATIST. GROESSEN DER W-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 64) .EQ. 1 .OR. ISELEP(1, 64) .EQ. 2) THEN
         NPROLD = ISELEA(2, 64)
         IF(ISELEP(2, 64) .GT. 0) THEN
            NPRRUS = ISELEP(2, 64)
            CALL STMNLI  (ARXWW  ,JJ2L,SRXWW  ,JJ1L,KKXL,ILIMXP,
     $                    2  ,ILINTX, 1 , 1 , 1 ,KKXL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2, 64) = NPROLD
         ISELEP(2, 64) = ISELEP(2, 64) + NPROLD
      ENDIF
      IF(ISELEP(1, 66) .EQ. 1 .OR. ISELEP(1, 66) .EQ. 2) THEN
         NPROLD = ISELEA(2, 66)
         IF(ISELEP(2, 66) .GT. 0) THEN
            NPRRUS = ISELEP(2, 66)
            CALL STMNLI  (ARYWW  ,JJ2L,SRYWW  ,JJ1L,KKYL,ILIMXP,
     $                    2  ,ILINTY, 1 , 1 , 1 ,KKYL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2, 66) = NPROLD
         ISELEP(2, 66) = ISELEP(2, 66) + NPROLD
      ENDIF
      IF(ISELEP(1, 68) .EQ. 1 .OR. ISELEP(1, 68) .EQ. 2) THEN
         NPROLD = ISELEA(2, 68)
         IF(ISELEP(2, 68) .GT. 0) THEN
            NPRRUS = ISELEP(2, 68)
            CALL STMNLI  (ARZWW  ,JJ2L,SRZWW  ,JJ1L,KKZL,ILIMXP,
     $                    2  ,ILINTZ, 1 , 1 , 1 ,KKZL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2, 68) = NPROLD
         ISELEP(2, 68) = ISELEP(2, 68) + NPROLD
      ENDIF
      IF(ISELEP(1, 70) .EQ. 1 .OR. ISELEP(1, 70) .EQ. 2) THEN
         NPROLD = ISELEA(2, 70)
         IF(ISELEP(2, 70) .GT. 0) THEN
            NPRRUS = ISELEP(2, 70)
            CALL STMNLI  (ASXWW  ,JJ1L,SSXWW  ,JJ1L,KSXL,ILIMXP,
     $                    2  ,ILINTX, 1 , 1 , 1 ,KSXL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2, 70) = NPROLD
         ISELEP(2, 70) = ISELEP(2, 70) + NPROLD
      ENDIF
      IF(ISELEP(1, 72) .EQ. 1 .OR. ISELEP(1, 72) .EQ. 2) THEN
         NPROLD = ISELEA(2, 72)
         IF(ISELEP(2, 72) .GT. 0) THEN
            NPRRUS = ISELEP(2, 72)
            CALL STMNLI  (ASYWW  ,JJ1L,SSYWW  ,JJ1L,KSYL,ILIMXP,
     $                    2  ,ILINTY, 1 , 1 , 1 ,KSYL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2, 72) = NPROLD
         ISELEP(2, 72) = ISELEP(2, 72) + NPROLD
      ENDIF
      IF(ISELEP(1, 74) .EQ. 1 .OR. ISELEP(1, 74) .EQ. 2) THEN
         NPROLD = ISELEA(2, 74)
         IF(ISELEP(2, 74) .GT. 0) THEN
            NPRRUS = ISELEP(2, 74)
            CALL STMNLI  (ASZWW  ,JJ1L,SSZWW  ,JJ1L,KSZL,ILIMXP,
     $                    2  ,ILINTZ, 1 , 1 , 1 ,KSZL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2, 74) = NPROLD
         ISELEP(2, 74) = ISELEP(2, 74) + NPROLD
      ENDIF
C
C                                  HIER: KREUZKORRELATIONSFUNKTIONEN
C                                  ---------------------------------
C
      IF(ISELEP(1,210) .EQ. 1 .OR. ISELEP(1,210) .EQ. 2) THEN
         NPROLD = ISELEA(2,210)
         IF(ISELEP(2,210) .GT. 0) THEN
            NPRRUS = ISELEP(2,210)
            CALL STMNLI  (ARXUV  ,JJ2L,SRXUV  ,JJ1L,KKXL,ILIMXP,
     $                    2  ,ILINTX, 1 , 1 , 1 ,KKXL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,210) = NPROLD
         ISELEP(2,210) = ISELEP(2,210) + NPROLD
      ENDIF
      IF(ISELEP(1,212) .EQ. 1 .OR. ISELEP(1,212) .EQ. 2) THEN
         NPROLD = ISELEA(2,212)
         IF(ISELEP(2,212) .GT. 0) THEN
            NPRRUS = ISELEP(2,212)
            CALL STMNLI  (ARYUV  ,JJ2L,SRYUV  ,JJ1L,KKYL,ILIMXP,
     $                    2  ,ILINTY, 1 , 1 , 1 ,KKYL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,212) = NPROLD
         ISELEP(2,212) = ISELEP(2,212) + NPROLD
      ENDIF
      IF(ISELEP(1,216) .EQ. 1 .OR. ISELEP(1,216) .EQ. 2) THEN
         NPROLD = ISELEA(2,216)
         IF(ISELEP(2,216) .GT. 0) THEN
            NPRRUS = ISELEP(2,216)
            CALL STMNLI  (ARXUW  ,JJ2L,SRXUW  ,JJ1L,KKXL,ILIMXP,
     $                    2  ,ILINTX, 1 , 1 , 1 ,KKXL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,216) = NPROLD
         ISELEP(2,216) = ISELEP(2,216) + NPROLD
      ENDIF
      IF(ISELEP(1,218) .EQ. 1 .OR. ISELEP(1,218) .EQ. 2) THEN
         NPROLD = ISELEA(2,218)
         IF(ISELEP(2,218) .GT. 0) THEN
            NPRRUS = ISELEP(2,218)
            CALL STMNLI  (ARYUW  ,JJ2L,SRYUW  ,JJ1L,KKYL,ILIMXP,
     $                    2  ,ILINTY, 1 , 1 , 1 ,KKYL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,218) = NPROLD
         ISELEP(2,218) = ISELEP(2,218) + NPROLD
      ENDIF
      IF(ISELEP(1,224) .EQ. 1 .OR. ISELEP(1,224) .EQ. 2) THEN
         NPROLD = ISELEA(2,224)
         IF(ISELEP(2,224) .GT. 0) THEN
            NPRRUS = ISELEP(2,224)
            CALL STMNLI  (ARYVW  ,JJ2L,SRYVW  ,JJ1L,KKYL,ILIMXP,
     $                    2  ,ILINTY, 1 , 1 , 1 ,KKYL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,224) = NPROLD
         ISELEP(2,224) = ISELEP(2,224) + NPROLD
      ENDIF
C
C                                  HIER: KORRELATIONSKOEFFIZIENTEN
C                                  -------------------------------
C
      IF(ISELEP(1,234) .EQ. 1 .OR. ISELEP(1,234) .EQ. 2) THEN
         NPROLD = ISELEA(2,234)
         IF(ISELEP(2,234) .GT. 0) THEN
            NPRRUS = ISELEP(2,234)
            CALL STMNLI  (ACXUW  ,JJ2L,SCXUW  ,JJ1L,KKXL,ILIMXP,
     $                    2  ,ILINTX, 1 , 1 , 1 ,KKXL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,234) = NPROLD
         ISELEP(2,234) = ISELEP(2,234) + NPROLD
      ENDIF
      IF(ISELEP(1,238) .EQ. 1 .OR. ISELEP(1,238) .EQ. 2) THEN
         NPROLD = ISELEA(2,238)
         IF(ISELEP(2,238) .GT. 0) THEN
            NPRRUS = ISELEP(2,238)
            CALL STMNLI  (ACZUW  ,JJ2L,SCZUW  ,JJ1L,KKZL,ILIMXP,
     $                    2  ,ILINTZ, 1 , 1 , 1 ,KKZL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,238) = NPROLD
         ISELEP(2,238) = ISELEP(2,238) + NPROLD
      ENDIF
C
C                                  OBERE DREIECKSMATRIX DES TENSORS
C                                  DER KORRELATIONSFUNKTIONEN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  --------------------------------
C
      IF(ISELEP(1,246) .EQ. 1 .OR. ISELEP(1,246) .EQ. 2) THEN
         NPROLD = ISELEA(2,246)
         IF(ISELEP(2,246) .GT. 0) THEN
            NPRRUS = ISELEP(2,246)
            CALL STMNLI  (ARXOXX ,JJ2L,SRXOXX ,JJ1L,KKXL,ILIMXP,
     $                    2  ,ILINTX, 1 , 1 , 1 ,KKXL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,246) = NPROLD
         ISELEP(2,246) = ISELEP(2,246) + NPROLD
      ENDIF
      IF(ISELEP(1,248) .EQ. 1 .OR. ISELEP(1,248) .EQ. 2) THEN
         NPROLD = ISELEA(2,248)
         IF(ISELEP(2,248) .GT. 0) THEN
            NPRRUS = ISELEP(2,248)
            CALL STMNLI  (ARXOXY ,JJ2L,SRXOXY ,JJ1L,KKXL,ILIMXP,
     $                    2  ,ILINTX, 1 , 1 , 1 ,KKXL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,248) = NPROLD
         ISELEP(2,248) = ISELEP(2,248) + NPROLD
      ENDIF
      IF(ISELEP(1,250) .EQ. 1 .OR. ISELEP(1,250) .EQ. 2) THEN
         NPROLD = ISELEA(2,250)
         IF(ISELEP(2,250) .GT. 0) THEN
            NPRRUS = ISELEP(2,250)
            CALL STMNLI  (ARXOXZ ,JJ2L,SRXOXZ ,JJ1L,KKXL,ILIMXP,
     $                    2  ,ILINTX, 1 , 1 , 1 ,KKXL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,250) = NPROLD
         ISELEP(2,250) = ISELEP(2,250) + NPROLD
      ENDIF
      IF(ISELEP(1,252) .EQ. 1 .OR. ISELEP(1,252) .EQ. 2) THEN
         NPROLD = ISELEA(2,252)
         IF(ISELEP(2,252) .GT. 0) THEN
            NPRRUS = ISELEP(2,252)
            CALL STMNLI  (ARXOYY ,JJ2L,SRXOYY ,JJ1L,KKXL,ILIMXP,
     $                    2  ,ILINTX, 1 , 1 , 1 ,KKXL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,252) = NPROLD
         ISELEP(2,252) = ISELEP(2,252) + NPROLD
      ENDIF
      IF(ISELEP(1,254) .EQ. 1 .OR. ISELEP(1,254) .EQ. 2) THEN
         NPROLD = ISELEA(2,254)
         IF(ISELEP(2,254) .GT. 0) THEN
            NPRRUS = ISELEP(2,254)
            CALL STMNLI  (ARXOYZ ,JJ2L,SRXOYZ ,JJ1L,KKXL,ILIMXP,
     $                    2  ,ILINTX, 1 , 1 , 1 ,KKXL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,254) = NPROLD
         ISELEP(2,254) = ISELEP(2,254) + NPROLD
      ENDIF
      IF(ISELEP(1,256) .EQ. 1 .OR. ISELEP(1,256) .EQ. 2) THEN
         NPROLD = ISELEA(2,256)
         IF(ISELEP(2,256) .GT. 0) THEN
            NPRRUS = ISELEP(2,256)
            CALL STMNLI  (ARXOZZ ,JJ2L,SRXOZZ ,JJ1L,KKXL,ILIMXP,
     $                    2  ,ILINTX, 1 , 1 , 1 ,KKXL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,256) = NPROLD
         ISELEP(2,256) = ISELEP(2,256) + NPROLD
      ENDIF
      IF(ISELEP(1,258) .EQ. 1 .OR. ISELEP(1,258) .EQ. 2) THEN
         NPROLD = ISELEA(2,258)
         IF(ISELEP(2,258) .GT. 0) THEN
            NPRRUS = ISELEP(2,258)
            CALL STMNLI  (ARYOXX ,JJ2L,SRYOXX ,JJ1L,KKYL,ILIMXP,
     $                    2  ,ILINTY, 1 , 1 , 1 ,KKYL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,258) = NPROLD
         ISELEP(2,258) = ISELEP(2,258) + NPROLD
      ENDIF
      IF(ISELEP(1,260) .EQ. 1 .OR. ISELEP(1,260) .EQ. 2) THEN
         NPROLD = ISELEA(2,260)
         IF(ISELEP(2,260) .GT. 0) THEN
            NPRRUS = ISELEP(2,260)
            CALL STMNLI  (ARYOXY ,JJ2L,SRYOXY ,JJ1L,KKYL,ILIMXP,
     $                    2  ,ILINTY, 1 , 1 , 1 ,KKYL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,260) = NPROLD
         ISELEP(2,260) = ISELEP(2,260) + NPROLD
      ENDIF
      IF(ISELEP(1,262) .EQ. 1 .OR. ISELEP(1,262) .EQ. 2) THEN
         NPROLD = ISELEA(2,262)
         IF(ISELEP(2,262) .GT. 0) THEN
            NPRRUS = ISELEP(2,262)
            CALL STMNLI  (ARYOXZ ,JJ2L,SRYOXZ ,JJ1L,KKYL,ILIMXP,
     $                    2  ,ILINTY, 1 , 1 , 1 ,KKYL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,262) = NPROLD
         ISELEP(2,262) = ISELEP(2,262) + NPROLD
      ENDIF
      IF(ISELEP(1,264) .EQ. 1 .OR. ISELEP(1,264) .EQ. 2) THEN
         NPROLD = ISELEA(2,264)
         IF(ISELEP(2,264) .GT. 0) THEN
            NPRRUS = ISELEP(2,264)
            CALL STMNLI  (ARYOYY ,JJ2L,SRYOYY ,JJ1L,KKYL,ILIMXP,
     $                    2  ,ILINTY, 1 , 1 , 1 ,KKYL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,264) = NPROLD
         ISELEP(2,264) = ISELEP(2,264) + NPROLD
      ENDIF
      IF(ISELEP(1,266) .EQ. 1 .OR. ISELEP(1,266) .EQ. 2) THEN
         NPROLD = ISELEA(2,266)
         IF(ISELEP(2,266) .GT. 0) THEN
            NPRRUS = ISELEP(2,266)
            CALL STMNLI  (ARYOYZ ,JJ2L,SRYOYZ ,JJ1L,KKYL,ILIMXP,
     $                    2  ,ILINTY, 1 , 1 , 1 ,KKYL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,266) = NPROLD
         ISELEP(2,266) = ISELEP(2,266) + NPROLD
      ENDIF
      IF(ISELEP(1,268) .EQ. 1 .OR. ISELEP(1,268) .EQ. 2) THEN
         NPROLD = ISELEA(2,268)
         IF(ISELEP(2,268) .GT. 0) THEN
            NPRRUS = ISELEP(2,268)
            CALL STMNLI  (ARYOZZ ,JJ2L,SRYOZZ ,JJ1L,KKYL,ILIMXP,
     $                    2  ,ILINTY, 1 , 1 , 1 ,KKYL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,268) = NPROLD
         ISELEP(2,268) = ISELEP(2,268) + NPROLD
      ENDIF
      IF(ISELEP(1,270) .EQ. 1 .OR. ISELEP(1,270) .EQ. 2) THEN
         NPROLD = ISELEA(2,270)
         IF(ISELEP(2,270) .GT. 0) THEN
            NPRRUS = ISELEP(2,270)
            CALL STMNLI  (ARZOXX ,JJ2L,SRZOXX ,JJ1L,KKZL,ILIMXP,
     $                    2  ,ILINTZ, 1 , 1 , 1 ,KKZL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,270) = NPROLD
         ISELEP(2,270) = ISELEP(2,270) + NPROLD
      ENDIF
      IF(ISELEP(1,272) .EQ. 1 .OR. ISELEP(1,272) .EQ. 2) THEN
         NPROLD = ISELEA(2,272)
         IF(ISELEP(2,272) .GT. 0) THEN
            NPRRUS = ISELEP(2,272)
            CALL STMNLI  (ARZOXY ,JJ2L,SRZOXY ,JJ1L,KKZL,ILIMXP,
     $                    2  ,ILINTZ, 1 , 1 , 1 ,KKZL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,272) = NPROLD
         ISELEP(2,272) = ISELEP(2,272) + NPROLD
      ENDIF
      IF(ISELEP(1,274) .EQ. 1 .OR. ISELEP(1,274) .EQ. 2) THEN
         NPROLD = ISELEA(2,274)
         IF(ISELEP(2,274) .GT. 0) THEN
            NPRRUS = ISELEP(2,274)
            CALL STMNLI  (ARZOXZ ,JJ2L,SRZOXZ ,JJ1L,KKZL,ILIMXP,
     $                    2  ,ILINTZ, 1 , 1 , 1 ,KKZL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,274) = NPROLD
         ISELEP(2,274) = ISELEP(2,274) + NPROLD
      ENDIF
      IF(ISELEP(1,276) .EQ. 1 .OR. ISELEP(1,276) .EQ. 2) THEN
         NPROLD = ISELEA(2,276)
         IF(ISELEP(2,276) .GT. 0) THEN
            NPRRUS = ISELEP(2,276)
            CALL STMNLI  (ARZOYY ,JJ2L,SRZOYY ,JJ1L,KKZL,ILIMXP,
     $                    2  ,ILINTZ, 1 , 1 , 1 ,KKZL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,276) = NPROLD
         ISELEP(2,276) = ISELEP(2,276) + NPROLD
      ENDIF
      IF(ISELEP(1,278) .EQ. 1 .OR. ISELEP(1,278) .EQ. 2) THEN
         NPROLD = ISELEA(2,278)
         IF(ISELEP(2,278) .GT. 0) THEN
            NPRRUS = ISELEP(2,278)
            CALL STMNLI  (ARZOYZ ,JJ2L,SRZOYZ ,JJ1L,KKZL,ILIMXP,
     $                    2  ,ILINTZ, 1 , 1 , 1 ,KKZL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,278) = NPROLD
         ISELEP(2,278) = ISELEP(2,278) + NPROLD
      ENDIF
      IF(ISELEP(1,280) .EQ. 1 .OR. ISELEP(1,280) .EQ. 2) THEN
         NPROLD = ISELEA(2,280)
         IF(ISELEP(2,280) .GT. 0) THEN
            NPRRUS = ISELEP(2,280)
            CALL STMNLI  (ARZOZZ ,JJ2L,SRZOZZ ,JJ1L,KKZL,ILIMXP,
     $                    2  ,ILINTZ, 1 , 1 , 1 ,KKZL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,280) = NPROLD
         ISELEP(2,280) = ISELEP(2,280) + NPROLD
      ENDIF
C
C                                  LEISTUNGSDICHTESPEKTREN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  ----------------------------
C
      IF(ISELEP(1,282) .EQ. 1 .OR. ISELEP(1,282) .EQ. 2) THEN
         NPROLD = ISELEA(2,282)
         IF(ISELEP(2,282) .GT. 0) THEN
            NPRRUS = ISELEP(2,282)
            CALL STMNLI  (ASXOXX ,JJ1L,SSXOXX ,JJ1L,KSXL,ILIMXP,
     $                    2  ,ILINTX, 1 , 1 , 1 ,KSXL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,282) = NPROLD
         ISELEP(2,282) = ISELEP(2,282) + NPROLD
      ENDIF
      IF(ISELEP(1,284) .EQ. 1 .OR. ISELEP(1,284) .EQ. 2) THEN
         NPROLD = ISELEA(2,284)
         IF(ISELEP(2,284) .GT. 0) THEN
            NPRRUS = ISELEP(2,284)
            CALL STMNLI  (ASXOYY ,JJ1L,SSXOYY ,JJ1L,KSXL,ILIMXP,
     $                    2  ,ILINTX, 1 , 1 , 1 ,KSXL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,284) = NPROLD
         ISELEP(2,284) = ISELEP(2,284) + NPROLD
      ENDIF
      IF(ISELEP(1,286) .EQ. 1 .OR. ISELEP(1,286) .EQ. 2) THEN
         NPROLD = ISELEA(2,286)
         IF(ISELEP(2,286) .GT. 0) THEN
            NPRRUS = ISELEP(2,286)
            CALL STMNLI  (ASXOZZ ,JJ1L,SSXOZZ ,JJ1L,KSXL,ILIMXP,
     $                    2  ,ILINTX, 1 , 1 , 1 ,KSXL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,286) = NPROLD
         ISELEP(2,286) = ISELEP(2,286) + NPROLD
      ENDIF
      IF(ISELEP(1,288) .EQ. 1 .OR. ISELEP(1,288) .EQ. 2) THEN
         NPROLD = ISELEA(2,288)
         IF(ISELEP(2,288) .GT. 0) THEN
            NPRRUS = ISELEP(2,288)
            CALL STMNLI  (ASYOXX ,JJ1L,SSYOXX ,JJ1L,KSYL,ILIMXP,
     $                    2  ,ILINTY, 1 , 1 , 1 ,KSYL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,288) = NPROLD
         ISELEP(2,288) = ISELEP(2,288) + NPROLD
      ENDIF
      IF(ISELEP(1,290) .EQ. 1 .OR. ISELEP(1,290) .EQ. 2) THEN
         NPROLD = ISELEA(2,290)
         IF(ISELEP(2,290) .GT. 0) THEN
            NPRRUS = ISELEP(2,290)
            CALL STMNLI  (ASYOYY ,JJ1L,SSYOYY ,JJ1L,KSYL,ILIMXP,
     $                    2  ,ILINTY, 1 , 1 , 1 ,KSYL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,290) = NPROLD
         ISELEP(2,290) = ISELEP(2,290) + NPROLD
      ENDIF
      IF(ISELEP(1,292) .EQ. 1 .OR. ISELEP(1,292) .EQ. 2) THEN
         NPROLD = ISELEA(2,292)
         IF(ISELEP(2,292) .GT. 0) THEN
            NPRRUS = ISELEP(2,292)
            CALL STMNLI  (ASYOZZ ,JJ1L,SSYOZZ ,JJ1L,KSYL,ILIMXP,
     $                    2  ,ILINTY, 1 , 1 , 1 ,KSYL-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,292) = NPROLD
         ISELEP(2,292) = ISELEP(2,292) + NPROLD
      ENDIF
C
C                                  HAEUFIGKEITSVERTEILUNG DER
C                                  INKLINATIONSWINKEL DER VORTICITY-
C                                  VEKTOREN
C                                  ---------------------------------
C
      IF(ISELEP(1,306) .EQ. 1 .OR. ISELEP(1,306) .EQ. 2) THEN
         NPROLD = ISELEA(2,306)
         IF(ISELEP(2,306) .GT. 0) THEN
            NPRRUS = ISELEP(2,306)
            CALL STMNLI  (AHOZOY ,JJ2L,SHOZOY ,JJ2L,KH0L,ILIMXP,
     $                    2  ,ILINT0, 1 , 2 , 1 ,KH0L-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,306) = NPROLD
         ISELEP(2,306) = ISELEP(2,306) + NPROLD
      ENDIF
      IF(ISELEP(1,308) .EQ. 1 .OR. ISELEP(1,308) .EQ. 2) THEN
         NPROLD = ISELEA(2,308)
         IF(ISELEP(2,308) .GT. 0) THEN
            NPRRUS = ISELEP(2,308)
            CALL STMNLI  (AHOZOX ,JJ2L,SHOZOX ,JJ2L,KH0L,ILIMXP,
     $                    2  ,ILINT0, 1 , 2 , 1 ,KH0L-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,308) = NPROLD
         ISELEP(2,308) = ISELEP(2,308) + NPROLD
      ENDIF
      IF(ISELEP(1,310) .EQ. 1 .OR. ISELEP(1,310) .EQ. 2) THEN
         NPROLD = ISELEA(2,310)
         IF(ISELEP(2,310) .GT. 0) THEN
            NPRRUS = ISELEP(2,310)
            CALL STMNLI  (AHOYOX ,JJ2L,SHOYOX ,JJ2L,KH0L,ILIMXP,
     $                    2  ,ILINT0, 1 , 2 , 1 ,KH0L-LINFB,
     $                    NPROLD,NPRRUS,DREAD,DCONT,NPRNEU)
         ENDIF
         ISELEA(2,310) = NPROLD
         ISELEP(2,310) = ISELEP(2,310) + NPROLD
      ENDIF
C
      RETURN
      END
