










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
      SUBROUTINE DOBCA1  (KANAL,MODUS,ISTALM,ISTPLM,JSTALM,JSTPLM,
     $                    ITTOT,ITFLUC,DT,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
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
C        D O B C A 1      AUSGABE VON FELDERN (WAHLWEISE FORMATIERT ODER
C                         UNFORMATIERT S. VAR. 'MODUS'),DIE UEBERWIEGEND
C                         ENSEMBLE-MITTELWERTE VON KORRELATIONSFUNKTIO-
C                         NEN, -KOEFFIZIENTEN, LEISTUNGSDICHTESPEKTREN
C                         UND HAEUFIGKEITSVERTEILUNGEN ENTHALTEN AUF
C                         KANAL 'KANAL'.
C        A C H T U N G    DIESES UNTERPROGRAMM MUSS ERWEITERT WERDEN,
C                         WENN DIE STATIST. AUSWERTUNG UMFANGREICHER
C                         WIRD !
C                         DIE NEU HINZUKOMMENDEN FELDER SIND GEMAESS
C                         DER ORDNUNG IM ISELEP-FELD EINZUREIHEN.
C
C                         BEI DER BERECHNUNG DER PHYS. ZEITEN, DIE
C                         SICH AUF DIE ENSEMBLE-MITTELUNG BEZIEHEN,
C                         WURDE DAVON AUSGEGANGEN, DASS SOWOHL DER
C                         ZEITSCHRITT 'DT' ALS AUCH DIE GROESSE
C                         'ITFLUC' WAEHREND ALLER FORTSETZUNGSLAEUFE
C                         KONSTANT BLEIBEN. DA DIE ZEITEN NUR IN-
C                         FORMATIVEN CHARAKTER HABEN, KOENNEN
C                         EVENTUELL AUFTRETENDE FEHLER TOLERIERT
C                         WERDEN.
C*STARLET***************************************************************
C
C PARAM: KANAL                    - DIE DATEN WERDEN UEBER DIESEN KANAL
C                                   AUSGEGEBEN
C        MODUS                    - CHARACTER-VARIABLE:
C                                   'BINAER  ' : DATEN WERDEN UNFORMA-
C                                                TIERT GESCHRIEBEN
C                                   'CODIERT ' : DATEN WERDEN FORMA-
C                                                TIERT GESCHRIEBEN
C        ISTALM                   - STARTINDEX FUER DIE LINIENMITTELUNG
C                                   IN X-RICHTUNG (FALLS X-RI. HOMOGEN)
C        ISTPLM                   - STOPINDEX FUER DIE LINIENMITTELUNG
C                                   IN X-RICHTUNG (FALLS X-RI. HOMOGEN)
C        JSTALM                   - STARTINDEX FUER DIE LINIENMITTELUNG
C                                   IN Y-RICHTUNG (FALLS Y-RI. HOMOGEN)
C        JSTPLM                   - STOPINDEX FUER DIE LINIENMITTELUNG
C                                   IN Y-RICHTUNG (FALLS Y-RI. HOMOGEN)
C        ITTOT                    - GESAMTANZAHL DER ZEITSCHRITTE (SUMME
C                                   MEHRERER LAEUFE)
C        ITFLUC                   - NACH JEWEILS ITFLUC ZEITSCHRITTEN
C                                   WIRD EINE STICHPROBE FUER DIE
C                                   BILDUNG DER ENSEMBLE-MITTELWERTE
C                                   ENTNOMMEN
C        DT                       - ZEITSCHRITT
C        X  (II), Y  (JJ), Z  (KK)- KOORDINATEN DER BASISZELLBEZUGS-
C                                   PUNKTE
C        DX (II), DY (JJ), DZ (KK)- ABSTAENDE DER BASISZELLBEZUGS-
C                                   PUNKTE
C        DDX(II), DDY(JJ), DDZ(KK)- ABMESSUNGEN DER BASISZELLEN
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
C UPROG                           : DOBC
C
C DEFINE-DIREKTIVEN               : XHOMOG, YHOMOG (IN COMMON-DECK
C                                   CDOBCA)
C
C        09.12.88 (HW)  : ORIGINAL AUS DOBCA ABGELEITET
C        09.12.88 (HW)  : CDOBCA EINGEFUEHRT
C        23.12.88 (HW)  : ILINT. STEUERT DIE AUSABE FUER DIE I-RICHTUNG
C
C*STARLET***************************************************************
C
      CHARACTER (LEN=8)   MODUS
      CHARACTER (LEN=80)  CT, CH1680
C
      REAL     X  (II),      Y  (JJ),      Z  (KK),
     $         DX (II),      DY (JJ),      DZ (KK),
     $         DDX(II),      DDY(JJ),      DDZ(KK)

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
C                                 START- UND STOP-KOORDINATEN DER
C                                 LINIEN- BZW. FLAECHENMITTELUNG
C
      XSLM   = 0.0
      XELM   = 0.0
C                                 FUER STAGGERED VARIABLE:
      XSLMST = 0.0
      XELMST = 0.0
      YSLM   = 0.0
      YELM   = 0.0
      YSLMST = 0.0
      YELMST = 0.0
C
C                                 IN Z-RICHTUNG IST KEINE LINIENMITTE-
C                                 LUNG VORGESEHEN, DAHER WIRD FORMAL EIN
C                                 GEFUEHRT:
C
      ZSLM   = 0.0
      ZELM   = 0.0
      ZSLMST = 0.0
      ZELMST = 0.0
C
C                                 PHYSIKALISCHE ZEIT ZU DER DIE LETZTE
C                                 STICHPROBE FUER DIE BILDUNG VON
C                                 ENSEMBLE-MITTELWERTEN GENOMMEN WURDE
C
      TEEM   = FLOAT(ITTOT - MOD(ITTOT,ITFLUC)) * DT
C
C                                 STICHPROBEN F. D. BILD. V. ENSEMBLE-M.
C                                 WERDEN JEWEILS NACH ABLAUF FOLGENDER
C                                 PHYSIKAL. ZEIT ENTNOMMEN:
C
      DTEM   = FLOAT(ITFLUC) * DT
C
C                                 AUSGABE DER FELDER, DIE ENSEMBLE-
C                                 MITTELWERTE ENTHALTEN.
C
C                                 STATIST. GROESSEN DER U-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 14) .EQ. 1  .OR.  ISELEP(1, 14) .EQ. 2) THEN
C
C                                 DAS FELD WURDE IM MOMENTANEN LAUF
C                                 AUSGEWERTET UND WIRD DAHER AUSGEGEBEN
C
C                                 PHYSIKAL. ZEIT ZU DER MIT DER BILDUNG
C                                 VON ENSEMBLE-MITTELWERTEN BEGONNEN
C                                 WURDE. (EXAKTER: ZU DIESEM ZEITPUNKT
C                                 WURDE DIE ERSTE STICHPROBE ENTNOMMEN)
C
         TSEM   = TEEM - FLOAT((ISELEP(2, 14)-1) * ITFLUC) * DT
C
         CT     = CH1680 (' R_U"U"(X)      ')
         CALL DOBC   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILINTX,KANAL,MODUS,
     $                      0      ,      0      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARXUU  ,CT)
      ENDIF
      IF(ISELEP(1, 16) .EQ. 1  .OR.  ISELEP(1, 16) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 16)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_U"U"(Y)      ')
         CALL DOBC   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILINTY,KANAL,MODUS,
     $                      0      ,      0      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARYUU  ,CT)
      ENDIF
      IF(ISELEP(1, 18) .EQ. 1  .OR.  ISELEP(1, 18) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 18)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_U"U"(Z)      ')
         CALL DOBC   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILINTZ,KANAL,MODUS,
     $                      0      ,      0      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARZUU  ,CT)
      ENDIF
      IF(ISELEP(1, 20) .EQ. 1  .OR.  ISELEP(1, 20) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 20)-1) * ITFLUC) * DT
         CT     = CH1680 (' S_U"U"(K_X)    ')
         CALL DOBC   (KSXL,JJ1L,ILIMXP,KSXL,JJ1L,ILINTX,KANAL,MODUS,
     $                      0      ,      0      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ASXUU  ,CT)
      ENDIF
      IF(ISELEP(1, 22) .EQ. 1  .OR.  ISELEP(1, 22) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 22)-1) * ITFLUC) * DT
         CT     = CH1680 (' S_U"U"(K_Y)    ')
         CALL DOBC   (KSYL,JJ1L,ILIMXP,KSYL,JJ1L,ILINTY,KANAL,MODUS,
     $                      0      ,      0      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ASYUU  ,CT)
      ENDIF
      IF(ISELEP(1, 24) .EQ. 1  .OR.  ISELEP(1, 24) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 24)-1) * ITFLUC) * DT
         CT     = CH1680 (' S_U"U"(K_Z)    ')
         CALL DOBC   (KSZL,JJ1L,ILIMXP,KSZL,JJ1L,ILINTZ,KANAL,MODUS,
     $                      0      ,      0      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ASZUU  ,CT)
      ENDIF
C
C                                 STATIST. GROESSEN DER V-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 39) .EQ. 1  .OR.  ISELEP(1, 39) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 39)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_V"V"(X)      ')
         CALL DOBC   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILINTX,KANAL,MODUS,
     $                      0      ,      1      ,      0      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,ARXVV  ,CT)
      ENDIF
      IF(ISELEP(1, 41) .EQ. 1  .OR.  ISELEP(1, 41) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 41)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_V"V"(Y)      ')
         CALL DOBC   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILINTY,KANAL,MODUS,
     $                      0      ,      1      ,      0      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,ARYVV  ,CT)
      ENDIF
      IF(ISELEP(1, 43) .EQ. 1  .OR.  ISELEP(1, 43) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 43)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_V"V"(Z)      ')
         CALL DOBC   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILINTZ,KANAL,MODUS,
     $                      0      ,      1      ,      0      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,ARZVV  ,CT)
      ENDIF
      IF(ISELEP(1, 45) .EQ. 1  .OR.  ISELEP(1, 45) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 45)-1) * ITFLUC) * DT
         CT     = CH1680 (' S_V"V"(K_X)    ')
         CALL DOBC   (KSXL,JJ1L,ILIMXP,KSXL,JJ1L,ILINTX,KANAL,MODUS,
     $                      0      ,      1      ,      0      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,ASXVV  ,CT)
      ENDIF
      IF(ISELEP(1, 47) .EQ. 1  .OR.  ISELEP(1, 47) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 47)-1) * ITFLUC) * DT
         CT     = CH1680 (' S_V"V"(K_Y)    ')
         CALL DOBC   (KSYL,JJ1L,ILIMXP,KSYL,JJ1L,ILINTY,KANAL,MODUS,
     $                      0      ,      1      ,      0      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,ASYVV  ,CT)
      ENDIF
      IF(ISELEP(1, 49) .EQ. 1  .OR.  ISELEP(1, 49) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 49)-1) * ITFLUC) * DT
         CT     = CH1680 (' S_V"V"(K_Z)    ')
         CALL DOBC   (KSZL,JJ1L,ILIMXP,KSZL,JJ1L,ILINTZ,KANAL,MODUS,
     $                      0      ,      1      ,      0      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,ASZVV  ,CT)
      ENDIF
C
C                                 STATIST. GROESSEN DER W-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 64) .EQ. 1  .OR.  ISELEP(1, 64) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 64)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_W"W"(X)      ')
         CALL DOBC   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILINTX,KANAL,MODUS,
     $                      1      ,      0      ,      0      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,ARXWW  ,CT)
      ENDIF
      IF(ISELEP(1, 66) .EQ. 1  .OR.  ISELEP(1, 66) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 66)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_W"W"(Y)      ')
         CALL DOBC   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILINTY,KANAL,MODUS,
     $                      1      ,      0      ,      0      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,ARYWW  ,CT)
      ENDIF
      IF(ISELEP(1, 68) .EQ. 1  .OR.  ISELEP(1, 68) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 68)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_W"W"(Z)      ')
         CALL DOBC   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILINTZ,KANAL,MODUS,
     $                      1      ,      0      ,      0      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,ARZWW  ,CT)
      ENDIF
      IF(ISELEP(1, 70) .EQ. 1  .OR.  ISELEP(1, 70) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 70)-1) * ITFLUC) * DT
         CT     = CH1680 (' S_W"W"(K_X)    ')
         CALL DOBC   (KSXL,JJ1L,ILIMXP,KSXL,JJ1L,ILINTX,KANAL,MODUS,
     $                      1      ,      0      ,      0      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,ASXWW  ,CT)
      ENDIF
      IF(ISELEP(1, 72) .EQ. 1  .OR.  ISELEP(1, 72) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 72)-1) * ITFLUC) * DT
         CT     = CH1680 (' S_W"W"(K_Y)    ')
         CALL DOBC   (KSYL,JJ1L,ILIMXP,KSYL,JJ1L,ILINTY,KANAL,MODUS,
     $                      1      ,      0      ,      0      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,ASYWW  ,CT)
      ENDIF
      IF(ISELEP(1, 74) .EQ. 1  .OR.  ISELEP(1, 74) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 74)-1) * ITFLUC) * DT
         CT     = CH1680 (' S_W"W"(K_Z)    ')
         CALL DOBC   (KSZL,JJ1L,ILIMXP,KSZL,JJ1L,ILINTZ,KANAL,MODUS,
     $                      1      ,      0      ,      0      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,ASZWW  ,CT)
      ENDIF
C
C                                  HIER: KREUZKORRELATIONSFUNKTIONEN
C                                  ---------------------------------
C
      IF(ISELEP(1,210) .EQ. 1  .OR.  ISELEP(1,210) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,210)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_U"V"(X)      ')
         CALL DOBC   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILINTX,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARXUV  ,CT)
      ENDIF
      IF(ISELEP(1,212) .EQ. 1  .OR.  ISELEP(1,212) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,212)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_U"V"(Y)      ')
         CALL DOBC   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILINTY,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARYUV  ,CT)
      ENDIF
      IF(ISELEP(1,216) .EQ. 1  .OR.  ISELEP(1,216) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,216)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_U"W"(X)      ')
         CALL DOBC   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILINTX,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARXUW  ,CT)
      ENDIF
      IF(ISELEP(1,218) .EQ. 1  .OR.  ISELEP(1,218) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,218)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_U"W"(Y)      ')
         CALL DOBC   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILINTY,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARYUW  ,CT)
      ENDIF
      IF(ISELEP(1,224) .EQ. 1  .OR.  ISELEP(1,224) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,224)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_V"W"(Y)      ')
         CALL DOBC   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILINTY,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,ARYVW  ,CT)
      ENDIF
C
C                                  HIER: KORRELATIONSKOEFFIZIENTEN
C                                  -------------------------------
C
      IF(ISELEP(1,234) .EQ. 1  .OR.  ISELEP(1,234) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,234)-1) * ITFLUC) * DT
         CT     = CH1680 (' RC_U"W"(X)     ')
         CALL DOBC   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILINTX,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ACXUW  ,CT)
      ENDIF
      IF(ISELEP(1,238) .EQ. 1  .OR.  ISELEP(1,238) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,238)-1) * ITFLUC) * DT
         CT     = CH1680 (' RC_U"W"(Z)     ')
         CALL DOBC   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILINTZ,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ACZUW  ,CT)
      ENDIF
C
C                                  OBERE DREIECKSMATRIX DES TENSORS
C                                  DER KORRELATIONSFUNKTIONEN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  --------------------------------
C
      IF(ISELEP(1,246) .EQ. 1  .OR.  ISELEP(1,246) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,246)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_OX"OX"(X)    ')
         CALL DOBC   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILINTX,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,ARXOXX ,CT)
      ENDIF
      IF(ISELEP(1,248) .EQ. 1  .OR.  ISELEP(1,248) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,248)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_OX"OY"(X)    ')
         CALL DOBC   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILINTX,KANAL,MODUS,
     $                      1      ,      1      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARXOXY ,CT)
      ENDIF
      IF(ISELEP(1,250) .EQ. 1  .OR.  ISELEP(1,250) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,250)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_OX"OZ"(X)    ')
         CALL DOBC   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILINTX,KANAL,MODUS,
     $                      1      ,      1      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARXOXZ ,CT)
      ENDIF
      IF(ISELEP(1,252) .EQ. 1  .OR.  ISELEP(1,252) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,252)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_OY"OY"(X)    ')
         CALL DOBC   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILINTX,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARXOYY ,CT)
      ENDIF
      IF(ISELEP(1,254) .EQ. 1  .OR.  ISELEP(1,254) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,254)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_OY"OZ"(X)    ')
         CALL DOBC   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILINTX,KANAL,MODUS,
     $                      1      ,      1      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARXOYZ ,CT)
      ENDIF
      IF(ISELEP(1,256) .EQ. 1  .OR.  ISELEP(1,256) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,256)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_OZ"OZ"(X)    ')
         CALL DOBC   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILINTX,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARXOZZ ,CT)
      ENDIF
      IF(ISELEP(1,258) .EQ. 1  .OR.  ISELEP(1,258) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,258)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_OX"OX"(Y)    ')
         CALL DOBC   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILINTY,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,ARYOXX ,CT)
      ENDIF
      IF(ISELEP(1,260) .EQ. 1  .OR.  ISELEP(1,260) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,260)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_OX"OY"(Y)    ')
         CALL DOBC   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILINTY,KANAL,MODUS,
     $                      1      ,      1      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARYOXY ,CT)
      ENDIF
      IF(ISELEP(1,262) .EQ. 1  .OR.  ISELEP(1,262) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,262)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_OX"OZ"(Y)    ')
         CALL DOBC   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILINTY,KANAL,MODUS,
     $                      1      ,      1      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARYOXZ ,CT)
      ENDIF
      IF(ISELEP(1,264) .EQ. 1  .OR.  ISELEP(1,264) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,264)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_OY"OY"(Y)    ')
         CALL DOBC   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILINTY,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARYOYY ,CT)
      ENDIF
      IF(ISELEP(1,266) .EQ. 1  .OR.  ISELEP(1,266) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,266)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_OY"OZ"(Y)    ')
         CALL DOBC   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILINTY,KANAL,MODUS,
     $                      1      ,      1      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARYOYZ ,CT)
      ENDIF
      IF(ISELEP(1,268) .EQ. 1  .OR.  ISELEP(1,268) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,268)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_OZ"OZ"(Y)    ')
         CALL DOBC   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILINTY,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARYOZZ ,CT)
      ENDIF
      IF(ISELEP(1,270) .EQ. 1  .OR.  ISELEP(1,270) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,270)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_OX"OX"(Z)    ')
         CALL DOBC   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILINTZ,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,ARZOXX ,CT)
      ENDIF
      IF(ISELEP(1,272) .EQ. 1  .OR.  ISELEP(1,272) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,272)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_OX"OY"(Z)    ')
         CALL DOBC   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILINTZ,KANAL,MODUS,
     $                      1      ,      1      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARZOXY ,CT)
      ENDIF
      IF(ISELEP(1,274) .EQ. 1  .OR.  ISELEP(1,274) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,274)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_OX"OZ"(Z)    ')
         CALL DOBC   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILINTZ,KANAL,MODUS,
     $                      1      ,      1      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARZOXZ ,CT)
      ENDIF
      IF(ISELEP(1,276) .EQ. 1  .OR.  ISELEP(1,276) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,276)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_OY"OY"(Z)    ')
         CALL DOBC   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILINTZ,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARZOYY ,CT)
      ENDIF
      IF(ISELEP(1,278) .EQ. 1  .OR.  ISELEP(1,278) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,278)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_OY"OZ"(Z)    ')
         CALL DOBC   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILINTZ,KANAL,MODUS,
     $                      1      ,      1      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARZOYZ ,CT)
      ENDIF
      IF(ISELEP(1,280) .EQ. 1  .OR.  ISELEP(1,280) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,280)-1) * ITFLUC) * DT
         CT     = CH1680 (' R_OZ"OZ"(Z)    ')
         CALL DOBC   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILINTZ,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ARZOZZ ,CT)
      ENDIF
C
C                                  LEISTUNGSDICHTESPEKTREN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  ----------------------------
C
      IF(ISELEP(1,282) .EQ. 1  .OR.  ISELEP(1,282) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,282)-1) * ITFLUC) * DT
         CT     = CH1680 (' S_OX"OX"(K_X)  ')
         CALL DOBC   (KSXL,JJ1L,ILIMXP,KSXL,JJ1L,ILINTX,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,ASXOXX ,CT)
      ENDIF
      IF(ISELEP(1,284) .EQ. 1  .OR.  ISELEP(1,284) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,284)-1) * ITFLUC) * DT
         CT     = CH1680 (' S_OY"OY"(K_X)  ')
         CALL DOBC   (KSXL,JJ1L,ILIMXP,KSXL,JJ1L,ILINTX,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ASXOYY ,CT)
      ENDIF
      IF(ISELEP(1,286) .EQ. 1  .OR.  ISELEP(1,286) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,286)-1) * ITFLUC) * DT
         CT     = CH1680 (' S_OZ"OZ"(K_X)  ')
         CALL DOBC   (KSXL,JJ1L,ILIMXP,KSXL,JJ1L,ILINTX,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ASXOZZ ,CT)
      ENDIF
      IF(ISELEP(1,288) .EQ. 1  .OR.  ISELEP(1,288) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,288)-1) * ITFLUC) * DT
         CT     = CH1680 (' S_OX"OX"(K_Y)  ')
         CALL DOBC   (KSYL,JJ1L,ILIMXP,KSYL,JJ1L,ILINTY,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,ASYOXX ,CT)
      ENDIF
      IF(ISELEP(1,290) .EQ. 1  .OR.  ISELEP(1,290) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,290)-1) * ITFLUC) * DT
         CT     = CH1680 (' S_OY"OY"(K_Y)  ')
         CALL DOBC   (KSYL,JJ1L,ILIMXP,KSYL,JJ1L,ILINTY,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ASYOYY ,CT)
      ENDIF
      IF(ISELEP(1,292) .EQ. 1  .OR.  ISELEP(1,292) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,292)-1) * ITFLUC) * DT
         CT     = CH1680 (' S_OZ"OZ"(K_Y)  ')
         CALL DOBC   (KSYL,JJ1L,ILIMXP,KSYL,JJ1L,ILINTY,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,ASYOZZ ,CT)
      ENDIF
C
C                                  HAEUFIGKEITSVERTEILUNG DER
C                                  INKLINATIONSWINKEL DER VORTICITY-
C                                  VEKTOREN
C                                  ---------------------------------
C
      IF(ISELEP(1,306) .EQ. 1  .OR.  ISELEP(1,306) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,306)-1) * ITFLUC) * DT
         CT     = CH1680 (' H (ATAN(OZ/OY))')
         CALL DOBC   (KH0L,JJ2L,ILIMXP,KH0L,JJ2L,ILINT0,KANAL,MODUS,
     $                      1      ,      1      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AHOZOY ,CT)
      ENDIF
      IF(ISELEP(1,308) .EQ. 1  .OR.  ISELEP(1,308) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,308)-1) * ITFLUC) * DT
         CT     = CH1680 (' H (ATAN(OZ/OX))')
         CALL DOBC   (KH0L,JJ2L,ILIMXP,KH0L,JJ2L,ILINT0,KANAL,MODUS,
     $                      1      ,      1      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AHOZOX ,CT)
      ENDIF
      IF(ISELEP(1,310) .EQ. 1  .OR.  ISELEP(1,310) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,310)-1) * ITFLUC) * DT
         CT     = CH1680 (' H (ATAN(OY/OX))')
         CALL DOBC   (KH0L,JJ2L,ILIMXP,KH0L,JJ2L,ILINT0,KANAL,MODUS,
     $                      1      ,      1      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AHOYOX ,CT)
      ENDIF
C
      RETURN
      END
