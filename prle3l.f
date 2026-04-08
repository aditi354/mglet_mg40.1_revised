










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
      SUBROUTINE PRLE3L  (ILINT0,ILINTX,ILINTY,ILINTZ,K1,KSLIN,
     $                    KANAL,NRRUN,LINPRN,

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
C        P R L E 3 L      KOMPAKTAUSGABE DER "LINIEN"-FELDER:
C                         - KORRELATIONSFUNKTIONEN
C                         - KORRELATIONSKOEFFIZIENTEN
C                         - LEISTUNGSDICHTESPEKTREN
C                         - HAEUFIGKEITSVERTEILUNGEN
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
C        K1             - STARTINDEX FUER DIE AUSZUGEBENDE "LINIE"
C        KSLIN          - INKREMENT FUER DIE AUSZUGEBENDE "LINIE"
C        NRRUN          - LAUFNUMMER
C        LINPRN         - LOGISCHE VARIABLE:
C                         .FALSE. : ES WERDEN NUR DIE INFORMATIONEN WIE
C                                   MAXIMA, MINIMA, INTEGRALE LAENGEN-
C                                   MASSE ETC. AUSGEGEBEN. DAS EIGENT-
C                                   LICHE FELD WIRD JEDOCH NICHT AUS-
C                                   GEGEBEN
C                         .TRUE.  : INFORMATIONEN + FELD WERDEN AUSGE-
C                                   GEBEN
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
C UPROG                 : PR1LIN
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        22.12.88 (HW)  : ORIGINAL
C        25.02.92 (MM)  : ACHTUNG VARIABLE "KANAL" IN UEBERGABE
C                         EIGEFUEGT
C
C*STARLET***************************************************************
C

      CHARACTER (LEN=32)  CKEINE, CEREF, CGREF, CPREF, CTAURE, CUREF,
     $               COMREF, CO2REF, CHEREF, CLREF, CRLREF, CASDUI,
     $               CASDOM


      LOGICAL   LINPRN
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

      CKEINE = '--------------------------------'
      CLREF  = 'L-REF                           '
      CRLREF = '1.0 / L-REF                     '
      CEREF  = '(U-REF)**2                      '
      CGREF  = 'RHO * L-REF * U-REF             '
      CPREF  = '1/2 * RHO * (U-REF)**2          '
      CTAURE = 'RHO * (U-REF)**2                '
      CUREF  = 'U-REF                           '
      COMREF = '(U-REF) / (L-REF)               '
      CO2REF = '(U-REF)**2 / (L-REF)**2         '
      CHEREF = '(U-REF)**2 / (L-REF)            '
      CASDUI = '1/2 * (U-REF)**2 * (L-REF)      '
      CASDOM = '1/2 * (U-REF)**2 / (L-REF)      '
C
C                                  HIER: AUTOKORRELATIONSFUNKTIONEN
C                                  --------------------------------
C
      WRITE       (KANAL,6110)
C
      IF(ISELEP(1, 14) .EQ. 1 .OR. ISELEP(1, 14) .EQ. 2) THEN
         CALL PR1LIN  (ARXUU  , 2 ,' R_U"U"(X)      ',KKXL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTX,
     $                 SRXUU  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1, 16) .EQ. 1 .OR. ISELEP(1, 16) .EQ. 2) THEN
         CALL PR1LIN  (ARYUU  , 2 ,' R_U"U"(Y)      ',KKYL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTY,
     $                 SRYUU  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1, 18) .EQ. 1 .OR. ISELEP(1, 18) .EQ. 2) THEN
         CALL PR1LIN  (ARZUU  , 2 ,' R_U"U"(Z)      ',KKZL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTZ,
     $                 SRZUU  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
C
      WRITE       (KANAL,6130)
C
      IF(ISELEP(1, 39) .EQ. 1 .OR. ISELEP(1, 39) .EQ. 2) THEN
         CALL PR1LIN  (ARXVV  , 2 ,' R_V"V"(X)      ',KKXL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTX,
     $                 SRXVV  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1, 41) .EQ. 1 .OR. ISELEP(1, 41) .EQ. 2) THEN
         CALL PR1LIN  (ARYVV  , 2 ,' R_V"V"(Y)      ',KKYL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTY,
     $                 SRYVV  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1, 43) .EQ. 1 .OR. ISELEP(1, 43) .EQ. 2) THEN
         CALL PR1LIN  (ARZVV  , 2 ,' R_V"V"(Z)      ',KKZL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTZ,
     $                 SRZVV  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
C
      WRITE       (KANAL,6150)
C
      IF(ISELEP(1, 64) .EQ. 1 .OR. ISELEP(1, 64) .EQ. 2) THEN
         CALL PR1LIN  (ARXWW  , 2 ,' R_W"W"(X)      ',KKXL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTX,
     $                 SRXWW  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1, 66) .EQ. 1 .OR. ISELEP(1, 66) .EQ. 2) THEN
         CALL PR1LIN  (ARYWW  , 2 ,' R_W"W"(Y)      ',KKYL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTY,
     $                 SRYWW  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1, 68) .EQ. 1 .OR. ISELEP(1, 68) .EQ. 2) THEN
         CALL PR1LIN  (ARZWW  , 2 ,' R_W"W"(Z)      ',KKZL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTZ,
     $                 SRZWW  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
C
C                                  HIER: KREUZKORRELATIONSFUNKTIONEN
C                                  ---------------------------------
C
      WRITE       (KANAL,6170)
C
      IF(ISELEP(1,210) .EQ. 1 .OR. ISELEP(1,210) .EQ. 2) THEN
         CALL PR1LIN  (ARXUV  , 2 ,' R_U"V"(X)      ',KKXL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTX,
     $                 SRXUV  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,212) .EQ. 1 .OR. ISELEP(1,212) .EQ. 2) THEN
         CALL PR1LIN  (ARYUV  , 2 ,' R_U"V"(Y)      ',KKYL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTY,
     $                 SRYUV  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,216) .EQ. 1 .OR. ISELEP(1,216) .EQ. 2) THEN
         CALL PR1LIN  (ARXUW  , 2 ,' R_U"W"(X)      ',KKXL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTX,
     $                 SRXUW  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,218) .EQ. 1 .OR. ISELEP(1,218) .EQ. 2) THEN
         CALL PR1LIN  (ARYUW  , 2 ,' R_U"W"(Y)      ',KKYL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTY,
     $                 SRYUW  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,224) .EQ. 1 .OR. ISELEP(1,224) .EQ. 2) THEN
         CALL PR1LIN  (ARYVW  , 2 ,' R_V"W"(Y)      ',KKYL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTY,
     $                 SRYVW  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
C
C                                  OBERE DREIECKSMATRIX DES TENSORS
C                                  DER KORRELATIONSFUNKTIONEN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  --------------------------------
C
      WRITE       (KANAL,6190)
C
      IF(ISELEP(1,246) .EQ. 1 .OR. ISELEP(1,246) .EQ. 2) THEN
         CALL PR1LIN  (ARXOXX , 2 ,' R_OX"OX"(X)    ',KKXL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTX,
     $                 SRXOXX ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,248) .EQ. 1 .OR. ISELEP(1,248) .EQ. 2) THEN
         CALL PR1LIN  (ARXOXY , 2 ,' R_OX"OY"(X)    ',KKXL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTX,
     $                 SRXOXY ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,250) .EQ. 1 .OR. ISELEP(1,250) .EQ. 2) THEN
         CALL PR1LIN  (ARXOXZ , 2 ,' R_OX"OZ"(X)    ',KKXL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTX,
     $                 SRXOXZ ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,252) .EQ. 1 .OR. ISELEP(1,252) .EQ. 2) THEN
         CALL PR1LIN  (ARXOYY , 2 ,' R_OY"OY"(X)    ',KKXL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTX,
     $                 SRXOYY ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,254) .EQ. 1 .OR. ISELEP(1,254) .EQ. 2) THEN
         CALL PR1LIN  (ARXOYZ , 2 ,' R_OY"OZ"(X)    ',KKXL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTX,
     $                 SRXOYZ ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,256) .EQ. 1 .OR. ISELEP(1,256) .EQ. 2) THEN
         CALL PR1LIN  (ARXOZZ , 2 ,' R_OZ"OZ"(X)    ',KKXL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTX,
     $                 SRXOZZ ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,258) .EQ. 1 .OR. ISELEP(1,258) .EQ. 2) THEN
         CALL PR1LIN  (ARYOXX , 2 ,' R_OX"OX"(Y)    ',KKYL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTY,
     $                 SRYOXX ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,260) .EQ. 1 .OR. ISELEP(1,260) .EQ. 2) THEN
         CALL PR1LIN  (ARYOXY , 2 ,' R_OX"OY"(Y)    ',KKYL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTY,
     $                 SRYOXY ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,262) .EQ. 1 .OR. ISELEP(1,262) .EQ. 2) THEN
         CALL PR1LIN  (ARYOXZ , 2 ,' R_OX"OZ"(Y)    ',KKYL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTY,
     $                 SRYOXZ ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,264) .EQ. 1 .OR. ISELEP(1,264) .EQ. 2) THEN
         CALL PR1LIN  (ARYOYY , 2 ,' R_OY"OY"(Y)    ',KKYL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTY,
     $                 SRYOYY ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,266) .EQ. 1 .OR. ISELEP(1,266) .EQ. 2) THEN
         CALL PR1LIN  (ARYOYZ , 2 ,' R_OY"OZ"(Y)    ',KKYL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTY,
     $                 SRYOYZ ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,268) .EQ. 1 .OR. ISELEP(1,268) .EQ. 2) THEN
         CALL PR1LIN  (ARYOZZ , 2 ,' R_OZ"OZ"(Y)    ',KKYL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTY,
     $                 SRYOZZ ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,270) .EQ. 1 .OR. ISELEP(1,270) .EQ. 2) THEN
         CALL PR1LIN  (ARZOXX , 2 ,' R_OX"OX"(Z)    ',KKZL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTZ,
     $                 SRZOXX ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,272) .EQ. 1 .OR. ISELEP(1,272) .EQ. 2) THEN
         CALL PR1LIN  (ARZOXY , 2 ,' R_OX"OY"(Z)    ',KKZL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTZ,
     $                 SRZOXY ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,274) .EQ. 1 .OR. ISELEP(1,274) .EQ. 2) THEN
         CALL PR1LIN  (ARZOXZ , 2 ,' R_OX"OZ"(Z)    ',KKZL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTZ,
     $                 SRZOXZ ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,276) .EQ. 1 .OR. ISELEP(1,276) .EQ. 2) THEN
         CALL PR1LIN  (ARZOYY , 2 ,' R_OY"OY"(Z)    ',KKZL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTZ,
     $                 SRZOYY ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,278) .EQ. 1 .OR. ISELEP(1,278) .EQ. 2) THEN
         CALL PR1LIN  (ARZOYZ , 2 ,' R_OY"OZ"(Z)    ',KKZL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTZ,
     $                 SRZOYZ ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,280) .EQ. 1 .OR. ISELEP(1,280) .EQ. 2) THEN
         CALL PR1LIN  (ARZOZZ , 2 ,' R_OZ"OZ"(Z)    ',KKZL,JJ2L,ILIMXP,
     $                 'KORRELATIONSFKT.','KORR.-RADIUS    ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTZ,
     $                 SRZOZZ ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
C
C                                  HIER: KORRELATIONSKOEFFIZIENTEN
C                                  -------------------------------
C
      WRITE       (KANAL,6210)
C
      IF(ISELEP(1,234) .EQ. 1 .OR. ISELEP(1,234) .EQ. 2) THEN
         CALL PR1LIN  (ACXUW  , 2 ,' RC_U"W"(X)     ',KKXL,JJ2L,ILIMXP,
     $                 'KORR.-KOEFF.    ','ZUGEORD. KOORD. ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTX,
     $                 SCXUW  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,238) .EQ. 1 .OR. ISELEP(1,238) .EQ. 2) THEN
         CALL PR1LIN  (ACZUW  , 2 ,' RC_U"W"(Z)     ',KKZL,JJ2L,ILIMXP,
     $                 'KORR.-KOEFF.    ','ZUGEORD. KOORD. ',
     $                 CKEINE ,CLREF  ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTZ,
     $                 SCZUW  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
C
C                                  LEISTUNGSDICHTESPEKTREN
C                                  ------------------------
C
      WRITE       (KANAL,6230)
C
      IF(ISELEP(1, 20) .EQ. 1 .OR. ISELEP(1, 20) .EQ. 2) THEN
         CALL PR1LIN  (ASXUU  , 1 ,' S_U"U"(K_X)    ',KSXL,JJ1L,ILIMXP,
     $                 'LEISTUNGSDICHTE ','WELLENZAHL      ',
     $                 CASDUI ,CRLREF ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTX,
     $                 SSXUU  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1, 22) .EQ. 1 .OR. ISELEP(1, 22) .EQ. 2) THEN
         CALL PR1LIN  (ASYUU  , 1 ,' S_U"U"(K_Y)    ',KSYL,JJ1L,ILIMXP,
     $                 'LEISTUNGSDICHTE ','WELLENZAHL      ',
     $                 CASDUI ,CRLREF ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTY,
     $                 SSYUU  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1, 24) .EQ. 1 .OR. ISELEP(1, 24) .EQ. 2) THEN
         CALL PR1LIN  (ASZUU  , 1 ,' S_U"U"(K_Z)    ',KSZL,JJ1L,ILIMXP,
     $                 'LEISTUNGSDICHTE ','WELLENZAHL      ',
     $                 CASDUI ,CRLREF ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTZ,
     $                 SSZUU  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
C
      WRITE       (KANAL,6250)
C
      IF(ISELEP(1, 45) .EQ. 1 .OR. ISELEP(1, 45) .EQ. 2) THEN
         CALL PR1LIN  (ASXVV  , 1 ,' S_V"V"(K_X)    ',KSXL,JJ1L,ILIMXP,
     $                 'LEISTUNGSDICHTE ','WELLENZAHL      ',
     $                 CASDUI ,CRLREF ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTX,
     $                 SSXVV  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1, 47) .EQ. 1 .OR. ISELEP(1, 47) .EQ. 2) THEN
         CALL PR1LIN  (ASYVV  , 1 ,' S_V"V"(K_Y)    ',KSYL,JJ1L,ILIMXP,
     $                 'LEISTUNGSDICHTE ','WELLENZAHL      ',
     $                 CASDUI ,CRLREF ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTY,
     $                 SSYVV  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1, 49) .EQ. 1 .OR. ISELEP(1, 49) .EQ. 2) THEN
         CALL PR1LIN  (ASZVV  , 1 ,' S_V"V"(K_Z)    ',KSZL,JJ1L,ILIMXP,
     $                 'LEISTUNGSDICHTE ','WELLENZAHL      ',
     $                 CASDUI ,CRLREF ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTZ,
     $                 SSZVV  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
C
      WRITE       (KANAL,6270)
C
      IF(ISELEP(1, 70) .EQ. 1 .OR. ISELEP(1, 70) .EQ. 2) THEN
         CALL PR1LIN  (ASXWW  , 1 ,' S_W"W"(K_X)    ',KSXL,JJ1L,ILIMXP,
     $                 'LEISTUNGSDICHTE ','WELLENZAHL      ',
     $                 CASDUI ,CRLREF ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTX,
     $                 SSXWW  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1, 72) .EQ. 1 .OR. ISELEP(1, 72) .EQ. 2) THEN
         CALL PR1LIN  (ASYWW  , 1 ,' S_W"W"(K_Y)    ',KSYL,JJ1L,ILIMXP,
     $                 'LEISTUNGSDICHTE ','WELLENZAHL      ',
     $                 CASDUI ,CRLREF ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTY,
     $                 SSYWW  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1, 74) .EQ. 1 .OR. ISELEP(1, 74) .EQ. 2) THEN
         CALL PR1LIN  (ASZWW  , 1 ,' S_W"W"(K_Z)    ',KSZL,JJ1L,ILIMXP,
     $                 'LEISTUNGSDICHTE ','WELLENZAHL      ',
     $                 CASDUI ,CRLREF ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTZ,
     $                 SSZWW  ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
C
      WRITE       (KANAL,6290)
C
      IF(ISELEP(1,282) .EQ. 1 .OR. ISELEP(1,282) .EQ. 2) THEN
         CALL PR1LIN  (ASXOXX , 1 ,' S_OX"OX"(K_X)  ',KSXL,JJ1L,ILIMXP,
     $                 'LEISTUNGSDICHTE ','WELLENZAHL      ',
     $                 CASDOM ,CRLREF ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTX,
     $                 SSXOXX ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,284) .EQ. 1 .OR. ISELEP(1,284) .EQ. 2) THEN
         CALL PR1LIN  (ASXOYY , 1 ,' S_OY"OY"(K_X)  ',KSXL,JJ1L,ILIMXP,
     $                 'LEISTUNGSDICHTE ','WELLENZAHL      ',
     $                 CASDOM ,CRLREF ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTX,
     $                 SSXOYY ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,286) .EQ. 1 .OR. ISELEP(1,286) .EQ. 2) THEN
         CALL PR1LIN  (ASXOZZ , 1 ,' S_OZ"OZ"(K_X)  ',KSXL,JJ1L,ILIMXP,
     $                 'LEISTUNGSDICHTE ','WELLENZAHL      ',
     $                 CASDOM ,CRLREF ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTX,
     $                 SSXOZZ ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,288) .EQ. 1 .OR. ISELEP(1,288) .EQ. 2) THEN
         CALL PR1LIN  (ASYOXX , 1 ,' S_OX"OX"(K_Y)  ',KSYL,JJ1L,ILIMXP,
     $                 'LEISTUNGSDICHTE ','WELLENZAHL      ',
     $                 CASDOM ,CRLREF ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTY,
     $                 SSYOXX ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,290) .EQ. 1 .OR. ISELEP(1,290) .EQ. 2) THEN
         CALL PR1LIN  (ASYOYY , 1 ,' S_OY"OY"(K_Y)  ',KSYL,JJ1L,ILIMXP,
     $                 'LEISTUNGSDICHTE ','WELLENZAHL      ',
     $                 CASDOM ,CRLREF ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTY,
     $                 SSYOYY ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,292) .EQ. 1 .OR. ISELEP(1,292) .EQ. 2) THEN
         CALL PR1LIN  (ASYOZZ , 1 ,' S_OZ"OZ"(K_Y)  ',KSYL,JJ1L,ILIMXP,
     $                 'LEISTUNGSDICHTE ','WELLENZAHL      ',
     $                 CASDOM ,CRLREF ,
     $                 K1,KSLIN,K2,  1  ,  1  ,  2  ,ILINTY,
     $                 SSYOZZ ,JJ1L,NRRUN,LINPRN,KANAL)
      ENDIF
C
C                                  HAEUFIGKEITSVERTEILUNG DER
C                                  INKLINATIONSWINKEL DER VORTICITY-
C                                  VEKTOREN
C                                  ---------------------------------
C
      WRITE       (KANAL,6310)
C
      IF(ISELEP(1,306) .EQ. 1 .OR. ISELEP(1,306) .EQ. 2) THEN
         CALL PR1LIN  (AHOZOY , 1 ,' H (ATAN(OZ/OY))',KH0L,JJ2L,ILIMXP,
     $                 'REL. HAEUFIGKEIT','ZUGEORD. WINKEL ',
     $                 CKEINE ,CKEINE ,
     $                 K1,KSLIN,K2,  1  ,  2  ,  2  ,ILINT0,
     $                 SHOZOY ,JJ2L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,308) .EQ. 1 .OR. ISELEP(1,308) .EQ. 2) THEN
         CALL PR1LIN  (AHOZOX , 1 ,' H (ATAN(OZ/OX))',KH0L,JJ2L,ILIMXP,
     $                 'REL. HAEUFIGKEIT','ZUGEORD. WINKEL ',
     $                 CKEINE ,CKEINE ,
     $                 K1,KSLIN,K2,  1  ,  2  ,  2  ,ILINT0,
     $                 SHOZOX ,JJ2L,NRRUN,LINPRN,KANAL)
      ENDIF
      IF(ISELEP(1,310) .EQ. 1 .OR. ISELEP(1,310) .EQ. 2) THEN
         CALL PR1LIN  (AHOYOX , 1 ,' H (ATAN(OY/OX))',KH0L,JJ2L,ILIMXP,
     $                 'REL. HAEUFIGKEIT','ZUGEORD. WINKEL ',
     $                 CKEINE ,CKEINE ,
     $                 K1,KSLIN,K2,  1  ,  2  ,  2  ,ILINT0,
     $                 SHOYOX ,JJ2L,NRRUN,LINPRN,KANAL)
      ENDIF
C
      RETURN
C
 6110 FORMAT (//,40X,'AUTOKORRELATIONSFUNKTION  (U-KOMP.)',
     $        /,40X,35(1H-),/)
 6130 FORMAT (//,40X,'AUTOKORRELATIONSFUNKTION  (V-KOMP.)',
     $        /,40X,35(1H-),/)
 6150 FORMAT (//,40X,'AUTOKORRELATIONSFUNKTION  (W-KOMP.)',
     $        /,40X,35(1H-),/)
 6170 FORMAT (//,40X,'KREUZKORRELATIONSFUNKTIONEN',
     $        /,40X,27(1H-),/)
 6190 FORMAT (//,40X,'OBERE DREIECKSMATRIX DES TENSORS DER',
     $        ' KORRELATIONSFUNKTIONEN',/,40X,'DER FLUKTUATIONEN',
     $        ' DER VORTICITY',/,40X,59(1H-),/)
 6210 FORMAT (//,40X,'KORRELATIONSKOEFFIZIENTEN',
     $        /,40X,25(1H-),/)
 6230 FORMAT (//,40X,'LEISTUNGSDICHTESPEKTREN  (U-KOMP.)',
     $        /,40X,34(1H-),/)
 6250 FORMAT (//,40X,'LEISTUNGSDICHTESPEKTREN  (V-KOMP.)',
     $        /,40X,34(1H-),/)
 6270 FORMAT (//,40X,'LEISTUNGSDICHTESPEKTREN  (W-KOMP.)',
     $        /,40X,34(1H-),/)
 6290 FORMAT (//,40X,'LEISTUNGSDICHTESPEKTREN  (OMEGA-KOMP.)',
     $        /,40X,38(1H-),/)
 6310 FORMAT (//,40X,'HAEUFIGKEITSVERTEILUNGEN DER INKLINATIONS',
     $        'WINKEL',/,40X,'DER VORTICITY-VEKTOREN',/,40X,47(1H-),/)
      END
