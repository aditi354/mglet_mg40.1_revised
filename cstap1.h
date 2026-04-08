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
