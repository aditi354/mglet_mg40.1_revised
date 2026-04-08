










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
      SUBROUTINE SETST1  (HILFL,KKML,RKOMXP,RKOMXA,

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
C        S E T S T 1      VORBELEGUNG DER FELDER, DIE IM WESENTLICHEN
C                         ENSEMBLE-MITTELWERTE VON KORRELATIONSFUNKTIO-
C                         NEN, -KOEFFIZIENTEN, LEISTUNGSDICHTESPEKTREN
C                         UND HAEUFIGKEITSVERTEILUNGEN ENTHALTEN.
C                                                      ALLE FELDER IN
C                         DENEN SUMMEN GEBILDET WERDEN, MUESSEN MIT 0.0
C                         VORBELEGT WERDEN. BEI EINEM VOELLIGEN NEUBE-
C                         GINN DER STATIST. AUSWERTUNG EINER GROESSE
C                         WIRD EINE "BESTMOEGLICHE" VORBELEGUNG DES-
C                         JENIGEN FELDES DURCHGEFUEHRT, DAS DEN ENSEM-
C                         BLE-MITTELWERT ENTHAELT.
C  A C H T U N G:         DIESE ROUTINE MUSS ERWEITERT WERDEN, FALLS
C                         DIE STATIST. AUSWERTUNG ERWEITERT WIRD !
C*STARLET***************************************************************
C
C PARAM: HILFL (KKML,JJ2L,ILIMXP) + HILFSFELD ZUM UMSORTIEREN
C        KKML                     - ARRAYDIMENSION (MAXIMALE ANZ. IN
C                                   K-RICHTUNG EINES "LINIEN-FELDES")
C        RKOMXP                   - MAXIMAL ZULAESSIGER KORRELATIONS-
C                                   RAD. WAEHREND DES MOMENTANEN LAUFES
C        RKOMXA                   - MAXIMAL ZULAESSIGER KORRELATIONS-
C                                   RADIUS WAEHREND DES VORANGEG. LAUFES
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
C UPROG                 : DPHI0, LINCTL
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        13.12.88 (HW)  : ORIGINAL AUS SETSTA ABGELEITET
C
C*STARLET***************************************************************
C
      CHARACTER (LEN=80)  CT,  CH1680
C
      REAL     HILFL (KKML,JJ2L,ILIMXP)
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
C                                 BELEGUNG DES HILFL-FELDES MIT 0.0
C
      CALL DPHI0   (KKML,JJ2L,ILIMXP,KKML,JJ2L,ILIMXP,HILFL  )
C
C                                 STATIST. GROESSEN DER U-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 14) .GE. 1) THEN
         IF(ISELEA(1, 14) .EQ. 0  .OR.  ISELEA(1, 14) .EQ. 3) THEN
            CALL DPHI0   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILIMXP,ARXUU  )
         ELSE
C
C                                  UMSORTIEREN UND NEUBELEGEN DES
C                                  ARXUU-FELDES
C
            CT = CH1680  (' R_U"U"(X)      ')
            CALL LINCTL  (ARXUU  ,  14 , 'X', HILFL, KKXL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1, 15) .GE. 1) THEN
            CALL DPHI0   (KKXL,JJ1L,ILIMXP,KKXL,JJ1L,ILIMXP,SRXUU  )
      ENDIF
      IF(ISELEP(1, 16) .GE. 1) THEN
         IF(ISELEA(1, 16) .EQ. 0  .OR.  ISELEA(1, 16) .EQ. 3) THEN
            CALL DPHI0   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILIMXP,ARYUU  )
         ELSE
            CT = CH1680  (' R_U"U"(Y)      ')
            CALL LINCTL  (ARYUU  ,  16 , 'Y', HILFL, KKYL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1, 17) .GE. 1) THEN
            CALL DPHI0   (KKYL,JJ1L,ILIMXP,KKYL,JJ1L,ILIMXP,SRYUU  )
      ENDIF
      IF(ISELEP(1, 18) .GE. 1) THEN
         IF(ISELEA(1, 18) .EQ. 0  .OR.  ISELEA(1, 18) .EQ. 3) THEN
            CALL DPHI0   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILIMXP,ARZUU  )
         ELSE
            CT = CH1680  (' R_U"U"(Z)      ')
            CALL LINCTL  (ARZUU  ,  18 , 'Z', HILFL, KKZL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1, 19) .GE. 1) THEN
            CALL DPHI0   (KKZL,JJ1L,ILIMXP,KKZL,JJ1L,ILIMXP,SRZUU  )
      ENDIF
      IF(ISELEP(1, 20) .GE. 1) THEN
         IF(ISELEA(1, 20) .EQ. 0  .OR.  ISELEA(1, 20) .EQ. 3) THEN
            CALL DPHI0   (KSXL,JJ1L,ILIMXP,KSXL,JJ1L,ILIMXP,ASXUU  )
         ELSE
            CT = CH1680  (' S_U"U"(K_X)    ')
            CALL LINCTL  (ASXUU  ,  20 , 'X', HILFL, KSXL, JJ1L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1, 21) .GE. 1) THEN
            CALL DPHI0   (KSXL,JJ1L,ILIMXP,KSXL,JJ1L,ILIMXP,SSXUU  )
      ENDIF
      IF(ISELEP(1, 22) .GE. 1) THEN
         IF(ISELEA(1, 22) .EQ. 0  .OR.  ISELEA(1, 22) .EQ. 3) THEN
            CALL DPHI0   (KSYL,JJ1L,ILIMXP,KSYL,JJ1L,ILIMXP,ASYUU  )
         ELSE
            CT = CH1680  (' S_U"U"(K_Y)    ')
            CALL LINCTL  (ASYUU  ,  22 , 'Y', HILFL, KSYL, JJ1L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1, 23) .GE. 1) THEN
            CALL DPHI0   (KSYL,JJ1L,ILIMXP,KSYL,JJ1L,ILIMXP,SSYUU  )
      ENDIF
      IF(ISELEP(1, 24) .GE. 1) THEN
         IF(ISELEA(1, 24) .EQ. 0  .OR.  ISELEA(1, 24) .EQ. 3) THEN
            CALL DPHI0   (KSZL,JJ1L,ILIMXP,KSZL,JJ1L,ILIMXP,ASZUU  )
         ELSE
            CT = CH1680  (' S_U"U"(K_Z)    ')
            CALL LINCTL  (ASZUU  ,  24 , 'Z', HILFL, KSZL, JJ1L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1, 25) .GE. 1) THEN
            CALL DPHI0   (KSZL,JJ1L,ILIMXP,KSZL,JJ1L,ILIMXP,SSZUU  )
      ENDIF
C
C                                 STATIST. GROESSEN DER V-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 39) .GE. 1) THEN
         IF(ISELEA(1, 39) .EQ. 0  .OR.  ISELEA(1, 39) .EQ. 3) THEN
            CALL DPHI0   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILIMXP,ARXVV  )
         ELSE
            CT = CH1680  (' R_V"V"(X)      ')
            CALL LINCTL  (ARXVV  ,  39 , 'X', HILFL, KKXL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1, 40) .GE. 1) THEN
            CALL DPHI0   (KKXL,JJ1L,ILIMXP,KKXL,JJ1L,ILIMXP,SRXVV  )
      ENDIF
      IF(ISELEP(1, 41) .GE. 1) THEN
         IF(ISELEA(1, 41) .EQ. 0  .OR.  ISELEA(1, 41) .EQ. 3) THEN
            CALL DPHI0   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILIMXP,ARYVV  )
         ELSE
            CT = CH1680  (' R_V"V"(Y)      ')
            CALL LINCTL  (ARYVV  ,  41 , 'Y', HILFL, KKYL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1, 42) .GE. 1) THEN
            CALL DPHI0   (KKYL,JJ1L,ILIMXP,KKYL,JJ1L,ILIMXP,SRYVV  )
      ENDIF
      IF(ISELEP(1, 43) .GE. 1) THEN
         IF(ISELEA(1, 43) .EQ. 0  .OR.  ISELEA(1, 43) .EQ. 3) THEN
            CALL DPHI0   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILIMXP,ARZVV  )
         ELSE
            CT = CH1680  (' R_V"V"(Z)      ')
            CALL LINCTL  (ARZVV  ,  43 , 'Z', HILFL, KKZL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1, 44) .GE. 1) THEN
            CALL DPHI0   (KKZL,JJ1L,ILIMXP,KKZL,JJ1L,ILIMXP,SRZVV  )
      ENDIF
      IF(ISELEP(1, 45) .GE. 1) THEN
         IF(ISELEA(1, 45) .EQ. 0  .OR.  ISELEA(1, 45) .EQ. 3) THEN
            CALL DPHI0   (KSXL,JJ1L,ILIMXP,KSXL,JJ1L,ILIMXP,ASXVV  )
         ELSE
            CT = CH1680  (' S_V"V"(K_X)    ')
            CALL LINCTL  (ASXVV  ,  45 , 'X', HILFL, KSXL, JJ1L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1, 46) .GE. 1) THEN
            CALL DPHI0   (KSXL,JJ1L,ILIMXP,KSXL,JJ1L,ILIMXP,SSXVV  )
      ENDIF
      IF(ISELEP(1, 47) .GE. 1) THEN
         IF(ISELEA(1, 47) .EQ. 0  .OR.  ISELEA(1, 47) .EQ. 3) THEN
            CALL DPHI0   (KSYL,JJ1L,ILIMXP,KSYL,JJ1L,ILIMXP,ASYVV  )
         ELSE
            CT = CH1680  (' S_V"V"(K_Y)    ')
            CALL LINCTL  (ASYVV  ,  47 , 'Y', HILFL, KSYL, JJ1L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1, 48) .GE. 1) THEN
            CALL DPHI0   (KSYL,JJ1L,ILIMXP,KSYL,JJ1L,ILIMXP,SSYVV  )
      ENDIF
      IF(ISELEP(1, 49) .GE. 1) THEN
         IF(ISELEA(1, 49) .EQ. 0  .OR.  ISELEA(1, 49) .EQ. 3) THEN
            CALL DPHI0   (KSZL,JJ1L,ILIMXP,KSZL,JJ1L,ILIMXP,ASZVV  )
         ELSE
            CT = CH1680  (' S_V"V"(K_Z)    ')
            CALL LINCTL  (ASZVV  ,  49 , 'Z', HILFL, KSZL, JJ1L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1, 50) .GE. 1) THEN
            CALL DPHI0   (KSZL,JJ1L,ILIMXP,KSZL,JJ1L,ILIMXP,SSZVV  )
      ENDIF
C
C                                 STATIST. GROESSEN DER W-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 64) .GE. 1) THEN
         IF(ISELEA(1, 64) .EQ. 0  .OR.  ISELEA(1, 64) .EQ. 3) THEN
            CALL DPHI0   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILIMXP,ARXWW  )
         ELSE
            CT = CH1680  (' R_W"W"(X)      ')
            CALL LINCTL  (ARXWW  ,  64 , 'X', HILFL, KKXL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1, 65) .GE. 1) THEN
            CALL DPHI0   (KKXL,JJ1L,ILIMXP,KKXL,JJ1L,ILIMXP,SRXWW  )
      ENDIF
      IF(ISELEP(1, 66) .GE. 1) THEN
         IF(ISELEA(1, 66) .EQ. 0  .OR.  ISELEA(1, 66) .EQ. 3) THEN
            CALL DPHI0   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILIMXP,ARYWW  )
         ELSE
            CT = CH1680  (' R_W"W"(Y)      ')
            CALL LINCTL  (ARYWW  ,  66 , 'Y', HILFL, KKYL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1, 67) .GE. 1) THEN
            CALL DPHI0   (KKYL,JJ1L,ILIMXP,KKYL,JJ1L,ILIMXP,SRYWW  )
      ENDIF
      IF(ISELEP(1, 68) .GE. 1) THEN
         IF(ISELEA(1, 68) .EQ. 0  .OR.  ISELEA(1, 68) .EQ. 3) THEN
            CALL DPHI0   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILIMXP,ARZWW  )
         ELSE
            CT = CH1680  (' R_W"W"(Z)      ')
            CALL LINCTL  (ARZWW  ,  68 , 'Z', HILFL, KKZL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1, 69) .GE. 1) THEN
            CALL DPHI0   (KKZL,JJ1L,ILIMXP,KKZL,JJ1L,ILIMXP,SRZWW  )
      ENDIF
      IF(ISELEP(1, 70) .GE. 1) THEN
         IF(ISELEA(1, 70) .EQ. 0  .OR.  ISELEA(1, 70) .EQ. 3) THEN
            CALL DPHI0   (KSXL,JJ1L,ILIMXP,KSXL,JJ1L,ILIMXP,ASXWW  )
         ELSE
            CT = CH1680  (' S_W"W"(K_X)    ')
            CALL LINCTL  (ASXWW  ,  70 , 'X', HILFL, KSXL, JJ1L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1, 71) .GE. 1) THEN
            CALL DPHI0   (KSXL,JJ1L,ILIMXP,KSXL,JJ1L,ILIMXP,SSXWW  )
      ENDIF
      IF(ISELEP(1, 72) .GE. 1) THEN
         IF(ISELEA(1, 72) .EQ. 0  .OR.  ISELEA(1, 72) .EQ. 3) THEN
            CALL DPHI0   (KSYL,JJ1L,ILIMXP,KSYL,JJ1L,ILIMXP,ASYWW  )
         ELSE
            CT = CH1680  (' S_W"W"(K_Y)    ')
            CALL LINCTL  (ASYWW  ,  72 , 'Y', HILFL, KSYL, JJ1L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1, 73) .GE. 1) THEN
            CALL DPHI0   (KSYL,JJ1L,ILIMXP,KSYL,JJ1L,ILIMXP,SSYWW  )
      ENDIF
      IF(ISELEP(1, 74) .GE. 1) THEN
         IF(ISELEA(1, 74) .EQ. 0  .OR.  ISELEA(1, 74) .EQ. 3) THEN
            CALL DPHI0   (KSZL,JJ1L,ILIMXP,KSZL,JJ1L,ILIMXP,ASZWW  )
         ELSE
            CT = CH1680  (' S_W"W"(K_Z)    ')
            CALL LINCTL  (ASZWW  ,  74 , 'Z', HILFL, KSZL, JJ1L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1, 75) .GE. 1) THEN
            CALL DPHI0   (KSZL,JJ1L,ILIMXP,KSZL,JJ1L,ILIMXP,SSZWW  )
      ENDIF
C
C                                  HIER: KREUZKORRELATIONSFUNKTIONEN
C                                  ---------------------------------
C
      IF(ISELEP(1,210) .GE. 1) THEN
         IF(ISELEA(1,210) .EQ. 0  .OR.  ISELEA(1,210) .EQ. 3) THEN
            CALL DPHI0   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILIMXP,ARXUV  )
         ELSE
            CT = CH1680  (' R_U"V"(X)      ')
            CALL LINCTL  (ARXUV  , 210 , 'X', HILFL, KKXL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,211) .GE. 1) THEN
            CALL DPHI0   (KKXL,JJ1L,ILIMXP,KKXL,JJ1L,ILIMXP,SRXUV  )
      ENDIF
      IF(ISELEP(1,212) .GE. 1) THEN
         IF(ISELEA(1,212) .EQ. 0  .OR.  ISELEA(1,212) .EQ. 3) THEN
            CALL DPHI0   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILIMXP,ARYUV  )
         ELSE
            CT = CH1680  (' R_U"V"(Y)      ')
            CALL LINCTL  (ARYUV  , 212 , 'Y', HILFL, KKYL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,213) .GE. 1) THEN
            CALL DPHI0   (KKYL,JJ1L,ILIMXP,KKYL,JJ1L,ILIMXP,SRYUV  )
      ENDIF
      IF(ISELEP(1,216) .GE. 1) THEN
         IF(ISELEA(1,216) .EQ. 0  .OR.  ISELEA(1,216) .EQ. 3) THEN
            CALL DPHI0   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILIMXP,ARXUW  )
         ELSE
            CT = CH1680  (' R_U"W"(X)      ')
            CALL LINCTL  (ARXUW  , 216 , 'X', HILFL, KKXL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,217) .GE. 1) THEN
            CALL DPHI0   (KKXL,JJ1L,ILIMXP,KKXL,JJ1L,ILIMXP,SRXUW  )
      ENDIF
      IF(ISELEP(1,218) .GE. 1) THEN
         IF(ISELEA(1,218) .EQ. 0  .OR.  ISELEA(1,218) .EQ. 3) THEN
            CALL DPHI0   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILIMXP,ARYUW  )
         ELSE
            CT = CH1680  (' R_U"W"(Y)      ')
            CALL LINCTL  (ARYUW  , 218 , 'Y', HILFL, KKYL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,219) .GE. 1) THEN
            CALL DPHI0   (KKYL,JJ1L,ILIMXP,KKYL,JJ1L,ILIMXP,SRYUW  )
      ENDIF
      IF(ISELEP(1,224) .GE. 1) THEN
         IF(ISELEA(1,224) .EQ. 0  .OR.  ISELEA(1,224) .EQ. 3) THEN
            CALL DPHI0   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILIMXP,ARYVW  )
         ELSE
            CT = CH1680  (' R_V"W"(Y)      ')
            CALL LINCTL  (ARYVW  , 224 , 'Y', HILFL, KKYL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,225) .GE. 1) THEN
            CALL DPHI0   (KKYL,JJ1L,ILIMXP,KKYL,JJ1L,ILIMXP,SRYVW  )
      ENDIF
C
C                                  HIER: KORRELATIONSKOEFFIZIENTEN
C                                  -------------------------------
C
      IF(ISELEP(1,234) .GE. 1) THEN
         IF(ISELEA(1,234) .EQ. 0  .OR.  ISELEA(1,234) .EQ. 3) THEN
            CALL DPHI0   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILIMXP,ACXUW  )
         ELSE
            CT = CH1680  (' RC_U"W"(X)     ')
            CALL LINCTL  (ACXUW  , 234 , 'X', HILFL, KKXL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,235) .GE. 1) THEN
            CALL DPHI0   (KKXL,JJ1L,ILIMXP,KKXL,JJ1L,ILIMXP,SCXUW  )
      ENDIF
      IF(ISELEP(1,238) .GE. 1) THEN
         IF(ISELEA(1,238) .EQ. 0  .OR.  ISELEA(1,238) .EQ. 3) THEN
            CALL DPHI0   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILIMXP,ACZUW  )
         ELSE
            CT = CH1680  (' RC_U"W"(Z)     ')
            CALL LINCTL  (ACZUW  , 238 , 'Z', HILFL, KKZL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,239) .GE. 1) THEN
            CALL DPHI0   (KKZL,JJ1L,ILIMXP,KKZL,JJ1L,ILIMXP,SCZUW  )
      ENDIF
C
C                                  OBERE DREIECKSMATRIX DES TENSORS
C                                  DER KORRELATIONSFUNKTIONEN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  --------------------------------
C
      IF(ISELEP(1,246) .GE. 1) THEN
         IF(ISELEA(1,246) .EQ. 0  .OR.  ISELEA(1,246) .EQ. 3) THEN
            CALL DPHI0   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILIMXP,ARXOXX )
         ELSE
            CT = CH1680  (' R_OX"OX"(X)    ')
            CALL LINCTL  (ARXOXX , 246 , 'X', HILFL, KKXL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,247) .GE. 1) THEN
            CALL DPHI0   (KKXL,JJ1L,ILIMXP,KKXL,JJ1L,ILIMXP,SRXOXX )
      ENDIF
      IF(ISELEP(1,248) .GE. 1) THEN
         IF(ISELEA(1,248) .EQ. 0  .OR.  ISELEA(1,248) .EQ. 3) THEN
            CALL DPHI0   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILIMXP,ARXOXY )
         ELSE
            CT = CH1680  (' R_OX"OY"(X)    ')
            CALL LINCTL  (ARXOXY , 248 , 'X', HILFL, KKXL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,249) .GE. 1) THEN
            CALL DPHI0   (KKXL,JJ1L,ILIMXP,KKXL,JJ1L,ILIMXP,SRXOXY )
      ENDIF
      IF(ISELEP(1,250) .GE. 1) THEN
         IF(ISELEA(1,250) .EQ. 0  .OR.  ISELEA(1,250) .EQ. 3) THEN
            CALL DPHI0   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILIMXP,ARXOXZ )
         ELSE
            CT = CH1680  (' R_OX"OZ"(X)    ')
            CALL LINCTL  (ARXOXZ , 250 , 'X', HILFL, KKXL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,251) .GE. 1) THEN
            CALL DPHI0   (KKXL,JJ1L,ILIMXP,KKXL,JJ1L,ILIMXP,SRXOXZ )
      ENDIF
      IF(ISELEP(1,252) .GE. 1) THEN
         IF(ISELEA(1,252) .EQ. 0  .OR.  ISELEA(1,252) .EQ. 3) THEN
            CALL DPHI0   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILIMXP,ARXOYY )
         ELSE
            CT = CH1680  (' R_OY"OY"(X)    ')
            CALL LINCTL  (ARXOYY , 252 , 'X', HILFL, KKXL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,253) .GE. 1) THEN
            CALL DPHI0   (KKXL,JJ1L,ILIMXP,KKXL,JJ1L,ILIMXP,SRXOYY )
      ENDIF
      IF(ISELEP(1,254) .GE. 1) THEN
         IF(ISELEA(1,254) .EQ. 0  .OR.  ISELEA(1,254) .EQ. 3) THEN
            CALL DPHI0   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILIMXP,ARXOYZ )
         ELSE
            CT = CH1680  (' R_OY"OZ"(X)    ')
            CALL LINCTL  (ARXOYZ , 254 , 'X', HILFL, KKXL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,255) .GE. 1) THEN
            CALL DPHI0   (KKXL,JJ1L,ILIMXP,KKXL,JJ1L,ILIMXP,SRXOYZ )
      ENDIF
      IF(ISELEP(1,256) .GE. 1) THEN
         IF(ISELEA(1,256) .EQ. 0  .OR.  ISELEA(1,256) .EQ. 3) THEN
            CALL DPHI0   (KKXL,JJ2L,ILIMXP,KKXL,JJ2L,ILIMXP,ARXOZZ )
         ELSE
            CT = CH1680  (' R_OZ"OZ"(X)    ')
            CALL LINCTL  (ARXOZZ , 256 , 'X', HILFL, KKXL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,257) .GE. 1) THEN
            CALL DPHI0   (KKXL,JJ1L,ILIMXP,KKXL,JJ1L,ILIMXP,SRXOZZ )
      ENDIF
      IF(ISELEP(1,258) .GE. 1) THEN
         IF(ISELEA(1,258) .EQ. 0  .OR.  ISELEA(1,258) .EQ. 3) THEN
            CALL DPHI0   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILIMXP,ARYOXX )
         ELSE
            CT = CH1680  (' R_OX"OX"(Y)    ')
            CALL LINCTL  (ARYOXX , 258 , 'Y', HILFL, KKYL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,259) .GE. 1) THEN
            CALL DPHI0   (KKYL,JJ1L,ILIMXP,KKYL,JJ1L,ILIMXP,SRYOXX )
      ENDIF
      IF(ISELEP(1,260) .GE. 1) THEN
         IF(ISELEA(1,260) .EQ. 0  .OR.  ISELEA(1,260) .EQ. 3) THEN
            CALL DPHI0   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILIMXP,ARYOXY )
         ELSE
            CT = CH1680  (' R_OX"OY"(Y)    ')
            CALL LINCTL  (ARYOXY , 260 , 'Y', HILFL, KKYL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,261) .GE. 1) THEN
            CALL DPHI0   (KKYL,JJ1L,ILIMXP,KKYL,JJ1L,ILIMXP,SRYOXY )
      ENDIF
      IF(ISELEP(1,262) .GE. 1) THEN
         IF(ISELEA(1,262) .EQ. 0  .OR.  ISELEA(1,262) .EQ. 3) THEN
            CALL DPHI0   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILIMXP,ARYOXZ )
         ELSE
            CT = CH1680  (' R_OX"OZ"(Y)    ')
            CALL LINCTL  (ARYOXZ , 262 , 'Y', HILFL, KKYL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,263) .GE. 1) THEN
            CALL DPHI0   (KKYL,JJ1L,ILIMXP,KKYL,JJ1L,ILIMXP,SRYOXZ )
      ENDIF
      IF(ISELEP(1,264) .GE. 1) THEN
         IF(ISELEA(1,264) .EQ. 0  .OR.  ISELEA(1,264) .EQ. 3) THEN
            CALL DPHI0   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILIMXP,ARYOYY )
         ELSE
            CT = CH1680  (' R_OY"OY"(Y)    ')
            CALL LINCTL  (ARYOYY , 264 , 'Y', HILFL, KKYL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,265) .GE. 1) THEN
            CALL DPHI0   (KKYL,JJ1L,ILIMXP,KKYL,JJ1L,ILIMXP,SRYOYY )
      ENDIF
      IF(ISELEP(1,266) .GE. 1) THEN
         IF(ISELEA(1,266) .EQ. 0  .OR.  ISELEA(1,266) .EQ. 3) THEN
            CALL DPHI0   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILIMXP,ARYOYZ )
         ELSE
            CT = CH1680  (' R_OY"OZ"(Y)    ')
            CALL LINCTL  (ARYOYZ , 266 , 'Y', HILFL, KKYL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,267) .GE. 1) THEN
            CALL DPHI0   (KKYL,JJ1L,ILIMXP,KKYL,JJ1L,ILIMXP,SRYOYZ )
      ENDIF
      IF(ISELEP(1,268) .GE. 1) THEN
         IF(ISELEA(1,268) .EQ. 0  .OR.  ISELEA(1,268) .EQ. 3) THEN
            CALL DPHI0   (KKYL,JJ2L,ILIMXP,KKYL,JJ2L,ILIMXP,ARYOZZ )
         ELSE
            CT = CH1680  (' R_OZ"OZ"(Y)    ')
            CALL LINCTL  (ARYOZZ , 268 , 'Y', HILFL, KKYL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,269) .GE. 1) THEN
            CALL DPHI0   (KKYL,JJ1L,ILIMXP,KKYL,JJ1L,ILIMXP,SRYOZZ )
      ENDIF
      IF(ISELEP(1,270) .GE. 1) THEN
         IF(ISELEA(1,270) .EQ. 0  .OR.  ISELEA(1,270) .EQ. 3) THEN
            CALL DPHI0   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILIMXP,ARZOXX )
         ELSE
            CT = CH1680  (' R_OX"OX"(Z)    ')
            CALL LINCTL  (ARZOXX , 270 , 'Z', HILFL, KKZL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,271) .GE. 1) THEN
            CALL DPHI0   (KKZL,JJ1L,ILIMXP,KKZL,JJ1L,ILIMXP,SRZOXX )
      ENDIF
      IF(ISELEP(1,272) .GE. 1) THEN
         IF(ISELEA(1,272) .EQ. 0  .OR.  ISELEA(1,272) .EQ. 3) THEN
            CALL DPHI0   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILIMXP,ARZOXY )
         ELSE
            CT = CH1680  (' R_OX"OY"(Z)    ')
            CALL LINCTL  (ARZOXY , 272 , 'Z', HILFL, KKZL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,273) .GE. 1) THEN
            CALL DPHI0   (KKZL,JJ1L,ILIMXP,KKZL,JJ1L,ILIMXP,SRZOXY )
      ENDIF
      IF(ISELEP(1,274) .GE. 1) THEN
         IF(ISELEA(1,274) .EQ. 0  .OR.  ISELEA(1,274) .EQ. 3) THEN
            CALL DPHI0   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILIMXP,ARZOXZ )
         ELSE
            CT = CH1680  (' R_OX"OZ"(Z)    ')
            CALL LINCTL  (ARZOXZ , 274 , 'Z', HILFL, KKZL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,275) .GE. 1) THEN
            CALL DPHI0   (KKZL,JJ1L,ILIMXP,KKZL,JJ1L,ILIMXP,SRZOXZ )
      ENDIF
      IF(ISELEP(1,276) .GE. 1) THEN
         IF(ISELEA(1,276) .EQ. 0  .OR.  ISELEA(1,276) .EQ. 3) THEN
            CALL DPHI0   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILIMXP,ARZOYY )
         ELSE
            CT = CH1680  (' R_OY"OY"(Z)    ')
            CALL LINCTL  (ARZOYY , 276 , 'Z', HILFL, KKZL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,277) .GE. 1) THEN
            CALL DPHI0   (KKZL,JJ1L,ILIMXP,KKZL,JJ1L,ILIMXP,SRZOYY )
      ENDIF
      IF(ISELEP(1,278) .GE. 1) THEN
         IF(ISELEA(1,278) .EQ. 0  .OR.  ISELEA(1,278) .EQ. 3) THEN
            CALL DPHI0   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILIMXP,ARZOYZ )
         ELSE
            CT = CH1680  (' R_OY"OZ"(Z)    ')
            CALL LINCTL  (ARZOYZ , 278 , 'Z', HILFL, KKZL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,279) .GE. 1) THEN
            CALL DPHI0   (KKZL,JJ1L,ILIMXP,KKZL,JJ1L,ILIMXP,SRZOYZ )
      ENDIF
      IF(ISELEP(1,280) .GE. 1) THEN
         IF(ISELEA(1,280) .EQ. 0  .OR.  ISELEA(1,280) .EQ. 3) THEN
            CALL DPHI0   (KKZL,JJ2L,ILIMXP,KKZL,JJ2L,ILIMXP,ARZOZZ )
         ELSE
            CT = CH1680  (' R_OZ"OZ"(Z)    ')
            CALL LINCTL  (ARZOZZ , 280 , 'Z', HILFL, KKZL, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,281) .GE. 1) THEN
            CALL DPHI0   (KKZL,JJ1L,ILIMXP,KKZL,JJ1L,ILIMXP,SRZOZZ )
      ENDIF
C
C                                  LEISTUNGSDICHTESPEKTREN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  ----------------------------
C
      IF(ISELEP(1,282) .GE. 1) THEN
         IF(ISELEA(1,282) .EQ. 0  .OR.  ISELEA(1,282) .EQ. 3) THEN
            CALL DPHI0   (KSXL,JJ1L,ILIMXP,KSXL,JJ1L,ILIMXP,ASXOXX )
         ELSE
            CT = CH1680  (' S_OX"OX"(K_X)  ')
            CALL LINCTL  (ASXOXX , 282 , 'X', HILFL, KSXL, JJ1L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,283) .GE. 1) THEN
            CALL DPHI0   (KSXL,JJ1L,ILIMXP,KSXL,JJ1L,ILIMXP,SSXOXX )
      ENDIF
      IF(ISELEP(1,284) .GE. 1) THEN
         IF(ISELEA(1,284) .EQ. 0  .OR.  ISELEA(1,284) .EQ. 3) THEN
            CALL DPHI0   (KSXL,JJ1L,ILIMXP,KSXL,JJ1L,ILIMXP,ASXOYY )
         ELSE
            CT = CH1680  (' S_OY"OY"(K_X)  ')
            CALL LINCTL  (ASXOYY , 284 , 'X', HILFL, KSXL, JJ1L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,285) .GE. 1) THEN
            CALL DPHI0   (KSXL,JJ1L,ILIMXP,KSXL,JJ1L,ILIMXP,SSXOYY )
      ENDIF
      IF(ISELEP(1,286) .GE. 1) THEN
         IF(ISELEA(1,286) .EQ. 0  .OR.  ISELEA(1,286) .EQ. 3) THEN
            CALL DPHI0   (KSXL,JJ1L,ILIMXP,KSXL,JJ1L,ILIMXP,ASXOZZ )
         ELSE
            CT = CH1680  (' S_OZ"OZ"(K_X)  ')
            CALL LINCTL  (ASXOZZ , 286 , 'X', HILFL, KSXL, JJ1L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,287) .GE. 1) THEN
            CALL DPHI0   (KSXL,JJ1L,ILIMXP,KSXL,JJ1L,ILIMXP,SSXOZZ )
      ENDIF
      IF(ISELEP(1,288) .GE. 1) THEN
         IF(ISELEA(1,288) .EQ. 0  .OR.  ISELEA(1,288) .EQ. 3) THEN
            CALL DPHI0   (KSYL,JJ1L,ILIMXP,KSYL,JJ1L,ILIMXP,ASYOXX )
         ELSE
            CT = CH1680  (' S_OX"OX"(K_Y)  ')
            CALL LINCTL  (ASYOXX , 288 , 'Y', HILFL, KSYL, JJ1L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,289) .GE. 1) THEN
            CALL DPHI0   (KSYL,JJ1L,ILIMXP,KSYL,JJ1L,ILIMXP,SSYOXX )
      ENDIF
      IF(ISELEP(1,290) .GE. 1) THEN
         IF(ISELEA(1,290) .EQ. 0  .OR.  ISELEA(1,290) .EQ. 3) THEN
            CALL DPHI0   (KSYL,JJ1L,ILIMXP,KSYL,JJ1L,ILIMXP,ASYOYY )
         ELSE
            CT = CH1680  (' S_OY"OY"(K_Y)  ')
            CALL LINCTL  (ASYOYY , 290 , 'Y', HILFL, KSYL, JJ1L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,291) .GE. 1) THEN
            CALL DPHI0   (KSYL,JJ1L,ILIMXP,KSYL,JJ1L,ILIMXP,SSYOYY )
      ENDIF
      IF(ISELEP(1,292) .GE. 1) THEN
         IF(ISELEA(1,292) .EQ. 0  .OR.  ISELEA(1,292) .EQ. 3) THEN
            CALL DPHI0   (KSYL,JJ1L,ILIMXP,KSYL,JJ1L,ILIMXP,ASYOZZ )
         ELSE
            CT = CH1680  (' S_OZ"OZ"(K_Y)  ')
            CALL LINCTL  (ASYOZZ , 292 , 'Y', HILFL, KSYL, JJ1L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,293) .GE. 1) THEN
            CALL DPHI0   (KSYL,JJ1L,ILIMXP,KSYL,JJ1L,ILIMXP,SSYOZZ )
      ENDIF
C
C                                  HAEUFIGKEITSVERTEILUNG DER
C                                  INKLINATIONSWINKEL DER VORTICITY-
C                                  VEKTOREN
C                                  ---------------------------------
C
      IF(ISELEP(1,306) .GE. 1) THEN
         IF(ISELEA(1,306) .EQ. 0  .OR.  ISELEA(1,306) .EQ. 3) THEN
            CALL DPHI0   (KH0L,JJ2L,ILIMXP,KH0L,JJ2L,ILIMXP,AHOZOY )
         ELSE
            CT = CH1680  (' H (ATAN(OZ/OY))')
            CALL LINCTL  (AHOZOY , 306 , '0', HILFL, KH0L, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,307) .GE. 1) THEN
            CALL DPHI0   (KH0L,JJ2L,ILIMXP,KH0L,JJ2L,ILIMXP,SHOZOY )
      ENDIF
      IF(ISELEP(1,308) .GE. 1) THEN
         IF(ISELEA(1,308) .EQ. 0  .OR.  ISELEA(1,308) .EQ. 3) THEN
            CALL DPHI0   (KH0L,JJ2L,ILIMXP,KH0L,JJ2L,ILIMXP,AHOZOX )
         ELSE
            CT = CH1680  (' H (ATAN(OZ/OX))')
            CALL LINCTL  (AHOZOX , 308 , '0', HILFL, KH0L, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,309) .GE. 1) THEN
            CALL DPHI0   (KH0L,JJ2L,ILIMXP,KH0L,JJ2L,ILIMXP,SHOZOX )
      ENDIF
      IF(ISELEP(1,310) .GE. 1) THEN
         IF(ISELEA(1,310) .EQ. 0  .OR.  ISELEA(1,310) .EQ. 3) THEN
            CALL DPHI0   (KH0L,JJ2L,ILIMXP,KH0L,JJ2L,ILIMXP,AHOYOX )
         ELSE
            CT = CH1680  (' H (ATAN(OY/OX))')
            CALL LINCTL  (AHOYOX , 310 , '0', HILFL, KH0L, JJ2L, ILIMXP,
     $                    ISLINP, ISLINA, ISLIDI, RKOMXP, RKOMXA, CT,
     $                    ISELEP, ISELEA)
         ENDIF
      ENDIF
      IF(ISELEP(1,311) .GE. 1) THEN
            CALL DPHI0   (KH0L,JJ2L,ILIMXP,KH0L,JJ2L,ILIMXP,SHOYOX )
      ENDIF
C
      RETURN
      END
