










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
      SUBROUTINE SUMST1  (ILIMX,RKOMXP,NBND,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     %                    NXGRAE,NYGRAE,NZGRAE,XHOMOG,YHOMOG,ZHOMOG,

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
C        S U M S T 1      IN SUMST1 WERDEN DIE MOMENTANWERTE VERSCHIE-
C                         DENER STATIST. GROESSEN AUS DEM AKTUELLEN GE-
C                         SCHWINDIGKEITSFELD GEBILDET. DIESE MOMENTAN-
C                         WERTE WERDEN (EVENTUELL UNTER ZWISCHENSCHAL-
C                         TUNG EINER LINIEN- ODER FLAECHENMITTELUNG)
C                         ANSCHLIESSEND IN DAS JEWEILIGE SUMMATIONS-
C                         FELD SUMMIERT, SO DASS AM ENDE EINES LAUFES
C                         DIE ENSEMBLE-MITTELWERTE NEU BERECHNET BZW.
C                         VERBESSERT WERDEN KOENNEN.
C  A C H T U N G:         DIESE ROUTINE MUSS ERWEITERT WERDEN, FALLS
C                         DIE STATIST. AUSWERTUNG ERWEITERT WIRD !
C*STARLET***************************************************************
C
C PARAM: ILIMX          - TATSAECHLICH BENOETIGTE I-LINIEN WAEHREND
C                         DES MOMENTANEN LAUFES
C        RKOMXP         - MAXIMAL ZULAESSIGE KORRELATIONSLAENGE WAEHREND
C                         DES MOMENTANEN LAUFES
C        NBND           - ANZAHL DER RANDSCHICHTEN DES BERECHNUNGS-
C                         GEBIETES
C        X,Y,Z          - KOORDINATEN DER ZELLMITTELPUNKTE
C        DX,DY,DZ       - ABSTAND DER BASISZELLMITTELPUNKTE
C        DDX,DDY,DDZ    - ABMESSUNGEN DER BASISZELLEN
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
C UPROG                 : KO2VAR, ASDFUN, HRELOM
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        14.12.88 (HW)  : ORIGINAL AUS SUMSTA ABGELEITET
C
C*STARLET***************************************************************
C
      REAL       X(II),     Y(JJ),     Z(KK),
     $          DX(II),    DY(JJ),    DZ(KK),
     $         DDX(II),   DDY(JJ),   DDZ(KK)
C
      LOGICAL  XHOMOG,YHOMOG,ZHOMOG
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
      IF(ISELEP(1, 14) .GE. 1) THEN
         CALL KO2VAR  (UFG   ,AURG  ,  0  ,  0  ,  1  ,
     $                 UFG   ,AURG  ,  0  ,  0  ,  1  ,
     $                 ARXUU ,SRXUU ,KKXL,JJ2L,JJ1L,ILIMXP,
     $                 'R X             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2, 14) = ISELEP(2, 14) + 1
      ENDIF
      IF(ISELEP(1, 16) .GE. 1) THEN
         CALL KO2VAR  (UFG   ,AURG  ,  0  ,  0  ,  1  ,
     $                 UFG   ,AURG  ,  0  ,  0  ,  1  ,
     $                 ARYUU ,SRYUU ,KKYL,JJ2L,JJ1L,ILIMXP,
     $                 'R Y             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2, 16) = ISELEP(2, 16) + 1
      ENDIF
      IF(ISELEP(1, 18) .GE. 1) THEN
         CALL KO2VAR  (UFG   ,AURG  ,  0  ,  0  ,  1  ,
     $                 UFG   ,AURG  ,  0  ,  0  ,  1  ,
     $                 ARZUU ,SRZUU ,KKZL,JJ2L,JJ1L,ILIMXP,
     $                 'R Z             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2, 18) = ISELEP(2, 18) + 1
      ENDIF
      IF(ISELEP(1, 20) .GE. 1) THEN
         CALL ASDFUN  (UFG  ,  0  ,  0  ,  1  ,
     $                 ASXUU ,SSXUU ,KSXL,JJ1L,ILIMXP,
     $                 'S X             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 ILIMX,ISLINP,ISLIDI,FVT,IDIMF,NBND,
     $                 LINFB,NZERON,NXGRAE,NYGRAE,NZGRAE,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ISELEP(2, 20) = ISELEP(2, 20) + NZERON
      ENDIF
      IF(ISELEP(1, 22) .GE. 1) THEN
         CALL ASDFUN  (UFG  ,  0  ,  0  ,  1  ,
     $                 ASYUU ,SSYUU ,KSYL,JJ1L,ILIMXP,
     $                 'S Y             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 ILIMX,ISLINP,ISLIDI,FVT,IDIMF,NBND,
     $                 LINFB,NZERON,NXGRAE,NYGRAE,NZGRAE,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ISELEP(2, 22) = ISELEP(2, 22) + NZERON
      ENDIF
      IF(ISELEP(1, 24) .GE. 1) THEN
         CALL ASDFUN  (UFG  ,  0  ,  0  ,  1  ,
     $                 ASZUU ,SSZUU ,KSZL,JJ1L,ILIMXP,
     $                 'S Z             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 ILIMX,ISLINP,ISLIDI,FVT,IDIMF,NBND,
     $                 LINFB,NZERON,NXGRAE,NYGRAE,NZGRAE,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ISELEP(2, 24) = ISELEP(2, 24) + NZERON
      ENDIF
C
C                                 STATIST. GROESSEN DER V-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 39) .GE. 1) THEN
         CALL KO2VAR  (VFG   ,AVRG  ,  0  ,  1  ,  0  ,
     $                 VFG   ,AVRG  ,  0  ,  1  ,  0  ,
     $                 ARXVV ,SRXVV ,KKXL,JJ2L,JJ1L,ILIMXP,
     $                 'R X             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2, 39) = ISELEP(2, 39) + 1
      ENDIF
      IF(ISELEP(1, 41) .GE. 1) THEN
         CALL KO2VAR  (VFG   ,AVRG  ,  0  ,  1  ,  0  ,
     $                 VFG   ,AVRG  ,  0  ,  1  ,  0  ,
     $                 ARYVV ,SRYVV ,KKYL,JJ2L,JJ1L,ILIMXP,
     $                 'R Y             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2, 41) = ISELEP(2, 41) + 1
      ENDIF
      IF(ISELEP(1, 43) .GE. 1) THEN
         CALL KO2VAR  (VFG   ,AVRG  ,  0  ,  1  ,  0  ,
     $                 VFG   ,AVRG  ,  0  ,  1  ,  0  ,
     $                 ARZVV ,SRZVV ,KKZL,JJ2L,JJ1L,ILIMXP,
     $                 'R Z             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2, 43) = ISELEP(2, 43) + 1
      ENDIF
      IF(ISELEP(1, 45) .GE. 1) THEN
         CALL ASDFUN  (VFG  ,  0  ,  1  ,  0  ,
     $                 ASXVV ,SSXVV ,KSXL,JJ1L,ILIMXP,
     $                 'S X             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 ILIMX,ISLINP,ISLIDI,FVT,IDIMF,NBND,
     $                 LINFB,NZERON,NXGRAE,NYGRAE,NZGRAE,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ISELEP(2, 45) = ISELEP(2, 45) + NZERON
      ENDIF
      IF(ISELEP(1, 47) .GE. 1) THEN
         CALL ASDFUN  (VFG  ,  0  ,  1  ,  0  ,
     $                 ASYVV ,SSYVV ,KSYL,JJ1L,ILIMXP,
     $                 'S Y             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 ILIMX,ISLINP,ISLIDI,FVT,IDIMF,NBND,
     $                 LINFB,NZERON,NXGRAE,NYGRAE,NZGRAE,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ISELEP(2, 47) = ISELEP(2, 47) + NZERON
      ENDIF
      IF(ISELEP(1, 49) .GE. 1) THEN
         CALL ASDFUN  (VFG  ,  0  ,  1  ,  0  ,
     $                 ASZVV ,SSZVV ,KSZL,JJ1L,ILIMXP,
     $                 'S Z             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 ILIMX,ISLINP,ISLIDI,FVT,IDIMF,NBND,
     $                 LINFB,NZERON,NXGRAE,NYGRAE,NZGRAE,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ISELEP(2, 49) = ISELEP(2, 49) + NZERON
      ENDIF
C
C                                 STATIST. GROESSEN DER W-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 64) .GE. 1) THEN
         CALL KO2VAR  (WFG   ,AWRG  ,  1  ,  0  ,  0  ,
     $                 WFG   ,AWRG  ,  1  ,  0  ,  0  ,
     $                 ARXWW ,SRXWW ,KKXL,JJ2L,JJ1L,ILIMXP,
     $                 'R X             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2, 64) = ISELEP(2, 64) + 1
      ENDIF
      IF(ISELEP(1, 66) .GE. 1) THEN
         CALL KO2VAR  (WFG   ,AWRG  ,  1  ,  0  ,  0  ,
     $                 WFG   ,AWRG  ,  1  ,  0  ,  0  ,
     $                 ARYWW ,SRYWW ,KKYL,JJ2L,JJ1L,ILIMXP,
     $                 'R Y             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2, 66) = ISELEP(2, 66) + 1
      ENDIF
      IF(ISELEP(1, 68) .GE. 1) THEN
         CALL KO2VAR  (WFG   ,AWRG  ,  1  ,  0  ,  0  ,
     $                 WFG   ,AWRG  ,  1  ,  0  ,  0  ,
     $                 ARZWW ,SRZWW ,KKZL,JJ2L,JJ1L,ILIMXP,
     $                 'R Z             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2, 68) = ISELEP(2, 68) + 1
      ENDIF
      IF(ISELEP(1, 70) .GE. 1) THEN
         CALL ASDFUN  (WFG  ,  1  ,  0  ,  0  ,
     $                 ASXWW ,SSXWW ,KSXL,JJ1L,ILIMXP,
     $                 'S X             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 ILIMX,ISLINP,ISLIDI,FVT,IDIMF,NBND,
     $                 LINFB,NZERON,NXGRAE,NYGRAE,NZGRAE,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ISELEP(2, 70) = ISELEP(2, 70) + NZERON
      ENDIF
      IF(ISELEP(1, 72) .GE. 1) THEN
         CALL ASDFUN  (WFG  ,  1  ,  0  ,  0  ,
     $                 ASYWW ,SSYWW ,KSYL,JJ1L,ILIMXP,
     $                 'S Y             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 ILIMX,ISLINP,ISLIDI,FVT,IDIMF,NBND,
     $                 LINFB,NZERON,NXGRAE,NYGRAE,NZGRAE,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ISELEP(2, 72) = ISELEP(2, 72) + NZERON
      ENDIF
      IF(ISELEP(1, 74) .GE. 1) THEN
         CALL ASDFUN  (WFG  ,  1  ,  0  ,  0  ,
     $                 ASZWW ,SSZWW ,KSZL,JJ1L,ILIMXP,
     $                 'S Z             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 ILIMX,ISLINP,ISLIDI,FVT,IDIMF,NBND,
     $                 LINFB,NZERON,NXGRAE,NYGRAE,NZGRAE,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ISELEP(2, 74) = ISELEP(2, 74) + NZERON
      ENDIF
C
C                                  HIER: KREUZKORRELATIONSFUNKTIONEN
C                                  ---------------------------------
C
      IF(ISELEP(1,210) .GE. 1) THEN
         CALL KO2VAR  (UFG   ,AURG  ,  0  ,  0  ,  1  ,
     $                 VFG   ,AVRG  ,  0  ,  1  ,  0  ,
     $                 ARXUV ,SRXUV ,KKXL,JJ2L,JJ1L,ILIMXP,
     $                 'R X             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,210) = ISELEP(2,210) + 1
      ENDIF
      IF(ISELEP(1,212) .GE. 1) THEN
         CALL KO2VAR  (UFG   ,AURG  ,  0  ,  0  ,  1  ,
     $                 VFG   ,AVRG  ,  0  ,  1  ,  0  ,
     $                 ARYUV ,SRYUV ,KKYL,JJ2L,JJ1L,ILIMXP,
     $                 'R Y             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,212) = ISELEP(2,212) + 1
      ENDIF
      IF(ISELEP(1,216) .GE. 1) THEN
         CALL KO2VAR  (UFG   ,AURG  ,  0  ,  0  ,  1  ,
     $                 WFG   ,AWRG  ,  1  ,  0  ,  0  ,
     $                 ARXUW ,SRXUW ,KKXL,JJ2L,JJ1L,ILIMXP,
     $                 'R X             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,216) = ISELEP(2,216) + 1
      ENDIF
      IF(ISELEP(1,218) .GE. 1) THEN
         CALL KO2VAR  (UFG   ,AURG  ,  0  ,  0  ,  1  ,
     $                 WFG   ,AWRG  ,  1  ,  0  ,  0  ,
     $                 ARYUW ,SRYUW ,KKYL,JJ2L,JJ1L,ILIMXP,
     $                 'R Y             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,218) = ISELEP(2,218) + 1
      ENDIF
      IF(ISELEP(1,224) .GE. 1) THEN
         CALL KO2VAR  (VFG   ,AVRG  ,  0  ,  1  ,  0  ,
     $                 WFG   ,AWRG  ,  1  ,  0  ,  0  ,
     $                 ARYVW ,SRYVW ,KKYL,JJ2L,JJ1L,ILIMXP,
     $                 'R Y             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,224) = ISELEP(2,224) + 1
      ENDIF
C
C                                  HIER: KORRELATIONSKOEFFIZIENTEN
C                                  -------------------------------
C
      IF(ISELEP(1,234) .GE. 1) THEN
         CALL KO2VAR  (UFG   ,AURG  ,  0  ,  0  ,  1  ,
     $                 WFG   ,AWRG  ,  1  ,  0  ,  0  ,
     $                 ACXUW ,SCXUW ,KKXL,JJ2L,JJ1L,ILIMXP,
     $                 'RCX             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,234) = ISELEP(2,234) + 1
      ENDIF
      IF(ISELEP(1,238) .GE. 1) THEN
         CALL KO2VAR  (UFG   ,AURG  ,  0  ,  0  ,  1  ,
     $                 WFG   ,AWRG  ,  1  ,  0  ,  0  ,
     $                 ACZUW ,SCZUW ,KKZL,JJ2L,JJ1L,ILIMXP,
     $                 'RCZ             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,238) = ISELEP(2,238) + 1
      ENDIF
C
C                                  OBERE DREIECKSMATRIX DES TENSORS
C                                  DER KORRELATIONSFUNKTIONEN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  --------------------------------
C
      IF(ISELEP(1,246) .GE. 1) THEN
         CALL KO2VAR  (OXFG  ,AOXRG ,  1  ,  1  ,  0  ,
     $                 OXFG  ,AOXRG ,  1  ,  1  ,  0  ,
     $                 ARXOXX,SRXOXX,KKXL,JJ2L,JJ1L,ILIMXP,
     $                 'R X             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,246) = ISELEP(2,246) + 1
      ENDIF
      IF(ISELEP(1,248) .GE. 1) THEN
         CALL KO2VAR  (OXFG  ,AOXRG ,  1  ,  1  ,  0  ,
     $                 OYFG  ,AOYRG ,  1  ,  0  ,  1  ,
     $                 ARXOXY,SRXOXY,KKXL,JJ2L,JJ1L,ILIMXP,
     $                 'R X             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,248) = ISELEP(2,248) + 1
      ENDIF
      IF(ISELEP(1,250) .GE. 1) THEN
         CALL KO2VAR  (OXFG  ,AOXRG ,  1  ,  1  ,  0  ,
     $                 OZFG  ,AOZRG ,  0  ,  1  ,  1  ,
     $                 ARXOXZ,SRXOXZ,KKXL,JJ2L,JJ1L,ILIMXP,
     $                 'R X             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,250) = ISELEP(2,250) + 1
      ENDIF
      IF(ISELEP(1,252) .GE. 1) THEN
         CALL KO2VAR  (OYFG  ,AOYRG ,  1  ,  0  ,  1  ,
     $                 OYFG  ,AOYRG ,  1  ,  0  ,  1  ,
     $                 ARXOYY,SRXOYY,KKXL,JJ2L,JJ1L,ILIMXP,
     $                 'R X             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,252) = ISELEP(2,252) + 1
      ENDIF
      IF(ISELEP(1,254) .GE. 1) THEN
         CALL KO2VAR  (OYFG  ,AOYRG ,  1  ,  0  ,  1  ,
     $                 OZFG  ,AOZRG ,  0  ,  1  ,  1  ,
     $                 ARXOYZ,SRXOYZ,KKXL,JJ2L,JJ1L,ILIMXP,
     $                 'R X             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,254) = ISELEP(2,254) + 1
      ENDIF
      IF(ISELEP(1,256) .GE. 1) THEN
         CALL KO2VAR  (OZFG  ,AOZRG ,  0  ,  1  ,  1  ,
     $                 OZFG  ,AOZRG ,  0  ,  1  ,  1  ,
     $                 ARXOZZ,SRXOZZ,KKXL,JJ2L,JJ1L,ILIMXP,
     $                 'R X             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,256) = ISELEP(2,256) + 1
      ENDIF
      IF(ISELEP(1,258) .GE. 1) THEN
         CALL KO2VAR  (OXFG  ,AOXRG ,  1  ,  1  ,  0  ,
     $                 OXFG  ,AOXRG ,  1  ,  1  ,  0  ,
     $                 ARYOXX,SRYOXX,KKYL,JJ2L,JJ1L,ILIMXP,
     $                 'R Y             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,258) = ISELEP(2,258) + 1
      ENDIF
      IF(ISELEP(1,260) .GE. 1) THEN
         CALL KO2VAR  (OXFG  ,AOXRG ,  1  ,  1  ,  0  ,
     $                 OYFG  ,AOYRG ,  1  ,  0  ,  1  ,
     $                 ARYOXY,SRYOXY,KKYL,JJ2L,JJ1L,ILIMXP,
     $                 'R Y             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,260) = ISELEP(2,260) + 1
      ENDIF
      IF(ISELEP(1,262) .GE. 1) THEN
         CALL KO2VAR  (OXFG  ,AOXRG ,  1  ,  1  ,  0  ,
     $                 OZFG  ,AOZRG ,  0  ,  1  ,  1  ,
     $                 ARYOXZ,SRYOXZ,KKYL,JJ2L,JJ1L,ILIMXP,
     $                 'R Y             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,262) = ISELEP(2,262) + 1
      ENDIF
      IF(ISELEP(1,264) .GE. 1) THEN
         CALL KO2VAR  (OYFG  ,AOYRG ,  1  ,  0  ,  1  ,
     $                 OYFG  ,AOYRG ,  1  ,  0  ,  1  ,
     $                 ARYOYY,SRYOYY,KKYL,JJ2L,JJ1L,ILIMXP,
     $                 'R Y             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,264) = ISELEP(2,264) + 1
      ENDIF
      IF(ISELEP(1,266) .GE. 1) THEN
         CALL KO2VAR  (OYFG  ,AOYRG ,  1  ,  0  ,  1  ,
     $                 OZFG  ,AOZRG ,  0  ,  1  ,  1  ,
     $                 ARYOYZ,SRYOYZ,KKYL,JJ2L,JJ1L,ILIMXP,
     $                 'R Y             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,266) = ISELEP(2,266) + 1
      ENDIF
      IF(ISELEP(1,268) .GE. 1) THEN
         CALL KO2VAR  (OZFG  ,AOZRG ,  0  ,  1  ,  1  ,
     $                 OZFG  ,AOZRG ,  0  ,  1  ,  1  ,
     $                 ARYOZZ,SRYOZZ,KKYL,JJ2L,JJ1L,ILIMXP,
     $                 'R Y             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,268) = ISELEP(2,268) + 1
      ENDIF
      IF(ISELEP(1,270) .GE. 1) THEN
         CALL KO2VAR  (OXFG  ,AOXRG ,  1  ,  1  ,  0  ,
     $                 OXFG  ,AOXRG ,  1  ,  1  ,  0  ,
     $                 ARZOXX,SRZOXX,KKZL,JJ2L,JJ1L,ILIMXP,
     $                 'R Z             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,270) = ISELEP(2,270) + 1
      ENDIF
      IF(ISELEP(1,272) .GE. 1) THEN
         CALL KO2VAR  (OXFG  ,AOXRG ,  1  ,  1  ,  0  ,
     $                 OYFG  ,AOYRG ,  1  ,  0  ,  1  ,
     $                 ARZOXY,SRZOXY,KKZL,JJ2L,JJ1L,ILIMXP,
     $                 'R Z             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,272) = ISELEP(2,272) + 1
      ENDIF
      IF(ISELEP(1,274) .GE. 1) THEN
         CALL KO2VAR  (OXFG  ,AOXRG ,  1  ,  1  ,  0  ,
     $                 OZFG  ,AOZRG ,  0  ,  1  ,  1  ,
     $                 ARZOXZ,SRZOXZ,KKZL,JJ2L,JJ1L,ILIMXP,
     $                 'R Z             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,274) = ISELEP(2,274) + 1
      ENDIF
      IF(ISELEP(1,276) .GE. 1) THEN
         CALL KO2VAR  (OYFG  ,AOYRG ,  1  ,  0  ,  1  ,
     $                 OYFG  ,AOYRG ,  1  ,  0  ,  1  ,
     $                 ARZOYY,SRZOYY,KKZL,JJ2L,JJ1L,ILIMXP,
     $                 'R Z             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,276) = ISELEP(2,276) + 1
      ENDIF
      IF(ISELEP(1,278) .GE. 1) THEN
         CALL KO2VAR  (OYFG  ,AOYRG ,  1  ,  0  ,  1  ,
     $                 OZFG  ,AOZRG ,  0  ,  1  ,  1  ,
     $                 ARZOYZ,SRZOYZ,KKZL,JJ2L,JJ1L,ILIMXP,
     $                 'R Z             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,278) = ISELEP(2,278) + 1
      ENDIF
      IF(ISELEP(1,280) .GE. 1) THEN
         CALL KO2VAR  (OZFG  ,AOZRG ,  0  ,  1  ,  1  ,
     $                 OZFG  ,AOZRG ,  0  ,  1  ,  1  ,
     $                 ARZOZZ,SRZOZZ,KKZL,JJ2L,JJ1L,ILIMXP,
     $                 'R Z             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                 NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
         ISELEP(2,280) = ISELEP(2,280) + 1
      ENDIF
C
C                                  LEISTUNGSDICHTESPEKTREN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  ----------------------------
C
      IF(ISELEP(1,282) .GE. 1) THEN
         CALL ASDFUN  (OXFG ,  1  ,  1  ,  0  ,
     $                 ASXOXX,SSXOXX,KSXL,JJ1L,ILIMXP,
     $                 'S X             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 ILIMX,ISLINP,ISLIDI,FVT,IDIMF,NBND,
     $                 LINFB,NZERON,NXGRAE,NYGRAE,NZGRAE,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ISELEP(2,282) = ISELEP(2,282) + NZERON
      ENDIF
      IF(ISELEP(1,284) .GE. 1) THEN
         CALL ASDFUN  (OYFG ,  1  ,  0  ,  1  ,
     $                 ASXOYY,SSXOYY,KSXL,JJ1L,ILIMXP,
     $                 'S X             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 ILIMX,ISLINP,ISLIDI,FVT,IDIMF,NBND,
     $                 LINFB,NZERON,NXGRAE,NYGRAE,NZGRAE,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ISELEP(2,284) = ISELEP(2,284) + NZERON
      ENDIF
      IF(ISELEP(1,286) .GE. 1) THEN
         CALL ASDFUN  (OZFG ,  0  ,  1  ,  1  ,
     $                 ASXOZZ,SSXOZZ,KSXL,JJ1L,ILIMXP,
     $                 'S X             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 ILIMX,ISLINP,ISLIDI,FVT,IDIMF,NBND,
     $                 LINFB,NZERON,NXGRAE,NYGRAE,NZGRAE,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ISELEP(2,286) = ISELEP(2,286) + NZERON
      ENDIF
      IF(ISELEP(1,288) .GE. 1) THEN
         CALL ASDFUN  (OXFG ,  1  ,  1  ,  0  ,
     $                 ASYOXX,SSYOXX,KSYL,JJ1L,ILIMXP,
     $                 'S Y             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 ILIMX,ISLINP,ISLIDI,FVT,IDIMF,NBND,
     $                 LINFB,NZERON,NXGRAE,NYGRAE,NZGRAE,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ISELEP(2,288) = ISELEP(2,288) + NZERON
      ENDIF
      IF(ISELEP(1,290) .GE. 1) THEN
         CALL ASDFUN  (OYFG ,  1  ,  0  ,  1  ,
     $                 ASYOYY,SSYOYY,KSYL,JJ1L,ILIMXP,
     $                 'S Y             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 ILIMX,ISLINP,ISLIDI,FVT,IDIMF,NBND,
     $                 LINFB,NZERON,NXGRAE,NYGRAE,NZGRAE,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ISELEP(2,290) = ISELEP(2,290) + NZERON
      ENDIF
      IF(ISELEP(1,292) .GE. 1) THEN
         CALL ASDFUN  (OZFG ,  0  ,  1  ,  1  ,
     $                 ASYOZZ,SSYOZZ,KSYL,JJ1L,ILIMXP,
     $                 'S Y             ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 ILIMX,ISLINP,ISLIDI,FVT,IDIMF,NBND,
     $                 LINFB,NZERON,NXGRAE,NYGRAE,NZGRAE,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ISELEP(2,292) = ISELEP(2,292) + NZERON
      ENDIF
C
C                                  HAEUFIGKEITSVERTEILUNG DER
C                                  INKLINATIONSWINKEL DER VORTICITY-
C                                  VEKTOREN
C                                  ---------------------------------
C
      IF(ISELEP(1,306) .GE. 1) THEN
         CALL HRELOM  (OY    ,AOY   ,  1  ,  0  ,  1  ,
     $                 OZ    ,AOZ   ,  0  ,  1  ,  1  ,
     $                 AHOZOY,SHOZOY,KH0L,JJ2L,ILIMXP,
     $                 'H               ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,HVK,KLASS,
     $                 NBND)
         ISELEP(2,306) = ISELEP(2,306) + 1
      ENDIF
      IF(ISELEP(1,308) .GE. 1) THEN
         CALL HRELOM  (OX    ,AOX   ,  1  ,  1  ,  0  ,
     $                 OZ    ,AOZ   ,  0  ,  1  ,  1  ,
     $                 AHOZOX,SHOZOX,KH0L,JJ2L,ILIMXP,
     $                 'H               ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,HVK,KLASS,
     $                 NBND)
         ISELEP(2,308) = ISELEP(2,308) + 1
      ENDIF
      IF(ISELEP(1,310) .GE. 1) THEN
         CALL HRELOM  (OX    ,AOX   ,  1  ,  1  ,  0  ,
     $                 OY    ,AOY   ,  1  ,  0  ,  1  ,
     $                 AHOYOX,SHOYOX,KH0L,JJ2L,ILIMXP,
     $                 'H               ',KK,JJ,II,
     $                 KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,HVK,KLASS,
     $                 NBND)
         ISELEP(2,310) = ISELEP(2,310) + 1
      ENDIF
C
      RETURN
      END
