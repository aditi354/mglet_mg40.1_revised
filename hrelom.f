










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
      SUBROUTINE HRELOM  (OMH   ,AOMH  ,KSTAGH,JSTAGH,ISTAGH,
     $                    OMV   ,AOMV  ,KSTAGV,JSTAGV,ISTAGV,
     $                    AHOVH ,SHOVH ,KH0L,JJ2L,ILIMXP,
     $                    IVAR            ,KK,JJ,II,
     $                    KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                    KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,HVK,KLASS,
     $                    NBND)
C*STARLET***************************************************************
C        H R E L O M      IN HRELOM WIRD DIE RELATIVE HAEUFIGKEITSVER-
C                         TEILUNG DER INKLINATIONSWINKEL DER WIRBEL-
C                         VEKTOREN IM BEREICH
C  -180 [GRAD] <= ATAN (OMEGA_VERTIKAL / OMEGA_HORIZONTAL) <= 180 [GRAD]
C                         BESTIMMT. ES WIRD SOWOHL DIE UNGEWICHTETE ALS
C                         AUCH DIE MIT DEM VERHAELTNIS  (OMEGA_VEKTOR)
C                         **2 / (<OMEGA_VEKTOR>)**2 GEWICHTETE HAEUFIG-
C                         KEITSVERTEILUNG BERECHNET.
C*STARLET***************************************************************
C
C PARAM: OMH  (KK,JJ,II)- ENTHAELT DEN FLAECHENMITTELWERT DER
C                         H_ORIZONTALEN OMEGA-KOMPONENTE
C        AOMH (KKA,JJA, - ENTHAELT DEN FLAECHENMITTELWERT (STATIST.
C              IIA)       GEMITTELT !) DER H_ORIZ. OMEGA-KOMPONENTE
C        K-, J-, ISTAGH - KENNZEICHNET DIE VERSCHIEBUNG DER H_ORIZ.
C                         OMEGA-KOMPONENTE
C                         NSTAGH = 0: KEINE VERSCHIEBUNG IN N-RICHTUNG
C                         NSTAGH = 1: OMEGA_HORIZONTAL IST IN POSITIVER
C                                     N-RICHTUNG IM MASCHENGITTER VER-
C                                     SCHOBEN
C        OMV  (KK,JJ,II)- ENTHAELT DEN FLAECHENMITTELWERT DER
C                         V_ERTIKALEN OMEGA-KOMPONENTE
C        AOMV (KKA,JJA, - ENTHAELT DEN FLAECHENMITTELWERT (STATIST.
C              IIA)       GEMITTELT !) DER V_ERT. OMEGA-KOMPONENTE
C        K-, J-, ISTAGV - KENNZEICHNET DIE VERSCHIEBUNG DER V_ERT.
C                         OMEGA-KOMPONENTE
C                         NSTAGV = 0: KEINE VERSCHIEBUNG IN N-RICHTUNG
C                         NSTAGV = 1: OMEGA_VERTIKAL IST IN POSITIVER
C                                     N-RICHTUNG IM MASCHENGITTER VER-
C                                     SCHOBEN
C        AHOVH (KH0L,   + ENTHAELT DIE STATISTISCHEN MITTELWERTE
C        JJ2L,ILIMXP)     DER HAEUFIGKEITSVERTEILUNGEN AN DEN AUF-
C                         PUNKTEN  2 ... ILIMX
C        SHOVH (KH0L,   + SUMMATIONSFELD FUER DIE HAEUFIGKEITS-
C        JJ2L,ILIMXP)     VERTEILUNGEN AN DEN AUFPUNKTEN  2 ... ILIMX
C        KH0L,JJ2L,     - ARRAYDIMENSIONEN
C        ILIMXP         - ARRAYDIMENSIONEN
C        ILIMXP         - MAXIMAL MOEGLICHE ANZAHL VON I-LINIEN WAEHREND
C                         DES MOMENTANEN LAUFES
C        IVAR           - CHARACTER-VARIABLE
C        KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        X, Y, Z        - KOORDINATEN DER ZELLDEFINITIONSPUNKTE
C        DX,DY,DZ       - ABSTAND DER GITTERPUNKTE
C        DDX,DDY,DDZ    - KANTENLAENGEN DER KONTROLLVOLUMINA
C        KKA,JJA,IIA    - ARRAYDIMENSIONEN DER AUSWERTEFELDER
C                         BEISPIEL: FALLS J-RICHTUNG NICHT HOMOGEN IST,
C                         GILT JJA = JJ. ANDERNFALLS IST JJA = 1.
C        ILIMX          - TATSAECHLICH BENOETIGTE I-LINIEN WAEHREND
C                         DES MOMENTANEN LAUFES
C        ISLINP(ISLIDI) - INHALTSVERZEICHNIS DER KODIERTEN AUFPUNKTE
C                         DES MOMENTANEN LAUFES (P FUER PRESENT)
C        ISLIDI         - DIMENSION DER ISLIN.-FELDER
C        HVK (KLASS,    + HILFSFELD FUER DIE HAEUFIGKEITSVERTEILUNG
C             JJ2L  )
C        KLASS          - ANZAHL DER KLASSEN DER HAEUFIGKEITSVERTEILUNG
C                         = ANZAHL DER BALKEN IM BALKENDIAGRAMM
C        NBND           - ANZAHL DER RANDSCHICHTEN DES BERECHNUNGS-
C                         GEBIETES
C
C UPROG                 : ERRR,  INTST,  INTSTA,  KJIDCO,  RKJIAE
C
C DEFINE-DIREKTIVEN     : XHOMOG,  YHOMOG,  ZHOMOG
C
C        30.11.88 (HW)  : ORIGINAL
C
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
      CHARACTER (LEN=16)  IVAR
      CHARACTER (LEN=1)  CDIRC
C
      INTEGER  ISLINP (ISLIDI)
C
      REAL     OMH   (KK,JJ,II),          AOMH  (KKA,JJA,IIA),
     $         OMV   (KK,JJ,II),          AOMV  (KKA,JJA,IIA),
     $         AHOVH (KH0L,JJ2L,ILIMXP),  SHOVH (KH0L,JJ2L,ILIMXP),
     $         HVK   (KLASS,JJ2L)
C
      REAL        X(II),         Y(JJ),         Z(KK),
     $           DX(II),        DY(JJ),        DZ(KK),
     $          DDX(II),       DDY(JJ),       DDZ(KK)
C
C                                 CHECK DER CHARACTER-VARIABLEN
C
      IF(IVAR(1:1) .NE. 'H')                CALL ERRR (501,' HRELOM   ')
C
C                                 BESTIMMUNG DER VERSCHIEBUNG, WELCHE
C                                 DIE INTERPOLIERTEN GROESSEN PHIA
C                                 UND PHIB ERHALTEN
C
      KSTAG  = MAX0 (KSTAGH, KSTAGV)
      JSTAG  = MAX0 (JSTAGH, JSTAGV)
      ISTAG  = MAX0 (ISTAGH, ISTAGV)
C
      IF(ILIMX .LT. 2) CALL ERRR (510,' HRELOM   ')
      IF(KLASS .LT. 1) CALL ERRR (511,' HRELOM   ')
C
      PI     = ACOS (-1.0)
      RPI    = 1.0 / PI
      AKM    = 0.5 * FLOAT (KLASS)
      BINWID = 360.0 / FLOAT (KLASS)
      BIN0   = -(180.0 + 0.5 * BINWID)
C
C                                 DO 100: DURCHBLAETTERN DES INHALTSVER-
C                                 ZEICHNISSES DER AUFPUNKTE
C
      ILI    = 1
      DO 100 ISLIN = 2,ISLIDI
         CALL KJIDCO (ISLINP(ISLIN), CDIRC, IDIRC, KAUFP, JAUFP, IAUFP)
C
C                                 FALLS DER LETZTE EINTRAG IM INHALTS-
C                                 VERZEICHNIS GELESEN WURDE --> RETURN
C
         IF(CDIRC .EQ. '0'  .AND.  KAUFP .EQ. 0  .AND.
     $      JAUFP .EQ.  0   .AND.  IAUFP .EQ. 0)       GOTO 9999
C
C                                 NUR WENN DIE RICHTUNGSINFORMATION
C                                 IDENTISCH '0' IST, WIRD DIE HAEUFIG-
C                                 KEITSVERTEILUNG BERECHNET. ANDERN-
C                                 FALLS WIRD IM INHALTSVERZEICHNIS
C                                 WEITERGEBLAETTERT
C
         IF(CDIRC .NE. '0') GOTO 95
C
         ILI    = ILI + 1
C
C
         IF(IABS ( IFIX (AHOVH(KH0L   , 1,ILI))) .EQ. 0) THEN
C
C                                 DER INFORMATIONSBLOCK IN DER
C                                 LINIE 'ILI' WAR BISHER NOCH NICHT
C                                 BELEGT WORDEN. DAHER WIRD JETZT
C                                 DIE BELEGUNG DURCHGEFUEHRT.
C
C                                 DIE NAUFP SIND VORZEICHENBEHAFTET
C                                 (BEI HOMOGENER N-RICHTUNG IST DAS
C                                 VORZEICHEN NEGATIV)
C
            IAKAUF = IABS (KAUFP)
            IAJAUF = IABS (JAUFP)
            IAIAUF = IABS (IAUFP)
C
            RKORR  = 10.0 * SMAONE
C
            CALL RKJIAE  (KK,JJ,II,KMX,JMX,IMX,X,Y,Z,IAKAUF,IAJAUF,
     $                    IAIAUF, 0  , 0  , 1  ,NBND,RKORR,
     $                    ISLIN,KBEGDU,KENDDU,KLAGDU)
C
C                                 BELEGUNG DES INFO-BLOCKES
C
            AHOVH (KH0L   , 1,ILI) =  1.1
            AHOVH (KH0L- 1, 1,ILI) = FLOAT (KAUFP)
            AHOVH (KH0L- 2, 1,ILI) = FLOAT (JAUFP)
            AHOVH (KH0L- 3, 1,ILI) = FLOAT (IAUFP)
            AHOVH (KH0L- 4, 1,ILI) = FLOAT (IDIRC)
*           AHOVH (KH0L- 5, 1,ILI) =
            AHOVH (KH0L- 6, 1,ILI) = Z(IAKAUF) + 0.5*DZ(IAKAUF)
     $                             * FLOAT(KSTAG)
            AHOVH (KH0L- 7, 1,ILI) = Y(IAJAUF) + 0.5*DY(IAJAUF)
     $                             * FLOAT(JSTAG)
            AHOVH (KH0L- 8, 1,ILI) = X(IAIAUF) + 0.5*DX(IAIAUF)
     $                             * FLOAT(ISTAG)
*           AHOVH (KH0L- 9, 1,ILI) =
*           AHOVH (KH0L-10, 1,ILI) =
*           AHOVH (KH0L-11, 1,ILI) =
            AHOVH (KH0L-12, 1,ILI) = FLOAT (KLASS)
            AHOVH (KH0L-13, 1,ILI) = FLOAT ( 0 )
         ELSE
C
C                                 HIER: DER INFORMATIONSBLOCK IST
C                                 BEREITS RICHTIG BELEGT
C
            IAKAUF = IABS (IFIX (AHOVH(KH0L- 1, 1,ILI)))
            IAJAUF = IABS (IFIX (AHOVH(KH0L- 2, 1,ILI)))
            IAIAUF = IABS (IFIX (AHOVH(KH0L- 3, 1,ILI)))
         END IF
C
C                                 NULLBELEGUNG DES HILFSFELDES FUER
C                                 DIE HAEUFIGKEITSVERTEILUNG
C
         DO 110 JH = 1,JJ2L
            DO 110 KH = 1,KLASS
  110          HVK (KH,JH) = 0.0
C
C                                 BEACHTE: X-RICHTUNG N I C H T HOMOGEN
C
         IBEG   = IAIAUF
         IEND   = IAIAUF
C
C                                 ANZAHL DER I-SCHLEIFEN
         IANZSL = IEND - IBEG + 1
C
         DO 200 I  = IBEG,IEND
C
C                                 BEACHTE: Y-RICHTUNG N I C H T HOMOGEN
C
            JBEG   = IAJAUF
            JEND   = IAJAUF
C
C                                 ANZAHL DER J-SCHLEIFEN
            JANZSL = JEND - JBEG + 1
C
            DO 300 J  = JBEG,JEND
C
C                                 BEACHTE: Z-RICHTUNG N I C H T HOMOGEN
C
               KBEG   = IAKAUF
               KEND   = IAKAUF
C
C                                 ANZAHL DER K-SCHLEIFEN
               KANZSL = KEND - KBEG + 1
C
C                                 IANZSL*JANZSL*KANZSL IST DIE
C                                 GESAMTE ZAHL VON PUNKTEN, UEBER DIE
C                                 GEMITTELT WIRD.
C
               RANZAH = 1.0 / FLOAT(IANZSL*JANZSL*KANZSL)
C
               DO 400 K  = KBEG,KEND
C
C                                 INTERPOLATION DER GROESSEN 'A' UND 'B'
C                                 AM AUFPUNKT  K,J,I
C
                  CALL INTST  (OMH  ,KSTAGH,JSTAGH,ISTAGH,
     $                         OMV  ,KSTAGV,JSTAGV,ISTAGV,
     $                         OMHIN ,OMVIN ,K ,J ,I ,KK, JJ, II)
                  CALL INTSTA (AOMH  ,KSTAGH,JSTAGH,ISTAGH,
     $                         AOMV  ,KSTAGV,JSTAGV,ISTAGV,
     $                         AOMHIN,AOMVIN,K ,J ,I ,KKA,JJA,IIA)
C
                  WRAD   = ATAN (         ABS ( OMVIN)
     $                   /       (AMAX1 ( ABS ( OMHIN), SMAONE)) )
                  WRAD   = (-  PI * (SIGN (0.5, OMHIN)-0.5)
     $                   +   WRAD * (SIGN (1.0, OMHIN)    ))
     $                   *           SIGN (1.0, OMVIN)
                  KLOC   = IFIX (AKM * (1.0 + WRAD*RPI)) + 1
                  KLOC   = MAX0 (MIN0 (KLASS, KLOC), 1)
                  OMIN2  =  OMHIN**2 +  OMVIN**2
                  AOMIN2 = AOMHIN**2 + AOMVIN**2
C
C                                 WENN DER BETRAG DES OMEGA-VEKTORS
C                                 KLEINER ALS 10.0*SMAONE IST, WIRD
C                                 DIE STICHPROBE IGNORIERT. EBENSO
C                                 FALLS <OMEGA> < 10.0*SMAONE .
C
                  ZERONE = SIGN (0.5, (SQRT ( OMIN2) - 10.0*SMAONE))+0.5
                  AZERON = SIGN (0.5, (SQRT (AOMIN2) - 10.0*SMAONE))+0.5
C
C                                 NICHT-GEWICHTETE HAEUFIGKEITSVERT.
C
                  HVK (KLOC, 1) = HVK (KLOC, 1) + RANZAH * ZERONE
C
C                                 GEWICHTETE HAEUFIGKEITSVERT.
C
                  WEIGHT        = OMIN2 / AMAX1 (AOMIN2, SMAONE)
                  HVK (KLOC, 2) = HVK (KLOC, 2) + RANZAH * WEIGHT
     $                          * ZERONE * AZERON
  400          CONTINUE
  300       CONTINUE
  200    CONTINUE
         IF(ILI .EQ. 2) THEN
C
C                        BESTIMMUNG DER BALKENMITTELPUNKTE
C                        (NUR EINMAL NOTWENDIG !)
C
            DO 500 KL = 1,KLASS
               AHOVH (KL,1,1) = BIN0 + FLOAT (KL) * BINWID
  500          AHOVH (KL,2,1) = BIN0 + FLOAT (KL) * BINWID
         END IF
C
C                        NORMIERUNG DER HAEUFIGKEITSVERTEILUNG
C
         SUMH1  = 0.0
         SUMH2  = 0.0
         DO 510 KL = 1,KLASS
            SUMH1          = SUMH1 + HVK (KL, 1)
  510       SUMH2          = SUMH2 + HVK (KL, 2)
         RSUMH1 = 1.0 / (AMAX1 (ABS (SUMH1), SMAONE)*SIGN (1.0, SUMH1))
         RSUMH2 = 1.0 / (AMAX1 (ABS (SUMH2), SMAONE)*SIGN (1.0, SUMH2))
         DO 520 KL = 1,KLASS
            SHOVH (KL, 1,ILI) = SHOVH (KL, 1,ILI) + HVK (KL, 1) * RSUMH1
  520       SHOVH (KL, 2,ILI) = SHOVH (KL, 2,ILI) + HVK (KL, 2) * RSUMH2
         SHOVH (KH0L-13, 1,ILI) = FLOAT( IFIX (SHOVH(KH0L-13, 1,ILI))+1)
   95    CONTINUE
  100 CONTINUE
C
 9999 RETURN
      END
