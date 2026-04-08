










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
      SUBROUTINE KO2PKT  (PHIFA ,APHIRA,KSTAGA,JSTAGA,ISTAGA,
     $                    PHIFB ,APHIRB,KSTAGB,JSTAGB,ISTAGB,
     $                    ARNAB ,SRNAB ,KKNL,JJ2L,JJ1L,ILIMXP,
     $                    IVAR            ,KK,JJ,II,
     $                    KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                    KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                    NBND,LINFB)
C*STARLET***************************************************************
C        K O 2 P K T      IN KO2PKT WERDEN RAEUMLICHE ZWEIPUNKT-
C                         KORRELATIONEN ZWISCHEN ZWEI ALLGEMEINEN GROES-
C                         SEN 'A' UND 'B' GEBILDET. DIE GROESSEN 'A' UND
C                         'B' DUERFEN BELIEBIG IM MASCHENGITTER VERSCHO-
C                         BEN SEIN, DA SIE MIT HILFE DER INTERPOLATIONS-
C                         ROUTINEN  INTST  BZW.  INTSTA  AUF EINEN
C                         (NEUEN) GEMEINSAMEN PUNKT INTERPOLIERT WERDEN.
C
C                         DIE KORRELATIONSRICHTUNG DARF EINE DER DREI
C                         KOORDINATENRICHTUNGEN X-, Y-, ODER Z  SEIN.
C                         KO2PKT SETZT VORAUS, DASS :
C                         A) DIE KORRELATIONSRICHTUNG  N I C H T  MIT
C                            EINER HOMOGENEN RICHTUNG DES STROEMUNGS-
C                            FELDES UEBEREINSTIMMT (KORKOE = 0)
C                                         ODER
C                         B) DER KORRELATIONSKOEFFIZIENT AUSGEWERTET
C                            WERDEN SOLL           (KORKOE = 1)
C*STARLET***************************************************************
C
C PARAM: PHIFA(KK,JJ,II)- ENTHAELT DIE FLUKTUATIONEN DER GROESSE PHIA
C        APHIRA(KKA,JJA,- ENTHAELT DIE STATISTISCHEN MITTELWERTE
C               IIA)      <PHIRMS> DER ROOT-MEAN-SQUARE-WERTE DER
C                         GROESSE 'A'
C        K-, J-, ISTAGA - KENNZEICHNET DIE VERSCHIEBUNG DER GROESSE 'A'
C                         NSTAGA = 0: KEINE VERSCHIEBUNG IN N-RICHTUNG
C                         NSTAGA = 1: DIE GROESSE 'A' IST IN POSITIVER
C                                     N-RICHTUNG IM MASCHENGITTER VER-
C                                     SCHOBEN
C        PHIFB(KK,JJ,II)- ENTHAELT DIE FLUKTUATIONEN DER GROESSE PHIB
C        APHIRB(KKA,JJA,- ENTHAELT DIE STATISTISCHEN MITTELWERTE
C               IIA)      <PHIRMS> DER ROOT-MEAN-SQUARE-WERTE DER
C                         GROESSE 'B'
C        K-, J-, ISTAGB - KENNZEICHNET DIE VERSCHIEBUNG DER GROESSE 'B'
C                         NSTAGB = 0: KEINE VERSCHIEBUNG IN N-RICHTUNG
C                         NSTAGB = 1: DIE GROESSE 'B' IST IN POSITIVER
C                                     N-RICHTUNG IM MASCHENGITTER VER-
C                                     SCHOBEN
C        ARNAB (KKNL,   + ENTHAELT DIE STATISTISCHEN MITTELWERTE
C        JJ2L,ILIMXP)     DER KORRELATIONSFUNKTIONEN AN DEN AUF-
C                         PUNKTEN  2 ... ILIMX
C        SRNAB (KKNL,   + SUMMATIONSFELD FUER DIE KORRELATIONS-
C        JJ1L,ILIMXP)     FUNKTIONEN AN DEN AUFPUNKTEN  2 ... ILIMX
C        KKNL,JJ2L,     - ARRAYDIMENSIONEN
C        JJ1L,ILIMXP    - ARRAYDIMENSIONEN
C        ILIMXP         - MAXIMAL MOEGLICHE ANZAHL VON I-LINIEN WAEHREND
C                         DES MOMENTANEN LAUFES
C        ILIMX          - TATSAECHLICH BENOETIGTE I-LINIEN WAEHREND
C                         DES MOMENTANEN LAUFES
C        IVAR           - CHARACTER-VARIABLE (DIENT HAUPTSAECHLICH ZUR
C                         IDENTIFIKATION DER KORRELATIONSRICHTUNG)
C        KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        X, Y, Z        - KOORDINATEN DER ZELLDEFINITIONSPUNKTE
C        DX,DY,DZ       - ABSTAND DER GITTERPUNKTE
C        DDX,DDY,DDZ    - KANTENLAENGEN DER KONTROLLVOLUMINA
C        KKA,JJA,IIA    - ARRAYDIMENSIONEN DER AUSWERTEFELDER
C                         BEISPIEL: FALLS J-RICHTUNG NICHT HOMOGEN IST,
C                         GILT JJA = JJ. ANDERNFALLS IST JJA = 1.
C        ISLIDI         - DIMENSION DER ISLIN.-FELDER
C        RKOMXP         - MAXIMAL ZULAESSIGE KORRELATIONSLAENGE WAEHREND
C                         DES MOMENTANEN LAUFES
C        ISLINP(ISLIDI) - INHALTSVERZEICHNIS DER KODIERTEN AUFPUNKTE
C                         DES MOMENTANEN LAUFES (P FUER PRESENT)
C        NBND           - ANZAHL DER RANDSCHICHTEN DES BERECHNUNGS-
C                         GEBIETES
C        LINFB          - ANZAHL DER INFORMATIONSBLOECKE (-FELDER)
C                         JEDER "LINIE"
C
C UPROG                 : ERRR,  INTST,  INTSTA,  KJIDCO,  RKJIAE
C
C DEFINE-DIREKTIVEN     : XHOMOG,  YHOMOG,  ZHOMOG
C
C        21.09.88 (HW)  : ORIGINAL
C        07.10.88 (HW)  : AENDERUNGEN
C        24.11.88 (HW)  : UEBERGABE GEAENDERT AN SRNAB (JJ2L --> JJ1L)
C        02.01.89 (HW)  : BEI KORRELATIONSKOEFFIZIENTEN WIRD BEI
C                         KJIBEG  DIE VERSCHIEBUNG IM MASCHENGITTER
C                         BERUECKSICHTIGT
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
      REAL     PHIFA (KK,JJ,II),          APHIRA(KKA,JJA,IIA),
     $         PHIFB (KK,JJ,II),          APHIRB(KKA,JJA,IIA),
     $         ARNAB (KKNL,JJ2L,ILIMXP),  SRNAB (KKNL,JJ1L,ILIMXP)
C
      REAL        X(II),         Y(JJ),         Z(KK),
     $           DX(II),        DY(JJ),        DZ(KK),
     $          DDX(II),       DDY(JJ),       DDZ(KK)
C
C                                 CHECK DER CHARACTER-VARIABLEN
C
      IF(IVAR(1:1) .NE. 'R')                CALL ERRR (501,' KO2PKT   ')
      IF(IVAR(2:2) .NE. 'C'  .AND.
     $   IVAR(2:2) .NE. ' ')                CALL ERRR (502,' KO2PKT   ')
      IF(RKOMXP    .LE. 0.0)                CALL ERRR (503,' KO2PKT   ')
C
C                                 WICHTIGE ENTSCHEIDUNG: WIRD DIE
C                                 KORRELATIONSFUNKTION ODER DER
C                                 KORRELATIONSKOEFFIZIENT AUSGEWERTET ?
C
      IF(IVAR(2:2) .EQ. 'C') THEN
C                                 KORRELATIONSKOEFFIZIENT WIRD AUSGEW.
         KORKOE = 1
         RKORR  = GREAT
      ELSE
C                                 KORRELATIONSFUNKTION WIRD AUSGEW.
         KORKOE = 0
         RKORR  = RKOMXP
      ENDIF
C
C                                 KORRELATIONSRICHTUNG FESTSTELLEN
C
      KDIR   = 0
      JDIR   = 0
      IDIR   = 0
      IF(IVAR(3:3) .NE. 'Z'  .AND.  IVAR(3:3) .NE. 'Y'  .AND.
     $   IVAR(3:3) .NE. 'X')                CALL ERRR (504,' KO2PKT   ')
      IF(IVAR(3:3) .EQ. 'Z') KDIR = 1
      IF(IVAR(3:3) .EQ. 'Y') JDIR = 1
      IF(IVAR(3:3) .EQ. 'X') IDIR = 1
C
C                                 BESTIMMUNG DER VERSCHIEBUNG, WELCHE
C                                 DIE INTERPOLIERTEN GROESSEN PHIA
C                                 UND PHIB ERHALTEN
C
      KSTAG  = MAX0 (KSTAGA, KSTAGB)
      JSTAG  = MAX0 (JSTAGA, JSTAGB)
      ISTAG  = MAX0 (ISTAGA, ISTAGB)
C
C                                 DO 100: DURCHBLAETTERN DES INHALTSVER-
C                                 ZEICHNISSES DER AUFPUNKTE
C
      IF(ILIMX .LT. 2) CALL ERRR (510,' KO2PKT   ')
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
C                                 NUR WENN DIE KORRELATIONSRICHTUNG
C                                 MIT DER DES GERADE BETRACHTETEN
C                                 AUFPUNKTES UEBEREINSTIMMT, WIRD
C                                 DIE KORRELATION BERECHNET. ANDERN-
C                                 FALLS WIRD IM INHALTSVERZEICHNIS
C                                 WEITERGEBLAETTERT
C
         IF(IVAR(3:3) .NE. CDIRC) GOTO 95
C
         ILI    = ILI + 1
C
C
         IF(IABS ( IFIX (ARNAB(KKNL   , 1,ILI))) .EQ. 0) THEN
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
            CALL RKJIAE  (KK,JJ,II,KMX,JMX,IMX,X,Y,Z,IAKAUF,IAJAUF,
     $                    IAIAUF,KDIR,JDIR,IDIR,NBND,RKORR ,
     $                    ISLIN,KJIBEG,KJIEND,KJILAG)
C
            IF(IVAR(1:2) .EQ. 'RC') THEN
               KJISTA = KSTAG*KDIR + JSTAG*JDIR + ISTAG*IDIR
               IF(KJISTA .NE. 0  .AND.  KJISTA .NE. 1) THEN
                  CALL ERRR (515,' KO2PKT   ')
               END IF
               KJIBEG = KJIBEG - KJISTA
            END IF
C
C                                 BELEGUNG DES INFO-BLOCKES
C
            ARNAB (KKNL   , 1,ILI) = -1.1
            ARNAB (KKNL- 1, 1,ILI) = FLOAT (KAUFP)
            ARNAB (KKNL- 2, 1,ILI) = FLOAT (JAUFP)
            ARNAB (KKNL- 3, 1,ILI) = FLOAT (IAUFP)
            ARNAB (KKNL- 4, 1,ILI) = FLOAT (IDIRC)
*           ARNAB (KKNL- 5, 1,ILI) =
            ARNAB (KKNL- 6, 1,ILI) = Z(IAKAUF) + 0.5*DZ(IAKAUF)
     $                             * FLOAT(KSTAG)
            ARNAB (KKNL- 7, 1,ILI) = Y(IAJAUF) + 0.5*DY(IAJAUF)
     $                             * FLOAT(JSTAG)
            ARNAB (KKNL- 8, 1,ILI) = X(IAIAUF) + 0.5*DX(IAIAUF)
     $                             * FLOAT(ISTAG)
            ARNAB (KKNL- 9, 1,ILI) = FLOAT (KJIEND)
            ARNAB (KKNL-10, 1,ILI) = FLOAT (KJIBEG)
            ARNAB (KKNL-11, 1,ILI) = FLOAT (KJILAG)
            ARNAB (KKNL-12, 1,ILI) = FLOAT (KJIEND - KJIBEG + 1)
C
C                                 HIER: ABPRUEFEN, OB DIE FELDDIMEN-
C                                 SIONIERUNG IN K-RICHTUNG AUSREICHT
C
            IF( (IFIX ( ARNAB (KKNL-12, 1,ILI)) + LINFB) .GT. KKNL)
     $         CALL ERRR (520,' KO2PKT   ')
            ARNAB (KKNL-13, 1,ILI) = FLOAT ( 0 )
         ELSE
C
C                                 HIER: DER INFORMATIONSBLOCK IST
C                                 BEREITS RICHTIG BELEGT
C
            IAKAUF = IABS (IFIX (ARNAB(KKNL- 1, 1,ILI)))
            IAJAUF = IABS (IFIX (ARNAB(KKNL- 2, 1,ILI)))
            IAIAUF = IABS (IFIX (ARNAB(KKNL- 3, 1,ILI)))
            KJIEND =       IFIX (ARNAB(KKNL- 9, 1,ILI))
            KJIBEG =       IFIX (ARNAB(KKNL-10, 1,ILI))
            KJILAG =       IFIX (ARNAB(KKNL-11, 1,ILI))
         END IF
C
C                                 LAGE DES AUFPUNKTES
C
         KJIAUF = IAKAUF*KDIR + IAJAUF*JDIR + IAIAUF*IDIR
C
         IBEG   = IAIAUF
         IEND   = IAIAUF
C
C                                 ANZAHL DER I-SCHLEIFEN
         IANZSL = IEND - IBEG + 1
C
         DO 200 I  = IBEG,IEND
            JBEG   = IAJAUF
            JEND   = IAJAUF
C
C                                 ANZAHL DER J-SCHLEIFEN
            JANZSL = JEND - JBEG + 1
C
            DO 300 J  = JBEG,JEND
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
C                                 LAUFSCHLEIFE FUER DEN AUFPUNKT
C
                  KJIAPB = KJIAUF * (1-KORKOE) + KJIBEG * KORKOE
                  KJIAPE = KJIAUF * (1-KORKOE) + KJIEND * KORKOE
C
                  KSTO   = 0
                  DO 500 KJIAP  = KJIAPB,KJIAPE
                     KJIDIB = (KJIBEG-KJIAUF) * (1-KORKOE) + 0 * KORKOE
                     KJIDIE = (KJIEND-KJIAUF) * (1-KORKOE) + 0 * KORKOE
C
C                                 LAUFSCHLEIFE FUER DIE BESTIMMUNG
C                                 DES ABSTANDES: AUFPUNKT D. GROESSE "A"
C                                 <---> AUFPUNKT DER GROESSE "B"
C
                     DO 600 KJIDI  = KJIDIB,KJIDIE
                        KA     = K * (1-KDIR) +  KJIAP        * KDIR
                        JA     = J * (1-JDIR) +  KJIAP        * JDIR
                        IA     = I * (1-IDIR) +  KJIAP        * IDIR
                        KB     = K * (1-KDIR) + (KJIAP+KJIDI) * KDIR
                        JB     = J * (1-JDIR) + (KJIAP+KJIDI) * JDIR
                        IB     = I * (1-IDIR) + (KJIAP+KJIDI) * IDIR
C
C                                 INTERPOLATION DER GROESSEN 'A' UND 'B'
C                                 AM AUFPUNKT  "NA" (N = K,J,I)
C
                        CALL INTST  (PHIFA,KSTAGA,JSTAGA,ISTAGA,
     $                               PHIFB,KSTAGB,JSTAGB,ISTAGB,
     $                               PHIA0 ,PHIB0 ,KA,JA,IA,KK, JJ, II)
                        CALL INTSTA (APHIRA,KSTAGA,JSTAGA,ISTAGA,
     $                               APHIRB,KSTAGB,JSTAGB,ISTAGB,
     $                               PHIRA0,PHIRB0,KA,JA,IA,KKA,JJA,IIA)
C
C                                 INTERPOLATION DER GROESSEN 'A' UND 'B'
C                                 AM AUFPUNKT  "NB" (N = K,J,I)
C
                        CALL INTST  (PHIFA,KSTAGA,JSTAGA,ISTAGA,
     $                               PHIFB,KSTAGB,JSTAGB,ISTAGB,
     $                               PHIA1 ,PHIB1 ,KB,JB,IB,KK, JJ, II)
                        CALL INTSTA (APHIRA,KSTAGA,JSTAGA,ISTAGA,
     $                               APHIRB,KSTAGB,JSTAGB,ISTAGB,
     $                               PHIRA1,PHIRB1,KB,JB,IB,KKA,JJA,IIA)
C
                        ANENN  = PHIRA0 * PHIRB1
                        ANENN  = AMAX1 (ABS (ANENN), SMAONE)
     $                         * SIGN  (1.0, ANENN)
                        AKORR  = (PHIA0 * PHIB1 / ANENN) * RANZAH
                        RADKOR = (((Z(KB)+0.5*DZ(KB)*FLOAT(KSTAG))
     $                         -   (Z(KA)+0.5*DZ(KA)*FLOAT(KSTAG))
     $                         *    FLOAT(1-KORKOE)       )*FLOAT(KDIR)
     $                         +  ((Y(JB)+0.5*DY(JB)*FLOAT(JSTAG))
     $                         -   (Y(JA)+0.5*DY(JA)*FLOAT(JSTAG))
     $                         *    FLOAT(1-KORKOE)       )*FLOAT(JDIR)
     $                         +  ((X(IB)+0.5*DX(IB)*FLOAT(ISTAG))
     $                         -   (X(IA)+0.5*DX(IA)*FLOAT(ISTAG))
     $                         *    FLOAT(1-KORKOE)       )*FLOAT(IDIR))
                        KSTO   = KSTO + 1
                        SRNAB (KSTO   , 1,ILI) = SRNAB (KSTO   , 1,ILI)
     $                                         + AKORR
                        ARNAB (KSTO   , 2,ILI) = RADKOR
  600                CONTINUE
  500             CONTINUE
  400          CONTINUE
  300       CONTINUE
  200    CONTINUE
         SRNAB (KKNL-13, 1,ILI) = FLOAT( IFIX (SRNAB(KKNL-13, 1,ILI))+1)
   95    CONTINUE
  100 CONTINUE
C
 9999 RETURN
      END
