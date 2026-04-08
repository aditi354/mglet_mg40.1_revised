










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
      SUBROUTINE KO2VAR  (PHIFA ,APHIRA,KSTAGA,JSTAGA,ISTAGA,
     $                    PHIFB ,APHIRB,KSTAGB,JSTAGB,ISTAGB,
     $                    ARNAB ,SRNAB ,KKNL,JJ2L,JJ1L,ILIMXP,
     $                    IVAR            ,KK,JJ,II,
     $                    KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                    KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $                    NBND,LINFB,NXGRAE,NYGRAE,NZGRAE)
C*STARLET***************************************************************
C        K O 2 V A R      IN KO2VAR WERDEN RAEUMLICHE ZWEIPUNKT-
C                         KORRELATIONEN ZWISCHEN ZWEI ALLGEMEINEN GROES-
C                         SEN 'A' UND 'B' GEBILDET. DIE GROESSEN 'A' UND
C                         'B' DUERFEN BELIEBIG IM MASCHENGITTER VERSCHO-
C                         BEN SEIN, DA SIE MIT HILFE DER INTERPOLATIONS-
C                         ROUTINEN  INTST  BZW.  INTSTA  AUF EINEN
C                         (NEUEN) GEMEINSAMEN PUNKT INTERPOLIERT WERDEN.
C
C                         DIE KORRELATIONSRICHTUNG DARF EINE DER DREI
C                         KOORDINATENRICHTUNGEN X-, Y-, ODER Z  SEIN.
C                         ES KOENNEN SOWOHL KORRELATIONSFUNKTIONEN ALS
C                         AUCH KORRELATIONSKOEFFIZIENTEN BESTIMMT WERDEN
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
C        ILIMX          - TATSAECHLICH BENOETIGTE I-LINIEN WAEHREND
C                         DES MOMENTANEN LAUFES
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
C UPROG                 : ERRR,   KO2HOM,   KO2PKT
C
C DEFINE-DIREKTIVEN     : XHOMOG,  YHOMOG,  ZHOMOG
C
C        23.09.88 (HW)  : ORIGINAL
C        24.11.88 (HW)  : UEBERGABE GEAENDERT AN SRNAB (JJ2L --> JJ1L)
C
C*STARLET***************************************************************
C
C
      CHARACTER (LEN=16)  IVAR
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
      IF(IVAR(1:1) .NE. 'R')                CALL ERRR (501,' KO2VAR   ')
      IF(IVAR(2:2) .NE. 'C'  .AND.
     $   IVAR(2:2) .NE. ' ')                CALL ERRR (502,' KO2VAR   ')
      IF(IVAR(3:3) .NE. 'Z'  .AND.  IVAR(3:3) .NE. 'Y'  .AND.
     $   IVAR(3:3) .NE. 'X')                CALL ERRR (503,' KO2VAR   ')
C
      IF(IVAR(2:2) .EQ. 'C')  GOTO 2200
C
C                                 KORRELATIONSRICHTUNG FESTSTELLEN
C
      KDIR   = 0
      JDIR   = 0
      IDIR   = 0
      IF(IVAR(3:3) .EQ. 'Z') KDIR = 1
      IF(IVAR(3:3) .EQ. 'Y') JDIR = 1
      IF(IVAR(3:3) .EQ. 'X') IDIR = 1
C
C                                 DIE X-RICHTUNG IST  N I C H T  HOMOGEN
      NXRIHO = 0
C
C                                 DIE Y-RICHTUNG IST  N I C H T  HOMOGEN
      NYRIHO = 0
C
C                                 DIE Z-RICHTUNG IST  N I C H T  HOMOGEN
      NZRIHO = 0
C
      IF(KDIR .EQ. 1) THEN
C
C                                 KORRELATIONSRICHTUNG IST DIE Z-RI.
C
         IF(NZRIHO .EQ. 1  .AND.  NZGRAE .EQ. 1) THEN
            GOTO 2100
         ELSE
            GOTO 2200
         END IF
      END IF
      IF(JDIR .EQ. 1) THEN
C
C                                 KORRELATIONSRICHTUNG IST DIE Y-RI.
C
         IF(NYRIHO .EQ. 1  .AND.  NYGRAE .EQ. 1) THEN
            GOTO 2100
         ELSE
            GOTO 2200
         END IF
      END IF
      IF(IDIR .EQ. 1) THEN
C
C                                 KORRELATIONSRICHTUNG IST DIE X-RI.
C
         IF(NXRIHO .EQ. 1  .AND.  NXGRAE .EQ. 1) THEN
            GOTO 2100
         ELSE
            GOTO 2200
         END IF
      END IF
C
      CALL ERRR (510,' KO2VAR   ')
C
C                                 DIE KORRELATIONSRICHTUNG IST EINE
C                                 HOMOGENE RICHTUNG, DAS GITTER IN DIESE
C                                 RICHTUNG IST AEQUIDISTANT, ES WIRD DIE
C                                 KORRELATIONS F U N K T I O N  BESTIMMT
C
 2100 CALL KO2HOM  (PHIFA ,APHIRA,KSTAGA,JSTAGA,ISTAGA,
     $              PHIFB ,APHIRB,KSTAGB,JSTAGB,ISTAGB,
     $              ARNAB ,SRNAB ,KKNL,JJ2L,JJ1L,ILIMXP,
     $              IVAR            ,KK,JJ,II,
     $              KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $              NBND,LINFB)
      GOTO 9999
C
C                                 HIER: ES WIRD ENTWEDER DER KORRELA-
C                                 TIONSKOEFFIZIENT AUSGEWERTET
C                                            ODER
C                                 DIE KORRELATIONSFUNKTION, WOBEI :
C                                 A) DIE KORRELATIONSRICHTUNG KEINE
C                                    HOMOGENE RICHTUNG IST
C                                            ODER
C                                 B) DIE KORRELATIONSRICHTUNG EINE
C                                    HOMOGENE RICHTUNG IST, DAS GITTER
C                                    JEDOCH NICHTAEQUIDISTANT IST
C
 2200 CALL KO2PKT  (PHIFA ,APHIRA,KSTAGA,JSTAGA,ISTAGA,
     $              PHIFB ,APHIRB,KSTAGB,JSTAGB,ISTAGB,
     $              ARNAB ,SRNAB ,KKNL,JJ2L,JJ1L,ILIMXP,
     $              IVAR            ,KK,JJ,II,
     $              KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              KKA,JJA,IIA,ILIMX,ISLINP,ISLIDI,RKOMXP,
     $              NBND,LINFB)
C
 9999 RETURN
      END
