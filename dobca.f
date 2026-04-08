










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
      SUBROUTINE DOBCA  (KANAL,MODUS,IIDENT,ISTALM,ISTPLM,JSTALM,JSTPLM,
     $                   ITTOT,ITFLUC,DT,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                   ILIMX ,IDOBD1,IDOBD2,IDOBD3,IDOBD4,RKOMXP,
     $                   RDOBD1,RDOBD2,RDOBD3,RDOBD4,IGRID,IDIM3D,
     $                   XHOMOG,YHOMOG,ZHOMOG,IDIMA,

     $ KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,KMXA,JMXA,IMXA,NBND,
     $ ISELEA,ISELEP,
     $ ISAMPA,ISAMPP,ISLINP,ISLINA,ISLIDI,
     $ RIDENT,HILF,B,U,UO,AU,SU,UFG,AURG,SURG,AUSKG,SUSKG,AUFLG,SUFLG,V,
     $ VO,AV,SV,VFG,AVRG,SVRG,AVSKG,SVSKG,AVFLG,SVFLG,W,WO,AW,SW,WFG,
     $ AWRG,SWRG,AWSKG,SWSKG,AWFLG,SWFLG,P,AP,SP,PFG,APRG,SPRG,
     $ APSKG,SPSKG,APFLG,SPFLG,
     $ G,AEFG,SEFG,AEFS,SEFS,ADFG,SDFG,
     $ ADUDX2,SDUDX2,ADUDY2,SDUDY2,ADUDZ2,SDUDZ2,ADVDX2,SDVDX2,ADVDY2,
     $ SDVDY2,ADVDZ2,SDVDZ2,ADWDX2,SDWDX2,ADWDY2,SDWDY2,ADWDZ2,SDWDZ2,
     $ AUFWFG,SUFWFG,AUFWFS,SUFWFS,AUFWFM,SUFWFM,AVFWFG,SVFWFG,AVFWFS,
     $ SVFWFS,AVFWFM,SVFWFM,AUFVFG,SUFVFG,AUFVFS,SUFVFS,AUFVFM,
     $ SUFVFM,OX,AOX,SOX,OXFG,AOXRG,SOXRG,OY,AOY,SOY,OYFG,
     $ AOYRG,SOYRG,OZ,AOZ,SOZ,OZFG,AOZRG,SOZRG,AO2,SO2,O2FG,AO2RG,
     $ SO2RG,AHE,SHE,HEFG,AHERG,SHERG
     $ ,BP,BU,BV,BW
     $ ,HILF3D1,HILF3D2,HILF3D3
     $ ,AUUM,SUUM,AVVM,SVVM,AWWM,SWWM,APPM,SPPM
     $ ,AUVM,SUVM,AUWM,SUWM,AVWM,SVWM
     $ ,AUXUXM,SUXUXM,AUYUYM,SUYUYM,AUZUZM,SUZUZM,AVXVXM,SVXVXM,AVYVYM
     $ ,SVYVYM,AVZVZM,SVZVZM,AWXWXM,SWXWXM,AWYWYM,SWYWYM,AWZWZM,SWZWZM
     $ )
C*STARLET***************************************************************
C        D O B C A        AUSGABE VON FELDERN (WAHLWEISE FORMATIERT ODER
C                         UNFORMATIERT S. VAR. 'MODUS'),DIE UEBERWIEGEND
C                         ENSEMBLE-MITTELWERTE ENTHALTEN AUF KANAL
C                         'KANAL'.
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
C        ILIMX                    - ANZAHL DER AUFPUNKTE + 1 DES MOMEN-
C                                   TANEN LAUFES
C        IDOBD1 ... IDOBD4        - NOCH UNBELEGTE DUMMY-VARIABLE,
C                                   WELCHE VERWENDET WERDEN KOENNEN,
C                                   UM WICHTIGE INFORMATIONEN DES
C                                   MOMENTANEN LAUFES HERAUSSCHREIBEN
C                                   ZU KOENNEN (INTEGER-VARIABLE !)
C        RKOMXP                   - MAXIMAL ZULAESSIGER KORRELATIONS-
C                                   RADIUS WAEHREND DES MOMENTANEN
C                                   LAUFES
C        RDOBD1 ... RDOBD4        - WIE IDOBD. (REAL-VARIABLE !)
C        KK ,JJ ,II                 ARRAYDIMENSIONEN
C        KKA,JJA,IIA                ARRAYDIMENSIONEN DER AUSWERTEFELDER
C        KMX,JMX,IMX                GRENZEN DES BERECHNUNGSGEBIETES
C                                   MIT BOUND.
C        KMXA,JMXA,IMXA             GRENZEN DER AUSWERTEFELDER
C        ISELEA(2,752)              STEUERFELD F. DIE STAT. AUSWERTUNG
C                                   (ENTHAELT DIE STEUERDATEN DES VORAN-
C                                   GEGANGENEN LAUFES)
C        ISELEP(2,752)              STEUERFELD F. DIE STAT. AUSWERTUNG
C                                   (ENTHAELT DIE STEUERDATEN DES MOMEN-
C                                   TANEN LAUFES)
C        ISLINP(ISLIDI)             ENTHAELT DIE AUFPUNKTE IN KODIERTER
C                                   FORM FUER DEN MOMENTANEN LAUF
C        ISLINA(ISLIDI)             ENTHAELT DIE AUFPUNKTE IN KODIERTER
C                                   FORM FUER DEN VORANGEGANGENEN LAUF
C        ISLIDI                     ARRAYDIMENSION
C        RIDENT (100)               IDENT-FELD FUER REAL-KONSTANTEN
C        HILF   (KK ,JJ ,II )       ALLGEMEIN VERWENDBARES HILFSFELD
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
C        APHISKG(KKA,JJA,IIA)       ENSEMBLE-MITTELWERT DER SCHIEFE DER
C                                   PHI-FLUKTUATIONEN
C        SPHISKG(KKA,JJA,IIA)       SUMMATIONSFELD FUER DIE MOMENTANEN
C                                   WERTE DER SCHIEFE
C        APHIFLG(KKA,JJA,IIA)       ENSEMBLE-MITTELWERT DES FLACHHEITS-
C                                   GRADES DER PHI-FLUKTUATIONEN
C        SPHIFLG(KKA,JJA,IIA)       SUMMATIONSFELD FUER DIE MOMENTANEN
C                                   WERTE DES FLACHHEITSGRADES
C        AEFG   (KKA,JJA,IIA)       ENSEMBLE-MITTELWERT DER KINETISCHEN
C                                   ENERGIE DER SCHWANKUNGSGESCHWINDIG-
C                                   KEITEN (ANTEIL DER GROBSTRUKTUR)
C        SEFG   (KKA,JJA,IIA)       SUMMATIONSFELD F. MOMENTANWERTE EFG
C        AEFS   (KKA,JJA,IIA)       ENSEMBLE-MITTELWERT DER KIN. EN. DER
C                                   SCHWANKUNGSGESCHW. (FEINSTRUKTURA.)
C        SEFS   (KKA,JJA,IIA)       SUMMATIONSFELD F. MOMENTANWERTE EFS
C        AUIUJFG(KKA,JJA,IIA)       ENSEMBLE-MITTELWERT DER SCHUBSPANNUN
C                                   (ANTEIL DER GROBSTRUKTUR)
C        SUIUJFG(KKA,JJA,IIA)       SUMMATIONSFELD F. MOMENTANWERTE
C                                   UIUJFG
C        AUIUJFS(KKA,JJA,IIA)       ENSEMBLE-MITTELWERT DER SCHUBSPANNUN
C                                   (ANTEIL DER FEINSTRUKTUR)
C        SUIUJFS(KKA,JJA,IIA)       SUMMATIONSFELD F. MOMENTANWERTE
C        AUIUJFM(KKA,JJA,IIA)       ENSEMBLE-MITTELWERT DER SCHUBSPANNUN
C                                   (MOLEKULARER ANTEIL)
C        SUIUJFM(KKA,JJA,IIA)       SUMMATIONSFELD F. MOMENTANWERTE
C        OI(KK ,JJ ,II )            I-KOMPONENTE DER VORTICITY
C                                   (MOMENTANWERT)
C        AOI(KKA,JJA,IIA)           ENSEMBLE-MITTELWERT DER I-KOMP.
C                                   DER VORTICITY
C        SOI(KKA,JJA,IIA)           SUMMATIONSFELD FUER MOMENTANWERTE
C                                   DER VORTICITY
C        OIFG(KK ,JJ ,II )          FLUKTUATIONEN DER I-KOMP.
C                                   DER VORTICITY
C        AOIRG(KKA,JJA,IIA)         ENSEMBLE-MITTELWERT DER RMS-WERTE
C                                   DER VORTICITY (I-KOMP.)
C        SOIRG(KKA,JJA,IIA)         SUMMATIONSFELD FUER DIE MOMENTANEN
C                                   RMS-WERTE DER VORTICITY (I-KOMP.)
C
C        ...  O2  ...               ENSTROPHIE
C
C        ...  HE  ...               HELIZITAET
C
C
C UPROG                           : DOBC, ERRR
C
C DEFINE-DIREKTIVEN               : XHOMOG, YHOMOG (IN COMMON-DECK
C                                   CDOBCA)
C
C        12.08.86 (HW)  : ORIGINAL
C        29.08.86 (HW)  : BEI DEN KOMPONENTEN DES SCHUBSPANNUNGSTENSORS
C                         WIRD DIE SUMME AUS GROB- UND FEINSTRUKTUR-
C                         ANTEIL (REYNOLDSSPANNUNG) AUSGEGEBEN
C        30.09.86 (HW)  : DIE MOMENTANWERTE DER GESCHW.KOMP. UND DER
C                         FLUKTUATIONEN WERDEN AUSGEGEBEN (WEGEN
C                         GRAPHISCHER DARSTELLUNG)
C        09.12.86 (HW)  : ES KANN WAHLWEISE FORMATIERT ODER UNFORMA-
C                         TIERT GESCHRIEBEN WERDEN. (VARIABLE 'MODUS'
C                         EINGEFUEHRT)
C        17.12.86 (HW)  : WEGEN NOS/VE MUSSTE DIE FUNCTION CH1680
C                         EINGEFUEHRT WERDEN
C        04.07.88 (HW)  : ERWEITERUNG: AUSWERTUNG DER DREI KOMPONENTEN
C                         DER VORTICITY, DER ENSTROPHIE UND DER
C                         HELIZITAET
C        19.08.88 (HW)  : VORTICITY IST JETZT IM DOPPELT VERSETZTEN
C                         MASCHENGITTER DEFINIERT
C        08.12.88 (HW)  : AENDERUNG DER UEBERGABE. DAS ERWEITERTE ISLEA-
C                         FELD UND DAS ISLINA-FELD WERDEN GESCHRIEBEN
C        09.12.88 (HW)  : CDOBCA EINGEFUEHRT
C        03.01.89 (HW)  : DIE HELIZITAET IST JETZT DREIFACH IM MASCHEN-
C                         GITTER VERSCHOBEN (IN X-, Y- UND Z-RICHTUNG)
C         6. 4.92 (MM)  : IIDENT EINGEFUEHRT 
C        26. 5.95 (MM)  : UEBERGABE NACH DOBC GEAENDERT
C         6.03.01 (SE)  : Neue Statistik eingebaut
C        11.02.03 (TB)  : SCALAR FIELD STATISTICS ADDED
C
C*STARLET***************************************************************
C
      CHARACTER (LEN=8)   MODUS
      CHARACTER (LEN=80)  CT, CH1680
      INTEGER  IIDENT(100)
C
      REAL     X  (II),      Y  (JJ),      Z  (KK),
     $         DX (II),      DY (JJ),      DZ (KK),
     $         DDX(II),      DDY(JJ),      DDZ(KK)

      REAL     AU(IDIMA),AV(IDIMA),AW(IDIMA),HILF3D1(IDIMA),B(IDIM3D)
      LOGICAL XHOMOG,YHOMOG,ZHOMOG


      INTEGER  ISELEA (2,752)      , ISELEP (2,752),
     $         ISAMPA (752)        , ISAMPP (752),
     $         ISLINP (ISLIDI)     , ISLINA (ISLIDI)
      REAL     RIDENT (100)
      REAL     HILF   (KK ,JJ ,II )
C
      REAL     U      (KK ,JJ ,II ), UO     (KK ,JJ ,II ),
     $         UFG    (KK ,JJ ,II ),
     $         AURG   (KKA,JJA,IIA), SURG   (KKA,JJA,IIA),
     $         AUSKG  (KKA,JJA,IIA), SUSKG  (KKA,JJA,IIA),
     $         AUFLG  (KKA,JJA,IIA), SUFLG  (KKA,JJA,IIA)
      REAL     V      (KK ,JJ ,II ), VO     (KK ,JJ ,II ),
     $         VFG    (KK ,JJ ,II ),
     $         AVRG   (KKA,JJA,IIA), SVRG   (KKA,JJA,IIA),
     $         AVSKG  (KKA,JJA,IIA), SVSKG  (KKA,JJA,IIA),
     $         AVFLG  (KKA,JJA,IIA), SVFLG  (KKA,JJA,IIA)
      REAL     W      (KK ,JJ ,II ), WO     (KK ,JJ ,II ),
     $         WFG    (KK ,JJ ,II ),
     $         AWRG   (KKA,JJA,IIA), SWRG   (KKA,JJA,IIA),
     $         AWSKG  (KKA,JJA,IIA), SWSKG  (KKA,JJA,IIA),
     $         AWFLG  (KKA,JJA,IIA), SWFLG  (KKA,JJA,IIA)
      REAL     P      (KK ,JJ ,II ),
     $         AP     (KKA,JJA,IIA), SP     (KKA,JJA,IIA),
     $         PFG    (KK ,JJ ,II ),
     $         APRG   (KKA,JJA,IIA), SPRG   (KKA,JJA,IIA),
     $         APSKG  (KKA,JJA,IIA), SPSKG  (KKA,JJA,IIA),
     $         APFLG  (KKA,JJA,IIA), SPFLG  (KKA,JJA,IIA)
      REAL     G      (KK ,JJ ,II ),
     $         AEFG   (KKA,JJA,IIA), SEFG   (KKA,JJA,IIA),
     $         AEFS   (KKA,JJA,IIA), SEFS   (KKA,JJA,IIA),
     $         ADFG   (KKA,JJA,IIA), SDFG   (KKA,JJA,IIA),
     $         ADUDX2 (KKA,JJA,IIA), SDUDX2 (KKA,JJA,IIA),
     $         ADUDY2 (KKA,JJA,IIA), SDUDY2 (KKA,JJA,IIA),
     $         ADUDZ2 (KKA,JJA,IIA), SDUDZ2 (KKA,JJA,IIA),
     $         ADVDX2 (KKA,JJA,IIA), SDVDX2 (KKA,JJA,IIA),
     $         ADVDY2 (KKA,JJA,IIA), SDVDY2 (KKA,JJA,IIA),
     $         ADVDZ2 (KKA,JJA,IIA), SDVDZ2 (KKA,JJA,IIA),
     $         ADWDX2 (KKA,JJA,IIA), SDWDX2 (KKA,JJA,IIA),
     $         ADWDY2 (KKA,JJA,IIA), SDWDY2 (KKA,JJA,IIA),
     $         ADWDZ2 (KKA,JJA,IIA), SDWDZ2 (KKA,JJA,IIA),
     $         AUFWFG (KKA,JJA,IIA), SUFWFG (KKA,JJA,IIA),
     $         AUFWFS (KKA,JJA,IIA), SUFWFS (KKA,JJA,IIA),
     $         AUFWFM (KKA,JJA,IIA), SUFWFM (KKA,JJA,IIA)
      REAL     AVFWFG (KKA,JJA,IIA), SVFWFG (KKA,JJA,IIA),
     $         AVFWFS (KKA,JJA,IIA), SVFWFS (KKA,JJA,IIA),
     $         AVFWFM (KKA,JJA,IIA), SVFWFM (KKA,JJA,IIA),
     $         AUFVFG (KKA,JJA,IIA), SUFVFG (KKA,JJA,IIA),
     $         AUFVFS (KKA,JJA,IIA), SUFVFS (KKA,JJA,IIA),
     $         AUFVFM (KKA,JJA,IIA), SUFVFM (KKA,JJA,IIA)
      REAL     OX     (KK ,JJ ,II ), OXFG   (KK ,JJ ,II ),
     $         AOX    (KKA,JJA,IIA), SOX    (KKA,JJA,IIA),
     $         AOXRG  (KKA,JJA,IIA), SOXRG  (KKA,JJA,IIA)
      REAL     OY     (KK ,JJ ,II ), OYFG   (KK ,JJ ,II ),
     $         AOY    (KKA,JJA,IIA), SOY    (KKA,JJA,IIA),
     $         AOYRG  (KKA,JJA,IIA), SOYRG  (KKA,JJA,IIA)
      REAL     OZ     (KK ,JJ ,II ), OZFG   (KK ,JJ ,II ),
     $         AOZ    (KKA,JJA,IIA), SOZ    (KKA,JJA,IIA),
     $         AOZRG  (KKA,JJA,IIA), SOZRG  (KKA,JJA,IIA)
      REAL                           O2FG   (KK ,JJ ,II ),
     $         AO2    (KKA,JJA,IIA), SO2    (KKA,JJA,IIA),
     $         AO2RG  (KKA,JJA,IIA), SO2RG  (KKA,JJA,IIA)
      REAL                           HEFG   (KK ,JJ ,II ),
     $         AHE    (KKA,JJA,IIA), SHE    (KKA,JJA,IIA),
     $         AHERG  (KKA,JJA,IIA), SHERG  (KKA,JJA,IIA)
      REAL     BP(KK,JJ,II),BU(KK,JJ,II),BV(KK,JJ,II),BW(KK,JJ,II)
      REAL     HILF3D2(KK ,JJ ,II ),
     $         HILF3D3(KK ,JJ ,II )
      REAL       AUUM(KKA,JJA,IIA),    SUUM(KKA,JJA,IIA),
     $           AVVM(KKA,JJA,IIA),    SVVM(KKA,JJA,IIA),
     $           AWWM(KKA,JJA,IIA),    SWWM(KKA,JJA,IIA),
     $           APPM(KKA,JJA,IIA),    SPPM(KKA,JJA,IIA)
      REAL       AUVM(KKA,JJA,IIA),    SUVM(KKA,JJA,IIA),
     $           AUWM(KKA,JJA,IIA),    SUWM(KKA,JJA,IIA),
     $           AVWM(KKA,JJA,IIA),    SVWM(KKA,JJA,IIA)
      REAL     AUXUXM(KKA,JJA,IIA),  SUXUXM(KKA,JJA,IIA),
     $         AUYUYM(KKA,JJA,IIA),  SUYUYM(KKA,JJA,IIA),
     $         AUZUZM(KKA,JJA,IIA),  SUZUZM(KKA,JJA,IIA),
     $         AVXVXM(KKA,JJA,IIA),  SVXVXM(KKA,JJA,IIA),
     $         AVYVYM(KKA,JJA,IIA),  SVYVYM(KKA,JJA,IIA),
     $         AVZVZM(KKA,JJA,IIA),  SVZVZM(KKA,JJA,IIA),
     $         AWXWXM(KKA,JJA,IIA),  SWXWXM(KKA,JJA,IIA),
     $         AWYWYM(KKA,JJA,IIA),  SWYWYM(KKA,JJA,IIA),
     $         AWZWZM(KKA,JJA,IIA),  SWZWZM(KKA,JJA,IIA)


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
C
C
C

C 2000 CONTINUE
C
C                                 AUSGABE DER FELDER, DIE ENSEMBLE-
C                                 MITTELWERTE ENTHALTEN. DIE REIHEN-
C                                 FOLGE ENTSPRICHT DER IM ISELEP-FELD
C
C                                 STATIST. GROESSEN DER U-KOMPONENTE
C                                 ----------------------------------
C
      IF((ISELEP(1,  1) .EQ. 1  .OR.  ISELEP(1,  1) .EQ. 2) 
     $   .AND.(IIDENT(6).LT.8)) THEN
         TSEM   = FLOAT(ITTOT) * DT
         CT     = CH1680 (' U              ')
         CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      1      ,
     $                 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     $                TSEM,TSEM, 0.0,U      ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,  4) .EQ. 1  .OR.  ISELEP(1,  4) .EQ. 2) THEN
C
C                                 DAS FELD WURDE IM MOMENTANEN LAUF
C                                 AUSGEWERTET UND WIRD DAHER AUSGEGEBEN
C
C                                 PHYSIKAL. ZEIT ZU DER MIT DER BILDUNG
C                                 VON ENSEMBLE-MITTELWERTEN BEGONNEN
C                                 WURDE. (EXAKTER: ZU DIESEM ZEITPUNKT
C                                 WURDE DIE ERSTE STICHPROBE ENTNOMMEN)
C
         TSEM   = TEEM - FLOAT((ISELEP(2,  4)-1) * ITFLUC) * DT
C
         CT     = CH1680 ('<U>             ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AU     ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,  6) .EQ. 1  .OR.  ISELEP(1,  6) .EQ. 2) THEN
         TSEM   = FLOAT(ITTOT) * DT
         CT     = CH1680 (' U-F G          ')
         CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      1      ,
     $                 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     $                TSEM,TSEM, 0.0,UFG    ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,  7) .EQ. 1  .OR.  ISELEP(1,  7) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,  7)-1) * ITFLUC) * DT
         CT     = CH1680 ('<U-RMS> G       ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AURG   ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,  9) .EQ. 1  .OR.  ISELEP(1,  9) .EQ. 2) THEN
C                                 2,  7 IST IN ORDNUNG !
         TSEM   = TEEM - FLOAT((ISELEP(2,  7)-1) * ITFLUC) * DT
         CT     = CH1680 ('<U-RMS> G+S     ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,SURG   ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1, 10) .EQ. 1  .OR.  ISELEP(1, 10) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 10)-1) * ITFLUC) * DT
         CT     = CH1680 ('<U-SKE> G       ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AUSKG  ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1, 12) .EQ. 1  .OR.  ISELEP(1, 12) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 12)-1) * ITFLUC) * DT
         CT     = CH1680 ('<U-FLA> G       ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AUFLG  ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
C                                 STATIST. GROESSEN DER V-KOMPONENTE
C                                 ----------------------------------
C
      IF((ISELEP(1, 26) .EQ. 1  .OR.  ISELEP(1, 26) .EQ. 2) 
     $   .AND.(IIDENT(6).LT.8)) THEN
         TSEM   = FLOAT(ITTOT) * DT
         CT     = CH1680 (' V              ')
         CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      0      ,
     $                 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     $                TSEM,TSEM, 0.0,V      ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1, 29) .EQ. 1  .OR.  ISELEP(1, 29) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 29)-1) * ITFLUC) * DT
         CT     = CH1680 ('<V>             ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      0      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,AV     ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1, 31) .EQ. 1  .OR.  ISELEP(1, 31) .EQ. 2) THEN
         TSEM   = FLOAT(ITTOT) * DT
         CT     = CH1680 (' V-F G          ')
         CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      0      ,
     $                 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     $                TSEM,TSEM, 0.0,VFG    ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1, 32) .EQ. 1  .OR.  ISELEP(1, 32) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 32)-1) * ITFLUC) * DT
         CT     = CH1680 ('<V-RMS> G       ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      0      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,AVRG   ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1, 34) .EQ. 1  .OR.  ISELEP(1, 34) .EQ. 2) THEN
C                                 2, 32 IST IN ORDNUNG !
         TSEM   = TEEM - FLOAT((ISELEP(2, 32)-1) * ITFLUC) * DT
         CT     = CH1680 ('<V-RMS> G+S     ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      0      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,SVRG   ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1, 35) .EQ. 1  .OR.  ISELEP(1, 35) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 35)-1) * ITFLUC) * DT
         CT     = CH1680 ('<V-SKE> G       ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      0      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,AVSKG  ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1, 37) .EQ. 1  .OR.  ISELEP(1, 37) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 37)-1) * ITFLUC) * DT
         CT     = CH1680 ('<V-FLA> G       ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      0      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,AVFLG  ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
C                                 STATIST. GROESSEN DER W-KOMPONENTE
C                                 ----------------------------------
C
      IF((ISELEP(1, 51) .EQ. 1  .OR.  ISELEP(1, 51) .EQ. 2) 
     $   .AND.(IIDENT(6).LT.8)) THEN
         TSEM   = FLOAT(ITTOT) * DT
         CT     = CH1680 (' W              ')
         CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      0      ,
     $                 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     $                TSEM,TSEM, 0.0,W      ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1, 54) .EQ. 1  .OR.  ISELEP(1, 54) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 54)-1) * ITFLUC) * DT
         CT     = CH1680 ('<W>             ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      0      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,AW     ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1, 56) .EQ. 1  .OR.  ISELEP(1, 56) .EQ. 2) THEN
         TSEM   = FLOAT(ITTOT) * DT
         CT     = CH1680 (' W-F G          ')
         CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      0      ,
     $                 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     $                TSEM,TSEM, 0.0,WFG    ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1, 57) .EQ. 1  .OR.  ISELEP(1, 57) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 57)-1) * ITFLUC) * DT
         CT     = CH1680 ('<W-RMS> G       ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      0      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,AWRG   ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1, 59) .EQ. 1  .OR.  ISELEP(1, 59) .EQ. 2) THEN
C                                 2, 57 IST IN ORDNUNG !
         TSEM   = TEEM - FLOAT((ISELEP(2, 57)-1) * ITFLUC) * DT
         CT     = CH1680 ('<W-RMS> G+S     ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      0      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,SWRG   ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1, 60) .EQ. 1  .OR.  ISELEP(1, 60) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 60)-1) * ITFLUC) * DT
         CT     = CH1680 ('<W-SKE> G       ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      0      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,AWSKG  ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1, 62) .EQ. 1  .OR.  ISELEP(1, 62) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 62)-1) * ITFLUC) * DT
         CT     = CH1680 ('<W-FLA> G       ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      0      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,AWFLG  ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
C                                 STATIST. GROESSEN DER P-KOMPONENTE
C                                 ----------------------------------
C
      IF((ISELEP(1, 76) .EQ. 1  .OR.  ISELEP(1, 76) .EQ. 2) 
     $   .AND.(IIDENT(6).LT.8)) THEN
         TSEM   = FLOAT(ITTOT) * DT
         CT     = CH1680 (' P              ')
         CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     $                TSEM,TSEM, 0.0,P      ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1, 79) .EQ. 1  .OR.  ISELEP(1, 79) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 79)-1) * ITFLUC) * DT
         CT     = CH1680 ('<P>             ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                ZSLM  ,ZELM  ,YSLM  ,YELM  ,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,AP     ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1, 81) .EQ. 1  .OR.  ISELEP(1, 81) .EQ. 2) THEN
         TSEM   = FLOAT(ITTOT) * DT
         CT     = CH1680 (' P-F G          ')
         CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     $                TSEM,TSEM, 0.0,PFG    ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1, 82) .EQ. 1  .OR.  ISELEP(1, 82) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 82)-1) * ITFLUC) * DT
         CT     = CH1680 ('<P-RMS> G       ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                ZSLM  ,ZELM  ,YSLM  ,YELM  ,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,APRG   ,CT,HILF,IGRID,IDIM3D)
      ENDIF
*     IF(ISELEP(1, 84) .EQ. 1  .OR.  ISELEP(1, 84) .EQ. 2) THEN
C                                 2, 82 IST IN ORDNUNG !
*        TSEM   = TEEM - FLOAT((ISELEP(2, 82)-1) * ITFLUC) * DT
*        CT     = CH1680 ('<P-RMS> G+S     ')
*        CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
*    $                      0      ,      0      ,      0      ,
*    $                ZSLM  ,ZELM  ,YSLM  ,YELM  ,XSLM  ,XELM  ,
*    $                TSEM,TEEM,DTEM,SPRG   ,CT,HILF,IGRID,IDIM3D)
*     ENDIF
      IF(ISELEP(1, 85) .EQ. 1  .OR.  ISELEP(1, 85) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 85)-1) * ITFLUC) * DT
         CT     = CH1680 ('<P-SKE> G       ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                ZSLM  ,ZELM  ,YSLM  ,YELM  ,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,APSKG  ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1, 87) .EQ. 1  .OR.  ISELEP(1, 87) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2, 87)-1) * ITFLUC) * DT
         CT     = CH1680 ('<P-FLA> G       ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                ZSLM  ,ZELM  ,YSLM  ,YELM  ,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,APFLG  ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
C                                 KINETISCHE ENERGIE DER SCHWANKUNGS-
C                                 KOMPONENTEN
C                                 -----------------------------------
C
      IF(ISELEP(1,106) .EQ. 1  .OR.  ISELEP(1,106) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,106)-1) * ITFLUC) * DT
         CT     = CH1680 ('<E-F> G         ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                ZSLM  ,ZELM  ,YSLM  ,YELM  ,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,AEFG   ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,108) .EQ. 1  .OR.  ISELEP(1,108) .EQ. 2) THEN
C                                 2,106 IST IN ORDNUNG !
         TSEM   = TEEM - FLOAT((ISELEP(2,106)-1) * ITFLUC) * DT
         CT     = CH1680 ('<E-F> G+S       ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                ZSLM  ,ZELM  ,YSLM  ,YELM  ,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,SEFG   ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,109) .EQ. 1  .OR.  ISELEP(1,109) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,109)-1) * ITFLUC) * DT
         CT     = CH1680 ('<E-F> S         ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                ZSLM  ,ZELM  ,YSLM  ,YELM  ,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,AEFS   ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
C                                 DISSIPATION
C                                 -----------
      IF(ISELEP(1,120) .EQ. 1  .OR.  ISELEP(1,120) .EQ. 2) THEN
            TSEM   = TEEM - FLOAT((ISELEP(2,120)-1) * ITFLUC) * DT
            CT = CH1680 ('<D-F> G         ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADFG   ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C                               FUER   TAYLOR''SCHES MIKROMASS
C                               QUADRATE DER ELEMENTE DES 
C                               DEFORMATIONSGESCHWINDIGKEITSTENSORS
C                               -----------------------------------
C
      IF(ISELEP(1,131) .EQ. 1  .OR.  ISELEP(1,131) .EQ. 2) THEN
            TSEM   = TEEM - FLOAT((ISELEP(2,131)-1) * ITFLUC) * DT
            CT = CH1680 ('<DUDX2-F> G         ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADUDX2 ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
      IF(ISELEP(1,132) .EQ. 1  .OR.  ISELEP(1,132) .EQ. 2) THEN
            TSEM   = TEEM - FLOAT((ISELEP(2,132)-1) * ITFLUC) * DT
            CT = CH1680 ('<DUDY2-F> G         ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADUDY2 ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
      IF(ISELEP(1,133) .EQ. 1  .OR.  ISELEP(1,133) .EQ. 2) THEN
            TSEM   = TEEM - FLOAT((ISELEP(2,133)-1) * ITFLUC) * DT
            CT = CH1680 ('<DUDZ2-F> G         ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADUDZ2 ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
      IF(ISELEP(1,134) .EQ. 1  .OR.  ISELEP(1,134) .EQ. 2) THEN
            TSEM   = TEEM - FLOAT((ISELEP(2,134)-1) * ITFLUC) * DT
            CT = CH1680 ('<DVDX2-F> G         ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADVDX2 ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
      IF(ISELEP(1,135) .EQ. 1  .OR.  ISELEP(1,135) .EQ. 2) THEN
            TSEM   = TEEM - FLOAT((ISELEP(2,135)-1) * ITFLUC) * DT
            CT = CH1680 ('<DVDY2-F> G         ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADVDY2 ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
      IF(ISELEP(1,136) .EQ. 1  .OR.  ISELEP(1,136) .EQ. 2) THEN
            TSEM   = TEEM - FLOAT((ISELEP(2,136)-1) * ITFLUC) * DT
            CT = CH1680 ('<DVDZ2-F> G         ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADVDZ2 ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
      IF(ISELEP(1,137) .EQ. 1  .OR.  ISELEP(1,137) .EQ. 2) THEN
            TSEM   = TEEM - FLOAT((ISELEP(2,137)-1) * ITFLUC) * DT
            CT = CH1680 ('<DWDX2-F> G         ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADWDX2 ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
      IF(ISELEP(1,138) .EQ. 1  .OR.  ISELEP(1,138) .EQ. 2) THEN
            TSEM   = TEEM - FLOAT((ISELEP(2,138)-1) * ITFLUC) * DT
            CT = CH1680 ('<DWDY2-F> G         ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADWDY2 ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
      IF(ISELEP(1,139) .EQ. 1  .OR.  ISELEP(1,139) .EQ. 2) THEN
            TSEM   = TEEM - FLOAT((ISELEP(2,139)-1) * ITFLUC) * DT
            CT = CH1680 ('<DWDZ2-F> G         ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADWDZ2 ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
C                                 ANTEILE DES SPANNUNGSTENSORS
C                                 ----------------------------
C
C                                 HIER:  U - W
C                                 ------------
C
      IF(ISELEP(1,146) .EQ. 1  .OR.  ISELEP(1,146) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,146)-1) * ITFLUC) * DT
         CT     = CH1680 ('<UW-F> G        ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AUFWFG ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,148) .EQ. 1  .OR.  ISELEP(1,148) .EQ. 2) THEN
C                                 2,146 IST IN ORDNUNG !
         TSEM   = TEEM - FLOAT((ISELEP(2,146)-1) * ITFLUC) * DT
         CT     = CH1680 ('<UW-F> G+S+M    ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,SUFWFG ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,149) .EQ. 1  .OR.  ISELEP(1,149) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,149)-1) * ITFLUC) * DT
         CT     = CH1680 ('<UW-F> S        ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AUFWFS ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,151) .EQ. 1  .OR.  ISELEP(1,151) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,151)-1) * ITFLUC) * DT
         CT     = CH1680 ('<UW-F> M        ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AUFWFM ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,153) .EQ. 1  .OR.  ISELEP(1,153) .EQ. 2) THEN
C                                 2,146 IST IN ORDNUNG !
         TSEM   = TEEM - FLOAT((ISELEP(2,146)-1) * ITFLUC) * DT
         CT     = CH1680 ('<UW-F> G+S      ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,SUFWFS ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
C                                 HIER:  V - W
C                                 ------------
C
      IF(ISELEP(1,156) .EQ. 1  .OR.  ISELEP(1,156) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,156)-1) * ITFLUC) * DT
         CT     = CH1680 ('<VW-F> G        ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,AVFWFG ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,158) .EQ. 1  .OR.  ISELEP(1,158) .EQ. 2) THEN
C                                 2,156 IST IN ORDNUNG !
         TSEM   = TEEM - FLOAT((ISELEP(2,156)-1) * ITFLUC) * DT
         CT     = CH1680 ('<VW-F> G+S+M    ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,SVFWFG ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,159) .EQ. 1  .OR.  ISELEP(1,159) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,159)-1) * ITFLUC) * DT
         CT     = CH1680 ('<VW-F> S        ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,AVFWFS ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,161) .EQ. 1  .OR.  ISELEP(1,161) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,161)-1) * ITFLUC) * DT
         CT     = CH1680 ('<VW-F> M        ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,AVFWFM ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,163) .EQ. 1  .OR.  ISELEP(1,163) .EQ. 2) THEN
C                                 2,156 IST IN ORDNUNG !
         TSEM   = TEEM - FLOAT((ISELEP(2,156)-1) * ITFLUC) * DT
         CT     = CH1680 ('<VW-F> G+S      ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,SVFWFS ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
C                                 HIER:  U - V
C                                 ------------
C
      IF(ISELEP(1,166) .EQ. 1  .OR.  ISELEP(1,166) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,166)-1) * ITFLUC) * DT
         CT     = CH1680 ('<UV-F> G        ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AUFVFG ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,168) .EQ. 1  .OR.  ISELEP(1,168) .EQ. 2) THEN
C                                 2,166 IST IN ORDNUNG !
         TSEM   = TEEM - FLOAT((ISELEP(2,166)-1) * ITFLUC) * DT
         CT     = CH1680 ('<UV-F> G+S+M    ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,SUFVFG ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,169) .EQ. 1  .OR.  ISELEP(1,169) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,169)-1) * ITFLUC) * DT
         CT     = CH1680 ('<UV-F> S        ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AUFVFS ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,171) .EQ. 1  .OR.  ISELEP(1,171) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,171)-1) * ITFLUC) * DT
         CT     = CH1680 ('<UV-F> M        ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AUFVFM ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,173) .EQ. 1  .OR.  ISELEP(1,173) .EQ. 2) THEN
C                                 2,166 IST IN ORDNUNG !
         TSEM   = TEEM - FLOAT((ISELEP(2,166)-1) * ITFLUC) * DT
         CT     = CH1680 ('<UV-F> G+S      ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,SUFVFS ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
         TSEM   = FLOAT(ITTOT) * DT
         CT     = CH1680 (' BP             ')
        CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     $                TSEM,TSEM, 0.0,BP     ,CT,HILF,IGRID,IDIM3D)

         CT     = CH1680 (' BU             ')
        CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     $                TSEM,TSEM, 0.0,BU     ,CT,HILF,IGRID,IDIM3D)

         CT     = CH1680 (' BV             ')
        CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     $                TSEM,TSEM, 0.0,BV     ,CT,HILF,IGRID,IDIM3D)

         CT     = CH1680 (' BW             ')
        CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,
     $                TSEM,TSEM, 0.0,BW     ,CT,HILF,IGRID,IDIM3D)


C
C
C
C                                 HIER:  X-KOMPONENTE DER VORTICITY
C                                 -------------------------------
C
      IF(ISELEP(1,176) .EQ. 1  .OR.  ISELEP(1,176) .EQ. 2) THEN
         TSEM   = FLOAT(ITTOT) * DT
         CT     = CH1680 (' OMEGA-X        ')
         CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,
     $                TSEM,TSEM,0.0 ,OX     ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,177) .EQ. 1  .OR.  ISELEP(1,177) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,177)-1) * ITFLUC) * DT
         CT     = CH1680 ('<OMEGA-X>       ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,AOX    ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,179) .EQ. 1  .OR.  ISELEP(1,179) .EQ. 2) THEN
         TSEM   = FLOAT(ITTOT) * DT
         CT     = CH1680 (' OMEGA-X F      ')
         CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,
     $                TSEM,TSEM,0.0 ,OXFG   ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,180) .EQ. 1  .OR.  ISELEP(1,180) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,180)-1) * ITFLUC) * DT
         CT     = CH1680 ('<OMEGA-X RMS> G ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLM  ,XELM  ,
     $                TSEM,TEEM,DTEM,AOXRG  ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
C                                 HIER:  Y-KOMPONENTE DER VORTICITY
C                                 -------------------------------
C
      IF(ISELEP(1,182) .EQ. 1  .OR.  ISELEP(1,182) .EQ. 2) THEN
         TSEM   = FLOAT(ITTOT) * DT
         CT     = CH1680 (' OMEGA-Y        ')
         CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,
     $                TSEM,TSEM,0.0 ,OY     ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,183) .EQ. 1  .OR.  ISELEP(1,183) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,183)-1) * ITFLUC) * DT
         CT     = CH1680 ('<OMEGA-Y>       ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AOY    ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,185) .EQ. 1  .OR.  ISELEP(1,185) .EQ. 2) THEN
         TSEM   = FLOAT(ITTOT) * DT
         CT     = CH1680 (' OMEGA-Y F      ')
         CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,
     $                TSEM,TSEM,0.0 ,OYFG   ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,186) .EQ. 1  .OR.  ISELEP(1,186) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,186)-1) * ITFLUC) * DT
         CT     = CH1680 ('<OMEGA-Y RMS> G ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLM  ,YELM  ,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AOYRG  ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
C                                 HIER:  Z-KOMPONENTE DER VORTICITY
C                                 -------------------------------
C
      IF(ISELEP(1,188) .EQ. 1  .OR.  ISELEP(1,188) .EQ. 2) THEN
         TSEM   = FLOAT(ITTOT) * DT
         CT     = CH1680 (' OMEGA-Z        ')
         CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,
     $                TSEM,TSEM,0.0 ,OZ     ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,189) .EQ. 1  .OR.  ISELEP(1,189) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,189)-1) * ITFLUC) * DT
         CT     = CH1680 ('<OMEGA-Z>       ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AOZ    ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,191) .EQ. 1  .OR.  ISELEP(1,191) .EQ. 2) THEN
         TSEM   = FLOAT(ITTOT) * DT
         CT     = CH1680 (' OMEGA-Z F      ')
         CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,
     $                TSEM,TSEM,0.0 ,OZFG   ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,192) .EQ. 1  .OR.  ISELEP(1,192) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,192)-1) * ITFLUC) * DT
         CT     = CH1680 ('<OMEGA-Z RMS> G ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                ZSLM  ,ZELM  ,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AOZRG  ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
C                                 HIER:  ENSTROPHIE
C                                 -------------------------------
C
      IF(ISELEP(1,195) .EQ. 1  .OR.  ISELEP(1,195) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,195)-1) * ITFLUC) * DT
         CT     = CH1680 ('<ENSTROPHY>     ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AO2    ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,197) .EQ. 1  .OR.  ISELEP(1,197) .EQ. 2) THEN
         TSEM   = FLOAT(ITTOT) * DT
         CT     = CH1680 (' ENSTROPHY-F G  ')
         CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      1      ,
     $                0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,
     $                TSEM,TSEM,0.0 ,O2FG   ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,198) .EQ. 1  .OR.  ISELEP(1,198) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,198)-1) * ITFLUC) * DT
         CT     = CH1680 ('<ENSTROPHY-RMS>G')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AO2RG  ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
C                                 HIER:  HELIZITAET
C                                 -------------------------------
C
      IF(ISELEP(1,200) .EQ. 1  .OR.  ISELEP(1,200) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,200)-1) * ITFLUC) * DT
         CT     = CH1680 ('<HELICITY>      ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AHE    ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,202) .EQ. 1  .OR.  ISELEP(1,202) .EQ. 2) THEN
         TSEM   = FLOAT(ITTOT) * DT
         CT     = CH1680 (' HELICITY-F G   ')
         CALL DOBC   (KK ,JJ ,II ,KMX ,JMX ,IMX ,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      1      ,
     $                0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,
     $                TSEM,TSEM,0.0 ,HEFG   ,CT,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEP(1,203) .EQ. 1  .OR.  ISELEP(1,203) .EQ. 2) THEN
         TSEM   = TEEM - FLOAT((ISELEP(2,203)-1) * ITFLUC) * DT
         CT     = CH1680 ('<HELICITY-RMS> G')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AHERG  ,CT,HILF,IGRID,IDIM3D)
      ENDIF
C
C
C
c     "Neue Statistik"
c
c

      IF(ISELEP(1,400) .EQ. 1  .OR.  ISELEP(1,400) .EQ. 2) THEN
         TSEM   = TEEM-FLOAT((ISELEP(2,400)-1) * ITFLUC) * DT
         CT     = CH1680 ('<UU> G          ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AUUM  ,CT,HILF,IGRID,IDIM3D)
      ENDIF


      IF(ISELEP(1,401) .EQ. 1  .OR.  ISELEP(1,401) .EQ. 2) THEN
         TSEM   = TEEM-FLOAT((ISELEP(2,401)-1) * ITFLUC) * DT
         CT     = CH1680 ('<VV> G          ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AVVM  ,CT,HILF,IGRID,IDIM3D)
      ENDIF


      IF(ISELEP(1,402) .EQ. 1  .OR.  ISELEP(1,402) .EQ. 2) THEN
         TSEM   = TEEM-FLOAT((ISELEP(2,402)-1) * ITFLUC) * DT
         CT     = CH1680 ('<WW> G          ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AWWM  ,CT,HILF,IGRID,IDIM3D)
      ENDIF


      IF(ISELEP(1,406) .EQ. 1  .OR.  ISELEP(1,406) .EQ. 2) THEN
         TSEM   = TEEM-FLOAT((ISELEP(2,406)-1) * ITFLUC) * DT
         CT     = CH1680 ('<PP> G          ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,APPM  ,CT,HILF,IGRID,IDIM3D)
      ENDIF




      IF(ISELEP(1,403) .EQ. 1  .OR.  ISELEP(1,403) .EQ. 2) THEN
         TSEM   = TEEM-FLOAT((ISELEP(2,403)-1) * ITFLUC) * DT
         CT     = CH1680 ('<UV> G          ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AUVM  ,CT,HILF,IGRID,IDIM3D)
      ENDIF


      IF(ISELEP(1,404) .EQ. 1  .OR.  ISELEP(1,404) .EQ. 2) THEN
         TSEM   = TEEM-FLOAT((ISELEP(2,404)-1) * ITFLUC) * DT
         CT     = CH1680 ('<UW> G          ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AUWM  ,CT,HILF,IGRID,IDIM3D)
      ENDIF


      IF(ISELEP(1,405) .EQ. 1  .OR.  ISELEP(1,405) .EQ. 2) THEN
         TSEM   = TEEM-FLOAT((ISELEP(2,405)-1) * ITFLUC) * DT
         CT     = CH1680 ('<VW> G          ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AVWM  ,CT,HILF,IGRID,IDIM3D)
      ENDIF




      IF(ISELEP(1,430) .EQ. 1  .OR.  ISELEP(1,430) .EQ. 2) THEN
         TSEM   = TEEM-FLOAT((ISELEP(2,430)-1) * ITFLUC) * DT
         CT     = CH1680 ('<(dU/dX)^2> G   ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AUXUXM,CT,HILF,IGRID,IDIM3D)
      ENDIF


      IF(ISELEP(1,431) .EQ. 1  .OR.  ISELEP(1,431) .EQ. 2) THEN
         TSEM   = TEEM-FLOAT((ISELEP(2,431)-1) * ITFLUC) * DT
         CT     = CH1680 ('<(dU/dY)^2> G   ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AUYUYM,CT,HILF,IGRID,IDIM3D)
      ENDIF


      IF(ISELEP(1,432) .EQ. 1  .OR.  ISELEP(1,432) .EQ. 2) THEN
         TSEM   = TEEM-FLOAT((ISELEP(2,432)-1) * ITFLUC) * DT
         CT     = CH1680 ('<(dU/dZ)^2> G   ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AUZUZM,CT,HILF,IGRID,IDIM3D)
      ENDIF


      IF(ISELEP(1,433) .EQ. 1  .OR.  ISELEP(1,433) .EQ. 2) THEN
         TSEM   = TEEM-FLOAT((ISELEP(2,433)-1) * ITFLUC) * DT
         CT     = CH1680 ('<(dV/dX)^2> G   ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      1      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AVXVXM,CT,HILF,IGRID,IDIM3D)
      ENDIF


      IF(ISELEP(1,434) .EQ. 1  .OR.  ISELEP(1,434) .EQ. 2) THEN
         TSEM   = TEEM-FLOAT((ISELEP(2,434)-1) * ITFLUC) * DT
         CT     = CH1680 ('<(dV/dY)^2> G   ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AVYVYM,CT,HILF,IGRID,IDIM3D)
      ENDIF


      IF(ISELEP(1,435) .EQ. 1  .OR.  ISELEP(1,435) .EQ. 2) THEN
         TSEM   = TEEM-FLOAT((ISELEP(2,435)-1) * ITFLUC) * DT
         CT     = CH1680 ('<(dV/dZ)^2> G   ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AVZVZM,CT,HILF,IGRID,IDIM3D)
      ENDIF


      IF(ISELEP(1,436) .EQ. 1  .OR.  ISELEP(1,436) .EQ. 2) THEN
         TSEM   = TEEM-FLOAT((ISELEP(2,436)-1) * ITFLUC) * DT
         CT     = CH1680 ('<(dW/dX)^2> G   ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      0      ,      1      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AWXWXM,CT,HILF,IGRID,IDIM3D)
      ENDIF


      IF(ISELEP(1,437) .EQ. 1  .OR.  ISELEP(1,437) .EQ. 2) THEN
         TSEM   = TEEM-FLOAT((ISELEP(2,437)-1) * ITFLUC) * DT
         CT     = CH1680 ('<(dW/dY)^2> G   ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      1      ,      1      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AWYWYM,CT,HILF,IGRID,IDIM3D)
      ENDIF


      IF(ISELEP(1,438) .EQ. 1  .OR.  ISELEP(1,438) .EQ. 2) THEN
         TSEM   = TEEM-FLOAT((ISELEP(2,438)-1) * ITFLUC) * DT
         CT     = CH1680 ('<(dW/dZ)^2> G   ')
         CALL DOBC   (KKA,JJA,IIA,KMXA,JMXA,IMXA,NBND,KANAL,MODUS,
     $                      0      ,      0      ,      0      ,
     $                ZSLMST,ZELMST,YSLMST,YELMST,XSLMST,XELMST,
     $                TSEM,TEEM,DTEM,AWZWZM,CT,HILF,IGRID,IDIM3D)
      ENDIF






C
C
      RETURN
 6010 FORMAT (5(I9,1X))
 6015 FORMAT (4(I19,1X))
 6016 FORMAT (  I19    )
 6020 FORMAT (4(E19.12E3,1X))
 6021 FORMAT (  E19.12E3    )
      END
