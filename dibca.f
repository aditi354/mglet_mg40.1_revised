










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
      SUBROUTINE DIBCA   (KANAL,MODUS,IIDENT,
     $                    ILIMXA,IDIBD1,IDIBD2,IDIBD3,IDIBD4,RKOMXA,
     $                    RDIBD1,RDIBD2,RDIBD3,RDIBD4,IGRID,IDIM3D,

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
C        D I B C A        EINLESEN VON FELDERN, DIE IM WESENTLICHEN
C                         ENSEMBLE-MITTELWERTE ENTHALTEN.
C        A C H T U N G    DIESES UNTERPROGRAMM MUSS ERWEITERT WERDEN,
C                         WENN DIE STATIST. AUSWERTUNG UMFANGREICHER
C                         WIRD !
C                         DIE NEU HINZUKOMMENDEN FELDER SIND GEMAESS
C                         DER ORDNUNG IM ISELEP-FELD EINZUREIHEN.
C*STARLET***************************************************************
C
C PARAM: KANAL                    - DIE DATEN WERDEN UEBER DIESEN KANAL
C                                   EINGELESEN
C        MODUS                    - CHARACTER-VARIABLE:
C                                   'BINAER  ' : DATEN WERDEN UNFORMA-
C                                                TIERT GESCHRIEBEN
C                                   'CODIERT ' : DATEN WERDEN FORMA-
C                                                TIERT GESCHRIEBEN
C        ILIMXA                   + ANZAHL DER AUFPUNKTE + 1 DES VORAN-
C                                   GEGANGENEN LAUFES
C        IDIBD1 ... IDIBD4        + NOCH UNBELEGTE DUMMY-VARIABLE,
C                                   WELCHE VERWENDET WERDEN KOENNEN,
C                                   UM WICHTIGE INFORMATIONEN DES
C                                   VORANGEGANGENEN LAUFES FUER DEN
C                                   MOMENTANEN LAUF BEREITSTELLEN ZU
C                                   KOENNEN (INTEGER-VARIABLE !)
C        RKOMXA                   + MAXIMAL ZULAESSIGER KORRELATIONS-
C                                   RADIUS WAEHREND DES VORANGEG. LAUFES
C        RDIBD1 ... RDIBD4        + WIE IDIBD. (REAL-VARIABLE !)
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
C UPROG                 : DIBC, ERRR
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        12.08.86 (HW)  : ORIGINAL
C        29.08.86 (HW)  : BEI DEN KOMPONENTEN DES SCHUBSPANNUNGSTENSORS
C                         WIRD DIE SUMME AUS GROB- UND FEINSTRUKTUR-
C                         ANTEIL (REYNOLDSSPANNUNG) EINGELESEN
C        30.09.86 (HW)  : MOMENTANWERTE DER GESCHW. SOWIE FLUKTUATIONEN
C                         WERDEN EINGELESEN (IN DAS DUMMY-FELD "HILF")
C        09.12.86 (HW)  : ES KANN WAHLWEISE FORMATIERT ODER UNFORMA-
C                         TIERT GELESEN WERDEN. (VARIABLE 'MODUS'
C                         EINGEFUEHRT)
C        17.12.86 (HW)  : WEGEN NOS/VE MUSSTE DIE FUNCTION CH1680
C                         EINGEFUEHRT WERDEN
C        05.07.88 (HW)  : ERWEITERUNG: AUSWERTUNG DER DREI KOMPONENTEN
C                         DER VORTICITY, DER ENSTROPHIE UND DER
C                         HELIZITAET
C        08.12.88 (HW)  : AENDERUNG DER UEBERGABE. DAS ERWEITERTE ISLEA-
C                         FELD UND DAS ISLINA-FELD WERDEN GELESEN
C         6. 4.92 (MM)  : IIDENT EINGEFUEHRT 
C        26. 5.95 (MM)  : UEBERGABE NACH DOBC GEAENDERT
C	  5. 3.01 (SE)	: Statistik fuer balances eingebaut
C        07.02.03 (TB)  : SCALAR FIELD STATISTICS ADDED
C
C*STARLET***************************************************************
C
      CHARACTER (LEN=8)   MODUS
      CHARACTER (LEN=80)  CT,  CH1680
      INTEGER  IIDENT(100)
C

      INTEGER  ISELEA (2,752)      , ISELEP (2,752),
     $         ISAMPA (752)        , ISAMPP (752),
     $         ISLINP (ISLIDI)     , ISLINA (ISLIDI)
      REAL     RIDENT (100)
      REAL     HILF   (KK ,JJ ,II ), B      (KK ,JJ ,II )
C
      REAL     U      (KK ,JJ ,II ), UO     (KK ,JJ ,II ),
     $         AU     (KKA,JJA,IIA), SU     (KKA,JJA,IIA),
     $         UFG    (KK ,JJ ,II ),
     $         AURG   (KKA,JJA,IIA), SURG   (KKA,JJA,IIA),
     $         AUSKG  (KKA,JJA,IIA), SUSKG  (KKA,JJA,IIA),
     $         AUFLG  (KKA,JJA,IIA), SUFLG  (KKA,JJA,IIA)
      REAL     V      (KK ,JJ ,II ), VO     (KK ,JJ ,II ),
     $         AV     (KKA,JJA,IIA), SV     (KKA,JJA,IIA),
     $         VFG    (KK ,JJ ,II ),
     $         AVRG   (KKA,JJA,IIA), SVRG   (KKA,JJA,IIA),
     $         AVSKG  (KKA,JJA,IIA), SVSKG  (KKA,JJA,IIA),
     $         AVFLG  (KKA,JJA,IIA), SVFLG  (KKA,JJA,IIA)
      REAL     W      (KK ,JJ ,II ), WO     (KK ,JJ ,II ),
     $         AW     (KKA,JJA,IIA), SW     (KKA,JJA,IIA),
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
      REAL BP(KK,JJ,II),BU(KK,JJ,II),BV(KK,JJ,II),BW(KK,JJ,II)
      REAL     HILF3D1(KK ,JJ ,II ),HILF3D2(KK ,JJ ,II ),
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
C
C
C                                 EINLESEN DER FELDER, DIE ENSEMBLE-
C                                 MITTELWERTE ENTHALTEN. DIE REIHEN-
C                                 FOLGE ENTSPRICHT DER IM ISELEP-FELD
C
C                                 STATIST. GROESSEN DER U-KOMPONENTE
C                                 ----------------------------------
C
10000 CONTINUE
C
      WRITE(6,*)'ISELEA: ',ISELEA(1,352)
      WRITE(6,*)'IIDENT 5 ',IIDENT(5)     
      IF((ISELEA(1,  1) .EQ. 1  .OR.  ISELEA(1,  1) .EQ. 2)
     $   .AND.(IIDENT(5).LT.8)) THEN
            CT = CH1680 (' U              ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,  4) .EQ. 1  .OR.  ISELEA(1,  4) .EQ. 2) THEN
C
C                                 DAS FELD WURDE IM VORANGEGANGENEN LAUF
C                                 AUSGEWERTET UND MUSS DAHER GELESEN
C                                 WERDEN
C
         IF(ISELEP(1,  4) .EQ. 1  .OR.  ISELEP(1,  4) .EQ. 2) THEN
C
C                                 DAS FELD WIRD AUCH WAEHREND DES MOMEN-
C                                 TANEN LAUFES BENOETIGT UND WIRD DAHER
C                                 IN DAS VORGESEHENE A-FELD EINGELESEN
C
            CT = CH1680 ('<U>             ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AU     ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
C
C                                 DAS FELD WIRD WAEHREND DES MOMEN-
C                                 TANEN LAUFES NICHT BENOETIGT UND
C                                 WIRD DAHER IN DAS HILF -FELD GE-
C                                 LESEN (TAPE1 WIRD KORREKT POSITI-
C                                 ONIERT !)
C
            CT = CH1680 ('<U>             ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,  6) .EQ. 1  .OR.  ISELEA(1,  6) .EQ. 2) THEN
            CT = CH1680 (' U-F G          ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,  7) .EQ. 1  .OR.  ISELEA(1,  7) .EQ. 2) THEN
         IF(ISELEP(1,  7) .EQ. 1  .OR.  ISELEP(1,  7) .EQ. 2) THEN
            CT = CH1680 ('<U-RMS> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AURG   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<U-RMS> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C
C                                 FUER <U-RMS> G+S  STEHT KEIN EIGENES
C                                 FELD ZUR VERFUEGUNG, DAHER MUSS BEIM
C                                 EINLESEN IN DAS FELD 'HILF' GE-
C                                 SCHRIEBEN WERDEN.
C
      IF(ISELEA(1,  9) .EQ. 1  .OR.  ISELEA(1,  9) .EQ. 2) THEN
            CT = CH1680 ('<U-RMS> G+S     ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1, 10) .EQ. 1  .OR.  ISELEA(1, 10) .EQ. 2) THEN
         IF(ISELEP(1, 10) .EQ. 1  .OR.  ISELEP(1, 10) .EQ. 2) THEN
            CT = CH1680 ('<U-SKE> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUSKG  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<U-SKE> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1, 12) .EQ. 1  .OR.  ISELEA(1, 12) .EQ. 2) THEN
         IF(ISELEP(1, 12) .EQ. 1  .OR.  ISELEP(1, 12) .EQ. 2) THEN
            CT = CH1680 ('<U-FLA> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUFLG  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<U-FLA> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C
C                                 STATIST. GROESSEN DER V-KOMPONENTE
C                                 ----------------------------------
C
      IF((ISELEA(1, 26) .EQ. 1  .OR.  ISELEA(1, 26) .EQ. 2) 
     $   .AND.(IIDENT(5).LT.8)) THEN
            CT = CH1680 (' V              ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1, 29) .EQ. 1  .OR.  ISELEA(1, 29) .EQ. 2) THEN
         IF(ISELEP(1, 29) .EQ. 1  .OR.  ISELEP(1, 29) .EQ. 2) THEN
            CT = CH1680 ('<V>             ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AV     ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<V>             ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1, 31) .EQ. 1  .OR.  ISELEA(1, 31) .EQ. 2) THEN
            CT = CH1680 (' V-F G          ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1, 32) .EQ. 1  .OR.  ISELEA(1, 32) .EQ. 2) THEN
         IF(ISELEP(1, 32) .EQ. 1  .OR.  ISELEP(1, 32) .EQ. 2) THEN
            CT = CH1680 ('<V-RMS> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVRG   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<V-RMS> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1, 34) .EQ. 1  .OR.  ISELEA(1, 34) .EQ. 2) THEN
            CT = CH1680 ('<V-RMS> G+S     ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1, 35) .EQ. 1  .OR.  ISELEA(1, 35) .EQ. 2) THEN
         IF(ISELEP(1, 35) .EQ. 1  .OR.  ISELEP(1, 35) .EQ. 2) THEN
            CT = CH1680 ('<V-SKE> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVSKG  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<V-SKE> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1, 37) .EQ. 1  .OR.  ISELEA(1, 37) .EQ. 2) THEN
         IF(ISELEP(1, 37) .EQ. 1  .OR.  ISELEP(1, 37) .EQ. 2) THEN
            CT = CH1680 ('<V-FLA> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVFLG  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<V-FLA> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C
C                                 STATIST. GROESSEN DER W-KOMPONENTE
C                                 ----------------------------------
C
      IF((ISELEA(1, 51) .EQ. 1  .OR.  ISELEA(1, 51) .EQ. 2) 
     $   .AND.(IIDENT(5).LT.8)) THEN
            CT = CH1680 (' W              ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1, 54) .EQ. 1  .OR.  ISELEA(1, 54) .EQ. 2) THEN
         IF(ISELEP(1, 54) .EQ. 1  .OR.  ISELEP(1, 54) .EQ. 2) THEN
            CT = CH1680 ('<W>             ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AW     ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<W>             ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1, 56) .EQ. 1  .OR.  ISELEA(1, 56) .EQ. 2) THEN
            CT = CH1680 (' W-F G          ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1, 57) .EQ. 1  .OR.  ISELEA(1, 57) .EQ. 2) THEN
         IF(ISELEP(1, 57) .EQ. 1  .OR.  ISELEP(1, 57) .EQ. 2) THEN
            CT = CH1680 ('<W-RMS> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AWRG   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<W-RMS> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1, 59) .EQ. 1  .OR.  ISELEA(1, 59) .EQ. 2) THEN
            CT = CH1680 ('<W-RMS> G+S     ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1, 60) .EQ. 1  .OR.  ISELEA(1, 60) .EQ. 2) THEN
         IF(ISELEP(1, 60) .EQ. 1  .OR.  ISELEP(1, 60) .EQ. 2) THEN
            CT = CH1680 ('<W-SKE> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AWSKG  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<W-SKE> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1, 62) .EQ. 1  .OR.  ISELEA(1, 62) .EQ. 2) THEN
         IF(ISELEP(1, 62) .EQ. 1  .OR.  ISELEP(1, 62) .EQ. 2) THEN
            CT = CH1680 ('<W-FLA> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AWFLG  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<W-FLA> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C
C                                 STATIST. GROESSEN DER P-KOMPONENTE
C                                 ----------------------------------
C
      IF((ISELEA(1, 76) .EQ. 1  .OR.  ISELEA(1, 76) .EQ. 2) 
     $   .AND.(IIDENT(5).LT.8)) THEN
            CT = CH1680 (' P              ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1, 79) .EQ. 1  .OR.  ISELEA(1, 79) .EQ. 2) THEN
         IF(ISELEP(1, 79) .EQ. 1  .OR.  ISELEP(1, 79) .EQ. 2) THEN
            CT = CH1680 ('<P>             ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AP     ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<P>             ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1, 81) .EQ. 1  .OR.  ISELEA(1, 81) .EQ. 2) THEN
            CT = CH1680 (' P-F G          ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1, 82) .EQ. 1  .OR.  ISELEA(1, 82) .EQ. 2) THEN
         IF(ISELEP(1, 82) .EQ. 1  .OR.  ISELEP(1, 82) .EQ. 2) THEN
            CT = CH1680 ('<P-RMS> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   APRG   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<P-RMS> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
*     IF(ISELEA(1, 84) .EQ. 1  .OR.  ISELEA(1, 84) .EQ. 2) THEN
*           CT = CH1680 ('<P-RMS> G+S     ')
*           CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
*    $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
*    $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
*     ENDIF
      IF(ISELEA(1, 85) .EQ. 1  .OR.  ISELEA(1, 85) .EQ. 2) THEN
         IF(ISELEP(1, 85) .EQ. 1  .OR.  ISELEP(1, 85) .EQ. 2) THEN
            CT = CH1680 ('<P-SKE> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   APSKG  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<P-SKE> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1, 87) .EQ. 1  .OR.  ISELEA(1, 87) .EQ. 2) THEN
         IF(ISELEP(1, 87) .EQ. 1  .OR.  ISELEP(1, 87) .EQ. 2) THEN
            CT = CH1680 ('<P-FLA> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   APFLG  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<P-FLA> G       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C
C                                 KINETISCHE ENERGIE DER SCHWANKUNGS-
C                                 KOMPONENTEN
C                                 -----------------------------------
C
      IF(ISELEA(1,106) .EQ. 1  .OR.  ISELEA(1,106) .EQ. 2) THEN
         IF(ISELEP(1,106) .EQ. 1  .OR.  ISELEP(1,106) .EQ. 2) THEN
            CT = CH1680 ('<E-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AEFG   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<E-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,108) .EQ. 1  .OR.  ISELEA(1,108) .EQ. 2) THEN
            CT = CH1680 ('<E-F> G+S       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,109) .EQ. 1  .OR.  ISELEA(1,109) .EQ. 2) THEN
         IF(ISELEP(1,109) .EQ. 1  .OR.  ISELEP(1,109) .EQ. 2) THEN
            CT = CH1680 ('<E-F> S         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AEFS   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<E-F> S         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C                                 DISSIPATION
C                                 -----------
      IF(ISELEA(1,120) .EQ. 1  .OR.  ISELEA(1,120) .EQ. 2) THEN
         IF(ISELEP(1,120) .EQ. 1  .OR.  ISELEP(1,120) .EQ. 2) THEN
            CT = CH1680 ('<D-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADFG   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<D-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C                               FUER   TAYLOR''SCHES MIKROMASS
C                               QUADRATE DER ELEMENTE DES
C                               DEFORMATIONSGESCHWINDIGKEITSTENSORS
C                               -----------------------------------
      IF(ISELEA(1,131) .EQ. 1  .OR.  ISELEA(1,131) .EQ. 2) THEN
         IF(ISELEP(1,131) .EQ. 1  .OR.  ISELEP(1,131) .EQ. 2) THEN
            CT = CH1680 ('<DUDX2-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADUDX2   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<DUDX2-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,132) .EQ. 1  .OR.  ISELEA(1,132) .EQ. 2) THEN
         IF(ISELEP(1,132) .EQ. 1  .OR.  ISELEP(1,132) .EQ. 2) THEN
            CT = CH1680 ('<DUDY2-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADUDY2   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<DUDY2-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,133) .EQ. 1  .OR.  ISELEA(1,133) .EQ. 2) THEN
         IF(ISELEP(1,133) .EQ. 1  .OR.  ISELEP(1,133) .EQ. 2) THEN
            CT = CH1680 ('<DUDZ2-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADUDZ2   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<DUDZ2-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,134) .EQ. 1  .OR.  ISELEA(1,134) .EQ. 2) THEN
         IF(ISELEP(1,134) .EQ. 1  .OR.  ISELEP(1,134) .EQ. 2) THEN
            CT = CH1680 ('<DVDX2-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADVDX2   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<DVDX2-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,135) .EQ. 1  .OR.  ISELEA(1,135) .EQ. 2) THEN
         IF(ISELEP(1,135) .EQ. 1  .OR.  ISELEP(1,135) .EQ. 2) THEN
            CT = CH1680 ('<DVDY2-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADVDY2   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<DVDY2-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,136) .EQ. 1  .OR.  ISELEA(1,136) .EQ. 2) THEN
         IF(ISELEP(1,136) .EQ. 1  .OR.  ISELEP(1,136) .EQ. 2) THEN
            CT = CH1680 ('<DVDZ2-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADVDZ2   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<DVDZ2-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,137) .EQ. 1  .OR.  ISELEA(1,137) .EQ. 2) THEN
         IF(ISELEP(1,137) .EQ. 1  .OR.  ISELEP(1,137) .EQ. 2) THEN
            CT = CH1680 ('<DWDX2-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADWDX2   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<DWDX2-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,138) .EQ. 1  .OR.  ISELEA(1,138) .EQ. 2) THEN
         IF(ISELEP(1,138) .EQ. 1  .OR.  ISELEP(1,138) .EQ. 2) THEN
            CT = CH1680 ('<DWDY2-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADWDY2   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<DWDY2-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,139) .EQ. 1  .OR.  ISELEA(1,139) .EQ. 2) THEN
         IF(ISELEP(1,139) .EQ. 1  .OR.  ISELEP(1,139) .EQ. 2) THEN
            CT = CH1680 ('<DWDZ2-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ADWDZ2   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<DWDZ2-F> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C
C                                 ANTEILE DES SPANNUNGSTENSORS
C                                 ----------------------------
C
C                                 HIER:  U - W
C                                 ------------
C
      IF(ISELEA(1,146) .EQ. 1  .OR.  ISELEA(1,146) .EQ. 2) THEN
         IF(ISELEP(1,146) .EQ. 1  .OR.  ISELEP(1,146) .EQ. 2) THEN
            CT = CH1680 ('<UW-F> G        ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUFWFG ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<UW-F> G        ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,148) .EQ. 1  .OR.  ISELEA(1,148) .EQ. 2) THEN
            CT = CH1680 ('<UW-F> G+S+M    ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,149) .EQ. 1  .OR.  ISELEA(1,149) .EQ. 2) THEN
         IF(ISELEP(1,149) .EQ. 1  .OR.  ISELEP(1,149) .EQ. 2) THEN
            CT = CH1680 ('<UW-F> S        ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUFWFS ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<UW-F> S        ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,151) .EQ. 1  .OR.  ISELEA(1,151) .EQ. 2) THEN
         IF(ISELEP(1,151) .EQ. 1  .OR.  ISELEP(1,151) .EQ. 2) THEN
            CT = CH1680 ('<UW-F> M        ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUFWFM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<UW-F> M        ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,153) .EQ. 1  .OR.  ISELEA(1,153) .EQ. 2) THEN
            CT = CH1680 ('<UW-F> G+S      ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
C
C                                 HIER:  V - W
C                                 ------------
C
      IF(ISELEA(1,156) .EQ. 1  .OR.  ISELEA(1,156) .EQ. 2) THEN
         IF(ISELEP(1,156) .EQ. 1  .OR.  ISELEP(1,156) .EQ. 2) THEN
            CT = CH1680 ('<VW-F> G        ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVFWFG ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<VW-F> G        ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,158) .EQ. 1  .OR.  ISELEA(1,158) .EQ. 2) THEN
            CT = CH1680 ('<VW-F> G+S+M    ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,159) .EQ. 1  .OR.  ISELEA(1,159) .EQ. 2) THEN
         IF(ISELEP(1,159) .EQ. 1  .OR.  ISELEP(1,159) .EQ. 2) THEN
            CT = CH1680 ('<VW-F> S        ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVFWFS ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<VW-F> S        ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,161) .EQ. 1  .OR.  ISELEA(1,161) .EQ. 2) THEN
         IF(ISELEP(1,161) .EQ. 1  .OR.  ISELEP(1,161) .EQ. 2) THEN
            CT = CH1680 ('<VW-F> M        ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVFWFM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<VW-F> M        ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,163) .EQ. 1  .OR.  ISELEA(1,163) .EQ. 2) THEN
            CT = CH1680 ('<VW-F> G+S      ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
C
C                                 HIER:  U - V
C                                 ------------
C
      IF(ISELEA(1,166) .EQ. 1  .OR.  ISELEA(1,166) .EQ. 2) THEN
         IF(ISELEP(1,166) .EQ. 1  .OR.  ISELEP(1,166) .EQ. 2) THEN
            CT = CH1680 ('<UV-F> G        ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUFVFG ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<UV-F> G        ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,168) .EQ. 1  .OR.  ISELEA(1,168) .EQ. 2) THEN
            CT = CH1680 ('<UV-F> G+S+M    ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,169) .EQ. 1  .OR.  ISELEA(1,169) .EQ. 2) THEN
         IF(ISELEP(1,169) .EQ. 1  .OR.  ISELEP(1,169) .EQ. 2) THEN
            CT = CH1680 ('<UV-F> S        ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUFVFS ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<UV-F> S        ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,171) .EQ. 1  .OR.  ISELEA(1,171) .EQ. 2) THEN
         IF(ISELEP(1,171) .EQ. 1  .OR.  ISELEP(1,171) .EQ. 2) THEN
            CT = CH1680 ('<UV-F> M        ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUFVFM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<UV-F> M        ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,173) .EQ. 1  .OR.  ISELEA(1,173) .EQ. 2) THEN
            CT = CH1680 ('<UV-F> G+S      ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,301) .EQ. 1) THEN
            CT = CH1680 (' BP             ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,302) .EQ. 1) THEN
            CT = CH1680 (' BU             ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,303) .EQ. 1) THEN
            CT = CH1680 (' BV             ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,304) .EQ. 1) THEN
            CT = CH1680 (' BW             ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF

C 
C                                 ANTEILE DES SPANNUNGSTENSORS
C                                 NICHT-NEWTONSCHE ANTEILE
C                                 ----------------------------
C
C                                 HIER:  U - U
C                                 ------------
      IF(ISELEA(1,352) .EQ. 1  .OR.  ISELEA(1,352) .EQ. 2) THEN
            CT = CH1680 (' TAU11          ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,351) .EQ. 1  .OR.  ISELEA(1,351) .EQ. 2) THEN
         IF(ISELEP(1,351) .EQ. 1  .OR.  ISELEP(1,351) .EQ. 2) THEN
            CT = CH1680 ('<TAU11>         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ATAU11 ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<TAU11>         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C                                 HIER:  U - V
C                                 ------------
C
      IF(ISELEA(1,357) .EQ. 1  .OR.  ISELEA(1,357) .EQ. 2) THEN
            CT = CH1680 (' TAU12          ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,356) .EQ. 1  .OR.  ISELEA(1,356) .EQ. 2) THEN
         IF(ISELEP(1,356) .EQ. 1  .OR.  ISELEP(1,356) .EQ. 2) THEN
            CT = CH1680 ('<TAU12>         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ATAU12 ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<TAU12>         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C                                 HIER:  U - W
C                                 ------------
C
      IF(ISELEA(1,362) .EQ. 1  .OR.  ISELEA(1,362) .EQ. 2) THEN
            CT = CH1680 (' TAU13          ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,361) .EQ. 1  .OR.  ISELEA(1,361) .EQ. 2) THEN
         IF(ISELEP(1,361) .EQ. 1  .OR.  ISELEP(1,361) .EQ. 2) THEN
            CT = CH1680 ('<TAU13>         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ATAU13 ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<TAU13>         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C                                 HIER:  V - V
C                                 ------------
C
      IF(ISELEA(1,367) .EQ. 1  .OR.  ISELEA(1,367) .EQ. 2) THEN
            CT = CH1680 (' TAU22          ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,366) .EQ. 1  .OR.  ISELEA(1,366) .EQ. 2) THEN
         IF(ISELEP(1,366) .EQ. 1  .OR.  ISELEP(1,366) .EQ. 2) THEN
            CT = CH1680 ('<TAU22>         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ATAU22 ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<TAU22>         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C                                 HIER:  V - W
C                                 ------------
C
      IF(ISELEA(1,372) .EQ. 1  .OR.  ISELEA(1,372) .EQ. 2) THEN
            CT = CH1680 (' TAU23          ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,371) .EQ. 1  .OR.  ISELEA(1,371) .EQ. 2) THEN
         IF(ISELEP(1,371) .EQ. 1  .OR.  ISELEP(1,371) .EQ. 2) THEN
            CT = CH1680 ('<TAU23>         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ATAU23 ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<TAU23>         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C                                 HIER:  W - W
C                                 ------------
      IF(ISELEA(1,377) .EQ. 1  .OR.  ISELEA(1,377) .EQ. 2) THEN
            CT = CH1680 (' TAU33          ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,376) .EQ. 1  .OR.  ISELEA(1,376) .EQ. 2) THEN
         IF(ISELEP(1,376) .EQ. 1  .OR.  ISELEP(1,376) .EQ. 2) THEN
            CT = CH1680 ('<TAU33>         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ATAU33 ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<TAU33>         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C
C                                 HIER: X-KOMPONENTE DER VORTICITY
C                                 ------------------------------
C
      IF(ISELEA(1,176) .EQ. 1  .OR.  ISELEA(1,176) .EQ. 2) THEN
            CT = CH1680 (' OMEGA-X        ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,177) .EQ. 1  .OR.  ISELEA(1,177) .EQ. 2) THEN
         IF(ISELEP(1,177) .EQ. 1  .OR.  ISELEP(1,177) .EQ. 2) THEN
            CT = CH1680 ('<OMEGA-X>       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AOX    ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<OMEGA-X>       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,179) .EQ. 1  .OR.  ISELEA(1,179) .EQ. 2) THEN
            CT = CH1680 (' OMEGA-X F      ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,180) .EQ. 1  .OR.  ISELEA(1,180) .EQ. 2) THEN
         IF(ISELEP(1,180) .EQ. 1  .OR.  ISELEP(1,180) .EQ. 2) THEN
            CT = CH1680 ('<OMEGA-X RMS> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AOXRG  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<OMEGA-X RMS> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C
C                                 HIER: Y-KOMPONENTE DER VORTICITY
C                                 ------------------------------
C
      IF(ISELEA(1,182) .EQ. 1  .OR.  ISELEA(1,182) .EQ. 2) THEN
            CT = CH1680 (' OMEGA-Y        ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,183) .EQ. 1  .OR.  ISELEA(1,183) .EQ. 2) THEN
         IF(ISELEP(1,183) .EQ. 1  .OR.  ISELEP(1,183) .EQ. 2) THEN
            CT = CH1680 ('<OMEGA-Y>       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AOY    ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<OMEGA-Y>       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,185) .EQ. 1  .OR.  ISELEA(1,185) .EQ. 2) THEN
            CT = CH1680 (' OMEGA-Y F      ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,186) .EQ. 1  .OR.  ISELEA(1,186) .EQ. 2) THEN
         IF(ISELEP(1,186) .EQ. 1  .OR.  ISELEP(1,186) .EQ. 2) THEN
            CT = CH1680 ('<OMEGA-Y RMS> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AOYRG  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<OMEGA-Y RMS> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C
C                                 HIER: Z-KOMPONENTE DER VORTICITY
C                                 ------------------------------
C
      IF(ISELEA(1,188) .EQ. 1  .OR.  ISELEA(1,188) .EQ. 2) THEN
            CT = CH1680 (' OMEGA-Z        ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,189) .EQ. 1  .OR.  ISELEA(1,189) .EQ. 2) THEN
         IF(ISELEP(1,189) .EQ. 1  .OR.  ISELEP(1,189) .EQ. 2) THEN
            CT = CH1680 ('<OMEGA-Z>       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AOZ    ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<OMEGA-Z>       ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,191) .EQ. 1  .OR.  ISELEA(1,191) .EQ. 2) THEN
            CT = CH1680 (' OMEGA-Z F      ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,192) .EQ. 1  .OR.  ISELEA(1,192) .EQ. 2) THEN
         IF(ISELEP(1,192) .EQ. 1  .OR.  ISELEP(1,192) .EQ. 2) THEN
            CT = CH1680 ('<OMEGA-Z RMS> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AOZRG  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<OMEGA-Z RMS> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C
C                                 HIER: ENSTROPHIE
C                                 ----------------
C
      IF(ISELEA(1,195) .EQ. 1  .OR.  ISELEA(1,195) .EQ. 2) THEN
         IF(ISELEP(1,195) .EQ. 1  .OR.  ISELEP(1,195) .EQ. 2) THEN
            CT = CH1680 ('<ENSTROPHY>     ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AO2    ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<ENSTROPHY>     ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,197) .EQ. 1  .OR.  ISELEA(1,197) .EQ. 2) THEN
            CT = CH1680 (' ENSTROPHY-F G  ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,198) .EQ. 1  .OR.  ISELEA(1,198) .EQ. 2) THEN
         IF(ISELEP(1,198) .EQ. 1  .OR.  ISELEP(1,198) .EQ. 2) THEN
            CT = CH1680 ('<ENSTROPHY-RMS>G')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AO2RG  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<ENSTROPHY-RMS>G')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C
C                                 HIER: HELIZITAET
C                                 ----------------
C
      IF(ISELEA(1,200) .EQ. 1  .OR.  ISELEA(1,200) .EQ. 2) THEN
         IF(ISELEP(1,200) .EQ. 1  .OR.  ISELEP(1,200) .EQ. 2) THEN
            CT = CH1680 ('<HELICITY>      ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AHE    ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<HELICITY>      ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,202) .EQ. 1  .OR.  ISELEA(1,202) .EQ. 2) THEN
            CT = CH1680 (' HELICITY-F G   ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF
      IF(ISELEA(1,203) .EQ. 1  .OR.  ISELEA(1,203) .EQ. 2) THEN
         IF(ISELEP(1,203) .EQ. 1  .OR.  ISELEP(1,203) .EQ. 2) THEN
            CT = CH1680 ('<HELICITY-RMS> G')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AHERG  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<HELICITY-RMS> G')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C
C      
c
c     "Neue Statistik"
c
c
ccccccif defined 

      IF(ISELEA(1,400) .EQ. 1  .OR.  ISELEA(1,400) .EQ. 2) THEN
         IF(ISELEP(1,400) .EQ. 1  .OR.  ISELEP(1,400) .EQ. 2) THEN
            CT = CH1680 ('<UU> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUUM   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<UU> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,401) .EQ. 1  .OR.  ISELEA(1,401) .EQ. 2) THEN
         IF(ISELEP(1,401) .EQ. 1  .OR.  ISELEP(1,401) .EQ. 2) THEN
            CT = CH1680 ('<VV> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVVM   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<VV> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,402) .EQ. 1  .OR.  ISELEA(1,402) .EQ. 2) THEN
         IF(ISELEP(1,402) .EQ. 1  .OR.  ISELEP(1,402) .EQ. 2) THEN
            CT = CH1680 ('<WW> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AWWM   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<WW> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,406) .EQ. 1  .OR.  ISELEA(1,406) .EQ. 2) THEN
         IF(ISELEP(1,406) .EQ. 1  .OR.  ISELEP(1,406) .EQ. 2) THEN
            CT = CH1680 ('<PP> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   APPM   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<PP> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
cccccendif

cccccif defined 

      IF(ISELEA(1,403) .EQ. 1  .OR.  ISELEA(1,403) .EQ. 2) THEN
         IF(ISELEP(1,403) .EQ. 1  .OR.  ISELEP(1,403) .EQ. 2) THEN
            CT = CH1680 ('<UV> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUVM   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<UV> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,404) .EQ. 1  .OR.  ISELEA(1,404) .EQ. 2) THEN
         IF(ISELEP(1,404) .EQ. 1  .OR.  ISELEP(1,404) .EQ. 2) THEN
            CT = CH1680 ('<UW> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUWM   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<UW> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,405) .EQ. 1  .OR.  ISELEA(1,405) .EQ. 2) THEN
         IF(ISELEP(1,405) .EQ. 1  .OR.  ISELEP(1,405) .EQ. 2) THEN
            CT = CH1680 ('<VW> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVWM   ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<VW> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      
ccccendif

ccccif defined 

      IF(ISELEA(1,430) .EQ. 1  .OR.  ISELEA(1,430) .EQ. 2) THEN
         IF(ISELEP(1,430) .EQ. 1  .OR.  ISELEP(1,430) .EQ. 2) THEN
            CT = CH1680 ('<(dU/dX)^2> G   ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUXUXM,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<(dU/dX)^2> G   ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,431) .EQ. 1  .OR.  ISELEA(1,431) .EQ. 2) THEN
         IF(ISELEP(1,431) .EQ. 1  .OR.  ISELEP(1,431) .EQ. 2) THEN
            CT = CH1680 ('<(dU/dY)^2> G   ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUYUYM,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<(dU/dY)^2> G   ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,432) .EQ. 1  .OR.  ISELEA(1,432) .EQ. 2) THEN
         IF(ISELEP(1,432) .EQ. 1  .OR.  ISELEP(1,432) .EQ. 2) THEN
            CT = CH1680 ('<(dU/dZ)^2> G   ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUZUZM,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<(dU/dZ)^2> G   ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,433) .EQ. 1  .OR.  ISELEA(1,433) .EQ. 2) THEN
         IF(ISELEP(1,433) .EQ. 1  .OR.  ISELEP(1,433) .EQ. 2) THEN
            CT = CH1680 ('<(dV/dX)^2> G   ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVXVXM,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<(dV/dX)^2> G   ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,434) .EQ. 1  .OR.  ISELEA(1,434) .EQ. 2) THEN
         IF(ISELEP(1,434) .EQ. 1  .OR.  ISELEP(1,434) .EQ. 2) THEN
            CT = CH1680 ('<(dV/dY)^2> G   ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVYVYM,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<(dV/dY)^2> G   ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,435) .EQ. 1  .OR.  ISELEA(1,435) .EQ. 2) THEN
         IF(ISELEP(1,435) .EQ. 1  .OR.  ISELEP(1,435) .EQ. 2) THEN
            CT = CH1680 ('<(dV/dZ)^2> G   ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVZVZM,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<(dV/dZ)^2> G   ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,436) .EQ. 1  .OR.  ISELEA(1,436) .EQ. 2) THEN
         IF(ISELEP(1,436) .EQ. 1  .OR.  ISELEP(1,436) .EQ. 2) THEN
            CT = CH1680 ('<(dW/dX)^2> G   ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AWXWXM,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<(dW/dX)^2> G   ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,437) .EQ. 1  .OR.  ISELEA(1,437) .EQ. 2) THEN
         IF(ISELEP(1,437) .EQ. 1  .OR.  ISELEP(1,437) .EQ. 2) THEN
            CT = CH1680 ('<(dW/dY)^2> G   ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AWYWYM,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<(dW/dY)^2> G   ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,438) .EQ. 1  .OR.  ISELEA(1,438) .EQ. 2) THEN
         IF(ISELEP(1,438) .EQ. 1  .OR.  ISELEP(1,438) .EQ. 2) THEN
            CT = CH1680 ('<(dW/dZ)^2> G   ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AWZWZM,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<(dW/dZ)^2> G   ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

ccccendif

ccccif defined _STAT14_

      IF(ISELEA(1,410) .EQ. 1  .OR.  ISELEA(1,410) .EQ. 2) THEN
         IF(ISELEP(1,410) .EQ. 1  .OR.  ISELEP(1,410) .EQ. 2) THEN
            CT = CH1680 ('<UUU> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUUUM  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<UUU> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,411) .EQ. 1  .OR.  ISELEA(1,411) .EQ. 2) THEN
         IF(ISELEP(1,411) .EQ. 1  .OR.  ISELEP(1,411) .EQ. 2) THEN
            CT = CH1680 ('<VVV> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVVVM  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<VVV> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,412) .EQ. 1  .OR.  ISELEA(1,412) .EQ. 2) THEN
         IF(ISELEP(1,412) .EQ. 1  .OR.  ISELEP(1,412) .EQ. 2) THEN
            CT = CH1680 ('<WWW> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AWWWM  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<WWW> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,413) .EQ. 1  .OR.  ISELEA(1,413) .EQ. 2) THEN
         IF(ISELEP(1,413) .EQ. 1  .OR.  ISELEP(1,413) .EQ. 2) THEN
            CT = CH1680 ('<UUV> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUUVM  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<UUV> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,414) .EQ. 1  .OR.  ISELEA(1,414) .EQ. 2) THEN
         IF(ISELEP(1,414) .EQ. 1  .OR.  ISELEP(1,414) .EQ. 2) THEN
            CT = CH1680 ('<UUW> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUUWM  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<UUW> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      
      IF(ISELEA(1,415) .EQ. 1  .OR.  ISELEA(1,415) .EQ. 2) THEN
         IF(ISELEP(1,415) .EQ. 1  .OR.  ISELEP(1,415) .EQ. 2) THEN
            CT = CH1680 ('<VVU> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVVUM  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<VVU> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,416) .EQ. 1  .OR.  ISELEA(1,416) .EQ. 2) THEN
         IF(ISELEP(1,416) .EQ. 1  .OR.  ISELEP(1,416) .EQ. 2) THEN
            CT = CH1680 ('<VVW> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVVWM  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<VVW> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,417) .EQ. 1  .OR.  ISELEA(1,417) .EQ. 2) THEN
         IF(ISELEP(1,417) .EQ. 1  .OR.  ISELEP(1,417) .EQ. 2) THEN
            CT = CH1680 ('<WWU> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AWWUM  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<WWU> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,418) .EQ. 1  .OR.  ISELEA(1,418) .EQ. 2) THEN
         IF(ISELEP(1,418) .EQ. 1  .OR.  ISELEP(1,418) .EQ. 2) THEN
            CT = CH1680 ('<WWV> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AWWVM  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<WWV> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,419) .EQ. 1  .OR.  ISELEA(1,419) .EQ. 2) THEN
         IF(ISELEP(1,419) .EQ. 1  .OR.  ISELEP(1,419) .EQ. 2) THEN
            CT = CH1680 ('<UVW> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUVWM  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<UVW> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

c     ***** Geschwindigkeits-Druck-Term

      IF(ISELEA(1,420) .EQ. 1  .OR.  ISELEA(1,420) .EQ. 2) THEN
         IF(ISELEP(1,420) .EQ. 1  .OR.  ISELEP(1,420) .EQ. 2) THEN
            CT = CH1680 ('<UP> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUPM  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<UP> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,421) .EQ. 1  .OR.  ISELEA(1,421) .EQ. 2) THEN
         IF(ISELEP(1,421) .EQ. 1  .OR.  ISELEP(1,421) .EQ. 2) THEN
            CT = CH1680 ('<VP> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVPM  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<VP> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,422) .EQ. 1  .OR.  ISELEA(1,422) .EQ. 2) THEN
         IF(ISELEP(1,422) .EQ. 1  .OR.  ISELEP(1,422) .EQ. 2) THEN
            CT = CH1680 ('<WP> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AWPM  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<WP> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

c     ***** Druck-Scher-Korrelationen *****

      IF(ISELEA(1,423) .EQ. 1  .OR.  ISELEA(1,423) .EQ. 2) THEN
         IF(ISELEP(1,423) .EQ. 1  .OR.  ISELEP(1,423) .EQ. 2) THEN
            CT = CH1680 ('<(dU/dX)P> G    ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUXPM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<(dU/dX)P> G    ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,424) .EQ. 1  .OR.  ISELEA(1,424) .EQ. 2) THEN
         IF(ISELEP(1,424) .EQ. 1  .OR.  ISELEP(1,424) .EQ. 2) THEN
            CT = CH1680 ('<(dV/dY)P> G    ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVYPM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<(dV/dY)P> G    ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,425) .EQ. 1  .OR.  ISELEA(1,425) .EQ. 2) THEN
         IF(ISELEP(1,425) .EQ. 1  .OR.  ISELEP(1,425) .EQ. 2) THEN
            CT = CH1680 ('<(dW/dZ)P> G    ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AWZPM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<(dW/dZ)P> G    ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,426) .EQ. 1  .OR.  ISELEA(1,426) .EQ. 2) THEN
         IF(ISELEP(1,426) .EQ. 1  .OR.  ISELEP(1,426) .EQ. 2) THEN
            CT = CH1680 ('<(dU/dY+dV/dX)P>')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUYVXP ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<(dU/dY+dV/dX)P>')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,427) .EQ. 1  .OR.  ISELEA(1,427) .EQ. 2) THEN
         IF(ISELEP(1,427) .EQ. 1  .OR.  ISELEP(1,427) .EQ. 2) THEN
            CT = CH1680 ('<(dU/dZ+dW/dX)P>')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUZWXP ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<(dU/dZ+dW/dX)P>')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,428) .EQ. 1  .OR.  ISELEA(1,428) .EQ. 2) THEN
         IF(ISELEP(1,428) .EQ. 1  .OR.  ISELEP(1,428) .EQ. 2) THEN
            CT = CH1680 ('<(dV/dZ+dW/dY)P>')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVZWYP ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<(dV/dZ+dW/dY)P>')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

c     ***** Scher-Scher-Korrelationen *****

      IF(ISELEA(1,440) .EQ. 1  .OR.  ISELEA(1,440) .EQ. 2) THEN
         IF(ISELEP(1,440) .EQ. 1  .OR.  ISELEP(1,440) .EQ. 2) THEN
            CT = CH1680 ('<dU/dX dV/dX> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUXVXM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<dU/dX dV/dX> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,441) .EQ. 1  .OR.  ISELEA(1,441) .EQ. 2) THEN
         IF(ISELEP(1,441) .EQ. 1  .OR.  ISELEP(1,441) .EQ. 2) THEN
            CT = CH1680 ('<dU/dY dV/dY> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUYVYM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<dU/dY dV/dY> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,442) .EQ. 1  .OR.  ISELEA(1,442) .EQ. 2) THEN
         IF(ISELEP(1,442) .EQ. 1  .OR.  ISELEP(1,442) .EQ. 2) THEN
            CT = CH1680 ('<dU/dZ dV/dZ> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUZVZM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<dU/dZ dV/dZ> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,443) .EQ. 1  .OR.  ISELEA(1,443) .EQ. 2) THEN
         IF(ISELEP(1,443) .EQ. 1  .OR.  ISELEP(1,443) .EQ. 2) THEN
            CT = CH1680 ('<dU/dX dW/dX> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUXWXM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<dU/dX dW/dX> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,444) .EQ. 1  .OR.  ISELEA(1,444) .EQ. 2) THEN
         IF(ISELEP(1,444) .EQ. 1  .OR.  ISELEP(1,444) .EQ. 2) THEN
            CT = CH1680 ('<dU/dY dW/dY> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUYWYM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<dU/dY dW/dY> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,445) .EQ. 1  .OR.  ISELEA(1,445) .EQ. 2) THEN
         IF(ISELEP(1,445) .EQ. 1  .OR.  ISELEP(1,445) .EQ. 2) THEN
            CT = CH1680 ('<dU/dZ dW/dZ> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUZWZM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<dU/dZ dW/dZ> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,446) .EQ. 1  .OR.  ISELEA(1,446) .EQ. 2) THEN
         IF(ISELEP(1,446) .EQ. 1  .OR.  ISELEP(1,446) .EQ. 2) THEN
            CT = CH1680 ('<dV/dX dW/dX> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVXWXM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<dV/dX dW/dX> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,447) .EQ. 1  .OR.  ISELEA(1,447) .EQ. 2) THEN
         IF(ISELEP(1,447) .EQ. 1  .OR.  ISELEP(1,447) .EQ. 2) THEN
            CT = CH1680 ('<dV/dY dW/dY> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVYWYM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<dV/dY dW/dY> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,448) .EQ. 1  .OR.  ISELEA(1,448) .EQ. 2) THEN
         IF(ISELEP(1,448) .EQ. 1  .OR.  ISELEP(1,448) .EQ. 2) THEN
            CT = CH1680 ('<dV/dZ dW/dZ> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVZWZM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<dV/dZ dW/dZ> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
ccccendif

ccccifcdefinedcstat15

      IF(ISELEA(1,601) .EQ. 1  .OR.  ISELEA(1,601) .EQ. 2) THEN
         IF(ISELEP(1,601) .EQ. 1  .OR.  ISELEP(1,601) .EQ. 2) THEN
            CT = CH1680 ('<UT> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUTM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<UT> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,602) .EQ. 1  .OR.  ISELEA(1,602) .EQ. 2) THEN
         IF(ISELEP(1,602) .EQ. 1  .OR.  ISELEP(1,602) .EQ. 2) THEN
            CT = CH1680 ('<VT> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVTM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<VT> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,603) .EQ. 1  .OR.  ISELEA(1,603) .EQ. 2) THEN
         IF(ISELEP(1,603) .EQ. 1  .OR.  ISELEP(1,603) .EQ. 2) THEN
            CT = CH1680 ('<WT> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AWTM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<WT> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,604) .EQ. 1  .OR.  ISELEA(1,604) .EQ. 2) THEN
         IF(ISELEP(1,604) .EQ. 1  .OR.  ISELEP(1,604) .EQ. 2) THEN
            CT = CH1680 ('<TT> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ATTM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<TT> G          ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,605) .EQ. 1  .OR.  ISELEA(1,605) .EQ. 2) THEN
         IF(ISELEP(1,604) .EQ. 1  .OR.  ISELEP(1,605) .EQ. 2) THEN
            CT = CH1680 ('<SIJSIJ> G      ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ASSM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<SIJSIJ> G      ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
C Dieses Feld wird nie eingelesen da nur aus mittleren Geschwindigkeiten Berenchnet
      IF(ISELEA(1,606) .EQ. 1  .OR.  ISELEA(1,606) .EQ. 2) THEN
            CT = CH1680 ('<SIJ><SIJ> G    ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
      ENDIF

ccccendif _STAT15_
ccccif STAT17
      IF(ISELEA(1,607) .EQ. 1  .OR.  ISELEA(1,607) .EQ. 2) THEN
         IF(ISELEP(1,607) .EQ. 1  .OR.  ISELEP(1,607) .EQ. 2) THEN
            CT = CH1680 ('<UTT> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUTTM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<UTT> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,608) .EQ. 1  .OR.  ISELEA(1,608) .EQ. 2) THEN
         IF(ISELEP(1,608) .EQ. 1  .OR.  ISELEP(1,608) .EQ. 2) THEN
            CT = CH1680 ('<VTT> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVTTM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<VTT> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,609) .EQ. 1  .OR.  ISELEA(1,609) .EQ. 2) THEN
         IF(ISELEP(1,609) .EQ. 1  .OR.  ISELEP(1,609) .EQ. 2) THEN
            CT = CH1680 ('<WTT> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AWTTM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<WTT> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,610) .EQ. 1  .OR.  ISELEA(1,610) .EQ. 2) THEN
         IF(ISELEP(1,610) .EQ. 1  .OR.  ISELEP(1,610) .EQ. 2) THEN
            CT = CH1680 ('<dT/dX dT/dX> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ATXTXM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<dT/dX dT/dX> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,611) .EQ. 1  .OR.  ISELEA(1,611) .EQ. 2) THEN
         IF(ISELEP(1,611) .EQ. 1  .OR.  ISELEP(1,611) .EQ. 2) THEN
            CT = CH1680 ('<dT/dY dT/dY> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ATYTYM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<dT/dY dT/dY> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
      IF(ISELEA(1,612) .EQ. 1  .OR.  ISELEA(1,612) .EQ. 2) THEN
         IF(ISELEP(1,612) .EQ. 1  .OR.  ISELEP(1,612) .EQ. 2) THEN
            CT = CH1680 ('<dT/dZ dT/dZ> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ATZTZM ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<dT/dZ dT/dZ> G ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF


C-------------------------------------------------------------------------
ccccifdef _STAT20_

      IF(ISELEA(1,490) .EQ. 1  .OR.  ISELEA(1,490) .EQ. 2) THEN
         IF(ISELEP(1,490) .EQ. 1  .OR.  ISELEP(1,490) .EQ. 2) THEN
            CT = CH1680 ('<UUUU> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AUUUUM  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<UUUU> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,491) .EQ. 1  .OR.  ISELEA(1,491) .EQ. 2) THEN
         IF(ISELEP(1,491) .EQ. 1  .OR.  ISELEP(1,491) .EQ. 2) THEN
            CT = CH1680 ('<VVVV> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AVVVVM  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<VVVV> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,492) .EQ. 1  .OR.  ISELEA(1,492) .EQ. 2) THEN
         IF(ISELEP(1,492) .EQ. 1  .OR.  ISELEP(1,492) .EQ. 2) THEN
            CT = CH1680 ('<WWWW> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   AWWWWM  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<WWWW> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,493) .EQ. 1  .OR.  ISELEA(1,493) .EQ. 2) THEN
         IF(ISELEP(1,493) .EQ. 1  .OR.  ISELEP(1,493) .EQ. 2) THEN
            CT = CH1680 ('<PPPP> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   APPPPM  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<PPPP> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF
Cccccifdef _TSCAL_
      IF(ISELEA(1,521) .EQ. 1  .OR.  ISELEA(1,521) .EQ. 2) THEN
         IF(ISELEP(1,521) .EQ. 1  .OR.  ISELEP(1,521) .EQ. 2) THEN
            CT = CH1680 ('<TTTT> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ATTTTM  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<TTTT> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

      IF(ISELEA(1,518) .EQ. 1  .OR.  ISELEA(1,518) .EQ. 2) THEN
         IF(ISELEP(1,518) .EQ. 1  .OR.  ISELEP(1,518) .EQ. 2) THEN
            CT = CH1680 ('<TTT> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   ATTTM  ,CT,NBND,HILF,IGRID,IDIM3D)
         ELSE
            CT = CH1680 ('<TTT> G         ')
            CALL DIBC   (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   HILF   ,CT,NBND,HILF,IGRID,IDIM3D)
         ENDIF
      ENDIF

Cccccendif
CCCCcendif


C
      RETURN
 6010 FORMAT (5(I9,1X))
 6015 FORMAT (4(I19,1X))
 6016 FORMAT (  I19    )
 6020 FORMAT (4(E19.12E3,1X))
 6021 FORMAT (  E19.12E3    )
      END
