










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
      SUBROUTINE STMIMP  (DREAD,DCONT,NPRNEU,FPRNEU
     $     ,DX,DY,DZ,DDX,DDY,DDZ,X,Y,Z,

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
C        S T M I M P      IN STMIMP WERDEN DIE ALTEN ENSEMBLE-
C                         MITTELWERTE VERBESSERT. DIE GEWICHTUNG DER
C                         ALTEN (WAEHREND DER VORANGEGANGENEN LAEUFE
C                         ERZEUGT) UND DER NEUEN ENSEMBLE-MITTELWERTE
C                         (WAEHREND DES MOMENTANEN LAUFES ERZEUGT)
C                         ERGIBT SICH AUS DEM VERHAELTNIS DER ANZAHL
C                         DER STICHPROBEN, DIE ZU DEN JEWEILIGEN
C                         MITTELWERTEN FUEHRTEN.
C
C                         STMIMP DIENT DAZU, DIE FELDER AUSZUWAEHLEN,
C                         DEREN ENSEMBLE-MITTELWERT VERBESSERT WERDEN
C                         SOLL. DIE EIGENTLICHE ARITHMET. OPERATION
C                         FINDET IN SUBR. STMNEW STATT.
C  A C H T U N G:         DIESE ROUTINE MUSS ERWEITERT WERDEN, FALLS
C                         DIE STATIST. AUSWERTUNG ERWEITERT WIRD !
C*STARLET***************************************************************
C
C PARAM: DREAD          - DREAD = .T. : DIE WAEHREND DES VORANGEGANGENEN
C                         LAUFES ERZEUGTEN ENSEMBLE-MITTELWERTE WERDEN
C                         EINGELESEN.
C                         DREAD = .F. : ES FINDET EINE VORBELEGUNG MIT
C                         "BESTMOEGLICHEN" ENSEMBLE-MITTELWERTEN STATT.
C        DCONT          - DCONT = .T. .AND. DREAD = .T. : DIE GEWICH-
C                         TUNG DER "ALTEN" ENSEMBLE-MITTELWERTE ERFOLGT
C                         ENTSPRECHEND DER ANZAHL DER "ALTEN" STICH-
C                         PROBEN.
C
C                         DCONT = .F. .AND. DREAD = .T. : DIE "ALTEN"
C                         ENSEMBLE-MITTELWERTE WERDEN MIT NPRNEU
C                         STICHPROBEN GEWICHTET.
C        NPRNEU         - SIEHE BESCHREIBUNG  "DCONT"
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
C UPROG                 : STMNEW, STMNRM
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        29.08.86 (HW)  : ORIGINAL
C        03.09.86 (HW)  : KORREKTUREN
C        10.09.86 (HW)  : KORREKTUREN
C        12.09.86 (HW)  : SONDERBEHANDLUNG FUER DIE RMS-WERTE
C                         (STMNRM EINGEFUEHRT)
C        01.07.88 (HW)  : ERWEITERUNG: AUSWERTUNG DER DREI KOMPONENTEN
C                         DER VORTICITY, DER ENSTROPHIE UND DER
C                         HELIZITAET
C        12.02.03 (TB)  : SCALAR FIELD STATISTICS IMPLEMENTED
C
C*STARLET***************************************************************
C
      LOGICAL  DREAD,  DCONT,XHOMOG,YHOMOG,ZHOMOG
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
C                                 STATIST. GROESSEN DER U-KOMPONENTE
C                                 ----------------------------------
C
C                                 VERBESSERUNG DER ALTEN <U>-WERTE
C
      IF(ISELEP(1,  4) .EQ. 1 .OR. ISELEP(1,  4) .EQ. 2) THEN
         NPROLD = ISAMPA(    4)
         IF(ISAMPP(    4) .GT. 0) THEN
            NPRRUN = ISAMPP(    4)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AU     ,SU     )
         ENDIF
         ISAMPA(    4) = NPROLD
         ISAMPP(    4) = ISAMPP(    4) + NPROLD
      ENDIF
C
C                                 VERBESSERUNG DER ALTEN <U-RMS>-WERTE
C
      IF(ISELEP(1,  7) .EQ. 1 .OR. ISELEP(1,  7) .EQ. 2) THEN
         NPROLD = ISAMPA(    7)
         IF(ISAMPP(    7) .GT. 0) THEN
            NPRRUN = ISAMPP(    7)
            CALL STMNRM  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AURG   ,SURG   )
         ENDIF
         ISAMPA(    7) = NPROLD
         ISAMPP(    7) = ISAMPP(    7) + NPROLD
      ENDIF
C
C                                 VERBESSERUNG DER ALTEN <U-SKE>-WERTE
C
      IF(ISELEP(1, 10) .EQ. 1 .OR. ISELEP(1, 10) .EQ. 2) THEN
         NPROLD = ISAMPA(   10)
         IF(ISAMPP(   10) .GT. 0) THEN
            NPRRUN = ISAMPP(   10)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AUSKG  ,SUSKG  )
         ENDIF
         ISAMPA(   10) = NPROLD
         ISAMPP(   10) = ISAMPP(   10) + NPROLD
      ENDIF
C
C                                 VERBESSERUNG DER ALTEN <U-FLA>-WERTE
C
      IF(ISELEP(1, 12) .EQ. 1 .OR. ISELEP(1, 12) .EQ. 2) THEN
         NPROLD = ISAMPA(   12)
         IF(ISAMPP(   12) .GT. 0) THEN
            NPRRUN = ISAMPP(   12)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AUFLG  ,SUFLG  )
         ENDIF
         ISAMPA(   12) = NPROLD
         ISAMPP(   12) = ISAMPP(   12) + NPROLD
      ENDIF
C
C                                 STATIST. GROESSEN DER V-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 29) .EQ. 1 .OR. ISELEP(1, 29) .EQ. 2) THEN
         NPROLD = ISAMPA(   29)
         IF(ISAMPP(   29) .GT. 0) THEN
            NPRRUN = ISAMPP(   29)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AV     ,SV     )
         ENDIF
         ISAMPA(   29) = NPROLD
         ISAMPP(   29) = ISAMPP(   29) + NPROLD
      ENDIF
      IF(ISELEP(1, 32) .EQ. 1 .OR. ISELEP(1, 32) .EQ. 2) THEN
         NPROLD = ISAMPA(   32)
         IF(ISAMPP(   32) .GT. 0) THEN
            NPRRUN = ISAMPP(   32)
            CALL STMNRM  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AVRG   ,SVRG   )
         ENDIF
         ISAMPA(   32) = NPROLD
         ISAMPP(   32) = ISAMPP(   32) + NPROLD
      ENDIF
      IF(ISELEP(1, 35) .EQ. 1 .OR. ISELEP(1, 35) .EQ. 2) THEN
         NPROLD = ISAMPA(   35)
         IF(ISAMPP(   35) .GT. 0) THEN
            NPRRUN = ISAMPP(   35)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AVSKG  ,SVSKG  )
         ENDIF
         ISAMPA(   35) = NPROLD
         ISAMPP(   35) = ISAMPP(   35) + NPROLD
      ENDIF
      IF(ISELEP(1, 37) .EQ. 1 .OR. ISELEP(1, 37) .EQ. 2) THEN
         NPROLD = ISAMPA(   37)
         IF(ISAMPP(   37) .GT. 0) THEN
            NPRRUN = ISAMPP(   37)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AVFLG  ,SVFLG  )
         ENDIF
         ISAMPA(   37) = NPROLD
         ISAMPP(   37) = ISAMPP(   37) + NPROLD
      ENDIF
C
C                                 STATIST. GROESSEN DER W-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 54) .EQ. 1 .OR. ISELEP(1, 54) .EQ. 2) THEN
         NPROLD = ISAMPA(   54)
         IF(ISAMPP(   54) .GT. 0) THEN
            NPRRUN = ISAMPP(   54)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AW     ,SW     )
         ENDIF
         ISAMPA(   54) = NPROLD
         ISAMPP(   54) = ISAMPP(   54) + NPROLD
      ENDIF
      IF(ISELEP(1, 57) .EQ. 1 .OR. ISELEP(1, 57) .EQ. 2) THEN
         NPROLD = ISAMPA(   57)
         IF(ISAMPP(   57) .GT. 0) THEN
            NPRRUN = ISAMPP(   57)
            CALL STMNRM  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AWRG   ,SWRG   )
         ENDIF
         ISAMPA(   57) = NPROLD
         ISAMPP(   57) = ISAMPP(   57) + NPROLD
      ENDIF
      IF(ISELEP(1, 60) .EQ. 1 .OR. ISELEP(1, 60) .EQ. 2) THEN
         NPROLD = ISAMPA(   60)
         IF(ISAMPP(   60) .GT. 0) THEN
            NPRRUN = ISAMPP(   60)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AWSKG  ,SWSKG  )
         ENDIF
         ISAMPA(   60) = NPROLD
         ISAMPP(   60) = ISAMPP(   60) + NPROLD
      ENDIF
      IF(ISELEP(1, 62) .EQ. 1 .OR. ISELEP(1, 62) .EQ. 2) THEN
         NPROLD = ISAMPA(   62)
         IF(ISAMPP(   62) .GT. 0) THEN
            NPRRUN = ISAMPP(   62)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AWFLG  ,SWFLG  )
         ENDIF
         ISAMPA(   62) = NPROLD
         ISAMPP(   62) = ISAMPP(   62) + NPROLD
      ENDIF
C
C                                 STATIST. GROESSEN DER P-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 79) .EQ. 1 .OR. ISELEP(1, 79) .EQ. 2) THEN
         NPROLD = ISAMPA(   79)
         IF(ISAMPP(   79) .GT. 0) THEN
            NPRRUN = ISAMPP(   79)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AP     ,SP     )
         ENDIF
         ISAMPA(   79) = NPROLD
         ISAMPP(   79) = ISAMPP(   79) + NPROLD
      ENDIF
      IF(ISELEP(1, 82) .EQ. 1 .OR. ISELEP(1, 82) .EQ. 2) THEN
         NPROLD = ISAMPA(   82)
         IF(ISAMPP(   82) .GT. 0) THEN
            NPRRUN = ISAMPP(   82)
            CALL STMNRM  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,APRG   ,SPRG   )
         ENDIF
         ISAMPA(   82) = NPROLD
         ISAMPP(   82) = ISAMPP(   82) + NPROLD
      ENDIF
      IF(ISELEP(1, 85) .EQ. 1 .OR. ISELEP(1, 85) .EQ. 2) THEN
         NPROLD = ISAMPA(   85)
         IF(ISAMPP(   85) .GT. 0) THEN
            NPRRUN = ISAMPP(   85)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,APSKG  ,SPSKG  )
         ENDIF
         ISAMPA(   85) = NPROLD
         ISAMPP(   85) = ISAMPP(   85) + NPROLD
      ENDIF
      IF(ISELEP(1, 87) .EQ. 1 .OR. ISELEP(1, 87) .EQ. 2) THEN
         NPROLD = ISAMPA(   87)
         IF(ISAMPP(   87) .GT. 0) THEN
            NPRRUN = ISAMPP(   87)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,APFLG  ,SPFLG  )
         ENDIF
         ISAMPA(   87) = NPROLD
         ISAMPP(   87) = ISAMPP(   87) + NPROLD
      ENDIF
C
C                                 KINETISCHE ENERGIE DER SCHWANKUNGS-
C                                 GESCHWINDIGKEITEN
C                                 -----------------------------------
C
      IF(ISELEP(1,106) .EQ. 1 .OR. ISELEP(1,106) .EQ. 2) THEN
         NPROLD = ISAMPA(  106)
         IF(ISAMPP(  106) .GT. 0) THEN
            NPRRUN = ISAMPP(  106)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AEFG   ,SEFG   )
         ENDIF
         ISAMPA(  106) = NPROLD
         ISAMPP(  106) = ISAMPP(  106) + NPROLD
      ENDIF
      IF(ISELEP(1,109) .EQ. 1 .OR. ISELEP(1,109) .EQ. 2) THEN
         NPROLD = ISAMPA(  109)
         IF(ISAMPP(  109) .GT. 0) THEN
            NPRRUN = ISAMPP(  109)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AEFS   ,SEFS   )
         ENDIF
         ISAMPA(  109) = NPROLD
         ISAMPP(  109) = ISAMPP(  109) + NPROLD
      ENDIF
C
C                                 DISSIPATION
C                                 -----------
      IF(ISELEP(1,120) .EQ. 1 .OR. ISELEP(1,120) .EQ. 2) THEN
         NPROLD = ISAMPA(  120)
         IF(ISAMPP(  120) .GT. 0) THEN
            NPRRUN = ISAMPP(  120)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,ADFG   ,SDFG   )
         ENDIF
         ISAMPA(  120) = NPROLD
         ISAMPP(  120) = ISAMPP(  120) + NPROLD
      ENDIF
C
C                               FUER   TAYLOR''SCHES MIKROMASS
C                               QUADRATE DER ELEMENTE DES 
C                               DEFORMATIONSGESCHWINDIGKEITSTENSORS
C                               -----------------------------------
C                                  U_KOMPONENTE IN X_RICHTUNG
      IF(ISELEP(1,131) .EQ. 1 .OR. ISELEP(1,131) .EQ. 2) THEN
         NPROLD = ISAMPA(  131)
         IF(ISAMPP(  131) .GT. 0) THEN
            NPRRUN = ISAMPP(  131)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,ADUDX2 ,SDUDX2 )
         ENDIF
         ISAMPA(  131) = NPROLD
         ISAMPP(  131) = ISAMPP(  131) + NPROLD
      ENDIF
C                                  U_KOMPONENTE IN Y_RICHTUNG
      IF(ISELEP(1,132) .EQ. 1 .OR. ISELEP(1,132) .EQ. 2) THEN
         NPROLD = ISAMPA(  132)
         IF(ISAMPP(  132) .GT. 0) THEN
            NPRRUN = ISAMPP(  132)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,ADUDY2 ,SDUDY2 )
         ENDIF
         ISAMPA(  132) = NPROLD
         ISAMPP(  132) = ISAMPP(  132) + NPROLD
      ENDIF
C                                  U_KOMPONENTE IN Z_RICHTUNG
      IF(ISELEP(1,133) .EQ. 1 .OR. ISELEP(1,133) .EQ. 2) THEN
         NPROLD = ISAMPA(  133)
         IF(ISAMPP(  133) .GT. 0) THEN
            NPRRUN = ISAMPP(  133)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,ADUDZ2 ,SDUDZ2 )
         ENDIF
         ISAMPA(  133) = NPROLD
         ISAMPP(  133) = ISAMPP(  133) + NPROLD
      ENDIF
C                                  V_KOMPONENTE IN X_RICHTUNG
      IF(ISELEP(1,134) .EQ. 1 .OR. ISELEP(1,134) .EQ. 2) THEN
         NPROLD = ISAMPA(  134)
         IF(ISAMPP(  134) .GT. 0) THEN
            NPRRUN = ISAMPP(  134)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,ADVDX2 ,SDVDX2 )
         ENDIF
         ISAMPA(  134) = NPROLD
         ISAMPP(  134) = ISAMPP(  134) + NPROLD
      ENDIF
C                                  V_KOMPONENTE IN Y_RICHTUNG
      IF(ISELEP(1,135) .EQ. 1 .OR. ISELEP(1,135) .EQ. 2) THEN
         NPROLD = ISAMPA(  135)
         IF(ISAMPP(  135) .GT. 0) THEN
            NPRRUN = ISAMPP(  135)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,ADVDY2 ,SDVDY2 )
         ENDIF
         ISAMPA(  135) = NPROLD
         ISAMPP(  135) = ISAMPP(  135) + NPROLD
      ENDIF
C                                  V_KOMPONENTE IN Z_RICHTUNG
      IF(ISELEP(1,136) .EQ. 1 .OR. ISELEP(1,136) .EQ. 2) THEN
         NPROLD = ISAMPA(  136)
         IF(ISAMPP(  136) .GT. 0) THEN
            NPRRUN = ISAMPP(  136)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,ADVDZ2 ,SDVDZ2 )
         ENDIF
         ISAMPA(  136) = NPROLD
         ISAMPP(  136) = ISAMPP(  136) + NPROLD
      ENDIF
C                                  W_KOMPONENTE IN X_RICHTUNG
      IF(ISELEP(1,137) .EQ. 1 .OR. ISELEP(1,137) .EQ. 2) THEN
         NPROLD = ISAMPA(  137)
         IF(ISAMPP(  137) .GT. 0) THEN
            NPRRUN = ISAMPP(  137)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,ADWDX2 ,SDWDX2 )
         ENDIF
         ISAMPA(  137) = NPROLD
         ISAMPP(  137) = ISAMPP(  137) + NPROLD
      ENDIF
C                                  W_KOMPONENTE IN Y_RICHTUNG
      IF(ISELEP(1,138) .EQ. 1 .OR. ISELEP(1,138) .EQ. 2) THEN
         NPROLD = ISAMPA(  138)
         IF(ISAMPP(  138) .GT. 0) THEN
            NPRRUN = ISAMPP(  138)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,ADWDY2 ,SDWDY2 )
         ENDIF
         ISAMPA(  138) = NPROLD
         ISAMPP(  138) = ISAMPP(  138) + NPROLD
      ENDIF
C                                  W_KOMPONENTE IN Z_RICHTUNG
      IF(ISELEP(1,139) .EQ. 1 .OR. ISELEP(1,139) .EQ. 2) THEN
         NPROLD = ISAMPA(  139)
         IF(ISAMPP(  139) .GT. 0) THEN
            NPRRUN = ISAMPP(  139)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,ADWDZ2 ,SDWDZ2 )
         ENDIF
         ISAMPA(  139) = NPROLD
         ISAMPP(  139) = ISAMPP(  139) + NPROLD
      ENDIF
C
C                                 KOMPONENTEN DES SPANNUNGSTENSORS
C                                 --------------------------------
C
C                                 HIER: U - W
C                                 -----------
C
      IF(ISELEP(1,146) .EQ. 1 .OR. ISELEP(1,146) .EQ. 2) THEN
         NPROLD = ISAMPA(  146)
         IF(ISAMPP(  146) .GT. 0) THEN
            NPRRUN = ISAMPP(  146)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AUFWFG ,SUFWFG )
         ENDIF
         ISAMPA(  146) = NPROLD
         ISAMPP(  146) = ISAMPP(  146) + NPROLD
      ENDIF
      IF(ISELEP(1,149) .EQ. 1 .OR. ISELEP(1,149) .EQ. 2) THEN
         NPROLD = ISAMPA(  149)
         IF(ISAMPP(  149) .GT. 0) THEN
            NPRRUN = ISAMPP(  149)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AUFWFS ,SUFWFS )
         ENDIF
         ISAMPA(  149) = NPROLD
         ISAMPP(  149) = ISAMPP(  149) + NPROLD
      ENDIF
      IF(ISELEP(1,151) .EQ. 1 .OR. ISELEP(1,151) .EQ. 2) THEN
         NPROLD = ISAMPA(  151)
         IF(ISAMPP(  151) .GT. 0) THEN
            NPRRUN = ISAMPP(  151)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AUFWFM ,SUFWFM )
         ENDIF
         ISAMPA(  151) = NPROLD
         ISAMPP(  151) = ISAMPP(  151) + NPROLD
      ENDIF
C
C                                 HIER: V - W
C                                 -----------
C
      IF(ISELEP(1,156) .EQ. 1 .OR. ISELEP(1,156) .EQ. 2) THEN
         NPROLD = ISAMPA(  156)
         IF(ISAMPP(  156) .GT. 0) THEN
            NPRRUN = ISAMPP(  156)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AVFWFG ,SVFWFG )
         ENDIF
         ISAMPA(  156) = NPROLD
         ISAMPP(  156) = ISAMPP(  156) + NPROLD
      ENDIF
      IF(ISELEP(1,159) .EQ. 1 .OR. ISELEP(1,159) .EQ. 2) THEN
         NPROLD = ISAMPA(  159)
         IF(ISAMPP(  159) .GT. 0) THEN
            NPRRUN = ISAMPP(  159)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AVFWFS ,SVFWFS )
         ENDIF
         ISAMPA(  159) = NPROLD
         ISAMPP(  159) = ISAMPP(  159) + NPROLD
      ENDIF
      IF(ISELEP(1,161) .EQ. 1 .OR. ISELEP(1,161) .EQ. 2) THEN
         NPROLD = ISAMPA(  161)
         IF(ISAMPP(  161) .GT. 0) THEN
            NPRRUN = ISAMPP(  161)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AVFWFM ,SVFWFM )
         ENDIF
         ISAMPA(  161) = NPROLD
         ISAMPP(  161) = ISAMPP(  161) + NPROLD
      ENDIF
C
C                                 HIER: U - V
C                                 -----------
C
      IF(ISELEP(1,166) .EQ. 1 .OR. ISELEP(1,166) .EQ. 2) THEN
         NPROLD = ISAMPA(  166)
         IF(ISAMPP(  166) .GT. 0) THEN
            NPRRUN = ISAMPP(  166)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AUFVFG ,SUFVFG )
         ENDIF
         ISAMPA(  166) = NPROLD
         ISAMPP(  166) = ISAMPP(  166) + NPROLD
      ENDIF
      IF(ISELEP(1,169) .EQ. 1 .OR. ISELEP(1,169) .EQ. 2) THEN
         NPROLD = ISAMPA(  169)
         IF(ISAMPP(  169) .GT. 0) THEN
            NPRRUN = ISAMPP(  169)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AUFVFS ,SUFVFS )
         ENDIF
         ISAMPA(  169) = NPROLD
         ISAMPP(  169) = ISAMPP(  169) + NPROLD
      ENDIF
      IF(ISELEP(1,171) .EQ. 1 .OR. ISELEP(1,171) .EQ. 2) THEN
         NPROLD = ISAMPA(  171)
         IF(ISAMPP(  171) .GT. 0) THEN
            NPRRUN = ISAMPP(  171)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AUFVFM ,SUFVFM )
         ENDIF
         ISAMPA(  171) = NPROLD
         ISAMPP(  171) = ISAMPP(  171) + NPROLD
      ENDIF
C
C
C                                 X-KOMPONENTE DER VORTICITY
C                                 ------------------------
C
      IF(ISELEP(1,177) .EQ. 1 .OR. ISELEP(1,177) .EQ. 2) THEN
         NPROLD = ISAMPA(  177)
         IF(ISAMPP(  177) .GT. 0) THEN
            NPRRUN = ISAMPP(  177)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AOX    ,SOX    )
         ENDIF
         ISAMPA(  177) = NPROLD
         ISAMPP(  177) = ISAMPP(  177) + NPROLD
      ENDIF
      IF(ISELEP(1,180) .EQ. 1 .OR. ISELEP(1,180) .EQ. 2) THEN
         NPROLD = ISAMPA(  180)
         IF(ISAMPP(  180) .GT. 0) THEN
            NPRRUN = ISAMPP(  180)
            CALL STMNRM  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AOXRG  ,SOXRG  )
         ENDIF
         ISAMPA(  180) = NPROLD
         ISAMPP(  180) = ISAMPP(  180) + NPROLD
      ENDIF
C
C                                 Y-KOMPONENTE DER VORTICITY
C                                 ------------------------
C
      IF(ISELEP(1,183) .EQ. 1 .OR. ISELEP(1,183) .EQ. 2) THEN
         NPROLD = ISAMPA(  183)
         IF(ISAMPP(  183) .GT. 0) THEN
            NPRRUN = ISAMPP(  183)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AOY    ,SOY    )
         ENDIF
         ISAMPA(  183) = NPROLD
         ISAMPP(  183) = ISAMPP(  183) + NPROLD
      ENDIF
      IF(ISELEP(1,186) .EQ. 1 .OR. ISELEP(1,186) .EQ. 2) THEN
         NPROLD = ISAMPA(  186)
         IF(ISAMPP(  186) .GT. 0) THEN
            NPRRUN = ISAMPP(  186)
            CALL STMNRM  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AOYRG  ,SOYRG  )
         ENDIF
         ISAMPA(  186) = NPROLD
         ISAMPP(  186) = ISAMPP(  186) + NPROLD
      ENDIF
C
C                                 Z-KOMPONENTE DER VORTICITY
C                                 ------------------------
C
      IF(ISELEP(1,189) .EQ. 1 .OR. ISELEP(1,189) .EQ. 2) THEN
         NPROLD = ISAMPA(  189)
         IF(ISAMPP(  189) .GT. 0) THEN
            NPRRUN = ISAMPP(  189)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AOZ    ,SOZ    )
         ENDIF
         ISAMPA(  189) = NPROLD
         ISAMPP(  189) = ISAMPP(  189) + NPROLD
      ENDIF
      IF(ISELEP(1,192) .EQ. 1 .OR. ISELEP(1,192) .EQ. 2) THEN
         NPROLD = ISAMPA(  192)
         IF(ISAMPP(  192) .GT. 0) THEN
            NPRRUN = ISAMPP(  192)
            CALL STMNRM  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AOZRG  ,SOZRG  )
         ENDIF
         ISAMPA(  192) = NPROLD
         ISAMPP(  192) = ISAMPP(  192) + NPROLD
      ENDIF
C
C                                 ENSTROPHIE
C                                 ----------
C
      IF(ISELEP(1,195) .EQ. 1 .OR. ISELEP(1,195) .EQ. 2) THEN
         NPROLD = ISAMPA(  195)
         IF(ISAMPP(  195) .GT. 0) THEN
            NPRRUN = ISAMPP(  195)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AO2    ,SO2    )
         ENDIF
         ISAMPA(  195) = NPROLD
         ISAMPP(  195) = ISAMPP(  195) + NPROLD
      ENDIF
      IF(ISELEP(1,198) .EQ. 1 .OR. ISELEP(1,198) .EQ. 2) THEN
         NPROLD = ISAMPA(  198)
         IF(ISAMPP(  198) .GT. 0) THEN
            NPRRUN = ISAMPP(  198)
            CALL STMNRM  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AO2RG  ,SO2RG  )
         ENDIF
         ISAMPA(  198) = NPROLD
         ISAMPP(  198) = ISAMPP(  198) + NPROLD
      ENDIF
C
C                                 HELIZITAET
C                                 ----------
C
      IF(ISELEP(1,200) .EQ. 1 .OR. ISELEP(1,200) .EQ. 2) THEN
         NPROLD = ISAMPA(  200)
         IF(ISAMPP(  200) .GT. 0) THEN
            NPRRUN = ISAMPP(  200)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AHE    ,SHE    )
         ENDIF
         ISAMPA(  200) = NPROLD
         ISAMPP(  200) = ISAMPP(  200) + NPROLD
      ENDIF
      IF(ISELEP(1,203) .EQ. 1 .OR. ISELEP(1,203) .EQ. 2) THEN
         NPROLD = ISAMPA(  203)
         IF(ISAMPP(  203) .GT. 0) THEN
            NPRRUN = ISAMPP(  203)
            CALL STMNRM  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AHERG  ,SHERG  )
         ENDIF
         ISAMPA(  203) = NPROLD
         ISAMPP(  203) = ISAMPP(  203) + NPROLD
      ENDIF
C
C

      IF(ISELEP(1,400) .GE. 1 ) THEN
         NPROLD = ISAMPA(  400)
         IF(ISAMPP(  400) .GT. 0) THEN
            NPRRUN = ISAMPP(  400)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AUUM ,SUUM )
         ENDIF
         ISAMPA(  400) = NPROLD
         ISAMPP(  400) = ISAMPP(  400) + NPROLD
      ENDIF

      IF(ISELEP(1,401) .GE. 1 ) THEN
         NPROLD = ISAMPA(  401)
         IF(ISAMPP(  401) .GT. 0) THEN
            NPRRUN = ISAMPP(  401)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AVVM ,SVVM )
         ENDIF
         ISAMPA(  401) = NPROLD
         ISAMPP(  401) = ISAMPP(  401) + NPROLD
      ENDIF

      IF(ISELEP(1,402) .GE. 1 ) THEN
         NPROLD = ISAMPA(  402)
         IF(ISAMPP(  402) .GT. 0) THEN
            NPRRUN = ISAMPP(  402)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AWWM ,SWWM )
         ENDIF
         ISAMPA(  402) = NPROLD
         ISAMPP(  402) = ISAMPP(  402) + NPROLD
      ENDIF

      IF(ISELEP(1,406) .GE. 1 ) THEN
         NPROLD = ISAMPA(  406)
         IF(ISAMPP(  406) .GT. 0) THEN
            NPRRUN = ISAMPP(  406)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,APPM ,SPPM )
         ENDIF
         ISAMPA(  406) = NPROLD
         ISAMPP(  406) = ISAMPP(  406) + NPROLD
      ENDIF


      IF(ISELEP(1,403) .GE. 1 ) THEN
         NPROLD = ISAMPA(  403)
         IF(ISAMPP(  403) .GT. 0) THEN
            NPRRUN = ISAMPP(  403)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AUVM ,SUVM )
         ENDIF
         ISAMPA(  403) = NPROLD
         ISAMPP(  403) = ISAMPP(  403) + NPROLD
      ENDIF

      IF(ISELEP(1,404) .GE. 1 ) THEN
         NPROLD = ISAMPA(  404)
         IF(ISAMPP(  404) .GT. 0) THEN
            NPRRUN = ISAMPP(  404)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AUWM ,SUWM )
         ENDIF
         ISAMPA(  404) = NPROLD
         ISAMPP(  404) = ISAMPP(  404) + NPROLD
      ENDIF

      IF(ISELEP(1,405) .GE. 1 ) THEN
         NPROLD = ISAMPA(  405)
         IF(ISAMPP(  405) .GT. 0) THEN
            NPRRUN = ISAMPP(  405)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $                    DREAD,DCONT,NPRNEU,FPRNEU,AVWM ,SVWM )
         ENDIF
         ISAMPA(  405) = NPROLD
         ISAMPP(  405) = ISAMPP(  405) + NPROLD
      ENDIF



C     QUADRATE DER ABLEITUNGEN VON U

      IF(ISELEP(1,430) .GE. 1 ) THEN
         NPROLD = ISAMPA(  430)
         IF(ISAMPP(  430) .GT. 0) THEN
            NPRRUN = ISAMPP(  430)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $           DREAD,DCONT,NPRNEU,FPRNEU,AUXUXM,SUXUXM)
         ENDIF
         ISAMPA(  430) = NPROLD
         ISAMPP(  430) = ISAMPP(  430) + NPROLD
      ENDIF

      IF(ISELEP(1,431) .GE. 1 ) THEN
         NPROLD = ISAMPA(  431)
         IF(ISAMPP(  431) .GT. 0) THEN
            NPRRUN = ISAMPP(  431)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $           DREAD,DCONT,NPRNEU,FPRNEU,AUYUYM,SUYUYM)
         ENDIF
         ISAMPA(  431) = NPROLD
         ISAMPP(  431) = ISAMPP(  431) + NPROLD
      ENDIF

      IF(ISELEP(1,432) .GE. 1 ) THEN
         NPROLD = ISAMPA(  432)
         IF(ISAMPP(  432) .GT. 0) THEN
            NPRRUN = ISAMPP(  432)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $           DREAD,DCONT,NPRNEU,FPRNEU,AUZUZM,SUZUZM)
         ENDIF
         ISAMPA(  432) = NPROLD
         ISAMPP(  432) = ISAMPP(  432) + NPROLD
      ENDIF

C     QUADRATE DER ABLEITUNGEN VON V

      IF(ISELEP(1,433) .GE. 1 ) THEN
         NPROLD = ISAMPA(  433)
         IF(ISAMPP(  433) .GT. 0) THEN
            NPRRUN = ISAMPP(  433)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $           DREAD,DCONT,NPRNEU,FPRNEU,AVXVXM,SVXVXM)
         ENDIF
         ISAMPA(  433) = NPROLD
         ISAMPP(  433) = ISAMPP(  433) + NPROLD
      ENDIF

      IF(ISELEP(1,434) .GE. 1 ) THEN
         NPROLD = ISAMPA(  434)
         IF(ISAMPP(  434) .GT. 0) THEN
            NPRRUN = ISAMPP(  434)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $           DREAD,DCONT,NPRNEU,FPRNEU,AVYVYM,SVYVYM)
         ENDIF
         ISAMPA(  434) = NPROLD
         ISAMPP(  434) = ISAMPP(  434) + NPROLD
      ENDIF

      IF(ISELEP(1,435) .GE. 1 ) THEN
         NPROLD = ISAMPA(  435)
         IF(ISAMPP(  435) .GT. 0) THEN
            NPRRUN = ISAMPP(  435)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $           DREAD,DCONT,NPRNEU,FPRNEU,AVZVZM,SVZVZM)
         ENDIF
         ISAMPA(  435) = NPROLD
         ISAMPP(  435) = ISAMPP(  435) + NPROLD
      ENDIF

C     QUADRATE DER ABLEITUNGEN VON W

      IF(ISELEP(1,436) .GE. 1 ) THEN
         NPROLD = ISAMPA(  436)
         IF(ISAMPP(  436) .GT. 0) THEN
            NPRRUN = ISAMPP(  436)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $           DREAD,DCONT,NPRNEU,FPRNEU,AWXWXM,SWXWXM)
         ENDIF
         ISAMPA(  436) = NPROLD
         ISAMPP(  436) = ISAMPP(  436) + NPROLD
      ENDIF

      IF(ISELEP(1,437) .GE. 1 ) THEN
         NPROLD = ISAMPA(  437)
         IF(ISAMPP(  437) .GT. 0) THEN
            NPRRUN = ISAMPP(  437)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $           DREAD,DCONT,NPRNEU,FPRNEU,AWYWYM,SWYWYM)
         ENDIF
         ISAMPA(  437) = NPROLD
         ISAMPP(  437) = ISAMPP(  437) + NPROLD
      ENDIF

      IF(ISELEP(1,438) .GE. 1 ) THEN
         NPROLD = ISAMPA(  438)
         IF(ISAMPP(  438) .GT. 0) THEN
            NPRRUN = ISAMPP(  438)
            CALL STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,NPRRUN,
     $           DREAD,DCONT,NPRNEU,FPRNEU,AWZWZM,SWZWZM)
         ENDIF
         ISAMPA(  438) = NPROLD
         ISAMPP(  438) = ISAMPP(  438) + NPROLD
      ENDIF






      RETURN
      END
