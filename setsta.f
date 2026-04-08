










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
      SUBROUTINE SETSTA (KPP,JPP,IPP,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,
     $                   XHOMOG,YHOMOG,ZHOMOG,
     $                   GMOL,RHO,UGRID,DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,
     $     RDX,RDY,RDZ,RDDX,RDDY,RDDZ,

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
C        S E T S T A      VORBELEGUNG ALLER FELDER, DIE FUER DIE STATIST
C                         AUSWERTUNG BENOETIGT WERDEN. ALLE FELDER IN
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
C PARAM: KPP,JPP,IPP    - ORT DES BEZUGSDRUCKES
C        IB1, IB2       - GRENZE DES KUBUSSES IN X-RI.
C        JB1, JB2       - GRENZE DES KUBUSSES IN Y-RI.
C        KB             - GRENZE DES KUBUSSES IN Z-RI. (TOP-FLAECHE)
C        GMOL           - MOLEKULARE DYNAMISCHE VISKOSITAET
C        RHO            - DICHTE DES FLUIDS (= CONST.)
C        PRTURB         - TURBULENT PRANDTL NUMBER
C        UGRID          - GESCHWINDIGKEIT DES KOORDINATENSYSTEMS
C                         IN X-RICHTUNG (GALILEI-TRANSFORMATION)
C        DDX,DDY,DDZ    - ABMESSUNGEN DER BASISZELLEN
C        DX,DY,DZ       - ABSTAND DER BASISZELLMITTELPUNKTE
C        X,Y,Z          - KOORDINATEN DER ZELLMITTELPUNKTE
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
C UPROG                 : AREAM, BF1ALL, BPHI0, DPHI0, ENERFG,
C                         ENERFS, GALIRT, PHIFLA, PHIFLU, PHIRMI,
C                         PHISKE, PLEVEL, TAUFG, TAUFM, TAUFS
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        28.08.86 (HW)  : ORIGINAL
C        12.09.86 (HW)  : INITIALISIERUNGSROUTINE F. RMS-WERTE
C        30.06.88 (HW)  : ERWEITERUNG: AUSWERTUNG DER DREI KOMPONENTEN
C                         DER VORTICITY, DER ENSTROPHIE UND DER
C                         HELIZITAET
C        21.08.88 (HW)  : VORTICITY IST JETZT IM DOPPELT VERSETZTEN
C                         MASCHENGITTER DEFINIERT
C        04.01.89 (HW)  : DIE HELIZITAET IST JETZT DREIFACH IM MASCHEN-
C                         GITTER VERSCHOBEN (IN X-, Y- UND Z-RICHTUNG)
C        19.12.89 (HW)  : BF1ALL SETZT FUER 'P' AN FESTEN WAENDEN
C                         DP/DN = 0.0
C        30.10.90 (HW)  : UEBERGABE AN PLEVEL GEAENDERT
C        12.02.03 (TB)  : SCALAR FIELD STATISTICS IMPLEMENTED
C
C*STARLET***************************************************************
C
C
      REAL     DDX(II),   DDY(JJ),   DDZ(KK),
     $          DX(II),    DY(JJ),    DZ(KK),
     $           X(II),     Y(JJ),     Z(KK)
C
      LOGICAL  XHOMOG,YHOMOG,ZHOMOG
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
C                                 VORBELEGUNG DES HILF-FELDES MIT 0.0
C
      CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,HILF   )
C
C                                 STATIST. GROESSEN DER U-KOMPONENTE
C                                 ----------------------------------
C
C                                 VORBELEGUNG DES <U> - FELDES MIT DEM
C                                 "BESTMOEGLICHEN" ENSEMBLE-MITTELWERT
C
      IF(ISELEP(1,  4) .GE. 1) THEN
         IF(ISAMPA(    4) .EQ. 0  ) THEN
C
C                                 ZUNAECHST WIRD DAS GESAMTE AU-FELD
C                                 MIT 0.0 VORBELEGT
C
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AU     )
C
C                                 DER RAND DES BERECHNUNGSGEBIETES
C                                 UND EINE SCHICHT IM INNEREN DES
C                                 KUBUSSES WERDEN AUF 0.0 GESETZT
C
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
C
C                                 GALILEI-RUECKTRANSFORMATION (NUR FUER
C                                 DAS U-FELD NOTWENDIG)
C
            CALL GALIRT  (KK,JJ,II,KMX,JMX,IMX,
     $                    UGRID,U,HILF,B,'   U',NFRO,NRGT,NBOT,
     $                                          NBAC,NLFT,NTOP)
C
C                                 SETZEN DER RANDBEDINGUNGEN
C
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'U     ')
C
C                                 BILDUNG EINES LINIEN/FLAECHENMITTEL-
C                                 WERTES UND SUMMATION IN DAS (MIT 0.0)
C                                 VORBELEGTE AU-FELD
C
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AU     ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
C
C                                 VORBELEGUNG DES SUMMATIONSFELDES
C                                 FUER DIE MOMENTANWERTE VON  U.
C
      IF(ISELEP(1,  5) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SU     )
      ENDIF
C
C                                 VORBELEGUNG DES U-FLUKTUATIONS FELDES
C
      IF(ISELEP(1,  6) .GE. 1) THEN
C
C                                 ZUNAECHST WIRD DAS GESAMTE UFG-FELD
C                                 MIT 0.0 VORBELEGT
C
            CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,UFG    )
C
C                                 BERECHNUNG DER U-FLUKTUATIONEN
C
            CALL PHIFLU  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,UGRID,
     $                    U      ,UFG    ,AU     ,'   U',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C
C                                 SETZEN DER RANDBEDINGUNGEN
C
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,UFG    ,IB1,IB2,JB1,JB2,KB,'UF G  ')
      ENDIF
C
C                                 VORBELEGUNG DES <U-RMS> - FELDES
C                                 MIT DEM "BESTMOEGLICHEN" ENSEMBLE-
C                                 MITTELWERT
C
      IF(ISELEP(1,  7) .GE. 1) THEN
         IF(ISAMPA(    7) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AURG   )
C
C                                 DER RAND DES BERECHNUNGSGEBIETES
C                                 UND EINE SCHICHT IM INNEREN DES
C                                 KUBUSSES WERDEN AUF 0.0 GESETZT
C
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
C
C                                 BERECHNUNG DER RMS-WERTE
C
            CALL PHIRMI  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,
     $                    UFG    ,HILF   ,'   U')
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'UR G  ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AURG   ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
C
C                                 VORBELEGUNG DES SUMMATIONSFELDES
C                                 FUER DIE MOMENTANEN RMS-WERTE
C
      IF(ISELEP(1,  8) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SURG   )
      ENDIF
C
C                                 VORBELEGUNG DES <U-SKE> - FELDES
C                                 MIT DEM "BESTMOEGLICHEN" ENSEMBLE-
C                                 MITTELWERT
C
      IF(ISELEP(1, 10) .GE. 1) THEN
         IF(ISAMPA(   10) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUSKG  )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
C
C                                 BERECHNUNG DER SKEWNESS-WERTE
C
            CALL PHISKE  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,
     $                    UFG    ,AURG   ,HILF   ,'   U',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'USKG  ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AUSKG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
C
C                                 VORBELEGUNG DES SUMMATIONSFELDES
C                                 FUER DIE MOMENTANEN SKEWNESS-WERTE
C
      IF(ISELEP(1, 11) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUSKG  )
      ENDIF
C
C                                 VORBELEGUNG DES <U-FLA> - FELDES
C                                 MIT DEM "BESTMOEGLICHEN" ENSEMBLE-
C                                 MITTELWERT
C
      IF(ISELEP(1, 12) .GE. 1) THEN
         IF(ISAMPA(   12) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUFLG  )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
C
C                                 BERECHNUNG DER FLATNESS-WERTE
C
            CALL PHIFLA  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,
     $                    UFG    ,AURG   ,HILF   ,'   U',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'UFLG  ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AUFLG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
C
C                                 VORBELEGUNG DES SUMMATIONSFELDES
C                                 FUER DIE MOMENTANEN FLATNESS-WERTE
C
      IF(ISELEP(1, 13) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFLG  )
      ENDIF
C
C                                 STATIST. GROESSEN DER V-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 29) .GE. 1) THEN
         IF(ISAMPA(   29) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AV     )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL GALIRT  (KK,JJ,II,KMX,JMX,IMX,
     $                    UGRID,V,HILF,B,'   V',NFRO,NRGT,NBOT,
     $                                          NBAC,NLFT,NTOP)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'V     ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AV     ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1, 30) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SV     )
      ENDIF
      IF(ISELEP(1, 31) .GE. 1) THEN
            CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,VFG    )
            CALL PHIFLU  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,UGRID,
     $                    V      ,VFG    ,AV     ,'   V',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,VFG    ,IB1,IB2,JB1,JB2,KB,'VF G  ')
      ENDIF
      IF(ISELEP(1, 32) .GE. 1) THEN
         IF(ISAMPA(   32) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AVRG   )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL PHIRMI  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,
     $                    VFG    ,HILF   ,'   V')
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'VR G  ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AVRG   ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1, 33) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVRG   )
      ENDIF
      IF(ISELEP(1, 35) .GE. 1) THEN
         IF(ISAMPA(   35) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AVSKG  )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL PHISKE  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,
     $                    VFG    ,AVRG   ,HILF   ,'   V',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'VSKG  ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AVSKG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1, 36) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVSKG  )
      ENDIF
      IF(ISELEP(1, 37) .GE. 1) THEN
         IF(ISAMPA(   37) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AVFLG  )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL PHIFLA  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,
     $                    VFG    ,AVRG   ,HILF   ,'   V',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'VFLG  ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AVFLG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1, 38) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVFLG  )
      ENDIF
C
C                                 STATIST. GROESSEN DER W-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 54) .GE. 1) THEN
         IF(ISAMPA(   54) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AW     )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL GALIRT  (KK,JJ,II,KMX,JMX,IMX,
     $                    UGRID,W,HILF,B,'   W',NFRO,NRGT,NBOT,
     $                                           NBAC,NLFT,NTOP)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'W     ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AW     ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1, 55) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SW     )
      ENDIF
      IF(ISELEP(1, 56) .GE. 1) THEN
            CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,WFG    )
            CALL PHIFLU  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,UGRID,
     $                    W      ,WFG    ,AW     ,'   W',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,WFG    ,IB1,IB2,JB1,JB2,KB,'WF G  ')
      ENDIF
      IF(ISELEP(1, 57) .GE. 1) THEN
         IF(ISAMPA(   57) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AWRG   )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL PHIRMI  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,
     $                    WFG    ,HILF   ,'   W')
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'WR G  ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AWRG   ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1, 58) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SWRG   )
      ENDIF
      IF(ISELEP(1, 60) .GE. 1) THEN
         IF(ISAMPA(   60) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AWSKG  )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL PHISKE  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,
     $                    WFG    ,AWRG   ,HILF   ,'   W',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'WSKG  ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AWSKG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1, 61) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SWSKG  )
      ENDIF
      IF(ISELEP(1, 62) .GE. 1) THEN
         IF(ISAMPA(   62) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AWFLG  )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL PHIFLA  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,
     $                    WFG    ,AWRG   ,HILF   ,'   W',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'WFLG  ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AWFLG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1, 63) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SWFLG  )
      ENDIF
C
C                                 STATIST. GROESSEN DER P-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 79) .GE. 1) THEN
         IF(ISAMPA(   79) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AP     )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
C
C                                 REDUZIERUNG DES DRUCKNIVEAUS
C                                 WIRD IM HAUPTPROGRAMM GEMACHT
C
C           CALL PLEVEL  (KK,JJ,II,KMX,JMX,IMX,KPP,JPP,IPP,P,
C    $                    DDX, DDY, DDZ )
            CALL GALIRT  (KK,JJ,II,KMX,JMX,IMX,
     $                    UGRID,P,HILF,B,'   P',NFRO,NRGT,NBOT
     $                                          ,NBAC,NLFT,NTOP)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'P     ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AP     ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1, 80) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SP     )
      ENDIF
      IF(ISELEP(1, 81) .GE. 1) THEN
            CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,PFG    )
            CALL PHIFLU  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,UGRID,
     $                    P      ,PFG    ,AP     ,'   P',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,PFG    ,IB1,IB2,JB1,JB2,KB,'PF G  ')
      ENDIF
      IF(ISELEP(1, 82) .GE. 1) THEN
         IF(ISAMPA(   82) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,APRG   )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL PHIRMI  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,
     $                    PFG    ,HILF   ,'   P')
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'PR G  ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,APRG   ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1, 83) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SPRG   )
      ENDIF
      IF(ISELEP(1, 85) .GE. 1) THEN
         IF(ISAMPA(   85) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,APSKG  )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL PHISKE  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,
     $                    PFG    ,APRG   ,HILF   ,'   P',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'PSKG  ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,APSKG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1, 86) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SPSKG  )
      ENDIF
      IF(ISELEP(1, 87) .GE. 1) THEN
         IF(ISAMPA(   87) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,APFLG  )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL PHIFLA  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,
     $                    PFG    ,APRG   ,HILF   ,'   P',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'PFLG  ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,APFLG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1, 88) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SPFLG  )
      ENDIF
C
C                                 KINETISCHE ENERGIE DER SCHWANKUNGS-
C                                 GESCHWINDIGKEITEN
C                                 -----------------------------------
C
      IF(ISELEP(1,106) .GE. 1) THEN
         IF(ISAMPA(  106) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AEFG   )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
C
C                                 BERECHNUNG DER KINETISCHEN
C                                 ENERGIE DER SCHWANKUNGSGESCHW.
C                                 (ANTEIL DER GROBSTRUKTUR)
C
            CALL ENERFG  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,
     $                    DX,DY,DZ,X,Y,Z,UFG,VFG,WFG,HILF,ESUM,
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,0)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'EF G  ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AEFG   ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,107) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SEFG   )
      ENDIF
      IF(ISELEP(1,109) .GE. 1) THEN
         IF(ISAMPA(  109) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AEFS   )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
C
C                                 BERECHNUNG DER KINETISCHEN
C                                 ENERGIE DER SCHWANKUNGSGESCHW.
C                                 (ANTEIL DER FEINSTRUKTUR)
C
            CALL ENERFS  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,
     $                    DX,DY,DZ,X,Y,Z,B,G,GMOL,RHO,HILF,ESUM,
     $                 IC1,IC2,JC1,JC2,KC1,KC2,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,0)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'EF S  ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AEFS   ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,110) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SEFS   )
      ENDIF
C
C                                 DISSIPATION
C                                 ___________
      IF(ISELEP(1,120) .GE. 1) THEN
         IF(ISAMPA(  120) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,ADFG   )
C
            CALL DISSIPG (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,
     $                    DX,DY,DZ,X,Y,Z,UFG,VFG,WFG,B,HILF,DSUM,
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,0)
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,ADFG   ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF

C                               FUER   TAYLOR''SCHES MIKROMASS
C                               QUADRATE DER ELEMENTE DES 
C                               DEFORMATIONSGESCHWINDIGKEITSTENSORS
C                               -----------------------------------
C
C                                  U_KOMPONENTE IN X_RICHTUNG
      IF(ISELEP(1,131) .GE. 1) THEN
         IF(ISAMPA(  131) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,ADUDX2 )
            CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   )
C
            CALL PHIDIDJ (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,HILF,
     $                    UFG,X,DX,DDX,II,'XX',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
             CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                     KMXA,JMXA,IMXA,
     $                     HILF   ,HILF   ,'  QP')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,ADUDX2 ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
C
C                                  U_KOMPONENTE IN Y_RICHTUNG
      IF(ISELEP(1,132) .GE. 1) THEN
         IF(ISAMPA(  132) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,ADUDY2 )
            CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   )
C
            CALL PHIDIDJ (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,HILF,
     $                    UFG,Y,DY,DDY,JJ,'XY',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
             CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                     KMXA,JMXA,IMXA,
     $                     HILF   ,HILF   ,'  VU')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,ADUDY2 ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF

C                                  U_KOMPONENTE IN Z_RICHTUNG
      IF(ISELEP(1,133) .GE. 1) THEN
         IF(ISAMPA(  133) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,ADUDZ2 )
            CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   )
C
            CALL PHIDIDJ (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,HILF,
     $                    UFG,Z,DZ,DDZ,KK,'XZ',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
             CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                     KMXA,JMXA,IMXA,
     $                     HILF   ,HILF   ,'  UW')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,ADUDZ2 ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
C
C                                  V_KOMPONENTE IN X_RICHTUNG
      IF(ISELEP(1,134) .GE. 1) THEN
         IF(ISAMPA(  134) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,ADVDX2 )
            CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   )
C
            CALL PHIDIDJ (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,HILF,
     $                    VFG,X,DX,DDX,II,'YX',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
             CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                     KMXA,JMXA,IMXA,
     $                     HILF   ,HILF   ,'  VU')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,ADVDX2 ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
C
C                                  V_KOMPONENTE IN Y_RICHTUNG
      IF(ISELEP(1,135) .GE. 1) THEN
         IF(ISAMPA(  135) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,ADVDY2 )
            CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   )
C
            CALL PHIDIDJ (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,HILF,
     $                    VFG,Y,DY,DDY,JJ,'YY',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
             CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                     KMXA,JMXA,IMXA,
     $                     HILF   ,HILF   ,'  QP')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,ADVDY2 ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
C
C                                  V_KOMPONENTE IN Z_RICHTUNG
      IF(ISELEP(1,136) .GE. 1) THEN
         IF(ISAMPA(  136) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,ADVDZ2 )
            CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   )
C
            CALL PHIDIDJ (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,HILF,
     $                    VFG,Z,DZ,DDZ,KK,'YZ',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
             CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                     KMXA,JMXA,IMXA,
     $                     HILF   ,HILF   ,'  WV')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,ADVDZ2 ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
C
C                                  W_KOMPONENTE IN X_RICHTUNG
      IF(ISELEP(1,137) .GE. 1) THEN
         IF(ISAMPA(  137) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,ADWDX2 )
            CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   )
C
            CALL PHIDIDJ (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,HILF,
     $                    WFG,X,DX,DDX,II,'ZX',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
             CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                     KMXA,JMXA,IMXA,
     $                     HILF   ,HILF   ,'  UW')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,ADWDX2 ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
C
C                                  W_KOMPONENTE IN Y_RICHTUNG
      IF(ISELEP(1,138) .GE. 1) THEN
         IF(ISAMPA(  138) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,ADWDY2 )
            CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   )
C
            CALL PHIDIDJ (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,HILF,
     $                    WFG,Y,DY,DDY,JJ,'ZY',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
             CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                     KMXA,JMXA,IMXA,
     $                     HILF   ,HILF   ,'  WV')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,ADWDY2 ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
C
C                                  W_KOMPONENTE IN Z_RICHTUNG
      IF(ISELEP(1,139) .GE. 1) THEN
         IF(ISAMPA(  139) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,ADWDZ2 )
            CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   )
C
            CALL PHIDIDJ (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,HILF,
     $                    WFG,Z,DZ,DDZ,KK,'ZZ',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
             CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                     KMXA,JMXA,IMXA,
     $                     HILF   ,HILF   ,'  QP')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,ADWDZ2 ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF

C
C                                 KOMPONENTEN DES SPANNUNGSTENSORS
C                                 --------------------------------
C
C                                 HIER: U - W
C                                 -----------
C
      IF(ISELEP(1,146) .GE. 1) THEN
         IF(ISAMPA(  146) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUFWFG )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
C
C                                 BERECHNUNG DES GROBSTRUKTUR-
C                                 ANTEILS AN DEN SCHUBSPANNUNGEN
C
            CALL TAUFG   (KK,JJ,II,KMX,JMX,IMX,GMOL,RHO,
     $                    HILF   ,UFG    ,WFG    ,'UW',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'UW    ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AUFWFG ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,147) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFWFG )
      ENDIF
      IF(ISELEP(1,149) .GE. 1) THEN
         IF(ISAMPA(  149) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUFWFS )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
C
C                                 BERECHNUNG DES FEINSTRUKTUR-
C                                 ANTEILS AN DEN SCHUBSPANNUNGEN
C
            CALL TAUFS   (KK,JJ,II,KMX,JMX,IMX,IB1,IB2,JB1,JB2,KB,
     $                    G,GMOL,RHO,UGRID,
     $                    HILF   ,U,X,DX,DDX,II,W,Z,DZ,DDZ,KK,'UW'
     $                    )     
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'UW    ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AUFWFS ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,150) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFWFS )
      ENDIF
      IF(ISELEP(1,151) .GE. 1) THEN
         IF(ISAMPA(  151) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUFWFM )
            CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,HILF   )
C
C                                 BERECHNUNG DES MOLEKULAREN
C                                 ANTEILS AN DEN SCHUBSPANNUNGEN
C                                 (EINSCHL. DER WANDSCHUBSPANNUNGEN)
C
            CALL TAUFM   (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,
     $                    HILF   ,U,X,DX,DDX,II,W,Z,DZ,DDZ,KK,'UW',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB
     $                    )      
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'UW    ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AUFWFM ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,152) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFWFM )
      ENDIF
C
C                                 HIER: V - W
C                                 -----------
C
      IF(ISELEP(1,156) .GE. 1) THEN
         IF(ISAMPA(  156) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AVFWFG )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL TAUFG   (KK,JJ,II,KMX,JMX,IMX,GMOL,RHO,
     $                    HILF   ,VFG    ,WFG    ,'VW',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'VW    ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AVFWFG ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,157) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVFWFG )
      ENDIF
      IF(ISELEP(1,159) .GE. 1) THEN
         IF(ISAMPA(  159) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AVFWFS )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL TAUFS   (KK,JJ,II,KMX,JMX,IMX,IB1,IB2,JB1,JB2,KB,
     $                    G,GMOL,RHO,UGRID,
     $                    HILF   ,V,Y,DY,DDY,JJ,W,Z,DZ,DDZ,KK,'VW'
     $                    )     
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'VW    ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AVFWFS ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,160) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVFWFS )
      ENDIF
      IF(ISELEP(1,161) .GE. 1) THEN
         IF(ISAMPA(  161) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AVFWFM )
            CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,HILF   )
            CALL TAUFM   (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,
     $                    HILF   ,V,Y,DY,DDY,JJ,W,Z,DZ,DDZ,KK,'VW',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB
     $                    )      
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'VW    ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AVFWFM ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,162) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVFWFM )
      ENDIF
C
C                                 HIER: U - V
C                                 -----------
C
      IF(ISELEP(1,166) .GE. 1) THEN
         IF(ISAMPA(  166) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUFVFG )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL TAUFG   (KK,JJ,II,KMX,JMX,IMX,GMOL,RHO,
     $                    HILF   ,UFG    ,VFG    ,'UV',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'UV    ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AUFVFG ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,167) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFVFG )
      ENDIF
      IF(ISELEP(1,169) .GE. 1) THEN
         IF(ISAMPA(  169) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUFVFS )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL TAUFS   (KK,JJ,II,KMX,JMX,IMX,IB1,IB2,JB1,JB2,KB,
     $                    G,GMOL,RHO,UGRID,
     $                    HILF   ,U,X,DX,DDX,II,V,Y,DY,DDY,JJ,'UV'
     $                    )     
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'UV    ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AUFVFS ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,170) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFVFS )
      ENDIF
      IF(ISELEP(1,171) .GE. 1) THEN
         IF(ISAMPA(  171) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUFVFM )
            CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,HILF   )
            CALL TAUFM   (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,
     $                    HILF   ,U,X,DX,DDX,II,V,Y,DY,DDY,JJ,'UV',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB
     $                    )      
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AUFVFM ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,172) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFVFM )
      ENDIF
C
C
C                                 STATIST. GROESSEN DER X-KOMPONENTE
C                                 DER VORTICITY
C                                 ----------------------------------
C
      IF(ISELEP(1,176) .GE. 1) THEN
         CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,OX     )
      ENDIF
C
C                                 VORBELEGUNG DES <OX> - FELDES MIT DEM
C                                 "BESTMOEGLICHEN" ENSEMBLE-MITTELWERT
C
      IF(ISELEP(1,177) .GE. 1) THEN
         IF(ISAMPA(  177) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AOX    )
C
C                                 BERECHNUNG DES MOMENTANWERTES DER
C                                 VORTICITY
C
            CALL OMEGA   (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    G, GMOL, RHO,    UGRID,
     $                    OX,W ,Z ,DZ ,DDZ ,KK,V ,Y ,DY ,DDY ,JJ,'OX',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,  OX   ,AOX    ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,178) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SOX    )
      ENDIF
      IF(ISELEP(1,179) .GE. 1) THEN
            CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,OXFG   )
            CALL PHIFLU  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,UGRID,
     $                    OX     ,OXFG   ,AOX    ,'  WV',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,OXFG   ,IB1,IB2,JB1,JB2,KB,'WV    ')
      ENDIF
      IF(ISELEP(1,180) .GE. 1) THEN
         IF(ISAMPA(  180) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AOXRG  )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL PHIRMI  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,
     $                    OXFG   ,HILF   ,'  WV')
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'WV    ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AOXRG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,181) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SOXRG  )
      ENDIF
C
C                                 STATIST. GROESSEN DER Y-KOMPONENTE
C                                 DER VORTICITY
C                                 ----------------------------------
C
      IF(ISELEP(1,182) .GE. 1) THEN
         CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,OY     )
      ENDIF
      IF(ISELEP(1,183) .GE. 1) THEN
         IF(ISAMPA(  183) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AOY    )
            CALL OMEGA   (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    G, GMOL, RHO,    UGRID,
     $                    OY,U ,X ,DX ,DDX ,II,W ,Z ,DZ ,DDZ ,KK,'OY',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,  OY   ,AOY    ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,184) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SOY    )
      ENDIF
      IF(ISELEP(1,185) .GE. 1) THEN
            CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,OYFG   )
            CALL PHIFLU  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,UGRID,
     $                    OY     ,OYFG   ,AOY    ,'  UW',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,OYFG   ,IB1,IB2,JB1,JB2,KB,'UW    ')
      ENDIF
      IF(ISELEP(1,186) .GE. 1) THEN
         IF(ISAMPA(  186) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AOYRG  )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL PHIRMI  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,
     $                    OYFG   ,HILF   ,'  UW')
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'UW    ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AOYRG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,187) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SOYRG  )
      ENDIF
C
C                                 STATIST. GROESSEN DER Z-KOMPONENTE
C                                 DER VORTICITY
C                                 ----------------------------------
C
      IF(ISELEP(1,188) .GE. 1) THEN
         CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,OZ     )
      ENDIF
      IF(ISELEP(1,189) .GE. 1) THEN
         IF(ISAMPA(  189) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AOZ    )
            CALL OMEGA   (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    G, GMOL, RHO,    UGRID,
     $                    OZ,V ,Y ,DY ,DDY ,JJ,U ,X ,DX ,DDX ,II,'OZ',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,  OZ   ,AOZ    ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,190) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SOZ    )
      ENDIF
      IF(ISELEP(1,191) .GE. 1) THEN
            CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,OZFG   )
            CALL PHIFLU  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,UGRID,
     $                    OZ     ,OZFG   ,AOZ    ,'  VU',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,OZFG   ,IB1,IB2,JB1,JB2,KB,'VU    ')
      ENDIF
      IF(ISELEP(1,192) .GE. 1) THEN
         IF(ISAMPA(  192) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AOZRG  )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL PHIRMI  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,
     $                    OZFG   ,HILF   ,'  VU')
C           CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                    X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'VU    ')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AOZRG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,193) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SOZRG  )
      ENDIF
C
C                                 STATIST. GROESSEN DER ENSTROPHIE
C                                 ----------------------------------
C
      IF(ISELEP(1,195) .GE. 1) THEN
         IF(ISAMPA(  195) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AO2    )
C
C                                 BERECHNUNG DES MOMENTANWERTES DER
C                                 ENSTROPHIE
C
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL ENSTRO  (KK,JJ,II,KMX,JMX,IMX,OX,OY,OZ,HILF  )
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA, HILF  ,AO2    ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,196) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SO2    )
      ENDIF
      IF(ISELEP(1,197) .GE. 1) THEN
            CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,O2FG   )
            CALL PHIFLU  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,UGRID,
     $                    HILF   ,O2FG   ,AO2    ,' UVW',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
      ENDIF
      IF(ISELEP(1,198) .GE. 1) THEN
         IF(ISAMPA(  198) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AO2RG  )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL PHIRMI  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,
     $                    O2FG   ,HILF   ,' UVW')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AO2RG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,199) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SO2RG  )
      ENDIF
C
C                                 STATIST. GROESSEN DER HELIZITAET
C                                 ----------------------------------
C
      IF(ISELEP(1,200) .GE. 1) THEN
         IF(ISAMPA(  200) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AHE    )
C
C                                 BERECHNUNG DES MOMENTANWERTES DER
C                                 HELIZITAET
C
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL HELICI  (KK,JJ,II,KMX,JMX,IMX,U,V,W,OX,OY,OZ,
     $                    HILF  )
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA, HILF  ,AHE    ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,201) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SHE    )
      ENDIF
      IF(ISELEP(1,202) .GE. 1) THEN
            CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,HEFG   )
            CALL PHIFLU  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,UGRID,
     $                    HILF   ,HEFG   ,AHE    ,' UVW',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
      ENDIF
      IF(ISELEP(1,203) .GE. 1) THEN
         IF(ISAMPA(  203) .EQ. 0  ) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AHERG  )
            CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                    IB1,IB2,JB1,JB2,KB)
            CALL PHIRMI  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,
     $                    HEFG   ,HILF   ,' UVW')
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,AHERG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF
      IF(ISELEP(1,204) .GE. 1) THEN
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SHERG  )
      ENDIF
C

C
C
C                                "NEUE" Statistik
C                                 ---------------
C


      IF (ISELEP(1,400) .GE. 1) THEN
         IF(ISAMPA(  400) .EQ. 0  ) THEN

         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUUM   )
         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUUM   )
         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $                 HILF3D1,'XX',U,U,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)

         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF3D1,AUUM,
     $                    XHOMOG,YHOMOG,ZHOMOG)

         ENDIF
      ENDIF


      IF (ISELEP(1,401) .GE. 1) THEN
         IF(ISAMPA(  401) .EQ. 0  ) THEN

         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AVVM   )
         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVVM   )
         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $                 HILF3D1,'XX',V,V,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF3D1,AVVM,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF

      IF (ISELEP(1,402) .GE. 1) THEN
         IF(ISAMPA(  402) .EQ. 0  ) THEN
c
         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AWWM   )
         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SWWM   )
         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $                 HILF3D1,'XX',W,W,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF3D1,AWWM,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF

      IF (ISELEP(1,406) .GE. 1) THEN
         IF(ISAMPA(  406) .EQ. 0  ) THEN

         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,APPM   )
         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SPPM   )
         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $                 HILF3D1,'XX',P,P,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)

         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF3D1,APPM,
     $                    XHOMOG,YHOMOG,ZHOMOG)

         ENDIF
      ENDIF



      IF (ISELEP(1,403) .GE. 1) THEN
         IF(ISAMPA(  403) .EQ. 0  ) THEN

         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUVM   )
         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUVM   )
         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $                 HILF3D1,'UV',U,V,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF3D1,AUVM,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF

      IF (ISELEP(1,404) .GE. 1) THEN
         IF(ISAMPA(  404) .EQ. 0  ) THEN

         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUWM   )
         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUWM   )
         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $                 HILF3D1,'UW',U,W,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF3D1,AUWM,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF

      IF (ISELEP(1,405) .GE. 1) THEN
         IF(ISAMPA(  405) .EQ. 0  ) THEN

         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AVWM   )
         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVWM   )
         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $                 HILF3D1,'VW',V,W,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF3D1,AVWM,
     $                    XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF



C     QUADRATE DER ABLEITUNGEN VON U

      IF (ISELEP(1,430) .GE. 1) THEN
         IF(ISAMPA(  430) .EQ. 0  ) THEN

            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUXUXM )
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUXUXM )

            CALL DFDX  (KK,JJ,II,RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $           'DDX',U,HILF3D1)

            CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $           HILF3D2,'XX',HILF3D1,HILF3D1,
     $           NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $           KMXA,JMXA,IMXA,HILF3D2,AUXUXM,
     $           XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF

      IF (ISELEP(1,431) .GE. 1) THEN
         IF(ISAMPA(  431) .EQ. 0  ) THEN

            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUYUYM )
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUYUYM )

            CALL DFDX  (KK,JJ,II,RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $           'DYS',U,HILF3D1)
            CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $           HILF3D2,'XX',HILF3D1,HILF3D1,
     $           NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $           KMXA,JMXA,IMXA,HILF3D2,AUYUYM,
     $           XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF

      IF (ISELEP(1,432) .GE. 1) THEN
         IF(ISAMPA(  432) .EQ. 0  ) THEN

            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUZUZM )
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUZUZM )

            CALL DFDX  (KK,JJ,II,RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $           'DZS',U,HILF3D1)
            CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $           HILF3D2,'XX',HILF3D1,HILF3D1,
     $           NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $           KMXA,JMXA,IMXA,HILF3D2,AUZUZM,
     $           XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF

C     QUADRATE DER ABLEITUNGEN VON V

      IF (ISELEP(1,433) .GE. 1) THEN
         IF(ISAMPA(  433) .EQ. 0  ) THEN

            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AVXVXM )
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVXVXM )
            CALL DFDX  (KK,JJ,II,RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $           'DXS',V,HILF3D1)
            CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $           HILF3D2,'XX',HILF3D1,HILF3D1,
     $           NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $           KMXA,JMXA,IMXA,HILF3D2,AVXVXM,
     $           XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF

      IF (ISELEP(1,434) .GE. 1) THEN
         IF(ISAMPA(  434) .EQ. 0  ) THEN

            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AVYVYM )
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVYVYM )

            CALL DFDX  (KK,JJ,II,RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $           'DDY',V,HILF3D1)

            CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $           HILF3D2,'XX',HILF3D1,HILF3D1,
     $           NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)

            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $           KMXA,JMXA,IMXA,HILF3D2,AVYVYM,
     $           XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF

      IF (ISELEP(1,435) .GE. 1) THEN
         IF(ISAMPA(  435) .EQ. 0  ) THEN

            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AVZVZM )
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVZVZM )

            CALL DFDX  (KK,JJ,II,RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $           'DZS',V,HILF3D1)

            CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $           HILF3D2,'XX',HILF3D1,HILF3D1,
     $           NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $           KMXA,JMXA,IMXA,HILF3D2,AVZVZM,
     $           XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF

C     QUADRATE DER ABLEITUNGEN VON W

      IF (ISELEP(1,436) .GE. 1) THEN
         IF(ISAMPA(  436) .EQ. 0  ) THEN

            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AWXWXM )
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SWXWXM )

            CALL DFDX  (KK,JJ,II,RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $           'DXS',W,HILF3D1)

            CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $           HILF3D2,'XX',HILF3D1,HILF3D1,
     $           NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)

            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $           KMXA,JMXA,IMXA,HILF3D2,AWXWXM,
     $           XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF

      IF (ISELEP(1,437) .GE. 1) THEN
         IF(ISAMPA(  437) .EQ. 0  ) THEN

            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AWYWYM )
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SWYWYM )

            CALL DFDX  (KK,JJ,II,RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $           'DYS',W,HILF3D1)

            CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $           HILF3D2,'XX',HILF3D1,HILF3D1,
     $           NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $           KMXA,JMXA,IMXA,HILF3D2,AWYWYM,
     $           XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF

      IF (ISELEP(1,438) .GE. 1) THEN
         IF(ISAMPA(  438) .EQ. 0  ) THEN

            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AWZWZM )
            CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SWZWZM )

            CALL DFDX  (KK,JJ,II,RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $           'DDZ',W,HILF3D1)

            CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $           HILF3D2,'XX',HILF3D1,HILF3D1,
     $           NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $           KMXA,JMXA,IMXA,HILF3D2,AWZWZM,
     $           XHOMOG,YHOMOG,ZHOMOG)
         ENDIF
      ENDIF






      RETURN
      END
