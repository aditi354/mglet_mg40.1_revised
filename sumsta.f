










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
      SUBROUTINE SUMSTA (KPP,JPP,IPP,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,
     $                   XHOMOG,YHOMOG,ZHOMOG,
     $                   GMOL,RHO,UGRID,
     $                   DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,
     $                 RDX,RDY,RDZ,RDDX,RDDY,RDDZ,


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
C        S U M S T A      IN SUMSTA WERDEN DIE MOMENTANWERTE VERSCHIE-
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
C PARAM: KPP,JPP,IPP    - ORT DES BEZUGSDRUCKES
C        IB1, IB2       - GRENZE DES KUBUSSES IN X-RI.
C        JB1, JB2       - GRENZE DES KUBUSSES IN Y-RI.
C        KB             - GRENZE DES KUBUSSES IN Z-RI. (TOP-FLAECHE)
C        GMOL           - MOLEKULARE DYNAMISCHE VISKOSITAET
C        RHO            - DICHTE DES FLUIDS (= CONST.)
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
C UPROG                 : AREAM, BF1ALL, (BF1CUB), BPHI0, DPHI0, ENERFG,
C                         ENERFS, GALIRT, PHIFLA, PHIFLU, PHIFL2,
C                         PHISKE, PLEVEL, TAUFG, TAUFM, TAUFS
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        28.08.86 (HW)  : ORIGINAL
C        12.09.86 (HW)  : GEAENDERTE BERECHNUNG DER RMS-WERTE
C        30.09.86 (HW)  : UEBERGABE DER BEZUGSGROESSEN IM RIDENT-FELD.
C                         (FUER DIE AUSGABE VON DIMENSIONSLOSEN
C                         MOMENTANWERTEN NOETIG)
C        01.07.88 (HW)  : ERWEITERUNG: AUSWERTUNG DER DREI KOMPONENTEN
C                         DER VORTICITY, DER ENSTROPHIE UND DER
C                         HELIZITAET
C        21.08.88 (HW)  : VORTICITY IST JETZT IM DOPPELT VERSETZTEN
C                         MASCHENGITTER DEFINIERT
C        15.12.88 (HW)  : BESTIMMUNG DER REFERENZGROESSEN ALS COMMON-
C                         DECK EINGEFUEHRT (CREFVA).
C        04.01.89 (HW)  : DIE HELIZITAET IST JETZT DREIFACH IM MASCHEN-
C                         GITTER VERSCHOBEN (IN X-, Y- UND Z-RICHTUNG)
C        19.12.89 (HW)  : BF1ALL SETZT FUER 'P' AN FESTEN WAENDEN
C                         DP/DN = 0.0
C        30.10.90 (HW)  : UEBERGABE AN PLEVEL GEAENDERT
C         1.11.93 (MM)  : RANDBEDINGUNGEN WERDEN UEBER KOPF WEITERGEGEBEN
C                         BF1ALL (SETZEN DER RANDBD. FUER STATISTISCHE
C                         GROESSEN ENTFAELLT)
C        12.02.03 (TB)  : SCALAR FIELD STATISTICS IMPLEMENTED
C
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
C
      REAL     DDX(II),   DDY(JJ),   DDZ(KK),
     $          DX(II),    DY(JJ),    DZ(KK),
     $           X(II),     Y(JJ),     Z(KK)
      REAL     RDX(II),   RDY(JJ),   RDZ(KK),
     $        RDDX(II),  RDDY(JJ),  RDDZ(KK)
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
Cinclude "crefva.h"
C
C                                 VORBELEGUNG DES HILF-FELDES MIT 0.0
C
      CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,HILF   )
C
C                                 STATIST. GROESSEN DER U-KOMPONENTE
C                                 ----------------------------------
C
C                                 SUMMATION DER U-MOMENTANWERTE IN DAS
C                                 SU-FELD
C
      IF(ISELEP(1,  4) .GE. 1) THEN
C
C                                 DER RAND DES BERECHNUNGSGEBIETES
C                                 UND EINE SCHICHT IM INNEREN DES
C                                 KUBUSSES WERDEN AUF 0.0 GESETZT
C
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
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
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'U     ')
C
C                                 BILDUNG EINES LINIEN/FLAECHENMITTEL-
C                                 WERTES UND SUMMATION IN DAS SU-FELD
C
         IF(ISELEP(1,  4) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SU     ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(    4) = ISAMPP(    4) + 1
         ENDIF
         IF(ISELEP(1,  4) .GE. 2) THEN
C
C                                 DIE MOMENTANWERTE DES U-FELDES WERDEN
C                                 FUER GRAPHISCHE ZWECKE RAUSGESCHRIEBEN
C
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
C
C                                 BERECHNUNG DER U-FLUKTUATIONEN
C
      IF(ISELEP(1,  6) .GE. 1) THEN
C
C                                 BPHI0 MUSS FUER UFG-FELD NICHT AUFGE-
C                                 RUFEN WERDEN, DA UFG (IM GEGENSATZ ZU
C                                 "HILF") FUER KEINE ANDEREN ZWECKE
C                                 VERWENDET WIRD.
C
         CALL PHIFLU  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,UGRID,
     $                 U      ,UFG    ,AU     ,'   U',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C
C                                 SETZEN DER RANDBEDINGUNGEN
C
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,UFG    ,IB1,IB2,JB1,JB2,KB,'UF G  ')
         IF(ISELEP(1,  6) .GE. 2) THEN
*           CALL DOBG    ( UFG,...
         ENDIF
      ENDIF
C
C                                 SUMMATION DER MOMENTANEN QUADRA-
C                                 TE DER FLUKTUATIONEN IN DAS SURG-FELD
C
      IF(ISELEP(1,  7) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
C
C                                 BERECHNUNG DER QUADRATE DER FLUK-
C                                 TUATIONEN
C
         CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,
     $                 UFG    ,HILF   ,'   U')
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'UR G  ')
         IF(ISELEP(1,  7) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SURG   ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(    7) = ISAMPP(    7) + 1
         ENDIF
         IF(ISELEP(1,  7) .GE. 2) THEN
C                                 VOR DER AUSGABE HILF = SQRT(HILF) !!!!
C           CALL WURZEL ( ....
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
C
C                                 SUMMATION DER MOMENTANEN U-SKE-
C                                 WERTE IN DAS SUSKG-FELD
C
      IF(ISELEP(1, 10) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
C
C                                 BERECHNUNG DER SKEWNESS-WERTE
C
         CALL PHISKE  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,
     $                 UFG    ,AURG   ,HILF   ,'   U',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'USKG  ')
         IF(ISELEP(1, 10) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SUSKG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(   10) = ISAMPP(   10) + 1
         ENDIF
         IF(ISELEP(1, 10) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
C
C                                 SUMMATION DER MOMENTANEN U-FLA-
C                                 WERTE IN DAS SUFLG-FELD
C
      IF(ISELEP(1, 12) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
C
C                                 BERECHNUNG DER FLATNESS-WERTE
C
         CALL PHIFLA  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,
     $                 UFG    ,AURG   ,HILF   ,'   U',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'UFLG  ')
         IF(ISELEP(1, 12) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SUFLG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(   12) = ISAMPP(   12) + 1
         ENDIF
         IF(ISELEP(1, 12) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
C
C                                 STATIST. GROESSEN DER V-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 29) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
            CALL GALIRT  (KK,JJ,II,KMX,JMX,IMX,
     $                    UGRID,V,HILF,B,'   V',NFRO,NRGT,NBOT,
     $                                          NBAC,NLFT,NTOP)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'V     ')
         IF(ISELEP(1, 29) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SV     ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(   29) = ISAMPP(   29) + 1
         ENDIF
         IF(ISELEP(1, 29) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
      IF(ISELEP(1, 31) .GE. 1) THEN
         CALL PHIFLU  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,UGRID,
     $                 V      ,VFG    ,AV     ,'   V',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,VFG    ,IB1,IB2,JB1,JB2,KB,'VF G  ')
         IF(ISELEP(1, 31) .GE. 2) THEN
*           CALL DOBG    ( VFG,...
         ENDIF
      ENDIF
      IF(ISELEP(1, 32) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,
     $                 VFG    ,HILF   ,'   V')
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'VR G  ')
         IF(ISELEP(1, 32) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SVRG   ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(   32) = ISAMPP(   32) + 1
         ENDIF
         IF(ISELEP(1, 32) .GE. 2) THEN
C                                 VOR DER AUSGABE HILF = SQRT(HILF) !!!!
C           CALL WURZEL ( ....
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
      IF(ISELEP(1, 35) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL PHISKE  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,
     $                 VFG    ,AVRG   ,HILF   ,'   V',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'VSKG  ')
         IF(ISELEP(1, 35) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SVSKG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(   35) = ISAMPP(   35) + 1
         ENDIF
         IF(ISELEP(1, 35) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
      IF(ISELEP(1, 37) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL PHIFLA  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,
     $                 VFG    ,AVRG   ,HILF   ,'   V',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'VFLG  ')
         IF(ISELEP(1, 37) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SVFLG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(   37) = ISAMPP(   37) + 1
         ENDIF
         IF(ISELEP(1, 37) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
C
C                                 STATIST. GROESSEN DER W-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 54) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
            CALL GALIRT  (KK,JJ,II,KMX,JMX,IMX,
     $                    UGRID,W,HILF,B,'   W',NFRO,NRGT,NBOT,
     $                                          NBAC,NLFT,NTOP)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'W     ')
         IF(ISELEP(1, 54) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SW     ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(   54) = ISAMPP(   54) + 1
         ENDIF
         IF(ISELEP(1, 54) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
      IF(ISELEP(1, 56) .GE. 1) THEN
         CALL PHIFLU  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,UGRID,
     $                 W      ,WFG    ,AW     ,'   W',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,WFG    ,IB1,IB2,JB1,JB2,KB,'WF G  ')
         IF(ISELEP(1, 56) .GE. 2) THEN
*           CALL DOBG    ( WFG,...
         ENDIF
      ENDIF
      IF(ISELEP(1, 57) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,
     $                 WFG    ,HILF   ,'   W')
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'WR G  ')
         IF(ISELEP(1, 57) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SWRG   ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(   57) = ISAMPP(   57) + 1
         ENDIF
         IF(ISELEP(1, 57) .GE. 2) THEN
C                                 VOR DER AUSGABE HILF = SQRT(HILF) !!!!
C           CALL WURZEL ( ....
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
      IF(ISELEP(1, 60) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL PHISKE  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,
     $                 WFG    ,AWRG   ,HILF   ,'   W',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'WSKG  ')
         IF(ISELEP(1, 60) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SWSKG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(   60) = ISAMPP(   60) + 1
         ENDIF
         IF(ISELEP(1, 60) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
      IF(ISELEP(1, 62) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL PHIFLA  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,
     $                 WFG    ,AWRG   ,HILF   ,'   W',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'WFLG  ')
         IF(ISELEP(1, 62) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SWFLG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(   62) = ISAMPP(   62) + 1
         ENDIF
         IF(ISELEP(1, 62) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
C
C                                 STATIST. GROESSEN DER P-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 79) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
C
C                                 REDUZIERUNG DES DRUCKNIVEAUS
C                                 WIRD IM HAUPTPROGRAMM GEMACHT
C
C        CALL PLEVEL  (KK,JJ,II,KMX,JMX,IMX,KPP,JPP,IPP,P,
C    $                 DDX, DDY, DDZ )
            CALL GALIRT  (KK,JJ,II,KMX,JMX,IMX,
     $                    UGRID,P,HILF,B,'   P',NFRO,NRGT,NBOT,
     $                                          NBAC,NLFT,NTOP)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'P     ')
         IF(ISELEP(1, 79) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SP     ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(   79) = ISAMPP(   79) + 1
         ENDIF
         IF(ISELEP(1, 79) .GE. 2) THEN
*           CALL DOBG    ( P,...
         ENDIF
      ENDIF
      IF(ISELEP(1, 81) .GE. 1) THEN
         CALL PHIFLU  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,UGRID,
     $                 P      ,PFG    ,AP     ,'   P',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,PFG    ,IB1,IB2,JB1,JB2,KB,'PF G  ')
         IF(ISELEP(1, 81) .GE. 2) THEN
*           CALL DOBG    ( PFG,...
         ENDIF
      ENDIF
      IF(ISELEP(1, 82) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,
     $                 PFG    ,HILF   ,'   P')
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'PR G  ')
         IF(ISELEP(1, 82) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SPRG   ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(   82) = ISAMPP(   82) + 1
         ENDIF
         IF(ISELEP(1, 82) .GE. 2) THEN
C                                 VOR DER AUSGABE HILF = SQRT(HILF) !!!!
C           CALL WURZEL ( ....
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
      IF(ISELEP(1, 85) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL PHISKE  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,
     $                 PFG    ,APRG   ,HILF   ,'   P',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'PSKG  ')
         IF(ISELEP(1, 85) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SPSKG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(   85) = ISAMPP(   85) + 1
         ENDIF
         IF(ISELEP(1, 85) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
      IF(ISELEP(1, 87) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL PHIFLA  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,
     $                 PFG    ,APRG   ,HILF   ,'   P',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'PFLG  ')
         IF(ISELEP(1, 87) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SPFLG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(   87) = ISAMPP(   87) + 1
         ENDIF
         IF(ISELEP(1, 87) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
C
C                                 KINETISCHE ENERGIE DER SCHWANKUNGS-
C                                 GESCHWINDIGKEITEN
C                                 -----------------------------------
C
      IF(ISELEP(1,106) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
C
C                                 BERECHNUNG DER KINETISCHEN
C                                 ENERGIE DER SCHWANKUNGSGESCHW.
C                                 (ANTEIL DER GROBSTRUKTUR)
C
         CALL ENERFG  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,
     $                 DX,DY,DZ,X,Y,Z,UFG,VFG,WFG,HILF,ESUM,
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,0)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'EF G  ')
         IF(ISELEP(1,106) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SEFG   ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  106) = ISAMPP(  106) + 1
         ENDIF
         IF(ISELEP(1,106) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
      IF(ISELEP(1,109) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
C
C                                 BERECHNUNG DER KINETISCHEN
C                                 ENERGIE DER SCHWANKUNGSGESCHW.
C                                 (ANTEIL DER FEINSTRUKTUR)
C
         CALL ENERFS  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,
     $                 DX,DY,DZ,X,Y,Z,B,G,GMOL,RHO,HILF,ESUM,
     $                 IC1,IC2,JC1,JC2,KC1,KC2,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,0)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'EF S  ')
         IF(ISELEP(1,109) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SEFS   ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  109) = ISAMPP(  109) + 1
         ENDIF
         IF(ISELEP(1,109) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
C                                 DISSIPATION
C                                 ___________
      IF(ISELEP(1,120) .GE. 1) THEN
C
C
            CALL DISSIPG (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,
     $                    DX,DY,DZ,X,Y,Z,UFG,VFG,WFG,B,HILF,DSUM,
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,0)
         IF(ISELEP(1,120) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SDFG   ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  120) = ISAMPP(  120) + 1
         ENDIF
      ENDIF
C                               FUER   TAYLOR''SCHES MIKROMASS
C                               QUADRATE DER ELEMENTE DES 
C                               DEFORMATIONSGESCHWINDIGKEITSTENSORS
C                               -----------------------------------
C
C                                  U_KOMPONENTE IN X_RICHTUNG
      IF(ISELEP(1,131) .GE. 1) THEN
C
            CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   )
            CALL PHIDIDJ (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,HILF,
     $                    UFG,X,DX,DDX,II,'XX',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
             CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                     KMXA,JMXA,IMXA,
     $                     HILF   ,HILF   ,'  QP')
         IF(ISELEP(1,131) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SDUDX2 ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  131) = ISAMPP(  131) + 1
         ENDIF
      ENDIF
C
C                                  U_KOMPONENTE IN Y_RICHTUNG
      IF(ISELEP(1,132) .GE. 1) THEN
C
            CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   )
            CALL PHIDIDJ (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,HILF,
     $                    UFG,Y,DY,DDY,JJ,'XY',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
             CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                     KMXA,JMXA,IMXA,
     $                     HILF   ,HILF   ,'  VU')
         IF(ISELEP(1,132) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SDUDY2 ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  132) = ISAMPP(  132) + 1
         ENDIF
      ENDIF
C
C                                  U_KOMPONENTE IN Z_RICHTUNG
      IF(ISELEP(1,133) .GE. 1) THEN
C
            CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   )
            CALL PHIDIDJ (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,HILF,
     $                    UFG,Z,DZ,DDZ,KK,'XZ',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
             CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                     KMXA,JMXA,IMXA,
     $                     HILF   ,HILF   ,'  UW')
         IF(ISELEP(1,133) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SDUDZ2 ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  133) = ISAMPP(  133) + 1
         ENDIF
      ENDIF
C
C                                  V_KOMPONENTE IN X_RICHTUNG
      IF(ISELEP(1,134) .GE. 1) THEN
C
            CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   )
            CALL PHIDIDJ (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,HILF,
     $                    VFG,X,DX,DDX,II,'YX',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
             CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                     KMXA,JMXA,IMXA,
     $                     HILF   ,HILF   ,'  VU')
         IF(ISELEP(1,134) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SDVDX2 ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  134) = ISAMPP(  134) + 1
         ENDIF
      ENDIF
C
C                                  V_KOMPONENTE IN Y_RICHTUNG
      IF(ISELEP(1,135) .GE. 1) THEN
C
            CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   )
            CALL PHIDIDJ (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,HILF,
     $                    VFG,Y,DY,DDY,JJ,'YY',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
             CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                     KMXA,JMXA,IMXA,
     $                     HILF   ,HILF   ,'  QP')
         IF(ISELEP(1,135) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SDVDY2 ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  135) = ISAMPP(  135) + 1
         ENDIF
      ENDIF
C
C                                  V_KOMPONENTE IN Z_RICHTUNG
      IF(ISELEP(1,136) .GE. 1) THEN
C
            CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   )
            CALL PHIDIDJ (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,HILF,
     $                    VFG,Z,DZ,DDZ,KK,'YZ',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
             CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                     KMXA,JMXA,IMXA,
     $                     HILF   ,HILF   ,'  WV')
         IF(ISELEP(1,136) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SDVDZ2 ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  136) = ISAMPP(  136) + 1
         ENDIF
      ENDIF
C
C                                  W_KOMPONENTE IN X_RICHTUNG
      IF(ISELEP(1,137) .GE. 1) THEN
C
            CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   )
            CALL PHIDIDJ (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,HILF,
     $                    WFG,X,DX,DDX,II,'ZX',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
             CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                     KMXA,JMXA,IMXA,
     $                     HILF   ,HILF   ,'  UW')
         IF(ISELEP(1,137) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SDWDX2 ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  137) = ISAMPP(  137) + 1
         ENDIF
      ENDIF
C
C                                  W_KOMPONENTE IN Y_RICHTUNG
      IF(ISELEP(1,138) .GE. 1) THEN
C
            CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   )
            CALL PHIDIDJ (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,HILF,
     $                    WFG,Y,DY,DDY,JJ,'ZY',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
             CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                     KMXA,JMXA,IMXA,
     $                     HILF   ,HILF   ,'  WV')
         IF(ISELEP(1,138) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SDWDY2 ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  138) = ISAMPP(  138) + 1
         ENDIF
      ENDIF
C
C                                  W_KOMPONENTE IN Z_RICHTUNG
      IF(ISELEP(1,139) .GE. 1) THEN
C
            CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   )
            CALL PHIDIDJ (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,HILF,
     $                    WFG,Z,DZ,DDZ,KK,'ZZ',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
             CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                     KMXA,JMXA,IMXA,
     $                     HILF   ,HILF   ,'  QP')
         IF(ISELEP(1,139) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SDWDZ2 ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  139) = ISAMPP(  139) + 1
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
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
C
C                                 BERECHNUNG DES GROBSTRUKTUR-
C                                 ANTEILS AN DEN SCHUBSPANNUNGEN
C
         CALL TAUFG   (KK,JJ,II,KMX,JMX,IMX,GMOL,RHO,
     $                 HILF   ,UFG    ,WFG    ,'UW',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'UW    ')
         IF(ISELEP(1,146) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SUFWFG ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  146) = ISAMPP(  146) + 1
         ENDIF
         IF(ISELEP(1,146) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
      IF(ISELEP(1,149) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
C
C                                 BERECHNUNG DES FEINSTRUKTUR-
C                                 ANTEILS AN DEN SCHUBSPANNUNGEN
C
         CALL TAUFS   (KK,JJ,II,KMX,JMX,IMX,IB1,IB2,JB1,JB2,KB,
     $                 G,GMOL,RHO,UGRID,
     $                 HILF   ,U,X,DX,DDX,II,W,Z,DZ,DDZ,KK,'UW'
     $                    )      
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'UW    ')
         IF(ISELEP(1,149) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SUFWFS ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  149) = ISAMPP(  149) + 1
         ENDIF
         IF(ISELEP(1,149) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
      IF(ISELEP(1,151) .GE. 1) THEN
         CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,HILF   )
C
C                                 BERECHNUNG DES MOLEKULAREN
C                                 ANTEILS AN DEN SCHUBSPANNUNGEN
C                                 (EINSCHL. DER WANDSCHUBSPANNUNGEN)
C
         CALL TAUFM   (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                 B,G,GMOL,RHO,UGRID,
     $                 HILF   ,U,X,DX,DDX,II,W,Z,DZ,DDZ,KK,'UW',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB
     $                    )     
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'UW    ')
         IF(ISELEP(1,151) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SUFWFM ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  151) = ISAMPP(  151) + 1
         ENDIF
         IF(ISELEP(1,151) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
C
C                                 HIER: V - W
C                                 -----------
C
      IF(ISELEP(1,156) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL TAUFG   (KK,JJ,II,KMX,JMX,IMX,GMOL,RHO,
     $                 HILF   ,VFG    ,WFG    ,'VW',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'VW    ')
         IF(ISELEP(1,156) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SVFWFG ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  156) = ISAMPP(  156) + 1
         ENDIF
         IF(ISELEP(1,156) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
      IF(ISELEP(1,159) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL TAUFS   (KK,JJ,II,KMX,JMX,IMX,IB1,IB2,JB1,JB2,KB,
     $                 G,GMOL,RHO,UGRID,
     $                 HILF   ,V,Y,DY,DDY,JJ,W,Z,DZ,DDZ,KK,'VW'
     $                    )      
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'VW    ')
         IF(ISELEP(1,159) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SVFWFS ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  159) = ISAMPP(  159) + 1
         ENDIF
         IF(ISELEP(1,159) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
      IF(ISELEP(1,161) .GE. 1) THEN
         CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,HILF   )
         CALL TAUFM   (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                 B,G,GMOL,RHO,UGRID,
     $                 HILF   ,V,Y,DY,DDY,JJ,W,Z,DZ,DDZ,KK,'VW',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB
     $                    )     
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'VW    ')
         IF(ISELEP(1,161) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SVFWFM ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  161) = ISAMPP(  161) + 1
         ENDIF
         IF(ISELEP(1,161) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
C
C                                 HIER: U - V
C                                 -----------
C
      IF(ISELEP(1,166) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL TAUFG   (KK,JJ,II,KMX,JMX,IMX,GMOL,RHO,
     $                 HILF   ,UFG    ,VFG    ,'UV',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'UV    ')
         IF(ISELEP(1,166) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SUFVFG ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  166) = ISAMPP(  166) + 1
         ENDIF
         IF(ISELEP(1,166) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
      IF(ISELEP(1,169) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL TAUFS   (KK,JJ,II,KMX,JMX,IMX,IB1,IB2,JB1,JB2,KB,
     $                 G,GMOL,RHO,UGRID,
     $                 HILF   ,U,X,DX,DDX,II,V,Y,DY,DDY,JJ,'UV'
     $                    )      
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'UV    ')
         IF(ISELEP(1,169) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SUFVFS ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  169) = ISAMPP(  169) + 1
         ENDIF
         IF(ISELEP(1,169) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
      IF(ISELEP(1,171) .GE. 1) THEN
         CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,HILF   )
         CALL TAUFM   (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                 B,G,GMOL,RHO,UGRID,
     $                 HILF   ,U,X,DX,DDX,II,V,Y,DY,DDY,JJ,'UV',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB
     $                    )     
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'UV    ')
         IF(ISELEP(1,171) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SUFVFM ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  171) = ISAMPP(  171) + 1
         ENDIF
         IF(ISELEP(1,171) .GE. 2) THEN
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
C
C                                 X-KOMPONENTE DER VORTICITY
C                                 ------------------------
C
      IF(ISELEP(1,176) .GE. 1) THEN
C
C                                 BERECHNUNG DER MOMENTANWERTE
C
         CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,OX     )
         CALL OMEGA   (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                 G, GMOL, RHO,    UGRID,
     $                 OX,W ,Z ,DZ ,DDZ ,KK,V ,Y ,DY ,DDY ,JJ,'OX',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,OX     ,IB1,IB2,JB1,JB2,KB,'WV    ')
         IF(ISELEP(1,176) .GE. 2) THEN
*           CALL DOBG    ( OX  ,...
         ENDIF
      ENDIF
C
      IF(ISELEP(1,177) .GE. 1) THEN
         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,OX     ,SOX    ,
     $                 XHOMOG,YHOMOG,ZHOMOG)
         ISAMPP(  177) = ISAMPP(  177) + 1
      ENDIF
C
      IF(ISELEP(1,179) .GE. 1) THEN
         CALL PHIFLU  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,UGRID,
     $                 OX     ,OXFG   ,AOX    ,'  WV',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,OXFG   ,IB1,IB2,JB1,JB2,KB,'WV    ')
         IF(ISELEP(1,179) .GE. 2) THEN
*           CALL DOBG    ( OXFG,...
         ENDIF
      ENDIF
      IF(ISELEP(1,180) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,
     $                 OXFG   ,HILF   ,'  WV')
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'WV    ')
         IF(ISELEP(1,180) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SOXRG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  180) = ISAMPP(  180) + 1
         ENDIF
         IF(ISELEP(1,180) .GE. 2) THEN
C                                 VOR DER AUSGABE HILF = SQRT(HILF) !!!!
C           CALL WURZEL ( ....
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
C
C                                 Y-KOMPONENTE DER VORTICITY
C                                 ------------------------
C
      IF(ISELEP(1,182) .GE. 1) THEN
         CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,OY     )
         CALL OMEGA   (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                 G, GMOL, RHO,    UGRID,
     $                 OY,U ,X ,DX ,DDX ,II,W ,Z ,DZ ,DDZ ,KK,'OY',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,OY     ,IB1,IB2,JB1,JB2,KB,'UW    ')
         IF(ISELEP(1,182) .GE. 2) THEN
*           CALL DOBG    ( OY  ,...
         ENDIF
      ENDIF
      IF(ISELEP(1,183) .GE. 1) THEN
         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,OY     ,SOY    ,
     $                 XHOMOG,YHOMOG,ZHOMOG)
         ISAMPP(  183) = ISAMPP(  183) + 1
      ENDIF
      IF(ISELEP(1,185) .GE. 1) THEN
         CALL PHIFLU  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,UGRID,
     $                 OY     ,OYFG   ,AOY    ,'  UW',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,OYFG   ,IB1,IB2,JB1,JB2,KB,'UW    ')
         IF(ISELEP(1,185) .GE. 2) THEN
*           CALL DOBG    ( OYFG,...
         ENDIF
      ENDIF
      IF(ISELEP(1,186) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,
     $                 OYFG   ,HILF   ,'  UW')
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'UW    ')
         IF(ISELEP(1,186) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SOYRG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  186) = ISAMPP(  186) + 1
         ENDIF
         IF(ISELEP(1,186) .GE. 2) THEN
C                                 VOR DER AUSGABE HILF = SQRT(HILF) !!!!
C           CALL WURZEL ( ....
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
C
C                                 Z-KOMPONENTE DER VORTICITY
C                                 ------------------------
C
      IF(ISELEP(1,188) .GE. 1) THEN
         CALL DPHI0   (KK ,JJ ,II ,KMX ,JMX ,IMX ,OZ     )
         CALL OMEGA   (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                 G, GMOL, RHO,    UGRID,
     $                 OZ,V ,Y ,DY ,DDY ,JJ,U ,X ,DX ,DDX ,II,'OZ',
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,OZ     ,IB1,IB2,JB1,JB2,KB,'VU    ')
         IF(ISELEP(1,188) .GE. 2) THEN
*           CALL DOBG    ( OZ  ,...
         ENDIF
      ENDIF
      IF(ISELEP(1,189) .GE. 1) THEN
         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,OZ     ,SOZ    ,
     $                 XHOMOG,YHOMOG,ZHOMOG)
         ISAMPP(  189) = ISAMPP(  189) + 1
      ENDIF
      IF(ISELEP(1,191) .GE. 1) THEN
         CALL PHIFLU  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,UGRID,
     $                 OZ     ,OZFG   ,AOZ    ,'  VU',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,OZFG   ,IB1,IB2,JB1,JB2,KB,'VU    ')
         IF(ISELEP(1,191) .GE. 2) THEN
*           CALL DOBG    ( OZFG,...
         ENDIF
      ENDIF
      IF(ISELEP(1,192) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,
     $                 OZFG   ,HILF   ,'  VU')
C        CALL BF1ALL  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,DX,DY,DZ,
C    $                 X,Y,Z,HILF   ,IB1,IB2,JB1,JB2,KB,'VU    ')
         IF(ISELEP(1,192) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SOZRG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  192) = ISAMPP(  192) + 1
         ENDIF
         IF(ISELEP(1,192) .GE. 2) THEN
C                                 VOR DER AUSGABE HILF = SQRT(HILF) !!!!
C           CALL WURZEL ( ....
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
C
C                                 ENSTROPHIE
C                                 ----------
C
      IF(ISELEP(1,195) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL ENSTRO  (KK,JJ,II,KMX,JMX,IMX,OX,OY,OZ,HILF  )
         IF(ISELEP(1,195) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SO2    ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  195) = ISAMPP(  195) + 1
         ENDIF
         IF(ISELEP(1,195) .GE. 2) THEN
*           CALL DOBG    ( HILF ,...
         ENDIF
      ENDIF
      IF(ISELEP(1,197) .GE. 1) THEN
         CALL PHIFLU  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,UGRID,
     $                 HILF   ,O2FG   ,AO2    ,' UVW',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
         IF(ISELEP(1,197) .GE. 2) THEN
*           CALL DOBG    ( O2FG,...
         ENDIF
      ENDIF
      IF(ISELEP(1,198) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,
     $                 O2FG   ,HILF   ,' UVW')
         IF(ISELEP(1,198) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SO2RG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  198) = ISAMPP(  198) + 1
         ENDIF
         IF(ISELEP(1,198) .GE. 2) THEN
C                                 VOR DER AUSGABE HILF = SQRT(HILF) !!!!
C           CALL WURZEL ( ....
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
C
C                                 HELIZITAET
C                                 ----------
C
      IF(ISELEP(1,200) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL HELICI  (KK,JJ,II,KMX,JMX,IMX,U,V,W,OX,OY,OZ,
     $                 HILF  )
         IF(ISELEP(1,200) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SHE    ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  200) = ISAMPP(  200) + 1
         ENDIF
         IF(ISELEP(1,200) .GE. 2) THEN
*           CALL DOBG    ( HILF ,...
         ENDIF
      ENDIF
      IF(ISELEP(1,202) .GE. 1) THEN
         CALL PHIFLU  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,UGRID,
     $                 HILF   ,HEFG   ,AHE    ,' UVW',
     $                 XHOMOG,YHOMOG,ZHOMOG,NFRO,NRGT)
         IF(ISELEP(1,202) .GE. 2) THEN
*           CALL DOBG    ( HEFG,...
         ENDIF
      ENDIF
      IF(ISELEP(1,203) .GE. 1) THEN
         CALL BPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF   ,
     $                 IB1,IB2,JB1,JB2,KB)
         CALL PHIFL2  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                 KMXA,JMXA,IMXA,
     $                 HEFG   ,HILF   ,' UVW')
         IF(ISELEP(1,203) .LE. 2) THEN
            CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF   ,SHERG  ,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  203) = ISAMPP(  203) + 1
         ENDIF
         IF(ISELEP(1,203) .GE. 2) THEN
C                                 VOR DER AUSGABE HILF = SQRT(HILF) !!!!
C           CALL WURZEL ( ....
*           CALL DOBG    ( HILF,...
         ENDIF
      ENDIF
C
C
C                                "NEUE" Statistik
C                                 ---------------
C

      IF (ISELEP(1,400) .GE. 1) THEN
         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D1)
         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $                 HILF3D1,'XX',U,U,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF3D1,SUUM,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  400) = ISAMPP(  400) + 1
      ENDIF

      IF (ISELEP(1,401) .GE. 1) THEN
         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D1)
         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $                 HILF3D1,'XX',V,V,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF3D1,SVVM,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  401) = ISAMPP(  401) + 1
      ENDIF

      IF (ISELEP(1,402) .GE. 1) THEN
         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D1)
         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $                 HILF3D1,'XX',W,W,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF3D1,SWWM,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  402) = ISAMPP(  402) + 1
      ENDIF

      IF (ISELEP(1,406) .GE. 1) THEN
         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D1)
         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $                 HILF3D1,'XX',P,P,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF3D1,SPPM,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  406) = ISAMPP(  406) + 1
      ENDIF


      IF (ISELEP(1,403) .GE. 1) THEN
         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D1)
         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $                 HILF3D1,'UV',U,V,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF3D1,SUVM,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  403) = ISAMPP(  403) + 1
      ENDIF

      IF (ISELEP(1,404) .GE. 1) THEN
         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D1)
         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $                 HILF3D1,'UW',U,W,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF3D1,SUWM,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  404) = ISAMPP(  404) + 1
      ENDIF

      IF (ISELEP(1,405) .GE. 1) THEN
         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D1)
         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $                 HILF3D1,'VW',V,W,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,HILF3D1,SVWM,
     $                    XHOMOG,YHOMOG,ZHOMOG)
            ISAMPP(  405) = ISAMPP(  405) + 1

      ENDIF



C     QUADRATE DER ABLEITUNGEN VON U

      IF (ISELEP(1,430) .GE. 1) THEN

         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D1)
         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D2)

         CALL DFDX  (KK,JJ,II,RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $        'DDX',U,HILF3D1)

         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $        HILF3D2,'XX',HILF3D1,HILF3D1,
     $        NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $        KMXA,JMXA,IMXA,HILF3D2,SUXUXM,
     $        XHOMOG,YHOMOG,ZHOMOG)

         ISAMPP(  430) = ISAMPP(  430) + 1

      ENDIF

      IF (ISELEP(1,431) .GE. 1) THEN

         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D1)
         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D2)

         CALL DFDX  (KK,JJ,II,RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $        'DYS',U,HILF3D1)

         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $        HILF3D2,'XX',HILF3D1,HILF3D1,
     $        NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)

         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $        KMXA,JMXA,IMXA,HILF3D2,SUYUYM,
     $        XHOMOG,YHOMOG,ZHOMOG)

         ISAMPP(  431) = ISAMPP(  431) + 1

      ENDIF

      IF (ISELEP(1,432) .GE. 1) THEN

         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D1)
         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D2)
         CALL DFDX  (KK,JJ,II,RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $        'DZS',U,HILF3D1)
         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $        HILF3D2,'XX',HILF3D1,HILF3D1,
     $        NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $        KMXA,JMXA,IMXA,HILF3D2,SUZUZM,
     $        XHOMOG,YHOMOG,ZHOMOG)
         ISAMPP(  432) = ISAMPP(  432) + 1
      ENDIF

C     QUADRATE DER ABLEITUNGEN VON V

      IF (ISELEP(1,433) .GE. 1) THEN

         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D1)
         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D2)

         CALL DFDX  (KK,JJ,II,RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $        'DXS',V,HILF3D1)

         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $        HILF3D2,'XX',HILF3D1,HILF3D1,
     $        NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)

         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $        KMXA,JMXA,IMXA,HILF3D2,SVXVXM,
     $        XHOMOG,YHOMOG,ZHOMOG)

         ISAMPP(  433) = ISAMPP(  433) + 1

      ENDIF

      IF (ISELEP(1,434) .GE. 1) THEN

         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D1)
         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D2)

         CALL DFDX  (KK,JJ,II,RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $        'DDY',V,HILF3D1)

         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $        HILF3D2,'XX',HILF3D1,HILF3D1,
     $        NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)

         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $        KMXA,JMXA,IMXA,HILF3D2,SVYVYM,
     $        XHOMOG,YHOMOG,ZHOMOG)

         ISAMPP(  434) = ISAMPP(  434) + 1

      ENDIF

      IF (ISELEP(1,435) .GE. 1) THEN

         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D1)
         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D2)

         CALL DFDX  (KK,JJ,II,RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $        'DZS',V,HILF3D1)

         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $        HILF3D2,'XX',HILF3D1,HILF3D1,
     $        NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)

         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $        KMXA,JMXA,IMXA,HILF3D2,SVZVZM,
     $        XHOMOG,YHOMOG,ZHOMOG)
         ISAMPP(  435) = ISAMPP(  435) + 1
      ENDIF

C     QUADRATE DER ABLEITUNGEN VON W

      IF (ISELEP(1,436) .GE. 1) THEN

         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D1)
         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D2)

         CALL DFDX  (KK,JJ,II,RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $        'DXS',W,HILF3D1)

         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $        HILF3D2,'XX',HILF3D1,HILF3D1,
     $        NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)

         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $        KMXA,JMXA,IMXA,HILF3D2,SWXWXM,
     $        XHOMOG,YHOMOG,ZHOMOG)

         ISAMPP(  436) = ISAMPP(  436) + 1

      ENDIF

      IF (ISELEP(1,437) .GE. 1) THEN

         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D1)
         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D2)

         CALL DFDX  (KK,JJ,II,RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $        'DYS',W,HILF3D1)

         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $        HILF3D2,'XX',HILF3D1,HILF3D1,
     $        NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)

         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $        KMXA,JMXA,IMXA,HILF3D2,SWYWYM,
     $        XHOMOG,YHOMOG,ZHOMOG)

         ISAMPP(  437) = ISAMPP(  437) + 1

      ENDIF

      IF (ISELEP(1,438) .GE. 1) THEN

         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D1)
         CALL DPHI0   (KK,JJ,II,KMX,JMX,IMX,HILF3D2)

            CALL DFDX  (KK,JJ,II,RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $           'DDZ',W,HILF3D1)

         CALL PHIMLT2 (KK,JJ,II,KMX,JMX,IMX,
     $        HILF3D2,'XX',HILF3D1,HILF3D1,
     $        NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)

         CALL AREAM   (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $        KMXA,JMXA,IMXA,HILF3D2,SWZWZM,
     $        XHOMOG,YHOMOG,ZHOMOG)
         ISAMPP(  438) = ISAMPP(  438) + 1
      ENDIF





      RETURN
      END
