










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
      SUBROUTINE PRLE3E  (LP1,LPS,LP2,EBE,MP1,MPS,MP2,NP1,NPS,NP2,
     $                    DX,DY,DZ,X,Y,Z,KANAL,NRRUN,ITFLUC,DT,TREF,
     $                   XHOMOG,YHOMOG,ZHOMOG,

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
C        P R L E 3 E      KOMPAKTAUSGABE DER ENSEMBLE-MITTELWERTE
C                         FOLGENDER FELDER:
C                         <PHIMEAN>        (U,V,W,T,P)
C                         <PHIRMS > G      (U,V,W,T,P)
C                         <PHIRMS > G + S  (U,V,W,T  )
C                         <E>       G
C                         <E>       G + S
C                         <E>       S
C                         <TAU>     G      (XZ-, YZ- UND XY-KOMP.)
C                         <TAU>     G + S  (XZ-, YZ- UND XY-KOMP.)
C                         <TAU>     G+S+M  (XZ-, YZ- UND XY-KOMP.)
C                         <PHISKEW> G      (U,V,W,P,T)
C                         <PHIFLAT> G      (U,V,W,P,T)
C  A C H T U N G:         DIESE ROUTINE MUSS ERWEITERT WERDEN, FALLS
C                         DIE STATIST. AUSWERTUNG ERWEITERT WIRD !
C*STARLET***************************************************************
C
C PARAM: LP1,LPS,LP2    - ANFANGSWERT, SCHRITTWEITE, ENDWERT DER EBENEN
C                         'EBE'=CONST, DIE AUSGEGEBEN WERDEN SOLLEN
C        EBE            - CHARACTER-VARIABLE DER FORM 'K', 'J' ODER 'I'
C        MP1,MPS,MP2        - STEUERPARAMETER F. PRINTAUSGABE
C        NP1,NPS,NP2        -        ""        "       "
C        DX(II)             - ABSTAND DER BASISZELLBEZUGSPUNKTE IN X
C        DY(JJ)             -   ""      "            "           IN Y
C        DZ(KK)             -   ""      "            "           IN Z
C        X(II),Y(JJ),Z(KK)  - KOORDINATEN D. BASISZELLBEZUGSPUNKTE
C        NRRUN          - KENNZIFFER DES LAUFES
C        ITFLUC         - NACH JEWEILS ITFLUC ZEITSCHRITTEN WIRD
C                         EINE STICHPROBE ENTNOMMEN
C        DT             - ZEITSCHRITT
C        TREF           - BEZUGSZEIT
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
C UPROG                 : PR3V, ERRR
C
C DEFINE-DIREKTIVEN     : XHOMOG, YHOMOG
C
C        02.09.86 (HW)  : ORIGINAL
C        30.06.88 (HW)  : ERWEITERUNG: AUSWERTUNG DER DREI KOMPONENTEN
C                         DER VORTICITY, DER ENSTROPHIE UND DER
C                         HELIZITAET
C        19.08.88 (HW)  : VORTICITY IST JETZT IM DOPPELT VERSETZTEN
C                         MASCHENGITTER DEFINIERT
C        17.12.88 (HW)  : COMMON-DECKS FUER DIE DIMENSIONIERUNG (CREFVD)
C                         UND BELEGUNG (CREFVC) DER REFERENZGROESSEN
C                         EINGEFUEHRT
C        03.01.89 (HW)  : DIE HELIZITAET IST JETZT DREIFACH IM MASCHEN-
C                         GITTER VERSCHOBEN (IN X-, Y- UND Z-RICHTUNG)
C        25.02.92 (MM)  : ACHTUNG VARIABLE "KANAL" IN UEBERGABE
C                         EIGEFUEGT
C        12.02.03 (TB)  : SCALAR FIELD STATISTICS IMPLEMENTED
C
C*STARLET***************************************************************
C
      CHARACTER (LEN=1)   EBE

      CHARACTER (LEN=32)  CKEINE, CEREF, CGREF, CPREF, CTAURE, CUREF,
     $               COMREF, CO2REF, CHEREF, CLREF, CRLREF, CASDUI,
     $               CASDOM

C
      LOGICAL  XHOMOG,YHOMOG,ZHOMOG

      REAL      DX(II),    DY(JJ),    DZ(KK),
     $           X(II),     Y(JJ),     Z(KK)
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
      TINTVA = FLOAT(ITFLUC) * DT
      TINTD  = TINTVA        / TREF

      CKEINE = '--------------------------------'
      CLREF  = 'L-REF                           '
      CRLREF = '1.0 / L-REF                     '
      CEREF  = '(U-REF)**2                      '
      CGREF  = 'RHO * L-REF * U-REF             '
      CPREF  = '1/2 * RHO * (U-REF)**2          '
      CTAURE = 'RHO * (U-REF)**2                '
      CUREF  = 'U-REF                           '
      COMREF = '(U-REF) / (L-REF)               '
      CO2REF = '(U-REF)**2 / (L-REF)**2         '
      CHEREF = '(U-REF)**2 / (L-REF)            '
      CASDUI = '1/2 * (U-REF)**2 * (L-REF)      '
      CASDOM = '1/2 * (U-REF)**2 / (L-REF)      '
C
C                                 VOREINSTELLUNG DER STEUERPARAMETER
C
      LH1    = LP1
      LHS    = LPS
      LH2    = LP2
      MH1    = MP1
      MHS    = MPS
      MH2    = MP2
      NH1    = NP1
      NHS    = NPS
      NH2    = NP2
C
      IF(EBE .NE. 'K') GOTO 2010

      IF (XHOMOG) THEN
C
C                                 X-RICHTUNG IST HOMOGEN !
C
         MH1 = 1
         MHS = 1
         MH2 = 1

      ENDIF

      IF (YHOMOG) THEN
C
C                                 Y-RICHTUNG IST HOMOGEN !
C
         NH1 = 1
         NHS = 1
         NH2 = 1

      ENDIF

         GOTO 2100
C
 2010 IF(EBE .NE. 'J') GOTO 2020

      IF (XHOMOG) THEN
C
C                                 X-RICHTUNG IST HOMOGEN !
C
         MH1 = 1
         MHS = 1
         MH2 = 1

      ENDIF

      IF (YHOMOG) THEN
C
C                                 Y-RICHTUNG IST HOMOGEN !
C
         LH1 = 1
         LHS = 1
         LH2 = 1

      ENDIF

         GOTO 2100
C
 2020 IF(EBE .NE. 'I') CALL ERRR (501,' PRLE3E   ')

      IF (XHOMOG) THEN
C
C                                 X-RICHTUNG IST HOMOGEN !
C
         LH1 = 1
         LHS = 1
         LH2 = 1

      ENDIF

      IF (YHOMOG) THEN
C
C                                 Y-RICHTUNG IST HOMOGEN !
C
         MH1 = 1
         MHS = 1
         MH2 = 1

      ENDIF
C
 2100 CONTINUE
C
      WRITE (KANAL,6010) ITFLUC,TINTVA,TINTD
      IF (XHOMOG) THEN
      WRITE (KANAL,6020)
      ENDIF
      IF (YHOMOG) THEN
      WRITE (KANAL,6030)
      ENDIF
C
C                                 HIER: AUSGABE DER <PHIMEAN>-WERTE
C
      WRITE       (KANAL,6110)
C
      IF(ISELEP(1,  4) .EQ. 1 .OR. ISELEP(1,  4) .EQ. 2) THEN
      WRITE (KANAL,6040) ISAMPP(    4),
     $                  NRRUN,ISAMPP(    4)-ISAMPA(    4),CUREF
         CALL PR3V  (KKA,JJA,IIA,AU     ,'<U>             ', 0 , 0 , 1 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1, 29) .EQ. 1 .OR. ISELEP(1, 29) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(   29),
     $                  NRRUN,ISAMPP(   29)-ISAMPA(   29),CUREF
         CALL PR3V  (KKA,JJA,IIA,AV     ,'<V>             ', 0 , 1 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1, 54) .EQ. 1 .OR. ISELEP(1, 54) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(   54),
     $                  NRRUN,ISAMPP(   54)-ISAMPA(   54),CUREF
         CALL PR3V  (KKA,JJA,IIA,AW     ,'<W>             ', 1 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1, 79) .EQ. 1 .OR. ISELEP(1, 79) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(   79),
     $                  NRRUN,ISAMPP(   79)-ISAMPA(   79),CPREF
         CALL PR3V  (KKA,JJA,IIA,AP     ,'<P>             ', 0 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
C                                 HIER: AUSGABE DER <PHIRMS>-WERTE
C                                 (ZUERST NUR GROBSTRUKTURANTEIL UND
C                                 ANSCHLIESSEND GROB+FEINSTRUKTURANTEIL)
C
      WRITE       (KANAL,6130)
C
      IF(ISELEP(1,  7) .EQ. 1 .OR. ISELEP(1,  7) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(    7),
     $                  NRRUN,ISAMPP(    7)-ISAMPA(    7),CUREF
         CALL PR3V  (KKA,JJA,IIA,AURG   ,'<U-RMS> G       ', 0 , 0 , 1 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1,  9) .EQ. 1 .OR. ISELEP(1,  9) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(    7),
     $                  NRRUN,ISAMPP(    7)-ISAMPA(    7),CUREF
         CALL PR3V  (KKA,JJA,IIA,SURG   ,'<U-RMS> G+S     ', 0 , 0 , 1 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1, 32) .EQ. 1 .OR. ISELEP(1, 32) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(   32),
     $                  NRRUN,ISAMPP(   32)-ISAMPA(   32),CUREF
         CALL PR3V  (KKA,JJA,IIA,AVRG   ,'<V-RMS> G       ', 0 , 1 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1, 34) .EQ. 1 .OR. ISELEP(1, 34) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(   32),
     $                  NRRUN,ISAMPP(   32)-ISAMPA(   32),CUREF
         CALL PR3V  (KKA,JJA,IIA,SVRG   ,'<V-RMS> G+S     ', 0 , 1 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1, 57) .EQ. 1 .OR. ISELEP(1, 57) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(   57),
     $                  NRRUN,ISAMPP(   57)-ISAMPA(   57),CUREF
         CALL PR3V  (KKA,JJA,IIA,AWRG   ,'<W-RMS> G       ', 1 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1, 59) .EQ. 1 .OR. ISELEP(1, 59) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(   57),
     $                  NRRUN,ISAMPP(   57)-ISAMPA(   57),CUREF
         CALL PR3V  (KKA,JJA,IIA,SWRG   ,'<W-RMS> G+S     ', 1 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1, 82) .EQ. 1 .OR. ISELEP(1, 82) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(   82),
     $                  NRRUN,ISAMPP(   82)-ISAMPA(   82),CPREF
         CALL PR3V  (KKA,JJA,IIA,APRG   ,'<P-RMS> G       ', 0 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
C                                 HIER: AUSGABE DER KINETISCHEN
C                                 ENERGIE DER SCHWANKUNGSGESCHW.
C
      WRITE       (KANAL,6150)
C
      IF(ISELEP(1,106) .EQ. 1 .OR. ISELEP(1,106) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  106),
     $                  NRRUN,ISAMPP(  106)-ISAMPA(  106),CEREF
         CALL PR3V  (KKA,JJA,IIA,AEFG   ,'<E-F> G         ', 0 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1,108) .EQ. 1 .OR. ISELEP(1,108) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  106),
     $                  NRRUN,ISAMPP(  106)-ISAMPA(  106),CEREF
         CALL PR3V  (KKA,JJA,IIA,SEFG   ,'<E-F> G+S       ', 0 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1,109) .EQ. 1 .OR. ISELEP(1,109) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  109),
     $                  NRRUN,ISAMPP(  109)-ISAMPA(  109),CEREF
         CALL PR3V  (KKA,JJA,IIA,AEFS   ,'<E-F> S         ', 0 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
C                                 DISSIPATION
C                                 -----------
C
      IF(ISELEP(1,120) .EQ. 1 .OR. ISELEP(1,120) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  120),
     $                  NRRUN,ISAMPP(  120)-ISAMPA(  120),CEREF
         CALL PR3V  (KKA,JJA,IIA,ADFG   ,'<D-F> G         ', 0 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
C                               FUER   TAYLOR''SCHES MIKROMASS
C                               QUADRATE DER ELEMENTE DES 
C                               DEFORMATIONSGESCHWINDIGKEITSTENSORS
C                               -----------------------------------
C
      IF(ISELEP(1,131) .EQ. 1 .OR. ISELEP(1,131) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  131),
     $                  NRRUN,ISAMPP(  131)-ISAMPA(  131),CEREF
         CALL PR3V  (KKA,JJA,IIA,ADUDX2 ,'<DUDX2-F> G     ', 0 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
      IF(ISELEP(1,132) .EQ. 1 .OR. ISELEP(1,132) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  132),
     $                  NRRUN,ISAMPP(  132)-ISAMPA(  132),CEREF
         CALL PR3V  (KKA,JJA,IIA,ADUDY2 ,'<DUDY2-F> G     ', 0 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
      IF(ISELEP(1,133) .EQ. 1 .OR. ISELEP(1,133) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  133),
     $                  NRRUN,ISAMPP(  133)-ISAMPA(  133),CEREF
         CALL PR3V  (KKA,JJA,IIA,ADUDZ2 ,'<DUDZ2-F> G     ', 0 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
      IF(ISELEP(1,134) .EQ. 1 .OR. ISELEP(1,134) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  134),
     $                  NRRUN,ISAMPP(  134)-ISAMPA(  134),CEREF
         CALL PR3V  (KKA,JJA,IIA,ADVDY2 ,'<DVDX2-F> G     ', 0 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
      IF(ISELEP(1,135) .EQ. 1 .OR. ISELEP(1,135) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  135),
     $                  NRRUN,ISAMPP(  135)-ISAMPA(  135),CEREF
         CALL PR3V  (KKA,JJA,IIA,ADVDY2 ,'<DVDY2-F> G     ', 0 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
      IF(ISELEP(1,136) .EQ. 1 .OR. ISELEP(1,136) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  136),
     $                  NRRUN,ISAMPP(  136)-ISAMPA(  136),CEREF
         CALL PR3V  (KKA,JJA,IIA,ADVDZ2 ,'<DVDZ2-F> G     ', 0 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
      IF(ISELEP(1,137) .EQ. 1 .OR. ISELEP(1,137) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  137),
     $                  NRRUN,ISAMPP(  137)-ISAMPA(  137),CEREF
         CALL PR3V  (KKA,JJA,IIA,ADWDX2 ,'<DWDX2-F> G     ', 0 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
      IF(ISELEP(1,138) .EQ. 1 .OR. ISELEP(1,138) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  138),
     $                  NRRUN,ISAMPP(  138)-ISAMPA(  138),CEREF
         CALL PR3V  (KKA,JJA,IIA,ADWDY2 ,'<DWDY2-F> G     ', 0 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
      IF(ISELEP(1,139) .EQ. 1 .OR. ISELEP(1,139) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  139),
     $                  NRRUN,ISAMPP(  139)-ISAMPA(  139),CEREF
         CALL PR3V  (KKA,JJA,IIA,ADWDZ2 ,'<DWDZ2-F> G     ', 0 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
C
C                                 HIER: AUSGABE DER KOMPONENTEN
C                                 DES SCHUBSPANNUNGSTENSORS IN FOL-
C                                 GENDER REIHENFOLGE: 1. ANTEIL DER
C                                 GROBSTRUKTUR, 2. GROBSTRUKTUR
C                                 + FEINSTR.-ANTEIL, 3. GROBSTR.+FEIN-
C                                 STRUKTUR + MOLEKULARER ANTEIL
C
C                                 HIER: U - W
C                                 -----------
C
      WRITE       (KANAL,6170)
C
      IF(ISELEP(1,146) .EQ. 1 .OR. ISELEP(1,146) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  146),
     $                  NRRUN,ISAMPP(  146)-ISAMPA(  146),CTAURE
         CALL PR3V  (KKA,JJA,IIA,AUFWFG ,'<UW-F> G        ', 1 , 0 , 1 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1,153) .EQ. 1 .OR. ISELEP(1,153) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  146),
     $                  NRRUN,ISAMPP(  146)-ISAMPA(  146),CTAURE
         CALL PR3V  (KKA,JJA,IIA,SUFWFS ,'<UW-F> G+S      ', 1 , 0 , 1 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1,148) .EQ. 1 .OR. ISELEP(1,148) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  146),
     $                  NRRUN,ISAMPP(  146)-ISAMPA(  146),CTAURE
         CALL PR3V  (KKA,JJA,IIA,SUFWFG ,'<UW-F> G+S+M    ', 1 , 0 , 1 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
C                                 HIER: V - W
C                                 -----------
C
      IF(ISELEP(1,156) .EQ. 1 .OR. ISELEP(1,156) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  156),
     $                  NRRUN,ISAMPP(  156)-ISAMPA(  156),CTAURE
         CALL PR3V  (KKA,JJA,IIA,AVFWFG ,'<VW-F> G        ', 1 , 1 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1,163) .EQ. 1 .OR. ISELEP(1,163) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  156),
     $                  NRRUN,ISAMPP(  156)-ISAMPA(  156),CTAURE
         CALL PR3V  (KKA,JJA,IIA,SVFWFS ,'<VW-F> G+S      ', 1 , 1 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1,158) .EQ. 1 .OR. ISELEP(1,158) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  156),
     $                  NRRUN,ISAMPP(  156)-ISAMPA(  156),CTAURE
         CALL PR3V  (KKA,JJA,IIA,SVFWFG ,'<VW-F> G+S+M    ', 1 , 1 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
C                                 HIER: U - V
C                                 -----------
C
      IF(ISELEP(1,166) .EQ. 1 .OR. ISELEP(1,166) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  166),
     $                  NRRUN,ISAMPP(  166)-ISAMPA(  166),CTAURE
         CALL PR3V  (KKA,JJA,IIA,AUFVFG ,'<UV-F> G        ', 0 , 1 , 1 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1,173) .EQ. 1 .OR. ISELEP(1,173) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  166),
     $                  NRRUN,ISAMPP(  166)-ISAMPA(  166),CTAURE
         CALL PR3V  (KKA,JJA,IIA,SUFVFS ,'<UV-F> G+S      ', 0 , 1 , 1 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1,168) .EQ. 1 .OR. ISELEP(1,168) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  166),
     $                  NRRUN,ISAMPP(  166)-ISAMPA(  166),CTAURE
         CALL PR3V  (KKA,JJA,IIA,SUFVFG ,'<UV-F> G+S+M    ', 0 , 1 , 1 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
C                                 HIER: AUSGABE DER SCHIEFE (SKEWNESS)
C                                 DER GESCHWINDIGKEITSFLUKTUATIONEN.
C                                 (NUR GROBSTRUKTURANTEIL)
C
      WRITE       (KANAL,6190)
C
      IF(ISELEP(1, 10) .EQ. 1 .OR. ISELEP(1, 10) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(   10),
     $                  NRRUN,ISAMPP(   10)-ISAMPA(   10),CKEINE
         CALL PR3V  (KKA,JJA,IIA,AUSKG  ,'<U-SKE> G       ', 0 , 0 , 1 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1, 35) .EQ. 1 .OR. ISELEP(1, 35) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(   35),
     $                  NRRUN,ISAMPP(   35)-ISAMPA(   35),CKEINE
         CALL PR3V  (KKA,JJA,IIA,AVSKG  ,'<V-SKE> G       ', 0 , 1 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1, 60) .EQ. 1 .OR. ISELEP(1, 60) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(   60),
     $                  NRRUN,ISAMPP(   60)-ISAMPA(   60),CKEINE
         CALL PR3V  (KKA,JJA,IIA,AWSKG  ,'<W-SKE> G       ', 1 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1, 85) .EQ. 1 .OR. ISELEP(1, 85) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(   85),
     $                  NRRUN,ISAMPP(   85)-ISAMPA(   85),CKEINE
         CALL PR3V  (KKA,JJA,IIA,APSKG  ,'<P-SKE> G       ', 0 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
C                                 HIER: AUSGABE DES FLACHHEITSGRADES
C                                 DER GESCHWINDIGKEITSFLUKTUATIONEN.
C                                 (NUR GROBSTRUKTURANTEIL)
C
      WRITE       (KANAL,6210)
C
      IF(ISELEP(1, 12) .EQ. 1 .OR. ISELEP(1, 12) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(   12),
     $                  NRRUN,ISAMPP(   12)-ISAMPA(   12),CKEINE
         CALL PR3V  (KKA,JJA,IIA,AUFLG  ,'<U-FLA> G       ', 0 , 0 , 1 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1, 37) .EQ. 1 .OR. ISELEP(1, 37) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(   37),
     $                  NRRUN,ISAMPP(   37)-ISAMPA(   37),CKEINE
         CALL PR3V  (KKA,JJA,IIA,AVFLG  ,'<V-FLA> G       ', 0 , 1 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1, 62) .EQ. 1 .OR. ISELEP(1, 62) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(   62),
     $                  NRRUN,ISAMPP(   62)-ISAMPA(   62),CKEINE
         CALL PR3V  (KKA,JJA,IIA,AWFLG  ,'<W-FLA> G       ', 1 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1, 87) .EQ. 1 .OR. ISELEP(1, 87) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(   87),
     $                  NRRUN,ISAMPP(   87)-ISAMPA(   87),CKEINE
         CALL PR3V  (KKA,JJA,IIA,APFLG  ,'<P-FLA> G       ', 0 , 0 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
C                                 HIER: X-KOMPONENTE DER VORTICITY
C
      WRITE       (KANAL,6230)
C
      IF(ISELEP(1,177) .EQ. 1 .OR. ISELEP(1,177) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  177),
     $                  NRRUN,ISAMPP(  177)-ISAMPA(  177),COMREF
         CALL PR3V  (KKA,JJA,IIA,AOX    ,'<OMEGA-X>       ', 1 , 1 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1,180) .EQ. 1 .OR. ISELEP(1,180) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  180),
     $                  NRRUN,ISAMPP(  180)-ISAMPA(  180),COMREF
         CALL PR3V  (KKA,JJA,IIA,AOXRG  ,'<OMEGA-X RMS> G ', 1 , 1 , 0 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
C                                 HIER: Y-KOMPONENTE DER VORTICITY
C
      WRITE       (KANAL,6250)
C
      IF(ISELEP(1,183) .EQ. 1 .OR. ISELEP(1,183) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  183),
     $                  NRRUN,ISAMPP(  183)-ISAMPA(  183),COMREF
         CALL PR3V  (KKA,JJA,IIA,AOY    ,'<OMEGA-Y>       ', 1 , 0 , 1 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1,186) .EQ. 1 .OR. ISELEP(1,186) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  186),
     $                  NRRUN,ISAMPP(  186)-ISAMPA(  186),COMREF
         CALL PR3V  (KKA,JJA,IIA,AOYRG  ,'<OMEGA-Y RMS> G ', 1 , 0 , 1 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
C                                 HIER: Z-KOMPONENTE DER VORTICITY
C
      WRITE       (KANAL,6270)
C
      IF(ISELEP(1,189) .EQ. 1 .OR. ISELEP(1,189) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  189),
     $                  NRRUN,ISAMPP(  189)-ISAMPA(  189),COMREF
         CALL PR3V  (KKA,JJA,IIA,AOZ    ,'<OMEGA-Z>       ', 0 , 1 , 1 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1,192) .EQ. 1 .OR. ISELEP(1,192) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  192),
     $                  NRRUN,ISAMPP(  192)-ISAMPA(  192),COMREF
         CALL PR3V  (KKA,JJA,IIA,AOZRG  ,'<OMEGA-Z RMS> G ', 0 , 1 , 1 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
C                                 HIER: ENSTROPHIE
C
      WRITE       (KANAL,6290)
C
      IF(ISELEP(1,195) .EQ. 1 .OR. ISELEP(1,195) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  195),
     $                  NRRUN,ISAMPP(  195)-ISAMPA(  195),CO2REF
         CALL PR3V  (KKA,JJA,IIA,AO2    ,'<ENSTROPHY>     ', 1 , 1 , 1 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1,198) .EQ. 1 .OR. ISELEP(1,198) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  198),
     $                  NRRUN,ISAMPP(  198)-ISAMPA(  198),CO2REF
         CALL PR3V  (KKA,JJA,IIA,AO2RG  ,'<ENSTROPHY-RMS>G', 1 , 1 , 1 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
C                                 HIER: HELIZITAET
C
      WRITE       (KANAL,6310)
C
      IF(ISELEP(1,200) .EQ. 1 .OR. ISELEP(1,200) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  200),
     $                  NRRUN,ISAMPP(  200)-ISAMPA(  200),CHEREF
         CALL PR3V  (KKA,JJA,IIA,AHE    ,'<HELICITY>      ', 1 , 1 , 1 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
      IF(ISELEP(1,203) .EQ. 1 .OR. ISELEP(1,203) .EQ. 2) THEN
         WRITE (KANAL,6040) ISAMPP(  203),
     $                  NRRUN,ISAMPP(  203)-ISAMPA(  203),CHEREF
         CALL PR3V  (KKA,JJA,IIA,AHERG  ,'<HELICITY-RMS> G', 1 , 1 , 1 ,
     $        EBE,LH1,LHS,LH2,MH1,MHS,MH2,NH1,NHS,NH2,DX,DY,DZ,X,Y,Z,
     $        KANAL)
      ENDIF
C
C
      RETURN
C
 6010 FORMAT (/,40X,47(1H*),/,40X,'S T A T I S T I S C H E     A',
     $        ' U S W E R T U N G',/,40X,47(1H*),/,1X,'DIE STICH',
     $        'PROBEN FUER DIE FOLGENDEN ENSEMBLE-MITTELWERTE',
     $        ' WURDEN NACH JEWEILS ',I6,' ZEITSCHRITTEN ENTNOMMEN.',
     $        /,1X,'DIES ENTSPRICHT EINEM ZEITINTERVALL  T = ',1PE11.4,
     $        '  (T/T-REF = ',1PE11.4,')')
 6020 FORMAT (1X,'DIE BEHANDELTE KONFIGURATION IST IN X-RICHTUNG',
     $        ' HOMOGEN')
 6030 FORMAT (1X,'DIE BEHANDELTE KONFIGURATION IST IN Y-RICHTUNG',
     $        ' HOMOGEN')
 6040 FORMAT (/,1X,'GESAMTANZAHL DER STICHPROBEN FUER DIE FOLGENDEN',
     $        ' ENSEMBLE-MITTELWERTE      : ',I8,/,1X,'DAVON WURDEN',
     $        ' WAEHREND DES VORLIEGENDEN LAUFES (NRRUN = ',I6,
     $        ') ERZEUGT',4X,': ',I8,
     $        /,1X,'ALS  B E Z U G S G R O E S S E  FUER DIE',
     $        ' FOLGENDE VARIABLE WURDE GEWAEHLT : ',A32)
 6110 FORMAT (/,40X,'STATISTISCHE MITTELWERTE  <PHIMEAN>/PHIREF',
     $        /,40X,42(1H-),/)
 6130 FORMAT (/,40X,'ROOT-MEAN-SQUARE WERTE DER SCHWANKUNGSGROESSEN',
     $        /,40X,'BEZOGEN AUF PHIRMS-REF',/,40X,46(1H-),/)
 6150 FORMAT (/,40X,'KINET. ENERGIE DER SCHWANKUNGSGESCHWINDIGKEITEN',
     $        /,40X,'BEZOGEN AUF (U-REF)**2',/,40X,47(1H-),/)
 6170 FORMAT (/,40X,'KOMPONENTEN DES SCHUBSPANNUNGSTENSORS',
     $        /,40X,'BEZOGEN AUF TAU-REF',/,40X,37(1H-),/)
 6190 FORMAT (/,40X,'SCHIEFE DER GESCHWINDIGKEITSFLUKTUATIONEN',
     $        /,40X,41(1H-),/)
 6210 FORMAT (/,40X,'FLACHHEITSGRAD DER GESCHWINDIGKEITSFLUKTUATIONEN',
     $        /,40X,48(1H-),/)
 6230 FORMAT (/,40X,'X-KOMPONENTE DER VORTICITY = <DW/DY - DV/DZ>',
     $        /,40X,'BEZOGEN AUF (U-REF)/(L-REF)',/,40X,44(1H-),/)
 6250 FORMAT (/,40X,'Y-KOMPONENTE DER VORTICITY = <DU/DZ - DW/DX>',
     $        /,40X,'BEZOGEN AUF (U-REF)/(L-REF)',/,40X,44(1H-),/)
 6270 FORMAT (/,40X,'Z-KOMPONENTE DER VORTICITY = <DV/DX - DU/DY>',
     $        /,40X,'BEZOGEN AUF (U-REF)/(L-REF)',/,40X,44(1H-),/)
 6290 FORMAT (/,40X,'ENSTROPHIE = 1/2*(OMEGA(I)*OMEGA(I))',
     $        /,40X,'BEZOGEN AUF (U-REF)**2/(L-REF)**2',/,40X,36(1H-),/)
 6310 FORMAT (/,40X,'HELIZITAET = 1/2*(U(I)*OMEGA(I))',
     $        /,40X,'BEZOGEN AUF (U-REF)**2/(L-REF)**2',/,40X,33(1H-),/)
      END
