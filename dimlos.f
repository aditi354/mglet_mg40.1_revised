










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
      SUBROUTINE DIMLOS  (KPP,JPP,IPP,GMOL,RHO,MTURB,
     $                    X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,

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
C        D I M L O S      IN DIMLOS WERDEN DIE DIMENSIONSBEHAFTETEN
C                         GROESSEN IN DIMENSIONSLOSE UMGEWANDELT.
C                         DIE FESTLEGUNG DER BASIS-BEZUGSGROESSEN
C                         ERFOLGT IN SUBR. SETREF.
C  A C H T U N G:         DIESE ROUTINE MUSS ERWEITERT WERDEN, FALLS
C                         DIE STATIST. AUSWERTUNG ERWEITERT WIRD !
C*STARLET***************************************************************
C
C PARAM: KPP,JPP,IPP    - ORT DES BEZUGSDRUCKES
C        GMOL           - MOLEKULARE DYNAMISCHE VISKOSITAET
C        RHO            - DICHTE DES FLUIDS (= CONST.)
C        MTURB          - SCHALTER: 0 = LAMINAR,  1 = TURBULENT
C        X, Y, Z        + KOORDINATEN DER ZELLDEFINITIONSPUNKTE
C        DX,DY,DZ       + ABSTAND DER GITTERPUNKTE
C        DDX,DDY,DDZ    + KANTENLAENGEN DER KONTROLLVOLUMINA
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
C UPROG                 : PLEVEL, AMULT
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        01.09.86 (HW)  : ORIGINAL
C        30.09.86 (HW)  : BEZUGSGROESSEN WERDEN IM RIDENT-FELD UEBER-
C                         GEBEN. UEBERGABE  G E A E N D E R T !!
C        01.07.88 (HW)  : ERWEITERUNG: AUSWERTUNG DER DREI KOMPONENTEN
C                         DER VORTICITY, DER ENSTROPHIE UND DER
C                         HELIZITAET
C        19.08.88 (HW)  : LAENGENABMESSUNGEN WERDEN EBENFALLS
C                         DIMENSIONSLOS GEMACHT.  A C H T U N G:
C                         UEBERGABE GEAENDERT !!
C        15.12.88 (HW)  : BESTIMMUNG DER REFERENZGROESSEN ALS COMMON-
C                         DECK EINGEFUEHRT (CREFVA).
C        30.10.90 (HW)  : UEBERGABE AN PLEVEL GEAENDERT
C        27.02.92 (MM)  : PRUEFUNG DER UNDEFINIERTEN WERTE MIT 
C                         UNDEF AUSGESCHALTET
C
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
C
      REAL        X(II),         Y(JJ),         Z(KK),
     $           DX(II),        DY(JJ),        DZ(KK),
     $          DDX(II),       DDY(JJ),       DDZ(KK)
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
C                                 REZIPROKE BEZUGSGROESSEN
C
      RUREF  = 1.0 / (AMAX1((10.0*SMALL),ABS(RIDENT(75))))
     $       *       SIGN(1.0,RIDENT(75))
      RLREF  = 1.0 / (AMAX1((10.0*SMALL),ABS(RIDENT(76))))
     $       *       SIGN(1.0,RIDENT(76))
      ROMREF =       (AMAX1((10.0*SMALL),ABS(RIDENT(77))))
     $       *       SIGN(1.0,RIDENT(77))
      REREF  = 1.0 / (AMAX1((10.0*SMALL),ABS(RIDENT(78))))
     $       *       SIGN(1.0,RIDENT(78))
      RGREF  = 1.0 / (AMAX1((10.0*SMALL),ABS(RIDENT(79))))
     $       *       SIGN(1.0,RIDENT(79))
      RPREF  = 1.0 / (AMAX1((10.0*SMALL),ABS(RIDENT(80))))
     $       *       SIGN(1.0,RIDENT(80))
      RTAURE = 1.0 / (AMAX1((10.0*SMALL),ABS(RIDENT(81))))
     $       *       SIGN(1.0,RIDENT(81))
      RO2REF = ROMREF**2
      RHEREF = RUREF**2 / RLREF
      RWAVEN = 1.0      / RLREF
      RASDUI = 2.0 * RUREF**2 * RLREF
      RASDOM = 2.0 * RUREF**2 / RLREF
C
C                                 ZUERST WERDEN DIE LAENGENABMESSUNGEN
C                                 DIMENSIONSLOS GEMACHT
C
      CALL AMULT   (II , 1 , 1 ,IMX , 1  , 1  ,X      ,RLREF  )
      CALL AMULT   (II , 1 , 1 ,IMX , 1  , 1  ,DX     ,RLREF  )
      CALL AMULT   (II , 1 , 1 ,IMX , 1  , 1  ,DDX    ,RLREF  )
C
      CALL AMULT   (JJ , 1 , 1 ,JMX , 1  , 1  ,Y      ,RLREF  )
      CALL AMULT   (JJ , 1 , 1 ,JMX , 1  , 1  ,DY     ,RLREF  )
      CALL AMULT   (JJ , 1 , 1 ,JMX , 1  , 1  ,DDY    ,RLREF  )
C
      CALL AMULT   (KK , 1 , 1 ,KMX , 1  , 1  ,Z      ,RLREF  )
      CALL AMULT   (KK , 1 , 1 ,KMX , 1  , 1  ,DZ     ,RLREF  )
      CALL AMULT   (KK , 1 , 1 ,KMX , 1  , 1  ,DDZ    ,RLREF  )
C
C                                 DIE MOMENTANEN GESCHWINDIGKEITS-,
C                                 DRUCK- UND G-FELDER WERDEN
C                                 DIMENSIONSLOS GEMACHT
C
      IF(ISELEP(1,  1) .EQ. 1 .OR. ISELEP(1,  1) .EQ. 2) THEN
C
C                                 UNDEFINIERTE WERTE IM FELD WERDEN
C                                 MIT 0.0 BELEGT WEGEN SUBR. AMULT !
C                                 SIND AB 27.02.92 IMMER MIT 0.0 BEL.
C
         CALL AMULT   (KK ,JJ ,II ,KMX ,JMX ,IMX ,U      ,RUREF  )
      ENDIF
      IF(ISELEP(1, 26) .EQ. 1 .OR. ISELEP(1, 26) .EQ. 2) THEN
         CALL AMULT   (KK ,JJ ,II ,KMX ,JMX ,IMX ,V      ,RUREF  )
      ENDIF
      IF(ISELEP(1, 51) .EQ. 1 .OR. ISELEP(1, 51) .EQ. 2) THEN
         CALL AMULT   (KK ,JJ ,II ,KMX ,JMX ,IMX ,W      ,RUREF  )
      ENDIF
      IF(ISELEP(1, 76) .EQ. 1 .OR. ISELEP(1, 76) .EQ. 2) THEN
         CALL PLEVEL  (KK ,JJ ,II ,KMX ,JMX ,IMX ,KPP,JPP,IPP,P,
     $                 DDX, DDY, DDZ )
         CALL AMULT   (KK ,JJ ,II ,KMX ,JMX ,IMX ,P      ,RPREF  )
      ENDIF
      IF(ISELEP(1,101) .EQ. 1 .OR. ISELEP(1,101) .EQ. 2) THEN
         CALL AMULT   (KK ,JJ ,II ,KMX ,JMX ,IMX ,G      ,RGREF  )
      ENDIF
C
      IF(MTURB .EQ. 0) RETURN
C
C                                 STATIST. GROESSEN DER U-KOMPONENTE
C                                 ----------------------------------
C
C                                 <U>-WERTE
C
      IF(ISELEP(1,  4) .EQ. 1 .OR. ISELEP(1,  4) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AU     ,RUREF  )
      ENDIF
C
C                                 UFG-WERTE
C
      IF(ISELEP(1,  6) .EQ. 1 .OR. ISELEP(1,  6) .EQ. 2) THEN
         CALL AMULT   (KK ,JJ ,II ,KMX ,JMX ,IMX ,UFG    ,RUREF  )
      ENDIF
C
C                                 <U-RMS>-WERTE
C
      IF(ISELEP(1,  7) .EQ. 1 .OR. ISELEP(1,  7) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AURG   ,RUREF  )
      ENDIF
C
C                                 <U-RMS> G+S -WERTE
C
      IF(ISELEP(1,  9) .EQ. 1 .OR. ISELEP(1,  9) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SURG   ,RUREF  )
      ENDIF
C
C                                 STATIST. GROESSEN DER V-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 29) .EQ. 1 .OR. ISELEP(1, 29) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AV     ,RUREF  )
      ENDIF
      IF(ISELEP(1, 31) .EQ. 1 .OR. ISELEP(1, 31) .EQ. 2) THEN
         CALL AMULT   (KK ,JJ ,II ,KMX ,JMX ,IMX ,VFG    ,RUREF  )
      ENDIF
      IF(ISELEP(1, 32) .EQ. 1 .OR. ISELEP(1, 32) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AVRG   ,RUREF  )
      ENDIF
      IF(ISELEP(1, 34) .EQ. 1 .OR. ISELEP(1, 34) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVRG   ,RUREF  )
      ENDIF
C
C                                 STATIST. GROESSEN DER W-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 54) .EQ. 1 .OR. ISELEP(1, 54) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AW     ,RUREF  )
      ENDIF
      IF(ISELEP(1, 56) .EQ. 1 .OR. ISELEP(1, 56) .EQ. 2) THEN
         CALL AMULT   (KK ,JJ ,II ,KMX ,JMX ,IMX ,WFG    ,RUREF  )
      ENDIF
      IF(ISELEP(1, 57) .EQ. 1 .OR. ISELEP(1, 57) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AWRG   ,RUREF  )
      ENDIF
      IF(ISELEP(1, 59) .EQ. 1 .OR. ISELEP(1, 59) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SWRG   ,RUREF  )
      ENDIF
C
C                                 STATIST. GROESSEN DER P-KOMPONENTE
C                                 ----------------------------------
C
      IF(ISELEP(1, 79) .EQ. 1 .OR. ISELEP(1, 79) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AP     ,RPREF  )
      ENDIF
      IF(ISELEP(1, 81) .EQ. 1 .OR. ISELEP(1, 81) .EQ. 2) THEN
         CALL AMULT   (KK ,JJ ,II ,KMX ,JMX ,IMX ,PFG    ,RPREF  )
      ENDIF
      IF(ISELEP(1, 82) .EQ. 1 .OR. ISELEP(1, 82) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,APRG   ,RPREF  )
      ENDIF
*     IF(ISELEP(1, 84) .EQ. 1 .OR. ISELEP(1, 84) .EQ. 2) THEN
*        CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,PHI    ,CMULT  )
*     ENDIF
C
C                                 KINETISCHE ENERGIE DER SCHWANKUNGS-
C                                 GESCHWINDIGKEITEN
C                                 -----------------------------------
C
      IF(ISELEP(1,106) .EQ. 1 .OR. ISELEP(1,106) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AEFG   ,REREF  )
      ENDIF
      IF(ISELEP(1,108) .EQ. 1 .OR. ISELEP(1,108) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SEFG   ,REREF  )
      ENDIF
      IF(ISELEP(1,109) .EQ. 1 .OR. ISELEP(1,109) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AEFS   ,REREF  )
      ENDIF
C
C                                 KOMPONENTEN DES SPANNUNGSTENSORS
C                                 --------------------------------
C
C                                 HIER: U - W
C                                 -----------
C
      IF(ISELEP(1,146) .EQ. 1 .OR. ISELEP(1,146) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUFWFG ,RTAURE )
      ENDIF
      IF(ISELEP(1,148) .EQ. 1 .OR. ISELEP(1,148) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFWFG ,RTAURE )
      ENDIF
      IF(ISELEP(1,149) .EQ. 1 .OR. ISELEP(1,149) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUFWFS ,RTAURE )
      ENDIF
      IF(ISELEP(1,151) .EQ. 1 .OR. ISELEP(1,151) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUFWFM ,RTAURE )
      ENDIF
      IF(ISELEP(1,153) .EQ. 1 .OR. ISELEP(1,153) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFWFS ,RTAURE )
      ENDIF
C
C                                 HIER: V - W
C                                 -----------
C
      IF(ISELEP(1,156) .EQ. 1 .OR. ISELEP(1,156) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AVFWFG ,RTAURE )
      ENDIF
      IF(ISELEP(1,158) .EQ. 1 .OR. ISELEP(1,158) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVFWFG ,RTAURE )
      ENDIF
      IF(ISELEP(1,159) .EQ. 1 .OR. ISELEP(1,159) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AVFWFS ,RTAURE )
      ENDIF
      IF(ISELEP(1,161) .EQ. 1 .OR. ISELEP(1,161) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AVFWFM ,RTAURE )
      ENDIF
      IF(ISELEP(1,163) .EQ. 1 .OR. ISELEP(1,163) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVFWFS ,RTAURE )
      ENDIF
C
C                                 HIER: U - V
C                                 -----------
C
      IF(ISELEP(1,166) .EQ. 1 .OR. ISELEP(1,166) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUFVFG ,RTAURE )
      ENDIF
      IF(ISELEP(1,168) .EQ. 1 .OR. ISELEP(1,168) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFVFG ,RTAURE )
      ENDIF
      IF(ISELEP(1,169) .EQ. 1 .OR. ISELEP(1,169) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUFVFS ,RTAURE )
      ENDIF
      IF(ISELEP(1,171) .EQ. 1 .OR. ISELEP(1,171) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AUFVFM ,RTAURE )
      ENDIF
      IF(ISELEP(1,173) .EQ. 1 .OR. ISELEP(1,173) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFVFS ,RTAURE )
      ENDIF
C
C                                 HIER: OMEGA-X
C                                 -------------
C
      IF(ISELEP(1,176) .EQ. 1 .OR. ISELEP(1,176) .EQ. 2) THEN
         CALL AMULT   (KK ,JJ ,II ,KMX ,JMX ,IMX ,OX     ,ROMREF )
      ENDIF
      IF(ISELEP(1,177) .EQ. 1 .OR. ISELEP(1,177) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AOX    ,ROMREF )
      ENDIF
      IF(ISELEP(1,179) .EQ. 1 .OR. ISELEP(1,179) .EQ. 2) THEN
         CALL AMULT   (KK ,JJ ,II ,KMX ,JMX ,IMX ,OXFG   ,ROMREF )
      ENDIF
      IF(ISELEP(1,180) .EQ. 1 .OR. ISELEP(1,180) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AOXRG  ,ROMREF )
      ENDIF
C
C                                 HIER: OMEGA-Y
C                                 -------------
C
      IF(ISELEP(1,182) .EQ. 1 .OR. ISELEP(1,182) .EQ. 2) THEN
         CALL AMULT   (KK ,JJ ,II ,KMX ,JMX ,IMX ,OY     ,ROMREF )
      ENDIF
      IF(ISELEP(1,183) .EQ. 1 .OR. ISELEP(1,183) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AOY    ,ROMREF )
      ENDIF
      IF(ISELEP(1,185) .EQ. 1 .OR. ISELEP(1,185) .EQ. 2) THEN
         CALL AMULT   (KK ,JJ ,II ,KMX ,JMX ,IMX ,OYFG   ,ROMREF )
      ENDIF
      IF(ISELEP(1,186) .EQ. 1 .OR. ISELEP(1,186) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AOYRG  ,ROMREF )
      ENDIF
C
C                                 HIER: OMEGA-Z
C                                 -------------
C
      IF(ISELEP(1,188) .EQ. 1 .OR. ISELEP(1,188) .EQ. 2) THEN
         CALL AMULT   (KK ,JJ ,II ,KMX ,JMX ,IMX ,OZ     ,ROMREF )
      ENDIF
      IF(ISELEP(1,189) .EQ. 1 .OR. ISELEP(1,189) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AOZ    ,ROMREF )
      ENDIF
      IF(ISELEP(1,191) .EQ. 1 .OR. ISELEP(1,191) .EQ. 2) THEN
         CALL AMULT   (KK ,JJ ,II ,KMX ,JMX ,IMX ,OZFG   ,ROMREF )
      ENDIF
      IF(ISELEP(1,192) .EQ. 1 .OR. ISELEP(1,192) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AOZRG  ,ROMREF )
      ENDIF
C
C                                 HIER: ENSTROPHIE
C                                 ----------------
C
      IF(ISELEP(1,195) .EQ. 1 .OR. ISELEP(1,195) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AO2    ,RO2REF )
      ENDIF
      IF(ISELEP(1,197) .EQ. 1 .OR. ISELEP(1,197) .EQ. 2) THEN
         CALL AMULT   (KK ,JJ ,II ,KMX ,JMX ,IMX ,O2FG   ,RO2REF )
      ENDIF
      IF(ISELEP(1,198) .EQ. 1 .OR. ISELEP(1,198) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AO2RG  ,RO2REF )
      ENDIF
C
C                                 HIER: HELIZITAET
C                                 ----------------
C
      IF(ISELEP(1,200) .EQ. 1 .OR. ISELEP(1,200) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AHE    ,RHEREF )
      ENDIF
      IF(ISELEP(1,202) .EQ. 1 .OR. ISELEP(1,202) .EQ. 2) THEN
         CALL AMULT   (KK ,JJ ,II ,KMX ,JMX ,IMX ,HEFG   ,RHEREF )
      ENDIF
      IF(ISELEP(1,203) .EQ. 1 .OR. ISELEP(1,203) .EQ. 2) THEN
         CALL AMULT   (KKA,JJA,IIA,KMXA,JMXA,IMXA,AHERG  ,RHEREF )
      ENDIF
C
      RETURN
      END
