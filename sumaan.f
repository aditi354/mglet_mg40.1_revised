










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
      SUBROUTINE SUMAAN  (

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
C        S U M A A N      IN SUMAAN WERDEN (NACH BILDUNG DER VER-
C                         BESSERTEN ENSEMBLE-MITTELWERTE IN SUBR.
C                         STMIMP) DIE ANTEILE BESTIMMTER GROESSEN AUS
C                         GROB- UND FEINSTRUKTUR ADDIERT.
C                         DIE RESULTIERENDEN WERTE SIND IN JEDEM FALL
C                         ENSEMBLE-MITTELWERTE.
C  A C H T U N G:         DIESE ROUTINE MUSS ERWEITERT WERDEN, FALLS
C                         DIE STATIST. AUSWERTUNG ERWEITERT WIRD !
C*STARLET***************************************************************
C
C PARAM:
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
C UPROG                 : DPHI0, SRMSGS, SUM2AR
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        01.09.86 (HW)  : ORIGINAL
C        12.02.03 (TB)  : SCALAR FIELD STATISTICS IMPLEMENTED
C
C*STARLET***************************************************************
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
C                                 ROOT-MEAN-SQUARE - WERTE
C                                 ------------------------
C
      IF(ISELEP(1,  9) .EQ. 1 .OR. ISELEP(1,  9) .EQ. 2) THEN
C
C                                 RMS-GESAMTWERT DER U-FLUKTUATIONEN
C                                 ZUM GROBSTRUKTURANTEIL WIRD DER
C                                 NAEHERUNGSWEISE BESTIMMTE FEIN-
C                                 STRUKTURANTEIL ADDIERT.
C
         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SURG   )
         CALL SRMSGS  (KKA,JJA,IIA,KMXA,JMXA,IMXA,SURG   ,
     $                 AURG   ,AEFS   ,'U')
      ENDIF
      IF(ISELEP(1, 34) .EQ. 1 .OR. ISELEP(1, 34) .EQ. 2) THEN
         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVRG   )
         CALL SRMSGS  (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVRG   ,
     $                 AVRG   ,AEFS   ,'V')
      ENDIF
      IF(ISELEP(1, 59) .EQ. 1 .OR. ISELEP(1, 59) .EQ. 2) THEN
         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SWRG   )
         CALL SRMSGS  (KKA,JJA,IIA,KMXA,JMXA,IMXA,SWRG   ,
     $                 AWRG   ,AEFS   ,'W')
      ENDIF
C
C                                 GESAMTE KINETISCHE ENERGIE DER
C                                 SCHWANKUNGSGESCHWINDIGKEITEN
C                                 ------------------------------
C
      IF(ISELEP(1,108) .EQ. 1 .OR. ISELEP(1,108) .EQ. 2) THEN
         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SEFG   )
         CALL SUM2AR  (KKA,JJA,IIA,KMXA,JMXA,IMXA,SEFG   ,
     $                 AEFG   ,AEFS   )
      ENDIF
C
C                                 SCHUBSPANNUNGEN
C                                 ---------------
C
C                                 HIER: U - W
C                                 -----------
C
C                                 GROBSTRUKTUR + FEINSTRUKTUR +
C                                 MOLEKULARER ANTEIL
C
      IF(ISELEP(1,148) .EQ. 1 .OR. ISELEP(1,148) .EQ. 2) THEN
         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFWFG )
         CALL SUM2AR  (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFWFG ,
     $                 AUFWFG ,AUFWFS )
         CALL SUM2AR  (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFWFG ,
     $                 SUFWFG ,AUFWFM )
      ENDIF
      IF(ISELEP(1,153) .EQ. 1 .OR. ISELEP(1,153) .EQ. 2) THEN
C
C                                 GROBSTRUKTUR + FEINSTRUKTURANTEIL
C                                 ERGIBT DEN REYNOLDSSTRESS
C
         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFWFS )
         CALL SUM2AR  (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFWFS ,
     $                 AUFWFG ,AUFWFS )
      ENDIF
C
C                                 HIER: V - W
C                                 -----------
C
      IF(ISELEP(1,158) .EQ. 1 .OR. ISELEP(1,158) .EQ. 2) THEN
         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVFWFG )
         CALL SUM2AR  (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVFWFG ,
     $                 AVFWFG ,AVFWFS )
         CALL SUM2AR  (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVFWFG ,
     $                 SVFWFG ,AVFWFM )
      ENDIF
      IF(ISELEP(1,163) .EQ. 1 .OR. ISELEP(1,163) .EQ. 2) THEN
         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVFWFS )
         CALL SUM2AR  (KKA,JJA,IIA,KMXA,JMXA,IMXA,SVFWFS ,
     $                 AVFWFG ,AVFWFS )
      ENDIF
C
C                                 HIER: U - V
C                                 -----------
C
      IF(ISELEP(1,168) .EQ. 1 .OR. ISELEP(1,168) .EQ. 2) THEN
         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFVFG )
         CALL SUM2AR  (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFVFG ,
     $                 AUFVFG ,AUFVFS )
         CALL SUM2AR  (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFVFG ,
     $                 SUFVFG ,AUFVFM )
      ENDIF
      IF(ISELEP(1,173) .EQ. 1 .OR. ISELEP(1,173) .EQ. 2) THEN
         CALL DPHI0   (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFVFS )
         CALL SUM2AR  (KKA,JJA,IIA,KMXA,JMXA,IMXA,SUFVFS ,
     $                 AUFVFG ,AUFVFS )
      ENDIF
C

      RETURN
      END
