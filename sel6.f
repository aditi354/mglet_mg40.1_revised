










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
      SUBROUTINE SEL6    (ISELEA,ISELEP,ISAMPA,ISAMPP,MAXGRIDS)
C*STARLET***************************************************************
C        S E L E C T      IN DIESER ROUTINE TRIFFT DER PROGRAMMBENUTZER
C                         DIE ENTSCHEIDUNG UEBER ART UND UMFANG DER
C                         STATISTISCHEN AUSWERTUNG BEI EINEM LARGE-
C                         EDDY-SIMULATIONSLAUF.
C                         ER MUSS DAZU DIE GEWUENSCHTE INTEGER-KENNZAHL
C                         (WIRD UNTEN ERLAEUTERT) AUF DER RECHTEN SEITE
C                         DES GLEICHHEITSZEICHENS EINTRAGEN. DIE LINKE
C                         SEITE (SPEICHERPLATZ IM ISELEP-FELD) DARF
C                         UNTER KEINEN UMSTAENDEN GEAENDERT WERDEN!!
C
C                         VOM BENUTZER MUESSEN NUR DIE TATSAECHLICH
C                         GEWUENSCHTEN GROESSEN BESTIMMT WERDEN; UEBER
C                         INTERNE ABHAENGIGKEITEN DER FELDER UNTEREIN-
C                         ANDER BRAUCHT ER SICH KEINE GEDANKEN ZU MACHEN
C                         (DIES GESCHIEHT IN SUBR. SELECA)
C
C                         FOLGENDE INTEGER-KENNZAHLEN SIND MOEGLICH:
C                         0 : KEINE AUSWERTUNG DER BETREFFENDEN GROESSE
C                         1 : ENSEMBLE-MITTELWERTE DER BETREFFENDEN
C                             GROESSE WERDEN DURCH SUMMATION VON (IN
C                             BESTIMMTEN ZEITLICHEN ABSTAENDEN GEBIL-
C                             DETEN) MOMENTANWERTEN ERMITTELT.
C                         2 : WIE 1, JEDOCH WIRD DER MOMENTANWERT FUER
C                             GRAPHISCHE ZWECKE HERAUSGESCHRIEBEN.
C                             N I C H T  MOEGLICH FUER KORRELATIONS-
C                             FUNKTIONEN, -KOEFFIZIENTEN, LEISTUNGS-
C                             DICHTESPEKTREN UND HAEUFIGKEITSVERTEILUN-
C                             GEN
C                         3 : ES WIRD KEINE ENSEMBLE-MITTELUNG DURCHGE-
C                             FUEHRT; ES ERFOLGT JEDOCH EINE AUSGABE
C                             F. GRAPHISCHE ZWECKE.
C                             N I C H T  MOEGLICH FUER KORRELATIONS-
C                             FUNKTIONEN, -KOEFFIZIENTEN, LEISTUNGS-
C                             DICHTESPEKTREN UND HAEUFIGKEITSVERTEILUN-
C                             GEN
C                         ANMERKUNG: DIE AUSGABE F. GRAPHIK IST NICHT
C                         BEI ALLEN GROESSEN MOEGLICH. DIE ROUTINE
C                         SELECA KORRIGIERT IN DIESEM FALL DIE EINGABE.
C  A C H T U N G:         DIESE ROUTINE MUSS ERWEITERT WERDEN, FALLS
C                         DIE STATIST. AUSWERTUNG ERWEITERT WIRD !
C                         DIE HIERARCHIE DER GROESSEN UNTEREINANDER
C                         SOLLTE OPTISCH DURCH ENTSPRECHENDE PLA-
C                         ZIERUNG DER ZUSAETZLICHEN PROGRAMMSTATEMENTS
C                         HERAUSGESTELLT WERDEN.
C
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KKA,JJA,IIA    - ARRAYDIMENSIONEN DER AUSWERTEFELDER
C        ISELEA(2,752)  + STEUERFELD FUER DIE AUSWERTUNG
C                         (ENTHAELT STEUERDATEN DES VORANGEGANGENEN
C                         LAUFES)
C        ISELEP(2,752)  + STEUERFELD FUER DIE AUSWERTUNG
C                         (ENTHAELT STEUERDATEN F. DEN MOMENTANEN LAUF)
C      .......          - ALLE UEBRIGEN VARIABLEN SIND ARRAYDIMENSIONEN
C
C UPROG                 : SELECA, SELECI
C
C DEFINE-DIREKTIVEN     : INFO
C
C VERS:  06.08.86 (HW)  : ORIGINAL
C        05.01.89 (HW)  : R03SEN AUS D02SEL FUER "RIPPE" ABGELEITET
C        02.06.89 (HW)  : R03SE1 NEUESTER STAND
C        26.06.89 (HW)  : D04SEL NEUESTER STAND
C        16.10.89 (HW)  : D05SEL AUS D04SEL ABGELEITET
C        17.10.89 (HW)  : R05SEL AUS D05SEL ABGELEITET
C        30.10.90 (HW)  : D06SEL AUS R05SEL ABGELEITET
C        27.10.93 (MM)  : NAME IST JETZT SEL1, HIER WERDEN DIE
C                         FELDER <U>, <V>, <W> UND <P> 
C                         UFG, VFG, WFG, PFG
C                         <U-RMS> G, <V.RMS> G, <W-RMS> G, <P-RMS> G
C                    <U-RMS> G+S, <V-RMS> G+S, <W-RMS> G+S,
C                    <UV-F> G,      <UW-F> G,     <VW-F> G,    <E-F> G
C                    <UV-F> S,      <UW-F> S,     <VW-F> S,    <E-F> S
C                    <UV-F> G+S,    <UW-F> G+S,   <VW-F> G+S,  <E-F> G+S
C                    <UV-F> M,      <UW-F> M,     <VW-F> M,
C                    <UV-F> G+S+M,  <UW-F> G+S+M, <VW-F> G+S+M,
C                    <U-SKE> G, <V-SKE> G, <W-SKE> G, <P-SKE> G
C                    <U-FLA> G, <V-FLA> G, <W-FLA> G, <P-FLA> G
C
C                         SELEKTIERT
C
C*STARLET***************************************************************
C
      INTEGER  ISELEA(2,752),  ISELEP(2,752),
     $         ISAMPA (752,MAXGRIDS),ISAMPP (752,MAXGRIDS)
C
      DO J=1,MAXGRIDS
      DO I=1,752
         ISAMPA(I,J) = 0
         ISAMPP(I,J) = 0
      ENDDO
      ENDDO
C                                  VORBELEGUNG (DARF NICHT GEAENDERT
C                                  WERDEN !)
C
      DO 100 I = 1,2
         DO 110 J = 1,752
            ISELEA(I,J) = 0
  110       ISELEP(I,J) = 0
  100 CONTINUE
C
C                                  STATIST. GROESSEN DER U-KOMPONENTE
C                                  ----------------------------------
C
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN Z-RI. DER U-FLUKT.
C
                  ISELEP(1, 24)      =  0
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN Y-RI. DER U-FLUKT.
C
                  ISELEP(1, 22)      =  0
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN X-RI. DER U-FLUKT.
C
                  ISELEP(1, 20)      =  0
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN Z-RI. DER U-FLUKT.
C
                  ISELEP(1, 18)      =  0
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN Y-RI. DER U-FLUKT.
C
                  ISELEP(1, 16)      =  0
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN X-RI. DER U-FLUKT.
C
                  ISELEP(1, 14)      =  0
C                                  ENSEMBLE-MITTELWERT DES FLACHHEITS-
C                                  GRADES DER U-FLUKTUATIONEN.
C
                  ISELEP(1, 12)      =  1
C
C                                  ENSEMBLE-MITTELWERT DER SCHIEFE
C                                  DER U-FLUKTUATIONEN.
                  ISELEP(1, 10)      =  1
C
C                                  ENSEMBLE-MITTELWERT DER ROOT-MEAN-
C                                  SQUARE - WERTE DER U-FLUKTUATIONEN
C                                  (ANTEIL AUS GROB- UND FEINSTRUKTUR)
C
               ISELEP(1,  9)         =  1
C
C                                  ENSEMBLE-MITTELWERT DER RMS-WERTE
C                                  DER U-FLUKT. (AUSSCHL. GROBSTRUKTUR)
C
            ISELEP(1,  7)            =  1
C
C                                  MOMENTANWERTE DER U-FLUKTUATIONEN
C                                  U-FLUK = U - <U>
         ISELEP(1,  6)               =  1
C
C                                  ENSEMBLE-MITTELWERT DER U-KOMPONENTE
C
      ISELEP(1,  4)                  =  1
C
C                                  STATIST. GROESSEN DER V-KOMPONENTE
C                                  ----------------------------------
C
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN Z-RI. DER V-FLUKT.
C
                  ISELEP(1, 49)      =  0
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN Y-RI. DER V-FLUKT.
C
                  ISELEP(1, 47)      =  0
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN X-RI. DER V-FLUKT.
C
                  ISELEP(1, 45)      =  0
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN Z-RI. DER V-FLUKT.
C
                  ISELEP(1, 43)      =  0
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN Y-RI. DER V-FLUKT.
C
                  ISELEP(1, 41)      =  0
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN X-RI. DER V-FLUKT.
C
                  ISELEP(1, 39)      =  0
C                                  ENSEMBLE-MITTELWERT DES FLACHHEITS-
C                                  GRADES DER V-FLUKTUATIONEN.
C
                  ISELEP(1, 37)      =  1
C
C                                  ENSEMBLE-MITTELWERT DER SCHIEFE
C                                  DER V-FLUKTUATIONEN.
C
                  ISELEP(1, 35)      =  1
C
C                                  ENSEMBLE-MITTELWERT DER ROOT-MEAN-
C                                  SQUARE - WERTE DER V-FLUKTUATIONEN
C                                  (ANTEIL AUS GROB- UND FEINSTRUKTUR)
C
               ISELEP(1, 34)         =  1
C
C                                  ENSEMBLE-MITTELWERT DER RMS-WERTE
C                                  DER V-FLUKT. (AUSSCHL. GROBSTRUKTUR)
C
            ISELEP(1, 32)            =  1
C
C                                  MOMENTANWERTE DER V-FLUKTUATIONEN
C                                  V-FLUK = V - <V>
         ISELEP(1, 31)               =  1
C
C                                  ENSEMBLE-MITTELWERT DER V-KOMPONENTE
C
      ISELEP(1, 29)                  =  1
C
C                                  STATIST. GROESSEN DER W-KOMPONENTE
C                                  ----------------------------------
C
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN Z-RI. DER W-FLUKT.
C
                  ISELEP(1, 74)      =  0
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN Y-RI. DER W-FLUKT.
C
                  ISELEP(1, 72)      =  0
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN X-RI. DER W-FLUKT.
C
                  ISELEP(1, 70)      =  0
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN Z-RI. DER W-FLUKT.
C
                  ISELEP(1, 68)      =  0
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN Y-RI. DER W-FLUKT.
C
                  ISELEP(1, 66)      =  0
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN X-RI. DER W-FLUKT.
C
                  ISELEP(1, 64)      =  0
C                                  ENSEMBLE-MITTELWERT DES FLACHHEITS-
C                                  GRADES DER W-FLUKTUATIONEN.
C
                  ISELEP(1, 62)      =  1
C
C                                  ENSEMBLE-MITTELWERT DER SCHIEFE
C                                  DER W-FLUKTUATIONEN.
C
                  ISELEP(1, 60)      =  1
C
C                                  ENSEMBLE-MITTELWERT DER ROOT-MEAN-
C                                  SQUARE - WERTE DER W-FLUKTUATIONEN
C                                  (ANTEIL AUS GROB- UND FEINSTRUKTUR)
C
               ISELEP(1, 59)         =  1
C
C                                  ENSEMBLE-MITTELWERT DER RMS-WERTE
C                                  DER W-FLUKT. (AUSSCHL. GROBSTRUKTUR)
C
            ISELEP(1, 57)            =  1
C
C                                  MOMENTANWERTE DER W-FLUKTUATIONEN
C                                  W-FLUK = W - <W>
         ISELEP(1, 56)               =  1
C
C                                  ENSEMBLE-MITTELWERT DER W-KOMPONENTE
C
      ISELEP(1, 54)                  =  1
C
C                                  STATIST. GROESSEN DER P-KOMPONENTE
C                                  ----------------------------------
C
C                                  ENSEMBLE-MITTELWERT DES FLACHHEITS-
C                                  GRADES DER P-FLUKTUATIONEN.
C
                  ISELEP(1, 87)      =  1
C
C                                  ENSEMBLE-MITTELWERT DER SCHIEFE
C                                  DER P-FLUKTUATIONEN.
C
                  ISELEP(1, 85)      =  1
C
C                                  ENSEMBLE-MITTELWERT DER ROOT-MEAN-
C                                  SQUARE - WERTE DER P-FLUKTUATIONEN
C                                  (ANTEIL AUS GROB- UND FEINSTRUKTUR)
C
*              ISELEP(1, 84)         =  0 EXISTIERT BEIM DRUCK NICHT !
C
C                                  ENSEMBLE-MITTELWERT DER RMS-WERTE
C                                  DER P-FLUKT. (AUSSCHL. GROBSTRUKTUR)
C
            ISELEP(1, 82)            =  1
C
C                                  MOMENTANWERTE DER P-FLUKTUATIONEN
C                                  P-FLUK = P - <P>
         ISELEP(1, 81)               =  1
C
C                                  ENSEMBLE-MITTELWERT DER P-KOMPONENTE
C
      ISELEP(1, 79)                  =  1
C
C                                  KINETISCHE ENERGIE DER SCHWANKUNGS-
C                                  GESCHWINDIGKEITEN
C                                  -----------------------------------
C
C                                  ENSEMBLE-MITTELWERT D. KINET. ENERGIE
C                                  (GROB- UND FEINSTRUKTURANTEIL)
C
         ISELEP(1,108)               =  1
C
C                                  ENSEMBLE-MITTELWERT D. KINET. ENERGIE
C                                  (AUSSCHLIESSLICH GROBSTRUKTURANTEIL)
C
      ISELEP(1,106)                  =  1
C
C                                  ENSEMBLE-MITTELWERT D. KINET. ENERGIE
C                                  (AUSSCHLIESSLICH FEINSTRUKTURANTEIL)
C
      ISELEP(1,109)                  =  1
C
C                                  ------------------------
C                                  DISSIPATION
C                                  --------------------------------
C
        ISELEP(1,115)                =  0
C
C                                  --------------------------------
C                                  TAYLOR''SCHES MIKROMASS
C                                  --------------------------------
C
C                                  U_KOMPONENTE IN X_RICHTUNG
        ISELEP(1,115)                =  0
C
C                                  U_KOMPONENTE IN Y_RICHTUNG
        ISELEP(1,115)                =  0

C                                  U_KOMPONENTE IN Z_RICHTUNG
        ISELEP(1,115)                =  0
C
C                                  V_KOMPONENTE IN X_RICHTUNG
        ISELEP(1,115)                =  0
C
C                                  V_KOMPONENTE IN Y_RICHTUNG
        ISELEP(1,115)                =  0
C
C                                  V_KOMPONENTE IN Z_RICHTUNG
        ISELEP(1,115)                =  0
C
C                                  W_KOMPONENTE IN X_RICHTUNG
        ISELEP(1,115)                =  0
C
C                                  W_KOMPONENTE IN Y_RICHTUNG
        ISELEP(1,115)                =  0
C
C                                  W_KOMPONENTE IN Z_RICHTUNG
        ISELEP(1,115)                =  0

C
C
C                                  KOMPONENTEN DES SPANNUNGSTENSORS
C                                  --------------------------------
C
C                                  HIER: U - W
C                                  -----------
C
C                                  ENSEMBLE-MITTELWERT DER GESAMTSCHUB-
C                                  SPANNUNG (SUMME AUS GROB- UND FEIN-
C                                  STRUKTUR SOWIE DEM MOLEKULAREN
C                                  ANTEIL)
C
         ISELEP(1,148)               =  1
C
C                                  ENSEMBLE-MITTELWERT DER REYNOLDS-
C                                  SPANNUNG (SUMME AUS GROB- UND FEIN-
C                                  STRUKTURANTEIL) 29.08.86
C
         ISELEP(1,153)               =  1
C
C                                  ENSEMBLE-M. (AUSSCHL. GROBSTRUKTURA.)
C
      ISELEP(1,146)                  =  1
C
C                                  ENSEMBLE-M. (AUSSCHL. FEINSTRUKTURA.)
C
      ISELEP(1,149)                  =  1
C
C                                  ENSEMBLE-M. (AUSSCHL. MOLEKULARER A.)
C
      ISELEP(1,151)                  =  1
C
C                                  HIER: V - W
C                                  -----------
C
C                                  ENSEMBLE-MITTELWERT DER GESAMTSCHUB-
C                                  SPANNUNG (SUMME AUS GROB- UND FEIN-
C                                  STRUKTUR SOWIE DEM MOLEKULAREN
C                                  ANTEIL)
C
         ISELEP(1,158)               =  1
C
C                                  ENSEMBLE-MITTELWERT DER REYNOLDS-
C                                  SPANNUNG (SUMME AUS GROB- UND FEIN-
C                                  STRUKTURANTEIL) 29.08.86
C
         ISELEP(1,163)               =  1
C
C                                  ENSEMBLE-M. (AUSSCHL. GROBSTRUKTURA.)
C
      ISELEP(1,156)                  =  1
C
C                                  ENSEMBLE-M. (AUSSCHL. FEINSTRUKTURA.)
C
      ISELEP(1,159)                  =  1
C
C                                  ENSEMBLE-M. (AUSSCHL. MOLEKULARER A.)
C
      ISELEP(1,161)                  =  1
C
C                                  HIER: U - V
C                                  -----------
C
C                                  ENSEMBLE-MITTELWERT DER GESAMTSCHUB-
C                                  SPANNUNG (SUMME AUS GROB- UND FEIN-
C                                  STRUKTUR SOWIE DEM MOLEKULAREN
C                                  ANTEIL)
C
         ISELEP(1,168)               =  1
C
C                                  ENSEMBLE-MITTELWERT DER REYNOLDS-
C                                  SPANNUNG (SUMME AUS GROB- UND FEIN-
C                                  STRUKTURANTEIL) 29.08.86
C
         ISELEP(1,173)               =  1
C
C                                  ENSEMBLE-M. (AUSSCHL. GROBSTRUKTURA.)
C
      ISELEP(1,166)                  =  1
C
C                                  ENSEMBLE-M. (AUSSCHL. FEINSTRUKTURA.)
C
      ISELEP(1,169)                  =  1
C
C                                  ENSEMBLE-M. (AUSSCHL. MOLEKULARER A.)
C
      ISELEP(1,171)                  =  1
C
C                                  RANDBEHANDLUNG VON FRED
C                                  -----------------------
          ISELEP(1,301)            =   1
          ISELEP(1,302)            =   1
          ISELEP(1,303)            =   1
          ISELEP(1,304)            =   1
C
C                                  HIER: STATISTISCHE GROESSEN DER
C                                  VORTICITY
C                                  -------------------------------
C
C                                  ENSEMBLE-MITTELWERT DER RMS-WERTE
C                                  DER OMEGA-X-FLUKTUATIONEN (GROBSTR.)
C
            ISELEP(1,180)            =  0
C
C                                  MOMENTANWERTE DER OMEGA-X-FLUKT.
C                                  OX-FLUK = OX - <OX>
         ISELEP(1,179)               =  0
C
C                                  ENSEMBLE-MITTELWERT DER OMEGA-X-KOMP.
C
      ISELEP(1,177)                  =  0
C
C                                  ENSEMBLE-MITTELWERT DER RMS-WERTE
C                                  DER OMEGA-Y-FLUKTUATIONEN (GROBSTR.)
C
            ISELEP(1,186)            =  0
C
C                                  MOMENTANWERTE DER OMEGA-Y-FLUKT.
C                                  OY-FLUK = OY - <OY>
         ISELEP(1,185)               =  0
C
C                                  ENSEMBLE-MITTELWERT DER OMEGA-Y-KOMP.
C
      ISELEP(1,183)                  =  0
C
C                                  ENSEMBLE-MITTELWERT DER RMS-WERTE
C                                  DER OMEGA-Z-FLUKTUATIONEN (GROBSTR.)
C
            ISELEP(1,192)            =  0
C
C                                  MOMENTANWERTE DER OMEGA-Z-FLUKT.
C                                  OZ-FLUK = OZ - <OZ>
         ISELEP(1,191)               =  0
C
C                                  ENSEMBLE-MITTELWERT DER OMEGA-Z-KOMP.
C
      ISELEP(1,189)                  =  0
C
C                                  HIER: ENSTROPHIE
C                                  ----------------
C
C                                  ENSEMBLE-MITTELWERT DER RMS-WERTE
C                                  DER ENSTROPHIE-FLUKT. (GROBSTR.)
C
            ISELEP(1,198)            =  0
C
C                                  MOMENTANWERTE DER ENSTROPHIE-FLUKT.
C                                  EN-FLUK = EN - <EN>
         ISELEP(1,197)               =  0
C
C                                  ENSEMBLE-MITTELWERT DER ENSTROPHIE
C
      ISELEP(1,195)                  =  0
C
C                                  HIER: HELIZITAET
C                                  ----------------
C
C                                  ENSEMBLE-MITTELWERT DER RMS-WERTE
C                                  DER HELIZITAET-FLUKT. (GROBSTR.)
C
            ISELEP(1,203)            =  0
C
C                                  MOMENTANWERTE DER HELIZITAET-FLUKT.
C                                  EN-FLUK = EN - <EN>
         ISELEP(1,202)               =  0
C
C                                  ENSEMBLE-MITTELWERT DER HELIZITAET
C
      ISELEP(1,200)                  =  0
C
C                                  KREUZKORRELATIONSFUNKTIONEN
C                                  ---------------------------
C
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN X-RI. (U"V")
C
                  ISELEP(1,210)      =  0
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN Y-RI. (U"V")
C
                  ISELEP(1,212)      =  0
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN X-RI. (U"W")
C
                  ISELEP(1,216)      =  0
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN Y-RI. (U"W")
C
                  ISELEP(1,218)      =  0
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN Y-RI. (V"W")
C
                  ISELEP(1,224)      =  0
C                                  KORRELATIONSKOEFFIZIENTEN
C                                  -------------------------
C
C                                  ENSEMBLE-MITTELWERT DES KORRELATIONS-
C                                  KOEFF. ZWISCHEN U" UND W" IN X-RI.
C
                  ISELEP(1,234)      =  0
C                                  ENSEMBLE-MITTELWERT DES KORRELATIONS-
C                                  KOEFF. ZWISCHEN U" UND W" IN Z-RI.
C
                  ISELEP(1,238)      =  0
C
C                                  OBERE DREIECKSMATRIX DES TENSORS
C                                  DER KORRELATIONSFUNKTIONEN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  --------------------------------
C
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN X-RI. (OX"OX")
C
                  ISELEP(1,246)      =  0
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN X-RI. (OX"OY")
C
                  ISELEP(1,248)      =  0
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN X-RI. (OX"OZ")
C
                  ISELEP(1,250)      =  0
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN X-RI. (OY"OY")
C
                  ISELEP(1,252)      =  0
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN X-RI. (OY"OZ")
C
                  ISELEP(1,254)      =  0
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN X-RI. (OZ"OZ")
C
                  ISELEP(1,256)      =  0
C
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN Y-RI. (OX"OX")
C
                  ISELEP(1,258)      =  0
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN Y-RI. (OX"OY")
C
                  ISELEP(1,260)      =  0
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN Y-RI. (OX"OZ")
C
                  ISELEP(1,262)      =  0
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN Y-RI. (OY"OY")
C
                  ISELEP(1,264)      =  0
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN Y-RI. (OY"OZ")
C
                  ISELEP(1,266)      =  0
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN Y-RI. (OZ"OZ")
C
                  ISELEP(1,268)      =  0
C
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN Z-RI. (OX"OX")
C
                  ISELEP(1,270)      =  0
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN Z-RI. (OX"OY")
C
                  ISELEP(1,272)      =  0
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN Z-RI. (OX"OZ")
C
                  ISELEP(1,274)      =  0
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN Z-RI. (OY"OY")
C
                  ISELEP(1,276)      =  0
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN Z-RI. (OY"OZ")
C
                  ISELEP(1,278)      =  0
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN Z-RI. (OZ"OZ")
C
                  ISELEP(1,280)      =  0
C
C                                  LEISTUNGSDICHTESPEKTREN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  ----------------------------
C
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN X-RI. (OX"OX")
C
                  ISELEP(1,282)      =  0
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN X-RI. (OY"OY")
C
                  ISELEP(1,284)      =  0
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN X-RI. (OZ"OZ")
C
                  ISELEP(1,286)      =  0
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN Y-RI. (OX"OX")
C
                  ISELEP(1,288)      =  0
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN Y-RI. (OY"OY")
C
                  ISELEP(1,290)      =  0
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN Y-RI. (OZ"OZ")
C
                  ISELEP(1,292)      =  0
C
C                                  HAEUFIGKEITSVERTEILUNG DER
C                                  INKLINATIONSWINKEL DER VORTICITY-
C                                  VEKTOREN
C                                  ---------------------------------
C
C                                  ENSEMBLE-MITTELWERT DER HAEUFIG-
C                                  KEITSVERTEILUNG VON  ATAN(OZ/OY)
C
                  ISELEP(1,306)      =  0
C                                  ENSEMBLE-MITTELWERT DER HAEUFIG-
C                                  KEITSVERTEILUNG VON  ATAN(OZ/OX)
C
                  ISELEP(1,308)      =  0
C                                  ENSEMBLE-MITTELWERT DER HAEUFIG-
C                                  KEITSVERTEILUNG VON  ATAN(OY/OX)
C
                  ISELEP(1,310)      =  0
C
C                                  JETZT UEBERPRUEFUNG AUF KONSISTENZ
C                                  UND VOLLSTAENDIGKEIT
C
      CALL SELECA (ISELEP)
C
C
      RETURN
      END
