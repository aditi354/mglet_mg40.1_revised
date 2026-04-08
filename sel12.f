










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
      SUBROUTINE SEL12    (ISELEP)
C*STARLET***************************************************************
C        S E L E C T      IN DIESER ROUTINE WERDEN ART UND UMFANG DER
C                         STATISTISCHEN AUSWERTUNG DEFINIERT.
C                         DAZU WERDEN DIE ENTSPRECHENDEN STELLEN DER
C                         VARIABLE ISELECT MIT NULLEN UND EINSEN BELEGT.
C
C                         FOLGENDE INTEGER-KENNZAHLEN SIND MOEGLICH:
C                         0 : KEINE AUSWERTUNG DER BETREFFENDEN GROESSE
C                         1 : ENSEMBLE-MITTELWERTE DER BETREFFENDEN
C                             GROESSE WERDEN DURCH SUMMATION VON (IN
C                             BESTIMMTEN ZEITLICHEN ABSTAENDEN GEBIL-
C                             DETEN) MOMENTANWERTEN ERMITTELT.
C
C  A C H T U N G:         DIESE ROUTINE MUSS ERWEITERT WERDEN, FALLS
C                         DIE STATIST. AUSWERTUNG ERWEITERT WIRD !
C                         DIE HIERARCHIE DER GROESSEN UNTEREINANDER
C                         SOLLTE OPTISCH DURCH ENTSPRECHENDE PLA-
C                         ZIERUNG DER ZUSAETZLICHEN PROGRAMMSTATEMENTS
C                         HERAUSGESTELLT WERDEN.
C
C
C                         Die neuen Versionen ab Maerz 2000 koennen nur
C                         zusammen mit _stat0_ aufgerufen werden. Dies
C                         ermöglicht, dass auch mehrere Levels aus-
C                         gewaehlt werden koennen.
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
C        28.03.99 (SE)  : NAME IST JETZT SEL12. FUER DIE BRECHNUNG DER
C                         BILANZTERME WIRD AUF EIN NEUES VERFAHREN ZUR
C                         BESTIMMUNG DER SCHWANKUNGSGROESSEN UMGESTELLT.
C                         ALLE NICHTMEHR FUNKTIONSFAEHIGEN OPTIONEN
C                         WERDEN ENTFERNT.
C                         UM KONFLIKTE MIT FRUEHEREN STATISTIKLEVELN ZU
C                         VERMEIDEN WERDEN FUER NEUE WERTE NUR DIE
C                         PALETZE VON ISELEP GROESSER 399 VERWENDET.
C*STARLET***************************************************************

      INTEGER  ISELEP(2,752)



C                                  ENSEMBLE-MITTELWERT DER U-KOMPONENTE
C
      ISELEP(1,  4)                  =  1
C
C
C                                  ENSEMBLE-MITTELWERT DER V-KOMPONENTE
C
      ISELEP(1, 29)                  =  1
C
C
C                                  ENSEMBLE-MITTELWERT DER W-KOMPONENTE
C
      ISELEP(1, 54)                  =  1
C
C
C                                  ENSEMBLE-MITTELWERT DER P-KOMPONENTE
C
      ISELEP(1, 79)                  =  1

c                                   Korrelationen der Momentan-
c                                   geschwindigkeiten
c                                   --------------------------------
c
c                                   <UU> G
c                  ISELEP(1,400)      =  1
c
c                                   <VV> G
c                  ISELEP(1,401)      =  1
c
c                                   <WW> G
c                  ISELEP(1,402)      =  1
c
c                                   <UV> G
                  ISELEP(1,403)      =  1
c
c                                   <UW> G
                  ISELEP(1,404)      =  1
c
c                                   <VW> G
                  ISELEP(1,405)      =  1
c
c
c                                   Trippel-Korrelationen der
c                                   Momentangeschwindigkeiten
c                                   --------------------------------
c
c                                   <UUU> G
c                  ISELEP(1,410)      =  1
c
c                                   <VVV> G
c                  ISELEP(1,411)      =  1
c
c                                   <WWW> G
c                  ISELEP(1,412)      =  1
c
c                                   <UUV> G
c                  ISELEP(1,413)      =  1
c
c                                   <UUW> G
c                  ISELEP(1,414)      =  1
c
c                                   <VVU> G
c                  ISELEP(1,415)      =  1
c
c                                   <VVW> G
c                  ISELEP(1,416)      =  1
c
c                                   <WWU> G
c                  ISELEP(1,417)      =  1
c
c                                   <WWV> G
c                  ISELEP(1,418)      =  1
c
c                                   <UVW> G
c                  ISELEP(1,419)      =  1
c
c
c                                   Druck-Momentangeschwindigkeits-
c                                   Korrelationen
c                                   -------------------------------
c
c                                  <UP> G
c        ISELEP(1,420)                =  0
c
c                                  <VP> G
c        ISELEP(1,421)                =  0
c
c                                  <WP> G
c        ISELEP(1,422)                =  0

c                                   Druck-Scher-Korrelationen
c                                   -------------------------------
c
c                                  <(dU/dX)P> G
c        ISELEP(1,423)                =  0
c
c                                  <(dV/dY)P> G
c        ISELEP(1,424)                =  0
c
c                                  <(dW/dZ)P> G
c        ISELEP(1,425)                =  0
c
c                                  <(dU/dY + dV/dX)P> G
c        ISELEP(1,426)                =  0
c
c                                  <(dU/dZ + dW/dX)P> G
c        ISELEP(1,427)                =  0
c
c                                  <(dV/dZ + dW/dY)P> G
c        ISELEP(1,428)                =  0

c
c
c                                   Elemente des Deformations-
c                                   geschwindigkeitstensors
c                                   -------------------------------
c
c                                   Quadratische Terme:
c
c                                  <(dU/dX)^2> G
c        ISELEP(1,430)                =  0
c
c                                  <(dU/dY)^2> G
c        ISELEP(1,431)                =  0
c
c                                  <(dU/dZ)^2> G
c        ISELEP(1,432)                =  0
c
c                                  <(dV/dX)^2> G
c        ISELEP(1,433)                =  0
c
c                                  <(dV/dY)^2> G
c        ISELEP(1,434)                =  0
c
c                                  <(dV/dZ)^2> G
c        ISELEP(1,435)                =  0
c
c                                  <(dW/dX)^2> G
c        ISELEP(1,436)                =  0
c
c                                  <(dW/dY)^2> G
c        ISELEP(1,437)                =  0
c
c                                  <(dW/dZ)^2> G
c        ISELEP(1,438)                =  0
c
c                                    Gemischte Terme
c
c                                  <dU/dX dV/dX > G
c        ISELEP(1,440)                =  0
c                                  <dU/dY dV/dY > G
c        ISELEP(1,441)                =  0
c                                  <dU/dZ dV/dZ > G
c        ISELEP(1,442)                =  0
c                                  <dU/dX dW/dX > G
c        ISELEP(1,443)                =  0
c                                  <dU/dY dW/dY > G
c        ISELEP(1,444)                =  0
c                                  <dU/dZ dW/dZ > G
c        ISELEP(1,445)                =  0
c                                  <dV/dX dW/dX > G
c        ISELEP(1,446)                =  0
c                                  <dV/dY dW/dY > G
c        ISELEP(1,447)                =  0
c                                  <dV/dZ dW/dz > G
c        ISELEP(1,448)                =  0 
C
C                                  RANDBEHANDLUNG VON FRED
C                                  -----------------------
          ISELEP(1,301)            =   1
          ISELEP(1,302)            =   1
          ISELEP(1,303)            =   1
          ISELEP(1,304)            =   1
C
      RETURN
      END
