










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
      SUBROUTINE SELECA  (ISELEP)
C*STARLET***************************************************************
C        S E L E C A      DIE VOM PROGRAMMBENUTZER IN SUBR. SELECT
C                         DURCHGEFUEHRTE AUSWAHL HINSICHTLICH DER AUS-
C                         WERTUNG STATIST. GROESSEN WIRD IN  'SELECA'
C                         AUF KONSISTENZ UND VOLLSTAENDIGKEIT HIN UEBER-
C                         PRUEFT.
C                         WURDE IN 'SELECT' UNABSICHTLICH DIE BELEGUNG
C                         VON FELDERN VERGESSEN, DIE ZUR AUSWERTUNG UN-
C                         BEDINGT NOETIG SIND, SO KORRIGIERT 'SELECA'
C                         (HOFFENTLICH !!) DIESE UNTERLASSUNG.
C  A C H T U N G:         DIESE ROUTINE MUSS ERWEITERT WERDEN, FALLS
C                         DIE STATIST. AUSWERTUNG ERWEITERT WIRD !
C                         DIE HIERARCHIE DER GROESSEN UNTEREINANDER
C                         SPIELT DIE ENTSCHEIDENDE ROLLE BEI DER
C                         EINORDNUNG DER ZUSAETZLICHEN PROGRAMM-
C                         STATEMENTS.
C*STARLET***************************************************************
C
C PARAM: ISELEP(2,752)  + STEUERFELD FUER DIE AUSWERTUNG
C
C UPROG                 : ERRR
C
C DEFINE-DIREKTIVEN     : INFO
C
C        05.08.86 (HW)  : ORIGINAL
C        29.08.86 (HW)  : ERWEITERUNG: BEI DEN KOMPONENTEN DES SCHUB-
C                         SPANNUNGSTENSORS KANN DIE SUMME AUS GROB-
C                         UND FEINSTRUKTURANTEIL (ENSEMBLE-MITTELWERT
C                         DER REYNOLDSSPANNUNG) AUSGEGEBEN WERDEN.
C        30.06.88 (HW)  : ERWEITERUNG: AUSWERTUNG DER DREI KOMPONENTEN
C                         DER VORTICITY, DER ENSTROPHIE UND DER
C                         HELIZITAET
C        02.12.88 (HW)  : DIMENSIONIERUNG DES ISELEA- UND
C                         ISELEP-FELDES VON 250 --> 752
C        05.12.88 (HW)  : ERWEITERUNG: AUSWERTUNG VON KORRELATIONS-
C                         FUNKTIONEN, -KOEFFIZIENTEN, LEISTUNGS-
C                         DICHTESPEKTREN UND HAEUFIGKEITSVERTEILUNGEN
C
C*STARLET***************************************************************
C
      INTEGER        ISELEP(2,752)
C
C                                  VERMERK IM ISELEP-FELD, DASS IN
C                                  DIBCA RICHTIG GELESEN WIRD (DAS ER-
C                                  WEITERTE ISELEP-FELD U. DAS ISLINP-
C                                  FELD)
C
      ISELEP (1,209) = 1
C
C                                  STATIST. GROESSEN DER U-KOMPONENTE
C                                  ----------------------------------
C
C                                  LEISTUNGSDICHTESPEKTRUM (Z-RI.)
                     IF(ISELEP(1, 24) .GE. 1) THEN
                                                 ISELEP(1, 24) = 1
                        IF(ISELEP(1, 24) .LE. 2) ISELEP(1, 25) = 1
                        IF(ISELEP(1,  6) .EQ. 0) ISELEP(1,  6) = 1
                     ENDIF
C
C                                  LEISTUNGSDICHTESPEKTRUM (Y-RI.)
                     IF(ISELEP(1, 22) .GE. 1) THEN
                                                 ISELEP(1, 22) = 1
                        IF(ISELEP(1, 22) .LE. 2) ISELEP(1, 23) = 1
                        IF(ISELEP(1,  6) .EQ. 0) ISELEP(1,  6) = 1
                     ENDIF
C
C                                  LEISTUNGSDICHTESPEKTRUM (X-RI.)
                     IF(ISELEP(1, 20) .GE. 1) THEN
                                                 ISELEP(1, 20) = 1
                        IF(ISELEP(1, 20) .LE. 2) ISELEP(1, 21) = 1
                        IF(ISELEP(1,  6) .EQ. 0) ISELEP(1,  6) = 1
                     ENDIF
C
C                                  AUTOKORRELATIONSFKT. (Z-RI.)
                     IF(ISELEP(1, 18) .GE. 1) THEN
                                                 ISELEP(1, 18) = 1
                        IF(ISELEP(1, 18) .LE. 2) ISELEP(1, 19) = 1
                        IF(ISELEP(1,  7) .EQ. 0) ISELEP(1,  7) = 1
                     ENDIF
C
C                                  AUTOKORRELATIONSFKT. (Y-RI.)
                     IF(ISELEP(1, 16) .GE. 1) THEN
                                                 ISELEP(1, 16) = 1
                        IF(ISELEP(1, 16) .LE. 2) ISELEP(1, 17) = 1
                        IF(ISELEP(1,  7) .EQ. 0) ISELEP(1,  7) = 1
                     ENDIF
C
C                                  AUTOKORRELATIONSFKT. (X-RI.)
                     IF(ISELEP(1, 14) .GE. 1) THEN
                                                 ISELEP(1, 14) = 1
                        IF(ISELEP(1, 14) .LE. 2) ISELEP(1, 15) = 1
                        IF(ISELEP(1,  7) .EQ. 0) ISELEP(1,  7) = 1
                     ENDIF
C
C                                  FLACHHEITSGRAD
                     IF(ISELEP(1, 12) .GE. 1) THEN
                        IF(ISELEP(1, 12) .LE. 2) ISELEP(1, 13) = 1
                        IF(ISELEP(1,  7) .EQ. 0) ISELEP(1,  7) = 1
                     ENDIF
C                                   SCHIEFE
                     IF(ISELEP(1, 10) .GE. 1) THEN
                        IF(ISELEP(1, 10) .LE. 2) ISELEP(1, 11) = 1
                        IF(ISELEP(1,  7) .EQ. 0) ISELEP(1,  7) = 1
                     ENDIF
C                                   <U-RMS> G+S WIRD AM ENDE DES LAUFES
C                                   NEU BESTIMMT
                  IF(ISELEP(1,  9) .GE. 1)    THEN
                     IF(ISELEP(1,  9) .GT. 1) THEN
                                                 ISELEP(1,  9) = 1
                     ENDIF
                     IF(ISELEP(1,  7) .EQ. 0)    ISELEP(1,  7) = 1
                     IF(ISELEP(1,  7) .EQ. 3)    ISELEP(1,  7) = 2
                     IF(ISELEP(1,109) .EQ. 0)    ISELEP(1,109) = 1
                     IF(ISELEP(1,109) .EQ. 3)    ISELEP(1,109) = 2
                  ENDIF
C                                   <U-RMS> GROBSTRUKTURANTEIL
 2007          IF(ISELEP(1,  7) .GE. 1)       THEN
                  IF(ISELEP(1,  7) .LE. 2)       ISELEP(1,  8) = 1
                  IF(ISELEP(1,  6) .EQ. 0)       ISELEP(1,  6) = 1
               ENDIF
C                                   U-FLUKTUATIONEN
 2006       IF(ISELEP(1,  6) .GE. 1)          THEN
               IF(ISELEP(1,  4) .EQ. 0)          ISELEP(1,  4) = 1
            ENDIF
C                                   <U> (ENSEMBLE-MITTELWERT)
         IF(ISELEP(1,  4) .GE. 1)             THEN
            IF(ISELEP(1,  4) .LE. 2)             ISELEP(1,  5) = 1
         ENDIF
C                                   U   (MOMENTANWERT)
      IF(ISELEP(1,  1) .EQ. 0)                   ISELEP(1,  1) = 1
      IF(ISELEP(1,  2) .EQ. 0)                   ISELEP(1,  2) = 1
C
C                                  STATIST. GROESSEN DER V-KOMPONENTE
C                                  ----------------------------------
C
C                                  LEISTUNGSDICHTESPEKTRUM (Z-RI.)
                     IF(ISELEP(1, 49) .GE. 1) THEN
                                                 ISELEP(1, 49) = 1
                        IF(ISELEP(1, 49) .LE. 2) ISELEP(1, 50) = 1
                        IF(ISELEP(1, 31) .EQ. 0) ISELEP(1, 31) = 1
                     ENDIF
C
C                                  LEISTUNGSDICHTESPEKTRUM (Y-RI.)
                     IF(ISELEP(1, 47) .GE. 1) THEN
                                                 ISELEP(1, 47) = 1
                        IF(ISELEP(1, 47) .LE. 2) ISELEP(1, 48) = 1
                        IF(ISELEP(1, 31) .EQ. 0) ISELEP(1, 31) = 1
                     ENDIF
C
C                                  LEISTUNGSDICHTESPEKTRUM (X-RI.)
                     IF(ISELEP(1, 45) .GE. 1) THEN
                                                 ISELEP(1, 45) = 1
                        IF(ISELEP(1, 45) .LE. 2) ISELEP(1, 46) = 1
                        IF(ISELEP(1, 31) .EQ. 0) ISELEP(1, 31) = 1
                     ENDIF
C
C                                  AUTOKORRELATIONSFKT. (Z-RI.)
                     IF(ISELEP(1, 43) .GE. 1) THEN
                                                 ISELEP(1, 43) = 1
                        IF(ISELEP(1, 43) .LE. 2) ISELEP(1, 44) = 1
                        IF(ISELEP(1, 32) .EQ. 0) ISELEP(1, 32) = 1
                     ENDIF
C
C                                  AUTOKORRELATIONSFKT. (Y-RI.)
                     IF(ISELEP(1, 41) .GE. 1) THEN
                                                 ISELEP(1, 41) = 1
                        IF(ISELEP(1, 41) .LE. 2) ISELEP(1, 42) = 1
                        IF(ISELEP(1, 32) .EQ. 0) ISELEP(1, 32) = 1
                     ENDIF
C
C                                  AUTOKORRELATIONSFKT. (X-RI.)
                     IF(ISELEP(1, 39) .GE. 1) THEN
                                                 ISELEP(1, 39) = 1
                        IF(ISELEP(1, 39) .LE. 2) ISELEP(1, 40) = 1
                        IF(ISELEP(1, 32) .EQ. 0) ISELEP(1, 32) = 1
                     ENDIF
C
C                                  FLACHHEITSGRAD
                     IF(ISELEP(1, 37) .GE. 1) THEN
                        IF(ISELEP(1, 37) .LE. 2) ISELEP(1, 38) = 1
                        IF(ISELEP(1, 32) .EQ. 0) ISELEP(1, 32) = 1
                     ENDIF
C                                   SCHIEFE
                     IF(ISELEP(1, 35) .GE. 1) THEN
                        IF(ISELEP(1, 35) .LE. 2) ISELEP(1, 36) = 1
                        IF(ISELEP(1, 32) .EQ. 0) ISELEP(1, 32) = 1
                     ENDIF
C                                   <V-RMS> G+S WIRD AM ENDE DES LAUFES
C                                   NEU BESTIMMT
                  IF(ISELEP(1, 34) .GE. 1)    THEN
                     IF(ISELEP(1, 34) .GT. 1) THEN
                                                 ISELEP(1, 34) = 1
                     ENDIF
                     IF(ISELEP(1, 32) .EQ. 0)    ISELEP(1, 32) = 1
                     IF(ISELEP(1, 32) .EQ. 3)    ISELEP(1, 32) = 2
                     IF(ISELEP(1,109) .EQ. 0)    ISELEP(1,109) = 1
                     IF(ISELEP(1,109) .EQ. 3)    ISELEP(1,109) = 2
                  ENDIF
C                                   <V-RMS> GROBSTRUKTURANTEIL
 2032          IF(ISELEP(1, 32) .GE. 1)       THEN
                  IF(ISELEP(1, 32) .LE. 2)       ISELEP(1, 33) = 1
                  IF(ISELEP(1, 31) .EQ. 0)       ISELEP(1, 31) = 1
               ENDIF
C                                   V-FLUKTUATIONEN
 2031       IF(ISELEP(1, 31) .GE. 1)          THEN
               IF(ISELEP(1, 29) .EQ. 0)          ISELEP(1, 29) = 1
            ENDIF
C                                   <V> (ENSEMBLE-MITTELWERT)
         IF(ISELEP(1, 29) .GE. 1)             THEN
            IF(ISELEP(1, 29) .LE. 2)             ISELEP(1, 30) = 1
         ENDIF
C                                   V   (MOMENTANWERT)
      IF(ISELEP(1, 26) .EQ. 0)                   ISELEP(1, 26) = 1
      IF(ISELEP(1, 27) .EQ. 0)                   ISELEP(1, 27) = 1
C
C                                  STATIST. GROESSEN DER W-KOMPONENTE
C                                  ----------------------------------
C
C                                  LEISTUNGSDICHTESPEKTRUM (Z-RI.)
                     IF(ISELEP(1, 74) .GE. 1) THEN
                                                 ISELEP(1, 74) = 1
                        IF(ISELEP(1, 74) .LE. 2) ISELEP(1, 75) = 1
                        IF(ISELEP(1, 56) .EQ. 0) ISELEP(1, 56) = 1
                     ENDIF
C
C                                  LEISTUNGSDICHTESPEKTRUM (Y-RI.)
                     IF(ISELEP(1, 72) .GE. 1) THEN
                                                 ISELEP(1, 72) = 1
                        IF(ISELEP(1, 72) .LE. 2) ISELEP(1, 73) = 1
                        IF(ISELEP(1, 56) .EQ. 0) ISELEP(1, 56) = 1
                     ENDIF
C
C                                  LEISTUNGSDICHTESPEKTRUM (X-RI.)
                     IF(ISELEP(1, 70) .GE. 1) THEN
                                                 ISELEP(1, 70) = 1
                        IF(ISELEP(1, 70) .LE. 2) ISELEP(1, 71) = 1
                        IF(ISELEP(1, 56) .EQ. 0) ISELEP(1, 56) = 1
                     ENDIF
C
C                                  AUTOKORRELATIONSFKT. (Z-RI.)
                     IF(ISELEP(1, 68) .GE. 1) THEN
                                                 ISELEP(1, 68) = 1
                        IF(ISELEP(1, 68) .LE. 2) ISELEP(1, 69) = 1
                        IF(ISELEP(1, 57) .EQ. 0) ISELEP(1, 57) = 1
                     ENDIF
C
C                                  AUTOKORRELATIONSFKT. (Y-RI.)
                     IF(ISELEP(1, 66) .GE. 1) THEN
                                                 ISELEP(1, 66) = 1
                        IF(ISELEP(1, 66) .LE. 2) ISELEP(1, 67) = 1
                        IF(ISELEP(1, 57) .EQ. 0) ISELEP(1, 57) = 1
                     ENDIF
C
C                                  AUTOKORRELATIONSFKT. (X-RI.)
                     IF(ISELEP(1, 64) .GE. 1) THEN
                                                 ISELEP(1, 64) = 1
                        IF(ISELEP(1, 64) .LE. 2) ISELEP(1, 65) = 1
                        IF(ISELEP(1, 57) .EQ. 0) ISELEP(1, 57) = 1
                     ENDIF
C
C                                  FLACHHEITSGRAD
                     IF(ISELEP(1, 62) .GE. 1) THEN
                        IF(ISELEP(1, 62) .LE. 2) ISELEP(1, 63) = 1
                        IF(ISELEP(1, 57) .EQ. 0) ISELEP(1, 57) = 1
                     ENDIF
C                                   SCHIEFE
                     IF(ISELEP(1, 60) .GE. 1) THEN
                        IF(ISELEP(1, 60) .LE. 2) ISELEP(1, 61) = 1
                        IF(ISELEP(1, 57) .EQ. 0) ISELEP(1, 57) = 1
                     ENDIF
C                                   <W-RMS> G+S WIRD AM ENDE DES LAUFES
C                                   NEU BESTIMMT
                  IF(ISELEP(1, 59) .GE. 1)    THEN
                     IF(ISELEP(1, 59) .GT. 1) THEN
                                                 ISELEP(1, 59) = 1
                     ENDIF
                     IF(ISELEP(1, 57) .EQ. 0)    ISELEP(1, 57) = 1
                     IF(ISELEP(1, 57) .EQ. 3)    ISELEP(1, 57) = 2
                     IF(ISELEP(1,109) .EQ. 0)    ISELEP(1,109) = 1
                     IF(ISELEP(1,109) .EQ. 3)    ISELEP(1,109) = 2
                  ENDIF
C                                   <W-RMS> GROBSTRUKTURANTEIL
 2057          IF(ISELEP(1, 57) .GE. 1)       THEN
                  IF(ISELEP(1, 57) .LE. 2)       ISELEP(1, 58) = 1
                  IF(ISELEP(1, 56) .EQ. 0)       ISELEP(1, 56) = 1
               ENDIF
C                                   W-FLUKTUATIONEN
 2056       IF(ISELEP(1, 56) .GE. 1)          THEN
               IF(ISELEP(1, 54) .EQ. 0)          ISELEP(1, 54) = 1
            ENDIF
C                                   <W> (ENSEMBLE-MITTELWERT)
         IF(ISELEP(1, 54) .GE. 1)             THEN
            IF(ISELEP(1, 54) .LE. 2)             ISELEP(1, 55) = 1
         ENDIF
C                                   W   (MOMENTANWERT)
      IF(ISELEP(1, 51) .EQ. 0)                   ISELEP(1, 51) = 1
      IF(ISELEP(1, 52) .EQ. 0)                   ISELEP(1, 52) = 1
C
C                                  STATIST. GROESSEN DER P-KOMPONENTE
C                                  ----------------------------------
C
C                                  FLACHHEITSGRAD
                     IF(ISELEP(1, 87) .GE. 1) THEN
                        IF(ISELEP(1, 87) .LE. 2) ISELEP(1, 88) = 1
                        IF(ISELEP(1, 82) .EQ. 0) ISELEP(1, 82) = 1
                     ENDIF
C                                   SCHIEFE
                     IF(ISELEP(1, 85) .GE. 1) THEN
                        IF(ISELEP(1, 85) .LE. 2) ISELEP(1, 86) = 1
                        IF(ISELEP(1, 82) .EQ. 0) ISELEP(1, 82) = 1
                     ENDIF
C                                   <P-RMS> G+S WIRD AM ENDE DES LAUFES
C                                   NEU BESTIMMT
*                 IF(ISELEP(1, 84) .GE. 1)    THEN
*                                                ISELEP(1, 84) = 1
*                    IF(ISELEP(1, 82) .EQ. 0)    ISELEP(1, 82) = 1
*                    IF(ISELEP(1, 82) .EQ. 3)    ISELEP(1, 82) = 2
*                    IF(ISELEP(1,109) .EQ. 0)    ISELEP(1,109) = 1
*                    IF(ISELEP(1,109) .EQ. 3)    ISELEP(1,109) = 2
*                 ENDIF
C                                   <P-RMS> GROBSTRUKTURANTEIL
               IF(ISELEP(1, 82) .GE. 1)       THEN
                  IF(ISELEP(1, 82) .LE. 2)       ISELEP(1, 83) = 1
                  IF(ISELEP(1, 81) .EQ. 0)       ISELEP(1, 81) = 1
               ENDIF
C                                   P-FLUKTUATIONEN
 2081       IF(ISELEP(1, 81) .GE. 1)          THEN
               IF(ISELEP(1, 79) .EQ. 0)          ISELEP(1, 79) = 1
            ENDIF
C                                   <P> (ENSEMBLE-MITTELWERT)
         IF(ISELEP(1, 79) .GE. 1)             THEN
            IF(ISELEP(1, 79) .LE. 2)             ISELEP(1, 80) = 1
         ENDIF
C                                   P   (MOMENTANWERT)
      IF(ISELEP(1, 76) .EQ. 0)                   ISELEP(1, 76) = 1
C
C                                  STATIST. GROESSEN DER G-KOMPONENTE
C                                  ----------------------------------
C
C                                   G   (MOMENTANWERT)
      IF(ISELEP(1,101) .EQ. 0)                   ISELEP(1,101) = 1
C
C                                  KINETISCHE ENERGIE DER SCHWANKUNGS-
C                                  KOMPONENTEN
C                                  -----------------------------------
C
         IF(ISELEP(1,108) .GE. 1)             THEN
            IF(ISELEP(1,108) .GT. 1) THEN
                                                 ISELEP(1,108) = 1
            ENDIF
C                                   GROBSTRUKTUR
            IF(ISELEP(1,106) .EQ. 0)             ISELEP(1,106) = 1
            IF(ISELEP(1,106) .EQ. 3)             ISELEP(1,106) = 2
C                                   FEINSTRUKTUR
            IF(ISELEP(1,109) .EQ. 0)             ISELEP(1,109) = 1
            IF(ISELEP(1,109) .EQ. 3)             ISELEP(1,109) = 2
         ENDIF
C
C                                   GROBSTRUKTUR
C
      IF(ISELEP(1,106) .GE. 1)                THEN
         IF(ISELEP(1,106) .LE. 2)                ISELEP(1,107) = 1
C                                   U-FLUKTUATION MUSS BEKANNT SEIN
         IF(ISELEP(1,  6) .EQ. 0) THEN
            ISELEP(1,  6) = 1
            GOTO 2006
         ENDIF
C                                   V-FLUKTUATION MUSS BEKANNT SEIN
         IF(ISELEP(1, 31) .EQ. 0) THEN
            ISELEP(1, 31) = 1
            GOTO 2031
         ENDIF
C                                   W-FLUKTUATION MUSS BEKANNT SEIN
         IF(ISELEP(1, 56) .EQ. 0) THEN
            ISELEP(1, 56) = 1
            GOTO 2056
         ENDIF
      ENDIF
C
C                                   FEINSTRUKTUR
C
      IF(ISELEP(1,109) .GE. 1)                THEN
         IF(ISELEP(1,109) .LE. 2)                ISELEP(1,110) = 1
         IF(ISELEP(1,101) .EQ. 0)                ISELEP(1,101) = 1
      ENDIF
C
C                                   AN DIESER STELLE WIRD SPAETER
C                                   DIE DISSIPATION BEHANDELT
C                                   -----------------------------
C
C****************
C
C                                   ANTEILE DES SPANNUNGSTENSORS
C                                   ----------------------------
C
C                                   HIER:  U - W
C                                   ------------
C
         IF(ISELEP(1,148) .GE. 1)             THEN
            IF(ISELEP(1,148) .GT. 1) THEN
                                                 ISELEP(1,148) = 1
            ENDIF
C                                   ANTEIL DER GROBSTRUKTUR
            IF(ISELEP(1,146) .EQ. 0)             ISELEP(1,146) = 1
            IF(ISELEP(1,146) .EQ. 3)             ISELEP(1,146) = 2
C                                   ANTEIL DER FEINSTRUKTUR
            IF(ISELEP(1,149) .EQ. 0)             ISELEP(1,149) = 1
            IF(ISELEP(1,149) .EQ. 3)             ISELEP(1,149) = 2
C                                   MOLEKULARER ANTEIL
            IF(ISELEP(1,151) .EQ. 0)             ISELEP(1,151) = 1
            IF(ISELEP(1,151) .EQ. 3)             ISELEP(1,151) = 2
         ENDIF
         IF(ISELEP(1,153) .GE. 1)             THEN
            IF(ISELEP(1,153) .GT. 1) THEN
                                                 ISELEP(1,153) = 1
            ENDIF
C                                   ANTEIL DER GROBSTRUKTUR
            IF(ISELEP(1,146) .EQ. 0)             ISELEP(1,146) = 1
            IF(ISELEP(1,146) .EQ. 3)             ISELEP(1,146) = 2
C                                   ANTEIL DER FEINSTRUKTUR
            IF(ISELEP(1,149) .EQ. 0)             ISELEP(1,149) = 1
            IF(ISELEP(1,149) .EQ. 3)             ISELEP(1,149) = 2
         ENDIF
C
C                                   ANTEIL DER GROBSTRUKTUR
C
      IF(ISELEP(1,146) .GE. 1)                THEN
         IF(ISELEP(1,146) .LE. 2)                ISELEP(1,147) = 1
C                                   U-FLUKTUATION MUSS BEKANNT SEIN
         IF(ISELEP(1,  6) .EQ. 0) THEN
            ISELEP(1,  6) = 1
            GOTO 2006
         ENDIF
C                                   W-FLUKTUATION MUSS BEKANNT SEIN
         IF(ISELEP(1, 56) .EQ. 0) THEN
            ISELEP(1, 56) = 1
            GOTO 2056
         ENDIF
      ENDIF
C
C                                   ANTEIL DER FEINSTRUKTUR
C
      IF(ISELEP(1,149) .GE. 1)                THEN
         IF(ISELEP(1,149) .LE. 2)                ISELEP(1,150) = 1
         IF(ISELEP(1,101) .EQ. 0)                ISELEP(1,101) = 1
      ENDIF
C
C                                   MOLEKULARER ANTEIL
C
      IF(ISELEP(1,151) .GE. 1)                THEN
         IF(ISELEP(1,151) .LE. 2)                ISELEP(1,152) = 1
      ENDIF
C
C                                   HIER:  V - W
C                                   ------------
C
         IF(ISELEP(1,158) .GE. 1)             THEN
            IF(ISELEP(1,158) .GT. 1) THEN
                                                 ISELEP(1,158) = 1
            ENDIF
C                                   ANTEIL DER GROBSTRUKTUR
            IF(ISELEP(1,156) .EQ. 0)             ISELEP(1,156) = 1
            IF(ISELEP(1,156) .EQ. 3)             ISELEP(1,156) = 2
C                                   ANTEIL DER FEINSTRUKTUR
            IF(ISELEP(1,159) .EQ. 0)             ISELEP(1,159) = 1
            IF(ISELEP(1,159) .EQ. 3)             ISELEP(1,159) = 2
C                                   MOLEKULARER ANTEIL
            IF(ISELEP(1,161) .EQ. 0)             ISELEP(1,161) = 1
            IF(ISELEP(1,161) .EQ. 3)             ISELEP(1,161) = 2
         ENDIF
         IF(ISELEP(1,163) .GE. 1)             THEN
            IF(ISELEP(1,163) .GT. 1) THEN
                                                 ISELEP(1,163) = 1
            ENDIF
C                                   ANTEIL DER GROBSTRUKTUR
            IF(ISELEP(1,156) .EQ. 0)             ISELEP(1,156) = 1
            IF(ISELEP(1,156) .EQ. 3)             ISELEP(1,156) = 2
C                                   ANTEIL DER FEINSTRUKTUR
            IF(ISELEP(1,159) .EQ. 0)             ISELEP(1,159) = 1
            IF(ISELEP(1,159) .EQ. 3)             ISELEP(1,159) = 2
         ENDIF
C
C                                   ANTEIL DER GROBSTRUKTUR
C
      IF(ISELEP(1,156) .GE. 1)                THEN
         IF(ISELEP(1,156) .LE. 2)                ISELEP(1,157) = 1
C                                   V-FLUKTUATION MUSS BEKANNT SEIN
         IF(ISELEP(1, 31) .EQ. 0) THEN
            ISELEP(1, 31) = 1
            GOTO 2031
         ENDIF
C                                   W-FLUKTUATION MUSS BEKANNT SEIN
         IF(ISELEP(1, 56) .EQ. 0) THEN
            ISELEP(1, 56) = 1
            GOTO 2056
         ENDIF
      ENDIF
C
C                                   ANTEIL DER FEINSTRUKTUR
C
      IF(ISELEP(1,159) .GE. 1)                THEN
         IF(ISELEP(1,159) .LE. 2)                ISELEP(1,160) = 1
         IF(ISELEP(1,101) .EQ. 0)                ISELEP(1,101) = 1
      ENDIF
C
C                                   MOLEKULARER ANTEIL
C
      IF(ISELEP(1,161) .GE. 1)                THEN
         IF(ISELEP(1,161) .LE. 2)                ISELEP(1,162) = 1
      ENDIF
C
C                                   HIER:  U - V
C                                   ------------
C
         IF(ISELEP(1,168) .GE. 1)             THEN
            IF(ISELEP(1,168) .GT. 1) THEN
                                                 ISELEP(1,168) = 1
            ENDIF
C                                   ANTEIL DER GROBSTRUKTUR
            IF(ISELEP(1,166) .EQ. 0)             ISELEP(1,166) = 1
            IF(ISELEP(1,166) .EQ. 3)             ISELEP(1,166) = 2
C                                   ANTEIL DER FEINSTRUKTUR
            IF(ISELEP(1,169) .EQ. 0)             ISELEP(1,169) = 1
            IF(ISELEP(1,169) .EQ. 3)             ISELEP(1,169) = 2
C                                   MOLEKULARER ANTEIL
            IF(ISELEP(1,171) .EQ. 0)             ISELEP(1,171) = 1
            IF(ISELEP(1,171) .EQ. 3)             ISELEP(1,171) = 2
         ENDIF
         IF(ISELEP(1,173) .GE. 1)             THEN
            IF(ISELEP(1,173) .GT. 1) THEN
                                                 ISELEP(1,173) = 1
            ENDIF
C                                   ANTEIL DER GROBSTRUKTUR
            IF(ISELEP(1,166) .EQ. 0)             ISELEP(1,166) = 1
            IF(ISELEP(1,166) .EQ. 3)             ISELEP(1,166) = 2
C                                   ANTEIL DER FEINSTRUKTUR
            IF(ISELEP(1,169) .EQ. 0)             ISELEP(1,169) = 1
            IF(ISELEP(1,169) .EQ. 3)             ISELEP(1,169) = 2
         ENDIF
C
C                                   ANTEIL DER GROBSTRUKTUR
C
      IF(ISELEP(1,166) .GE. 1)                THEN
         IF(ISELEP(1,166) .LE. 2)                ISELEP(1,167) = 1
C                                   U-FLUKTUATION MUSS BEKANNT SEIN
         IF(ISELEP(1,  6) .EQ. 0) THEN
            ISELEP(1,  6) = 1
            GOTO 2006
         ENDIF
C                                   V-FLUKTUATION MUSS BEKANNT SEIN
         IF(ISELEP(1, 31) .EQ. 0) THEN
            ISELEP(1, 31) = 1
            GOTO 2031
         ENDIF
      ENDIF
C
C                                   ANTEIL DER FEINSTRUKTUR
C
      IF(ISELEP(1,169) .GE. 1)                THEN
         IF(ISELEP(1,169) .LE. 2)                ISELEP(1,170) = 1
         IF(ISELEP(1,101) .EQ. 0)                ISELEP(1,101) = 1
      ENDIF
C
C                                   MOLEKULARER ANTEIL
C
      IF(ISELEP(1,171) .GE. 1)                THEN
         IF(ISELEP(1,171) .LE. 2)                ISELEP(1,172) = 1
      ENDIF
C
C                                  STATIST. GROESSEN DER HELIZITAET
C                                  ----------------------------------
C
C                                   <HELIZITAET-RMS> GROBSTRUKTURANTEIL
            IF(ISELEP(1,203) .GE. 1)       THEN
               IF(ISELEP(1,203) .LE. 2)       ISELEP(1,204) = 1
               IF(ISELEP(1,202) .EQ. 0)       ISELEP(1,202) = 1
            ENDIF
C                                   HELIZITAET-FLUKTUATIONEN
         IF(ISELEP(1,202) .GE. 1)          THEN
            IF(ISELEP(1,200) .EQ. 0)          ISELEP(1,200) = 1
         ENDIF
C                                   <HELIZITAET> (ENSEMBLE-MITTELWERT)
      IF(ISELEP(1,200) .GE. 1)             THEN
         IF(ISELEP(1,200) .LE. 2)             ISELEP(1,201) = 1
         IF(ISELEP(1,176) .EQ. 0)             ISELEP(1,176) =-1
         IF(ISELEP(1,182) .EQ. 0)             ISELEP(1,182) =-1
         IF(ISELEP(1,188) .EQ. 0)             ISELEP(1,188) =-1
      ENDIF
C
C                                  STATIST. GROESSEN DER ENSTROPHIE
C                                  ----------------------------------
C
C                                   <ENSTROPHIE-RMS> GROBSTRUKTURANTEIL
            IF(ISELEP(1,198) .GE. 1)       THEN
               IF(ISELEP(1,198) .LE. 2)       ISELEP(1,199) = 1
               IF(ISELEP(1,197) .EQ. 0)       ISELEP(1,197) = 1
            ENDIF
C                                   ENSTROPHIE-FLUKTUATIONEN
         IF(ISELEP(1,197) .GE. 1)          THEN
            IF(ISELEP(1,195) .EQ. 0)          ISELEP(1,195) = 1
         ENDIF
C                                   <ENSTROPHIE> (ENSEMBLE-MITTELWERT)
      IF(ISELEP(1,195) .GE. 1)             THEN
         IF(ISELEP(1,195) .LE. 2)             ISELEP(1,196) = 1
         IF(ISELEP(1,176) .EQ. 0)             ISELEP(1,176) =-1
         IF(ISELEP(1,182) .EQ. 0)             ISELEP(1,182) =-1
         IF(ISELEP(1,188) .EQ. 0)             ISELEP(1,188) =-1
      ENDIF
C
C                                  STATIST. GROESSEN DER X-KOMPONENTE
C                                  DER VORTICITY
C                                  ----------------------------------
C
C                                   <OMEGA-X-RMS> GROBSTRUKTURANTEIL
 2180          IF(ISELEP(1,180) .GE. 1)       THEN
                  IF(ISELEP(1,180) .LE. 2)       ISELEP(1,181) = 1
                  IF(ISELEP(1,179) .EQ. 0)       ISELEP(1,179) = 1
               ENDIF
C                                   OMEGA-X-FLUKTUATIONEN
 2179       IF(ISELEP(1,179) .GE. 1)          THEN
               IF(ISELEP(1,177) .EQ. 0)          ISELEP(1,177) = 1
            ENDIF
C                                   <OMEGA-X> (ENSEMBLE-MITTELWERT)
 2177    IF(ISELEP(1,177) .GE. 1)             THEN
            IF(ISELEP(1,177) .LE. 2)             ISELEP(1,178) = 1
            IF(ISELEP(1,176) .LE. 0)             ISELEP(1,176) = 1
         ENDIF
C                                   OMEGA-X (MOMENTANWERT)
      IF(ISELEP(1,176) .LT. 0)                   ISELEP(1,176) = 1
C
C                                  STATIST. GROESSEN DER Y-KOMPONENTE
C                                  DER VORTICITY
C                                  ----------------------------------
C
C                                   <OMEGA-Y-RMS> GROBSTRUKTURANTEIL
 2186          IF(ISELEP(1,186) .GE. 1)       THEN
                  IF(ISELEP(1,186) .LE. 2)       ISELEP(1,187) = 1
                  IF(ISELEP(1,185) .EQ. 0)       ISELEP(1,185) = 1
               ENDIF
C                                   OMEGA-Y-FLUKTUATIONEN
 2185       IF(ISELEP(1,185) .GE. 1)          THEN
               IF(ISELEP(1,183) .EQ. 0)          ISELEP(1,183) = 1
            ENDIF
C                                   <OMEGA-Y> (ENSEMBLE-MITTELWERT)
 2183    IF(ISELEP(1,183) .GE. 1)             THEN
            IF(ISELEP(1,183) .LE. 2)             ISELEP(1,184) = 1
            IF(ISELEP(1,182) .LE. 0)             ISELEP(1,182) = 1
         ENDIF
C                                   OMEGA-Y (MOMENTANWERT)
      IF(ISELEP(1,182) .LT. 0)                   ISELEP(1,182) = 1
C
C                                  STATIST. GROESSEN DER Z-KOMPONENTE
C                                  DER VORTICITY
C                                  ----------------------------------
C
C                                   <OMEGA-Z-RMS> GROBSTRUKTURANTEIL
 2192          IF(ISELEP(1,192) .GE. 1)       THEN
                  IF(ISELEP(1,192) .LE. 2)       ISELEP(1,193) = 1
                  IF(ISELEP(1,191) .EQ. 0)       ISELEP(1,191) = 1
               ENDIF
C                                   OMEGA-Z-FLUKTUATIONEN
 2191       IF(ISELEP(1,191) .GE. 1)          THEN
               IF(ISELEP(1,189) .EQ. 0)          ISELEP(1,189) = 1
            ENDIF
C                                   <OMEGA-Z> (ENSEMBLE-MITTELWERT)
 2189    IF(ISELEP(1,189) .GE. 1)             THEN
            IF(ISELEP(1,189) .LE. 2)             ISELEP(1,190) = 1
            IF(ISELEP(1,188) .LE. 0)             ISELEP(1,188) = 1
         ENDIF
C                                   OMEGA-Z (MOMENTANWERT)
      IF(ISELEP(1,188) .LT. 0)                   ISELEP(1,188) = 1
C
C                                  KREUZKORRELATIONSFUNKTIONEN
C                                  ---------------------------
C
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN X-RI. (U"V")
         IF(ISELEP(1,210) .GE. 1)             THEN
                                                 ISELEP(1,210) = 1
            IF(ISELEP(1,210) .LE. 2)             ISELEP(1,211) = 1
C                                   U-RMS MUSS BEKANNT SEIN
            IF(ISELEP(1,  7) .EQ. 0) THEN
               ISELEP(1,  7) = 1
               GOTO 2007
            ENDIF
C                                   V-RMS MUSS BEKANNT SEIN
            IF(ISELEP(1, 32) .EQ. 0) THEN
               ISELEP(1, 32) = 1
               GOTO 2032
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN Y-RI. (U"V")
         IF(ISELEP(1,212) .GE. 1)             THEN
                                                 ISELEP(1,212) = 1
            IF(ISELEP(1,212) .LE. 2)             ISELEP(1,213) = 1
C                                   U-RMS MUSS BEKANNT SEIN
            IF(ISELEP(1,  7) .EQ. 0) THEN
               ISELEP(1,  7) = 1
               GOTO 2007
            ENDIF
C                                   V-RMS MUSS BEKANNT SEIN
            IF(ISELEP(1, 32) .EQ. 0) THEN
               ISELEP(1, 32) = 1
               GOTO 2032
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN X-RI. (U"W")
         IF(ISELEP(1,216) .GE. 1)             THEN
                                                 ISELEP(1,216) = 1
            IF(ISELEP(1,216) .LE. 2)             ISELEP(1,217) = 1
C                                   U-RMS MUSS BEKANNT SEIN
            IF(ISELEP(1,  7) .EQ. 0) THEN
               ISELEP(1,  7) = 1
               GOTO 2007
            ENDIF
C                                   W-RMS MUSS BEKANNT SEIN
            IF(ISELEP(1, 57) .EQ. 0) THEN
               ISELEP(1, 57) = 1
               GOTO 2057
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN Y-RI. (U"W")
         IF(ISELEP(1,218) .GE. 1)             THEN
                                                 ISELEP(1,218) = 1
            IF(ISELEP(1,218) .LE. 2)             ISELEP(1,219) = 1
C                                   U-RMS MUSS BEKANNT SEIN
            IF(ISELEP(1,  7) .EQ. 0) THEN
               ISELEP(1,  7) = 1
               GOTO 2007
            ENDIF
C                                   W-RMS MUSS BEKANNT SEIN
            IF(ISELEP(1, 57) .EQ. 0) THEN
               ISELEP(1, 57) = 1
               GOTO 2057
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN Y-RI. (V"W")
         IF(ISELEP(1,224) .GE. 1)             THEN
                                                 ISELEP(1,224) = 1
            IF(ISELEP(1,224) .LE. 2)             ISELEP(1,225) = 1
C                                   V-RMS MUSS BEKANNT SEIN
            IF(ISELEP(1, 32) .EQ. 0) THEN
               ISELEP(1, 32) = 1
               GOTO 2032
            ENDIF
C                                   W-RMS MUSS BEKANNT SEIN
            IF(ISELEP(1, 57) .EQ. 0) THEN
               ISELEP(1, 57) = 1
               GOTO 2057
            ENDIF
         ENDIF
C
C                                  KORRELATIONSKOEFFIZIENTEN
C                                  -------------------------
C
C                                  ENSEMBLE-MITTELWERT DES KORRELATIONS-
C                                  KOEFF. ZWISCHEN U" UND W" IN X-RI.
         IF(ISELEP(1,234) .GE. 1)             THEN
                                                 ISELEP(1,234) = 1
            IF(ISELEP(1,234) .LE. 2)             ISELEP(1,235) = 1
C                                   U-RMS MUSS BEKANNT SEIN
            IF(ISELEP(1,  7) .EQ. 0) THEN
               ISELEP(1,  7) = 1
               GOTO 2007
            ENDIF
C                                   W-RMS MUSS BEKANNT SEIN
            IF(ISELEP(1, 57) .EQ. 0) THEN
               ISELEP(1, 57) = 1
               GOTO 2057
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DES KORRELATIONS-
C                                  KOEFF. ZWISCHEN U" UND W" IN Z-RI.
         IF(ISELEP(1,238) .GE. 1)             THEN
                                                 ISELEP(1,238) = 1
            IF(ISELEP(1,238) .LE. 2)             ISELEP(1,239) = 1
C                                   U-RMS MUSS BEKANNT SEIN
            IF(ISELEP(1,  7) .EQ. 0) THEN
               ISELEP(1,  7) = 1
               GOTO 2007
            ENDIF
C                                   W-RMS MUSS BEKANNT SEIN
            IF(ISELEP(1, 57) .EQ. 0) THEN
               ISELEP(1, 57) = 1
               GOTO 2057
            ENDIF
         ENDIF
C
C                                  OBERE DREIECKSMATRIX DES TENSORS
C                                  DER KORRELATIONSFUNKTIONEN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  --------------------------------
C
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN X-RI. (OX"OX")
         IF(ISELEP(1,246) .GE. 1)             THEN
                                                 ISELEP(1,246) = 1
            IF(ISELEP(1,246) .LE. 2)             ISELEP(1,247) = 1
C                                   OMEGA-X RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,180) .EQ. 0) THEN
               ISELEP(1,180) = 1
               GOTO 2180
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN X-RI. (OX"OY")
         IF(ISELEP(1,248) .GE. 1)             THEN
                                                 ISELEP(1,248) = 1
            IF(ISELEP(1,248) .LE. 2)             ISELEP(1,249) = 1
C                                   OMEGA-X RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,180) .EQ. 0) THEN
               ISELEP(1,180) = 1
               GOTO 2180
            ENDIF
C                                   OMEGA-Y RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,186) .EQ. 0) THEN
               ISELEP(1,186) = 1
               GOTO 2186
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN X-RI. (OX"OZ")
         IF(ISELEP(1,250) .GE. 1)             THEN
                                                 ISELEP(1,250) = 1
            IF(ISELEP(1,250) .LE. 2)             ISELEP(1,251) = 1
C                                   OMEGA-X RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,180) .EQ. 0) THEN
               ISELEP(1,180) = 1
               GOTO 2180
            ENDIF
C                                   OMEGA-Z RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,192) .EQ. 0) THEN
               ISELEP(1,192) = 1
               GOTO 2192
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN X-RI. (OY"OY")
         IF(ISELEP(1,252) .GE. 1)             THEN
                                                 ISELEP(1,252) = 1
            IF(ISELEP(1,252) .LE. 2)             ISELEP(1,253) = 1
C                                   OMEGA-Y RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,186) .EQ. 0) THEN
               ISELEP(1,186) = 1
               GOTO 2186
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN X-RI. (OY"OZ")
         IF(ISELEP(1,254) .GE. 1)             THEN
                                                 ISELEP(1,254) = 1
            IF(ISELEP(1,254) .LE. 2)             ISELEP(1,255) = 1
C                                   OMEGA-Y RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,186) .EQ. 0) THEN
               ISELEP(1,186) = 1
               GOTO 2186
            ENDIF
C                                   OMEGA-Z RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,192) .EQ. 0) THEN
               ISELEP(1,192) = 1
               GOTO 2192
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN X-RI. (OZ"OZ")
         IF(ISELEP(1,256) .GE. 1)             THEN
                                                 ISELEP(1,256) = 1
            IF(ISELEP(1,256) .LE. 2)             ISELEP(1,257) = 1
C                                   OMEGA-Z RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,192) .EQ. 0) THEN
               ISELEP(1,192) = 1
               GOTO 2192
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN Y-RI. (OX"OX")
         IF(ISELEP(1,258) .GE. 1)             THEN
                                                 ISELEP(1,258) = 1
            IF(ISELEP(1,258) .LE. 2)             ISELEP(1,259) = 1
C                                   OMEGA-X RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,180) .EQ. 0) THEN
               ISELEP(1,180) = 1
               GOTO 2180
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN Y-RI. (OX"OY")
         IF(ISELEP(1,260) .GE. 1)             THEN
                                                 ISELEP(1,260) = 1
            IF(ISELEP(1,260) .LE. 2)             ISELEP(1,261) = 1
C                                   OMEGA-X RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,180) .EQ. 0) THEN
               ISELEP(1,180) = 1
               GOTO 2180
            ENDIF
C                                   OMEGA-Y RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,186) .EQ. 0) THEN
               ISELEP(1,186) = 1
               GOTO 2186
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN Y-RI. (OX"OZ")
         IF(ISELEP(1,262) .GE. 1)             THEN
                                                 ISELEP(1,262) = 1
            IF(ISELEP(1,262) .LE. 2)             ISELEP(1,263) = 1
C                                   OMEGA-X RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,180) .EQ. 0) THEN
               ISELEP(1,180) = 1
               GOTO 2180
            ENDIF
C                                   OMEGA-Z RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,192) .EQ. 0) THEN
               ISELEP(1,192) = 1
               GOTO 2192
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN Y-RI. (OY"OY")
         IF(ISELEP(1,264) .GE. 1)             THEN
                                                 ISELEP(1,264) = 1
            IF(ISELEP(1,264) .LE. 2)             ISELEP(1,265) = 1
C                                   OMEGA-Y RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,186) .EQ. 0) THEN
               ISELEP(1,186) = 1
               GOTO 2186
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN Y-RI. (OY"OZ")
         IF(ISELEP(1,266) .GE. 1)             THEN
                                                 ISELEP(1,266) = 1
            IF(ISELEP(1,266) .LE. 2)             ISELEP(1,267) = 1
C                                   OMEGA-Y RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,186) .EQ. 0) THEN
               ISELEP(1,186) = 1
               GOTO 2186
            ENDIF
C                                   OMEGA-Z RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,192) .EQ. 0) THEN
               ISELEP(1,192) = 1
               GOTO 2192
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN Y-RI. (OZ"OZ")
         IF(ISELEP(1,268) .GE. 1)             THEN
                                                 ISELEP(1,268) = 1
            IF(ISELEP(1,268) .LE. 2)             ISELEP(1,269) = 1
C                                   OMEGA-Z RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,192) .EQ. 0) THEN
               ISELEP(1,192) = 1
               GOTO 2192
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN Z-RI. (OX"OX")
         IF(ISELEP(1,270) .GE. 1)             THEN
                                                 ISELEP(1,270) = 1
            IF(ISELEP(1,270) .LE. 2)             ISELEP(1,271) = 1
C                                   OMEGA-X RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,180) .EQ. 0) THEN
               ISELEP(1,180) = 1
               GOTO 2180
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN Z-RI. (OX"OY")
         IF(ISELEP(1,272) .GE. 1)             THEN
                                                 ISELEP(1,272) = 1
            IF(ISELEP(1,272) .LE. 2)             ISELEP(1,273) = 1
C                                   OMEGA-X RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,180) .EQ. 0) THEN
               ISELEP(1,180) = 1
               GOTO 2180
            ENDIF
C                                   OMEGA-Y RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,186) .EQ. 0) THEN
               ISELEP(1,186) = 1
               GOTO 2186
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN Z-RI. (OX"OZ")
         IF(ISELEP(1,274) .GE. 1)             THEN
                                                 ISELEP(1,274) = 1
            IF(ISELEP(1,274) .LE. 2)             ISELEP(1,275) = 1
C                                   OMEGA-X RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,180) .EQ. 0) THEN
               ISELEP(1,180) = 1
               GOTO 2180
            ENDIF
C                                   OMEGA-Z RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,192) .EQ. 0) THEN
               ISELEP(1,192) = 1
               GOTO 2192
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN Z-RI. (OY"OY")
         IF(ISELEP(1,276) .GE. 1)             THEN
                                                 ISELEP(1,276) = 1
            IF(ISELEP(1,276) .LE. 2)             ISELEP(1,277) = 1
C                                   OMEGA-Y RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,186) .EQ. 0) THEN
               ISELEP(1,186) = 1
               GOTO 2186
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER KREUZKORRE-
C                                  LATIONSFUNKTION IN Z-RI. (OY"OZ")
         IF(ISELEP(1,278) .GE. 1)             THEN
                                                 ISELEP(1,278) = 1
            IF(ISELEP(1,278) .LE. 2)             ISELEP(1,279) = 1
C                                   OMEGA-Y RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,186) .EQ. 0) THEN
               ISELEP(1,186) = 1
               GOTO 2186
            ENDIF
C                                   OMEGA-Z RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,192) .EQ. 0) THEN
               ISELEP(1,192) = 1
               GOTO 2192
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER AUTOKORRE-
C                                  LATIONSFUNKTION IN Z-RI. (OZ"OZ")
         IF(ISELEP(1,280) .GE. 1)             THEN
                                                 ISELEP(1,280) = 1
            IF(ISELEP(1,280) .LE. 2)             ISELEP(1,281) = 1
C                                   OMEGA-Z RMS  MUSS BEKANNT SEIN
            IF(ISELEP(1,192) .EQ. 0) THEN
               ISELEP(1,192) = 1
               GOTO 2192
            ENDIF
         ENDIF
C
C                                  LEISTUNGSDICHTESPEKTREN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  ----------------------------
C
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN X-RI. (OX"OX")
         IF(ISELEP(1,282) .GE. 1)             THEN
                                                 ISELEP(1,282) = 1
            IF(ISELEP(1,282) .LE. 2)             ISELEP(1,283) = 1
C                                   OMEGA-X FLUKT. MUSS BEKANNT SEIN
            IF(ISELEP(1,179) .EQ. 0) THEN
               ISELEP(1,179) = 1
               GOTO 2179
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN X-RI. (OY"OY")
         IF(ISELEP(1,284) .GE. 1)             THEN
                                                 ISELEP(1,284) = 1
            IF(ISELEP(1,284) .LE. 2)             ISELEP(1,285) = 1
C                                   OMEGA-Y FLUKT. MUSS BEKANNT SEIN
            IF(ISELEP(1,185) .EQ. 0) THEN
               ISELEP(1,185) = 1
               GOTO 2185
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN X-RI. (OZ"OZ")
         IF(ISELEP(1,286) .GE. 1)             THEN
                                                 ISELEP(1,286) = 1
            IF(ISELEP(1,286) .LE. 2)             ISELEP(1,287) = 1
C                                   OMEGA-Z FLUKT. MUSS BEKANNT SEIN
            IF(ISELEP(1,191) .EQ. 0) THEN
               ISELEP(1,191) = 1
               GOTO 2191
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN Y-RI. (OX"OX")
         IF(ISELEP(1,288) .GE. 1)             THEN
                                                 ISELEP(1,288) = 1
            IF(ISELEP(1,288) .LE. 2)             ISELEP(1,289) = 1
C                                   OMEGA-X FLUKT. MUSS BEKANNT SEIN
            IF(ISELEP(1,179) .EQ. 0) THEN
               ISELEP(1,179) = 1
               GOTO 2179
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN Y-RI. (OY"OY")
         IF(ISELEP(1,290) .GE. 1)             THEN
                                                 ISELEP(1,290) = 1
            IF(ISELEP(1,290) .LE. 2)             ISELEP(1,291) = 1
C                                   OMEGA-Y FLUKT. MUSS BEKANNT SEIN
            IF(ISELEP(1,185) .EQ. 0) THEN
               ISELEP(1,185) = 1
               GOTO 2185
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DES LEISTUNGS-
C                                  DICHTESPEKTRUMS IN Y-RI. (OZ"OZ")
         IF(ISELEP(1,292) .GE. 1)             THEN
                                                 ISELEP(1,292) = 1
            IF(ISELEP(1,292) .LE. 2)             ISELEP(1,293) = 1
C                                   OMEGA-Z FLUKT. MUSS BEKANNT SEIN
            IF(ISELEP(1,191) .EQ. 0) THEN
               ISELEP(1,191) = 1
               GOTO 2191
            ENDIF
         ENDIF
C
C                                  HAEUFIGKEITSVERTEILUNG DER
C                                  INKLINATIONSWINKEL DER VORTICITY-
C                                  VEKTOREN
C                                  ---------------------------------
C
C                                  ENSEMBLE-MITTELWERT DER HAEUFIG-
C                                  KEITSVERTEILUNG VON  ATAN(OZ/OY)
         IF(ISELEP(1,306) .GE. 1)             THEN
                                                 ISELEP(1,306) = 1
            IF(ISELEP(1,306) .LE. 2)             ISELEP(1,307) = 1
C                                   <OMEGA-Z> MUSS BEKANNT SEIN
            IF(ISELEP(1,189) .EQ. 0) THEN
               ISELEP(1,189) = 1
               GOTO 2189
            ENDIF
C                                   <OMEGA-Y> MUSS BEKANNT SEIN
            IF(ISELEP(1,183) .EQ. 0) THEN
               ISELEP(1,183) = 1
               GOTO 2183
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER HAEUFIG-
C                                  KEITSVERTEILUNG VON  ATAN(OZ/OX)
         IF(ISELEP(1,308) .GE. 1)             THEN
                                                 ISELEP(1,308) = 1
            IF(ISELEP(1,308) .LE. 2)             ISELEP(1,309) = 1
C                                   <OMEGA-Z> MUSS BEKANNT SEIN
            IF(ISELEP(1,189) .EQ. 0) THEN
               ISELEP(1,189) = 1
               GOTO 2189
            ENDIF
C                                   <OMEGA-X> MUSS BEKANNT SEIN
            IF(ISELEP(1,177) .EQ. 0) THEN
               ISELEP(1,177) = 1
               GOTO 2177
            ENDIF
         ENDIF
C
C                                  ENSEMBLE-MITTELWERT DER HAEUFIG-
C                                  KEITSVERTEILUNG VON  ATAN(OY/OX)
         IF(ISELEP(1,310) .GE. 1)             THEN
                                                 ISELEP(1,310) = 1
            IF(ISELEP(1,310) .LE. 2)             ISELEP(1,311) = 1
C                                   <OMEGA-Y> MUSS BEKANNT SEIN
            IF(ISELEP(1,183) .EQ. 0) THEN
               ISELEP(1,183) = 1
               GOTO 2183
            ENDIF
C                                   <OMEGA-X> MUSS BEKANNT SEIN
            IF(ISELEP(1,177) .EQ. 0) THEN
               ISELEP(1,177) = 1
               GOTO 2177
            ENDIF
         ENDIF
C
C
C                                  SICHERHEITS-CHECK DES FELD-INHALTES
C
      DO 100 N = 1,752
         IF(ISELEP(1,N) .LT. 0  .OR. ISELEP(1,N) .GT. 3)
     $      CALL ERRR (501,' SELECA   ')
  100 CONTINUE
C
C
      RETURN
      END
