










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
      SUBROUTINE SEL17    (ISELEP)
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
C                         FELDER <UT>, <VT>, <WT> UND <TT> 
C
C                         SELEKTIERT
C        04.05.04 (FS)  : VON SEL15 UEBERNOMMEN. STATISTISCHE AUSWERTUNG
C                         SKALARER GROESSEN
C*STARLET***************************************************************

      INTEGER  ISELEP(2,752)

      WRITE(6,*)'MELDUNG AUS SEL17: '
      WRITE(6,*)'STAT17 FUNKTIONIERT NUR MIT SCALARTRANSPORT'
      STOP
c
c                                   Momentangeschwindigkeits-Skalar
c                                   Korrelationen
c                                   Einpunktkorellation Skalar
c                                   -------------------------------
c
c                                  <UTT> G
        ISELEP(1,607)                =  1
c
c                                  <VTT> G
        ISELEP(1,608)                =  1
c
c                                  <WTT> G
        ISELEP(1,609)                =  1
c
c                                  <dT/dX^2> G
        ISELEP(1,610)                =  1
c
c                                  <dT/dY^2> G
        ISELEP(1,611)                =  1
C
C                                  <dT/dZ^2> G 
        ISELEP(1,612)                =  1


      RETURN
      END
