










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
C#ifdef _TSCAL_
      SUBROUTINE WIFAKSCA   (ITSTEP,ITMIT,WPHISCA,WCONSCA,WDIFSCA,
     $                       WSORSCA,IDUZSCA)
C*STARLET***************************************************************
C        W I F A K S C A  WEIGHTING FACTOR CALCULATION
C*STARLET***************************************************************
C
C PARAM: ITSTEP         - TIME STEP COUNTER
C        ITMIT          - LEAPFROG METHOD: EACH ITMIT TIMESTEP 
C                         ONE CLOSING STEP IS CARRIED OUT
C        WPHI           + WEIGHTING FACTOR FOR POINT VALUE TERM
C        WCON           + WEIGHTING FACTOR FOR CONVECTIVE TERM
C        WDIF           + WEIGHTING FACTOR FOR DIFFUSIVE TERM
C        WSOR           + WEIGHTING FACTOR FOR SOURCE TERM
C        IDUZSCA        - TIME LEVEL
C                         IDUZSCA = 1: OLD OLD 
C                         IDUZSCA = 2:     OLD 
C
C UPROG                 : ERRR
C
C VERS:  27.11.02 (TB)  : ORIGINAL BASED ON WIFAK FOR MOMENTUM EQUATION
C
C*STARLET***************************************************************
      END
