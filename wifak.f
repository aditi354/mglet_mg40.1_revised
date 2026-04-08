










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
      SUBROUTINE WIFAK   (ITSTEP,ITMIT,WPHI,WKON,WDIF,WSOR,IDUZ)
C*STARLET***************************************************************
C        W I F A K        DIE WICHTUNGSFAKTOREN FUER DIE KONVEKTIVEN,
C                         DIFFUSIVEN UND PHI-TERME DES ZEITSCHRITTES
C                         WERDEN BERECHNET.
C*STARLET***************************************************************
C
C PARAM: ITSTEP         - ZEITSCHRITTZAEHLER
C        ITMIT          - NACH JEWEILS ITMIT ZEITSCHRITTEN WIRD EIN AB-
C                         SCHLUSSCHRITT FUER DAS LEAPFROG-VERFAHREN GE-
C                         MACHT
C        WPHI           + WICHTUNGSFAKTOR FUER DEN PUNKTWERT
C        WKON           + WICHTUNGSFAKTOR FUER DIE KONVEKTIVEN TERME
C        WDIF           + WICHTUNGSFAKTOR FUER DIE DIFFUSIVEN TERME
C        WSOR           + WICHTUNGSFAKTOR FUER DEN QUELLTERM
C        IDUZ           - ENTHAELT DIE INFORMATION UEBER DAS ZEITNIVEAU
C                         IDUZ = 1: OLD OLD -WERTE IM KERNSPEICHER
C                         IDUZ = 2:     OLD -WERTE IM KERNSPEICHER
C
C UPROG                 : ERRR
C
C VERS:  23.08.85 (HW)  : ORIGINAL
C        27.08.85 (HW)  : WSOR ZUSAETZLICH EINGEFUEHRT
C
C*STARLET***************************************************************
      IF(IDUZ .NE. 2) CALL ERRR(506,' WIFAK    ')
         WPHI = 0.0
         WKON = 1.0
         WDIF = 1.0
         WSOR = 1.0
         RETURN
      END
