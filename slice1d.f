










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
      SUBROUTINE SLICE1D (IIP,XP,IIS,XS,IPOS,CS,CSI)
C
C*MGLET*****************************************************************
C  S L I C E 1 D      GEBIETSZERLEGUNG VON 1D-FELDERN
C
C   ORIGINAL:   11. 2. 94 (MM)
C
C*MGLET*****************************************************************
C
C  PARAMETER
C             IIP            - DIMENSIONIERUNG DES PARENT
C              XP            - 1D-FELD DES PARENT
C
C             IIS            - DIMENSIONIERUNG DES SUBGITTERS
C              XS            - 1D-FELD DES SUBGITTERS
C
C            IPOS            - POSITION DES SUBGITTERS IM PARENTGITTER
C
C************************************************************************
C
      REAL XP(IIP)
      REAL XS(IIS)

      CHARACTER (LEN=1) CS,CSI
C
C-------------------------------------------------- CHECK AUF KONSISTENZ
C
      IF ( IIP .LT. IIS ) CALL ERRR (501,'SLICE1D')
C
C----------------------------------------------- SETZEN DER GITTERPUNKTE
C
      DO IS = 1,IIS

         XS(IS) = XP(IS+IPOS-3)

      ENDDO

C
C------------------------------------------------------------
C
C
C----------------------------------------------------------------
C
      RETURN
      END


