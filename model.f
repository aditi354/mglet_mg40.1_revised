










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
      SUBROUTINE MODEL (KK,JJ,II,FINC,DUI,DUJ,TAU)
C
C**********************************************************************
C     M O D E L
C----------------------------------------------------------------------
C     SGS STRESS TENSOR: INC MODEL
C----------------------------------------------------------------------
C     C. BRUN
C     LEHRSTUHL FUER FLUIDMECHANIK, TU MUENCHEN
C----------------------------------------------------------------------
C     PROGRAM DEVELOPING:                          C. BRUN / 04.03.1999
C     LAST MODIFICATION :                                  /
C**********************************************************************
C
      IMPLICIT NONE
C
C
      integer KK,JJ,II
      integer K,J,I
C
      REAL    DUI(KK,JJ,II),   DUJ(KK,JJ,II),   TAU(KK,JJ,II)
      REAL    FINC(KK,JJ,II)
C
C---- CALCULATE (Tau_ij)
C
C     2D HOMOGENEOUS FILTERING ONLY !!!!
C     BRUN (ERCOF 99)
C
C
      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1
            TAU(K,J,I) = FINC(K,J,I)*DUI(K,J,I)*DUJ(K,J,I)
            ENDDO
         ENDDO
      ENDDO
C
C
      RETURN
      END
