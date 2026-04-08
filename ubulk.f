










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
      SUBROUTINE CAL_UBULK(KK,JJ,II,NBND,DDX,
     $                    DDY,DDZ,U,W,B,
     $                    UGRID,UBULK,BP)
C***********************************************************************
C  U B U L K     BERECHNUNG VON UBULK BEI CHANNEL
C*************************************************** T.HUETTL 26.01.98 *
C
C***********************************************************************
C
C     UBULK = Sum ( U*Flaeche) / Sum (Flaeche)
C   
      IMPLICIT NONE
C
      INTEGER KK,JJ,II,NBND,IP1,IP3,I,J,K
      REAL    LVRFLU,LVRFLA,DDX(II),DDY(JJ),DDZ(KK),U(KK,JJ,II),
     $     UGRID,UBULK
      REAL    B(KK,JJ,II),W(KK,JJ,II),BP(KK,JJ,II)
C
      LVRFLA = 0.
      LVRFLU = 0. 
C

      I=II-NBND
C

      DO 1100 J=NBND+1,JJ-NBND
         DO 1000 K=NBND+1,KK-NBND
            LVRFLA = LVRFLA + (DDY(J) * DDZ (K))
     $                       *(SIGN(0.5,B(K,J,I)) + 0.5)
     $                       *BP(K,J,I)
            LVRFLU = LVRFLU + (U(K,J,I)*(DDY(J) * DDZ(K)))
     $                       *(SIGN(0.5,B(K,J,I)) + 0.5)
     $                       *BP(K,J,I)
 1000     CONTINUE
 1100  CONTINUE

C
      UBULK = LVRFLU / LVRFLA

C
      RETURN
      END
