










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
      SUBROUTINE INTERUX_PER(KK,JJ,II,KSTART,JSTART,ISTART,
     $     KSTOP,JSTOP,ISTOP,RSGS,U,UINI,LCOL,DIAG,RCOL,UZ,RSP)
C*MGLET***************************************************************
C        I N T E R U X P E R
C                  = =
C        INTERPOLIERT U IN X-RICHTUNG (UINI) PERIODISCHE RB
C*MGLET***************************************************************
C
C PARAM: U(K,J,I)         - ZU INTERPOLIERENDE GROESSE(U)
C      : UINI(K,J,I)      - INTERPOLIERTE U IN I(X)-RICHTUNG
C
C DEFINE DIREKTIVEN     : KEINE 
C
C UPROG                 : TRIZYK3D
C
C VERS:  07.06.97 (AM)  : ORIGINAL        (KOMPAKT 4. ORDNUNG)
C
C*MGLET***************************************************************
C
      IMPLICIT NONE

      INTEGER KK,JJ,II,K,J,I,KSTART,JSTART,ISTART,KSTOP,JSTOP,ISTOP

      REAL       LCOL(II),DIAG(II),RCOL(II),UZ(II),RSP(II),
     $           U(KK,JJ,II),UINI(KK,JJ,II),RSGS(KK,JJ,II)
C
C
C


C
C
C*MGLET***************************************************************
C
       RETURN
       END

