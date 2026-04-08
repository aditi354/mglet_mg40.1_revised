










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
      SUBROUTINE INTERVY_PER(KK,JJ,II,KSTART,JSTART,ISTART,
     $     KSTOP,JSTOP,ISTOP,RSGS,V,VINJ,LCOL,DIAG,RCOL,UZ,RSP) 
C*MGLET***************************************************************
C        I N T E R V Y P E R
C                  = =
C        INTERPOLIERT V IN Y-RICHTUNG (VINJ) PERIODISCHE RB
C*MGLET***************************************************************
C
C PARAM: V(K,J,I)         - ZU INTERPOLIERENDE GROESSE(V)
C      : VINJ(K,J,I)      - INTERPOLIERTE V IN J(Y)-RICHTUNG
C
C DEFINE DIREKTIVEN     : KEINE 
C
C UPROG                 : TRIZYK3DJ
C
C VERS:  07.10.97 (AM)  : ORIGINAL        (KOMPAKT 4. ORDNUNG)
C
C*MGLET***************************************************************
C
      IMPLICIT NONE

      INTEGER KK,JJ,II,K,J,I,KSTART,JSTART,ISTART,KSTOP,JSTOP,ISTOP

      REAL    LCOL(JJ),DIAG(JJ),RCOL(JJ),UZ(JJ),RSP(JJ),
     $         V(KK,JJ,II), VINJ(KK,JJ,II), RSGS(KK,JJ,II)
C
C
C




       WRITE(6,*)'NOT IMPLEMENTED IN INTERVYPER.SRC'
       STOP
C
C                                ACHTER  ORDNUNG !!!!!!!!!!
C
C*MGLET***************************************************************
C
       RETURN
       END

