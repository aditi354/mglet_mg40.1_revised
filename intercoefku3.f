










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
C
      SUBROUTINE INTERCOEFKU3(I,II,DX,LCOL,DIAG,RCOL,RSIDE1,S)
C*MGLET*****************************************************************
C              I N T E R C O E F K U 3
C*MGLET*****************************************************************
C
C     BERECHNUNG DER KOEFFIZIENTEN FUER DIE KOMPAKTE INTER-
C     INTERPOLATION dritter ORDNUNG
C
C   ! ACHTUNG : DIESE FORMULIERUNG IST NUR FUER M G L E T                
C   !           (GESTAEGERT -------> NICHTGESTAEGERT)
C
C*MGLET*****************************************************************
C
C DEFINE DIREKTIVEN     : KEINE
C
C UPROG                 : KEINE
C
C 29.11.1999            : ADNAN MERI (ORIGINAL)
C
C*MGLET*****************************************************************
C
      REAL   b,DX (II),S,RSIDE1,RCOL,LCOL,DIAG
C
      b = (DX(I)/DX(I-1))
C
      IF (S .EQ. 1.) THEN
C
      DIAG  = 1.0
      RSIDE1  = (4.0*b)/((2.0*b) +1.0)
      RCOL = (-1.0 +(0.5 * RSIDE1)) / (b +1.0)
      LCOL = -1.0 +RSIDE1 -RCOL
C    
      ELSE
C
      DIAG  = 1.0
      RSIDE1  = (4.0)/(b +2.0)
      RCOL = ( -1.0 +((0.5 *b +1.0)*RSIDE1) ) / (b +1.0)
      LCOL = -1.0 +RSIDE1 -RCOL
C
      ENDIF

C     write (41,*) I, b, LCOL,DIAG,RCOL,RSIDE1
C
c*******************************************************************
c 
      RETURN
      END
C

      SUBROUTINE INTERCOEFKU32(I,II,DX,DDX,LCOL,DIAG,RCOL,RSIDE1,S)
C*MGLET*****************************************************************
C              I N T E R C O E F K U 3 2
C*MGLET*****************************************************************
C
C     BERECHNUNG DER KOEFFIZIENTEN FUER DIE KOMPAKTE INTER-
C     INTERPOLATION dritter ORDNUNG
C
C   ! ACHTUNG : DIESE FORMULIERUNG IST NUR FUER M G L E T                
C   !           (NICHT-GESTAEGERT -------> GESTAEGERT)
C
C*MGLET*****************************************************************
C
C DEFINE DIREKTIVEN     : KEINE
C
C UPROG                 : KEINE
C
C 29.11.1999            : ADNAN MERI (ORIGINAL)
C
C*MGLET*****************************************************************
C
      REAL   b,k,DX (II),DDX(II),S,RSIDE1,RCOL,LCOL,DIAG
C
      b = (DDX(I)/DDX(I-1))
C
      IF (S .EQ. 1.) THEN
C
      k = (DDX(I-1)/DX(I-1))
C
      DIAG  = 1.0
      RSIDE1  = ((k**2.0) *b)/( (k**2.0) *b - 0.5*k*b 
     $          + 0.5*k -0.25 )
      RCOL = ( -k +((k-0.5) * RSIDE1) ) / (k*(b +1.0))
      LCOL = -1.0 +RSIDE1 -RCOL
C    
      ELSE
C
      k = (DDX(I-1)/DX(I))
C
      DIAG  = 1.0
      RSIDE1  = ((k**2.0) *b)/( (k**2.0) *b + 0.5*k*b 
     $          - 0.5*k -0.25 )
      RCOL = ( -k +((k+0.5) * RSIDE1) ) / (k*(b +1.0))
      LCOL = -1.0 +RSIDE1 -RCOL
C
      ENDIF

C     write (41,*) I, b, LCOL,DIAG,RCOL,RSIDE1
C
c*******************************************************************
c 
      RETURN
      END

