










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
      SUBROUTINE BOUNDSTARTINTERSCA(II,I,DX,LCOL,DIAG,RCOL,
     $                              RSIDE1,RSIDE2,RSIDE3)
C---------------------------------------------------------------
C    BERECHNUNG DER KOEFFIZIENTEN FUER RANDBEDINGUNG SCALAR
C    INTERPOLATION  FIXED GRADIENT
C
C    26.11.03 (FS)
C---------------------------------------------------------------

      REAL LCOL,RCOL,DIAG,RSIDE1, RSIDE2,RSIDE3,
     $     alpha,beta,gamma

C     Fi + alpha*f(i+1/2) + beta*F(i+1) + gamma*f(i+3/2) = O**3
C     DIAG      RSIDE1       RCOL           RSIDE2 
C          h1=DX(i-1)/2     h2=DX(i-1)/2    h3=DX(i-1)/2
C                           +DX(i)/2        +DX(i)
      h1 = DX(I-1)/2.0
      h2 = DX(I-1)/2.0+DX(I)/2.0
      h3 = DX(I-1)/2.0+DX(I)

C--- KOEFFIZIENTEN DER 3x3 MATRIX
      a1 = h1
      a2 = h2
      a3 = h3
      b1 = 0.5*h1**2 
      b2=0.5*h2**2 
      b3=0.5*h3**2

      a= a1*b2 - a2*b1 - a1*b3 + b1*a3 + a2*b3 - a3*b2
      b= a3*b2 - a2*b3
      C= a1*b3 - b1*a3
      d= a2*b1 - a1*b2
      
      alpha = -b/a
      beta = c/a
      gamma = -d/a
      
      LCOL = 0.0
      DIAG = 1.0
      RCOL = beta
      RSIDE1 = alpha
      RSIDE2 = gamma
      RSIDE3 = 0.0
      
      
      RETURN
      END

      SUBROUTINE BOUNDSTOPINTERSCA(II,I,DX,LCOL,DIAG,RCOL,
     $                              RSIDE1,RSIDE2,RSIDE3)
C---------------------------------------------------------------
C    BERECHNUNG DER KOEFFIZIENTEN FUER RANDBEDINGUNG SCALAR
C    INTERPOLATION  FIXED GRADIENT
C
C    26.11.03 (FS)
C---------------------------------------------------------------

      REAL LCOL,RCOL,DIAG,RSIDE1, RSIDE2,RSIDE3,
     $     alpha,beta,gamma

C     Fi - alpha*f(i-1/2) - beta*F(i-1) - gamma*f(i-3/2) = O**3
C     DIAG      RSIDE1       RCOL           RSIDE2 
C          h1=DX(i-1)/2     h2=DX(i-1)/2    h3=DX(i-1)/2
C                           +DX(i-2)/2        +DX(i-2)
      h1 = DX(I-1)/2.0
      h2 = DX(I-1)/2.0+DX(I-2)/2.0
      h3 = DX(I-1)/2.0+DX(I-2)

C--- KOEFFIZIENTEN DER 3x3 MATRIX
      a1 = -h1 
      a2 = -h2
      a3 = -h3
      b1 = 0.5*h1**2 
      b2=0.5*h2**2 
      b3=0.5*h3**2

      a= a1*b2 - a2*b1 - a1*b3 + b1*a3 + a2*b3 - a3*b2
      b= a3*b2 - a2*b3
      C= a1*b3 - b1*a3
      d= a2*b1 - a1*b2
      
      alpha = -b/a
      beta = c/a
      gamma = -d/a
      
      LCOL = beta
      DIAG = 1.0
      RCOL = 0.0
      RSIDE1 = 0.0
      RSIDE2 = alpha
      RSIDE3 = gamma
      
      
      RETURN
      END
