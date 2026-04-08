










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
      SUBROUTINE INTERCOEF3(II,I,DX,DDX,LCOL,DIAG,RCOL,
     $                      RSIDE1,RSIDE2,RSIDE3)         
C*MGLET*****************************************************************
C              I N T E R C O E F 3
C*MGLET*****************************************************************
C
C     BERECHNUNG DER KOEFFIZIENTEN FUER DIE KOMPAKTE BERECHNUNG
C     DER     E R S T E N     A B L E I T U N G 
C             V I E R T E R   O R D N U N G ! !
C
C     DIESER ANSATZ ENTHAELT  D R E I  UNBEKANNTE AUF DER RECHTEN
C     SEITE (Modifizierter Ansatz nach A D A M )
C
C   ! ACHTUNG : DIESE FORMULIERUNG IST NUR FUER M G L E T                
C   !           (GESTAEGERT -------> NICHTGESTAEGERT)
C
C*MGLET*****************************************************************
C
C
C DEFINE DIREKTIVEN     : KEINE
C
C UPROG                 : KEINE
C
C 03.10.1996            : ADNAN MERI (ORIGINAL)
C
C*MGLET*****************************************************************
C
      IMPLICIT NONE
      REAL LCOL,RCOL,DIAG,RSIDE1, RSIDE2,RSIDE3,ALFA,BETA,
     $     GAMA,LAMDA,PHI,
     $     s,k,a,b,c,d,e,f,g,p,h,q,l,m,n

      INTEGER I,II
      REAL DX(II),DDX(II)



      s = (DX(I)/DX(I-1))
      k = (DX(I-1)/DDX(I))
      
      a = s + 1.0
      c = s**2 - 1.0           
      e = s**3 + 1.0    
      b = - k * s * ((1./8.)*s +0.5)
      d = - k * s * ((1./24.)*s**2 - 0.5)
      f = - k * s * ((1./64.)*s**3 + 0.5)

      g = a*d - b*c 
      p = a*f - b*e                   
      h = -(3./8.)*k*c - (11./24.)*k*a	
      q = -(3./8.)*k*e + (31./64.)*k*a	

      PHI  = ((e-a)*g -(c+a)*p)/(q*g-p*h)
      GAMA = (c+a -h*PHI)/(a*d -b*c) 
      LAMDA= -GAMA - PHI
      BETA = (-1.0 -(3./8.)*k*PHI -b*GAMA)/a
      ALFA = -1.0 -BETA -0.5*k*(PHI - s*GAMA)

c
c
c*************BILDUNG DER MATRIX************************************
c
      LCOL=ALFA
      DIAG=1.0
      RCOL=BETA
      RSIDE1=GAMA
      RSIDE2=LAMDA
      RSIDE3=PHI
c 
c*******************************************************************
c 




      RETURN
      END

