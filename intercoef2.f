










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
C$    LIST(ALL=0)
      SUBROUTINE INTERCOEF2(K,DX,DDX,IMX,LCOL,DIAG,RCOL,RSIDE1,RSIDE2)         
C*MGLET*****************************************************************
C              I N T E R C O E F 2
C*MGLET*****************************************************************
C
C     BERECHNUNG DER KOEFFIZIENTEN FUER DIE ERSTE ABLEITUNG
C     ANSATZ D R I T T E R  O R D N U N G   NACH LELE  
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
C 17.10.1996            : ADNAN MERI (ORIGINAL)
C
C*MGLET*****************************************************************
C

      REAL LCOL,RCOL,DIAG,RSIDE1,RSIDE2,ALFA,GAMA,BETA,PHI,
     $     l,a,b,c

      REAL DX(K),DDX(K)

      l = (DX(K)/DX(K-1))
      
      a = (1.0/2.0 ) * (DX(K)/DDX(K)) * ( 1.0 + 1.0/l)
      b = (1.0/8.0 ) * (DX(K)/DDX(K)) * ( l - 1.0/l)
      c = (1.0/24.0) * (DX(K)/DDX(K)) * ( l + 1.0/l)

      PHI  = l/((b+a)*(l-1) + (a-c))
      GAMA = l/((b+a)*(l-1) + (a-c))
      BETA = (-1.0 + (a+b)*GAMA)/(l + 1.0)
      ALFA = -1.0-BETA + a*GAMA
      WRITE (45,*) ALFA,BETA,GAMA,PHI
c
c
c*************BILDUNG DER MATRIX************************************
c
      LCOL=ALFA/PHI 
      DIAG=1.0/PHI 
      RCOL=BETA/PHI 
      RSIDE1=1.0       
      RSIDE2=1.0       
c 
c*******************************************************************
c 
      RETURN
      END

