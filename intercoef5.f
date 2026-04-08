










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
      SUBROUTINE INTERCOEF5(II,I,DX,DDX,LCOL,DIAG,RCOL,
     $                      RSIDE1,RSIDE2,RSIDE3)         
C*MGLET*****************************************************************
C              I N T E R C O E F 5
C*MGLET*****************************************************************
C
C     BERECHNUNG DER KOEFFIZIENTEN FUER DIE KOMPAKTE BERECHNUNG
C     DER     E R S T E N     A B L E I T U N G 
C
C     DIESER ANSATZ ENTHAELT  D R E I  UNBEKANNTE AUF DER RECHTEN
C     SEITE (Modifizierter Ansatz nach A D A M )
C
C   ! ACHTUNG : DIESE FORMULIERUNG IST NUR FUER M G L E T                
C   !           (NICHT-GESTAEGERT -------> GESTAEGERT)
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

      REAL LCOL,RCOL,DIAG,RSIDE1,RSIDE2,RSIDE3,ALFA,BETA,
     $     GAMA,LAMDA,PHI,
     $     a,b,c,d,e,f,g,h,l,m,n

      REAL DX(II),DDX(II)


      a = DDX(I+1)/DDX(I)
      b = DX(I)/DDX(I)
      
      c = -3.0*(b+4.)*(a-1.) + (b**2 -12.)
      d = -3.0*(b-4.)*(a-1.) - (b**2 -12.)
      e = -24.*a             
      f = -8.*(b**3 +32.)*(a+1.) + 64.*(a**3 +1.)*(b+4.)
      g = -8.*(b**3 -32.)*(a+1.) + 64.*(a**3 +1.)*(b-4.)
      h = 64.*8.*(-a + a**3 )


      PHI  = ((h*c)-(f*e)) /((g*c)-(f*d))
      GAMA = (e -d*PHI)/c
      LAMDA= -GAMA - PHI
      BETA = (-8. + (b-4.)*PHI + (b+4.)*GAMA)/(8.*(a+1.))
      ALFA = -1. + 0.5*(GAMA-PHI) - BETA              

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

