










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
      SUBROUTINE INTERCOEFKOBS0(KK,K,DX,DDX,LCOL,DIAG,RCOL,
     $			    RSIDE1,RSIDE2,RSIDE3)         
C*MGLET*****************************************************************
C              I N T E R C O E F 1
C*MGLET*****************************************************************
C
C     BERECHNUNG DER KOEFFIZIENTEN FUER DIE KOMPAKTE INTER-
C     INTERPOLATION VIERTER ORDNUNG
C
C     DIESER ANSATZ ENTHAELT  Z W E I  UNBEKANNTE AUF DER RECHTEN
C     SEITE (Modifizierter Ansatz)
C
C   ! ACHTUNG : DIESE FORMULIERUNG IST NUR FUER M G L E T             
C   !           (NICHTGESTAEGERT -------> GESTAEGERT)
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

      REAL LCOL,RCOL,DIAG,RSIDE1,RSIDE2,RSIDE3,
     $     ALFA,GAMA,BETA,PHI,DX(KK),DDX(KK),
     $     a,b,c,d,e,f



      a = DDX(K-1) /DX(K-1)
      b = DDX(K)/DX(K-1)
      c = -4.0 * a**2 + 1.0
      d = -4.0 * b**2 + 1.0
      e = -8.0 * a**3 + 2.0 * a
      f =  8.0 * b**3 - 2.0 * b

      ALFA = f/(-c*f + e*d)         
      BETA = (- 1.0 -ALFA*c)/d    
      GAMA = 0.5 * ( (-2.0 * a +1)*ALFA + (2.0*b +1)*BETA + 1.0)
      PHI  = 1.0 - GAMA + BETA +ALFA

C     write (45,*) ALFA, BETA,GAMA,PHI
c
c*************BILDUNG DER MATRIX************************************
c
c
      LCOL=ALFA
      DIAG=1.0
      RCOL=BETA
      RSIDE1 =GAMA
      RSIDE2 =PHI        
      RSIDE3 =0.0       



C------------- Sith-order is defined---------------





c 
c
c*******************************************************************
c
c 
      return
      end

