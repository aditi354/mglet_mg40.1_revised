










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
      SUBROUTINE INTERCOEF1(II,I,DX,DDX,LCOL,DIAG,RCOL,
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

      REAL LCOL,RCOL,DIAG,RSIDE1,RSIDE2,RSIDE3,
     $     ALFA,GAMA,TETA,PHI,DX(II),
     $     l,a,b,c,d,e,f,g,p,h,r,t,q


      l = (DX(I)/DX(I-1))
      
      IF (l.ge.0.9999 .and. l.le.1.0001) THEN
 
      ALFA=(1./6.)
      GAMA=(1./6.)
      TETA=(2./3.)
      PHI =(2./3.)


      ELSE

      a = l + 1.0 
      b = - (0.5 *l) -1.0
      c = (l**2) - 1.0
      d = -(0.25 * l**2) + 1.0
      e = (l**3) + 1.0
      f = -(0.125 * l**3) -1.0 

      g = -e*d + c*f
      h = - (0.75 * e) - ((7./8.)*c)
      r = -b*c + a*d
      t = 0.5*c + 0.75*a
      q = -t*g + h*r
      p = (-g*(c+a)) +(r*(-e-c))

      TETA = p/q
      PHI  = ((c+a) - t*TETA)/r
      GAMA = (-1.0 + 0.5*TETA -b*PHI)/a
      ALFA = -1.0+PHI+TETA-GAMA

      ENDIF
c
c
c*************BILDUNG DER MATRIX************************************
c
      LCOL=ALFA
      DIAG=1.0
      RCOL=GAMA
      RSIDE1=PHI
      RSIDE2=TETA            
      RSIDE3=0.0             
c 
c*******************************************************************
c 


      RETURN
      END

