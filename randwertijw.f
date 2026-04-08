










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
      SUBROUTINE RANDWERTIJW(KK,K,DX,DDX,LCOL,DIAG,RCOL,
     $			      RSIDE1,RSIDE2,RSIDE3)
C*MGLET***************************************************************
C
C     UNTERPROGRAMM ZU BERECHNUNG DER KOEFFIZIENTEN AM RAND
C     HIER DIE WERTE FUER DEN AUSSTROEMRAND FALLS NO-SLIP
C     WAND VORHANDEN IST
C
C   ! ACHTUNG : DIESE FORMULIERUNG IST NUR FUER M G L E T        
C   !           (NICHTGESTAEGERT -------> GESTAEGERT)
C*MGLET***************************************************************
C
C
C DEFINE DIREKTIVEN     : KEINE
C
C UPROG                 : KEINE                                      
C                                                            
C 06.08.1996            : ADNAN MERI (ORIGINAL)
C
C*MGLET***************************************************************
C

      REAL LCOL,RCOL,DIAG,RSIDE1, RSIDE2,RSIDE3,
     $     alfa,gama,teta,phi
      Real a,b,c,d,e,f,g,h,DX(KK),DDX(KK)

      b = DX(K)/DDX(K)
      a = (DX(K-1)+ (0.5 * DX(K)))/DDX(K)

      c = (-b**2 +4.)*(b+2.) - (-b**2 +4.)*(2.-b)       
      d = 4.*(-a**2 +1.)*(b+2.) - 2.*(1.-a)*(-b**2 +4.)
      e = 4.*(b +2.) -2.*(-b**2 +4.)                                  
      f = (8. -b**3)*(b +2.) - (b**3 +8.)*(2.-b)
      g = 8.*(-a**3 +1.)*(b +2.) -2.*(1. -a)*(b**3 +8.)
      h = 8.*(b+2.) - 2.*(b**3 +8.)



      teta = (h*c -e*f)/(g*c - d*f)               
      phi  = (e - d*teta)/c                 
      gama = (2. -2.*(1. -a)*teta - (2.- b)*phi)/(b+2.)  
      alfa = -1.0 + phi + teta + gama

      LCOL = alfa 
      DIAG = 1.0
      RCOL = 0.0          
      RSIDE1 = gama 
      RSIDE2 = phi 
      RSIDE3 = teta

c
c*******************************************************************
c
      RETURN
      END



