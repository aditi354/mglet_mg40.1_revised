










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
      SUBROUTINE RANDWERTIJAN(KK,K,DX,DDX,LCOL,DIAG,RCOL,
     $			    RSIDE1,RSIDE2,RSIDE3,NFRO)
C*MGLET***************************************************************
C
C     UNTERPROGRAMM ZU BERECHNUNG DER KOEFFIZIENTEN AM RAND
C     HIER DIE WERTE FUER DEN EINSTROEMRAND !
C     MODEFIZIERTER ANSATZ !!!!!!!!!!!!!!!!!!
C
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
      REAL a,b,c,d,e,f,g,h,q,r,DX(KK),DDX(KK)

C     IF (NFRO .NE.6) THEN
C                                       EINSTROEMRAND
      a = 0.5 * DX(K)
      b = DX(K+1)+ (0.5 * DX(K))
      c = DX(K+2) + DX(K+1) + (0.5*DX(K))
      d = DDX(K)                                     
C
      e = ((b-d)*(d+a)) + ((d-b)*(d+b))                  
      f = ((c-d)*(d+a)) + ((d-c)*(d+c))                  
      g = (b-d)*((a**3) - (d**3)) - (a-d)*((b**3) - (d**3)) 

      h = (c-d)*((a**3) - (d**3)) - (a-d)*((c**3) - (d**3))             
      r = -d*((a**3) - (d**3)) + ((d**3) * (a-d))  
      q = -d*a               
C
      teta = (r*e - q*g)/(h*e - g*f)               
      phi  = (q - f*teta)/e                 
      gama = (-d - (c-d)*teta -(b-d)*phi) / (a-d)         
      alfa = -1.0 + phi + teta + gama
C
C     ELSE 
C                                       SLIP-WAND
C     b = DX(K-1)/DDX(K)
C     a = (DX(K)+ (0.5 * DX(K-1)))/DDX(K)
C
C     c = (-b**2 +4.)*(b+2.) - (-b**2 +4.)*(2.-b)       
C     d = 4.*(-a**2 +1.)*(b+2.) - 2.*(1.-a)*(-b**2 +4.)
C     e = 4.*(b +2.) -2.*(-b**2 +4.)                                  
C     f = (8. -b**3)*(b +2.) - (b**3 +8.)*(2.-b)
C     g = 8.*(-a**3 +1.)*(b +2.) -2.*(1. -a)*(b**3 +8.)
C     h = 8.*(b+2.) - 2.*(b**3 +8.)
C
C     teta = (h*c -e*f)/(g*c - d*f)               
C     phi  = (e - d*teta)/c                 
C     gama = (2. -2.*(1. -a)*teta - (2.- b)*phi)/(b+2.)  
C     alfa = -1.0 + phi + teta + gama
c
C     ENDIF
c
c*******************************************************************
c
      LCOL = 0.0
      DIAG = 1.0
      RCOL = alfa
      RSIDE1 = gama 
      RSIDE2 = phi 
      RSIDE3 = teta
c
c*******************************************************************
c
      RETURN
      END



