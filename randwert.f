










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
      SUBROUTINE RANDWERT (KK,K,DX,STAG,LCOL,DIAG,RCOL,RSIDE1,RSIDE2,
     $                     RSIDE3)
C*MGLET***************************************************************
C
C     UNTERPROGRAMM ZU BERECHNUNG DER KOEFFIZIENTEN AM RAND
C     HIER DIE WERTE FUER DEN AUSSTROEMRAND !
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
      IMPLICIT NONE
      LOGICAL STAG
      INTEGER K,KK
      REAL LCOL,RCOL,DIAG,RSIDE1, RSIDE2,RSIDE3,
     $     alpha,beta,gamma,phi
      Real a1,a2,a3,a4,b1,b2,b3,b4,
     $     c1,c2,c3,c4,d1,d2,d3,d4,
     $     h1,h2,h3,h4,
     $     a,b,c,d,e,DX(KK)

C Fi + alpha*f(i-0.5)[h1]+beta*F(i-1)[h2]+gamma*f(i-3/2)[h3]
C    + phi*f(i-5/2)[h4] = O**4
      IF(STAG) THEN
         h1 = DX(K-1)/2.0
         h2 = DX(K-1)
         h3 = DX(K-1)+DX(K-2)/2.0
         h4 = DX(K-1)+DX(K-2)+DX(K-3)/2.0
      ELSE
         h1 = DX(K-1)/2.0
         h2 = DX(K-1)/2.0 + DX(K-2)/2.0
         h3 = DX(K-1)/2.0 + DX(K-2)
         h4 = DX(K-1)/2.0 + DX(K-2) + DX(K-3)
      ENDIF

      a1 = 1.0
      a2 = 1.0 
      a3 = 1.0 
      a4 = 1.0
      b1 = -h1 
      b2 = -h2
      b3 = -h3
      b4 = -h4
      c1 = 0.5*h1**2
      c2 = 0.5*h2**2
      c3 = 0.5*h3**2
      c4 = 0.5*h4**2
      d1 = -1.0/6.0*h1**3 
      d2 = -1.0/6.0*h2**3
      d3 = -1.0/6.0*h3**3 
      d4 = -1.0/6.0*h4**3

      a= a1*b2*c3*d4-a1*b2*c4*d3-a1*b3*c2*d4+a1*b3*d2*c4+ a1*c2*b4*d3 -
     $   a1*b4*c3*d2-a2*b1*c3*d4+ a2*b1*c4*d3+a2*c1*b3*d4-a2*c1*b4*d3 -
     $   a2*b3*d1*c4+a2*d1*b4*c3+b1*a3*c2*d4-b1*a3*d2*c4-b1*a4*c2*d3 +
     $   b1*a4*c3*d2-a3*b2*c1*d4+a3*b2*d1*c4+a3*c1*b4*d2-a3*c2*d1*b4 +
     $   b2*c1*a4*d3-b2*a4*d1*c3-c1*a4*b3*d2+a4*b3*c2*d1

      b= b2*c4*d3 - b2*c3*d4 + b3*c2*d4 - b3*d2*c4 - c2*b4*d3 + b4*c3*d2
      c= b1*c3*d4 - b1*c4*d3 - c1*b3*d4 + c1*b4*d3 + b3*d1*c4 - d1*b4*c3
      d= b1*d2*c4 - b1*c2*d4 + b2*c1*d4 - b2*d1*c4 - c1*b4*d2 + c2*d1*b4
      e= b1*c2*d3 - b1*c3*d2 - b2*c1*d3 + b2*d1*c3 + c1*b3*d2 - b3*c2*d1
      
      alpha = -b/a
      beta = c/a
      gamma = -d/a
      phi = -e/a

C      a = DX(K-1)+ (0.5 * DX(K-2))
C      b = DX(K-1) + DX(K-2) + (0.5*DX(K-3))
C      c = DX(K-1)              
C      d = (6.0 * c *(-c+a)) + (4.0 * (c**2 - a**2))
C      e = (6.0 * c *(-c+b)) + (4.0 * (c**2 - b**2))
C      f = (7.0 * c *(c**2 - a**2)) + (6.0 * (-c**3 + a**3))
C      g = (7.0 * c *(c**2 - b**2)) + (6.0 * (-c**3 + b**3))
C      h = (f*(-2.0 * c**2)) - (d*(c**3))
C
C
C
C      teta = h / (e*f - d*g)                       
C      phi  = ((-2.0 * c**2)- (e*teta))/d
C      gama = (c + ((b-c)*teta) + ((a-c)*phi))/(0.5*c)
C      alfa = -1.0 + phi + teta + gama
c
c*******************************************************************
c
      LCOL = beta 
      DIAG = 1.0
      RCOL = 0.0
      RSIDE1 = alpha
      RSIDE2 = gamma
      RSIDE3 = phi

c
c*******************************************************************
c
      RETURN
      END



