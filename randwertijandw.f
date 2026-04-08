










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
      SUBROUTINE RANDWERTIJANDW(II,K,DX,DDX,LCOL,DIAG,RCOL,RSIDE1,
     $                         RSIDE2,RSIDE3)
C*MGLET***************************************************************
C     R A N D W E R T I J D W
C     UNTERPROGRAMM ZU BERECHNUNG DER KOEFFIZIENTEN AM RAND
C     FALLS
C     EINE NO-SLIP-WAND VORHANDEN IST
C     BERECHNUNG DER E R S T E N  A B L E I T U N G ! ! ! ! !
C
C      ! ACHTUNG : DIESE FORMULIERUNG IST NUR FUER M G L E T
C      !           (NICHT-GESTAEGERT -------> GESTAEGERT)
C
C*MGLET***************************************************************
C
C
C DEFINE DIREKTIVEN     : KEINE
C
C UPROG                 : KEINE                                      
C                                                            
C 21.10.1996            : ADNAN MERI (ORIGINAL)
C 15.01.2003            : FLORIAN SCHWERTFIRM
C*MGLET***************************************************************
C

      REAL LCOL,RCOL,DIAG,RSIDE1, RSIDE2,RSIDE3,
     $     alpha,beta,gamma,phi

      REAL h1,h2,h3,a0,a,b,c,d,e,DX(II),DDX(II),
     $     a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4

C
C F(i) + alpha*f(i)+beta*f(i+1/2)[h1]+gamma*F(i+1)[h2]+phi*f(i+3/2)[h3] = O(3)
C F = Ableitung unbekannt f = Funktionswert bekannt
C Beachte: f(i) in obiger Gleichung ist die Interpolierte Geschwindigkeit an
C der Wand

       h1 = DX(K-1)/2.0
       h2 = DX(K-1)/2.0+DX(K)/2.0
       h3 = DX(K-1)/2.0+DX(K)
       
       a1 = 0.0
       a2 = h1
       a3 = 1.0
       a4 = h3
       b1 = 0.0
       b2 = 0.5*h1**2
       b3 = h2
       b4 = 0.5*h3**2
       c1 = 0
       c2 = h1**3/6.0
       c3 = 0.5*h2**2
       c4 = h3**3/6

       a = (a2*(b3*c4-c3*b4)-a3*(b2*c4-c2*b4)+a4*(b2*c3-c2*b3))
     $    -(a1*(b3*c4-c3*b4)-a3*(b1*c4-c1*b4)+a4*(b1*c3-c1*b3))
     $    -(a1*(b2*c3-c2*b3)-a2*(b1*c3-c1*b3)+a3*(b1*c2-c1*b2))
       b = (b3*c4-c3*b4)+(b2*c3-c2*b3)
       c = -(b3*c4-c3*b4)-(b1*c3-c1*b3)
       d = (b2*c4-c2*b4)-(b1*c4-c1*b4)+(b1*c2-c1*b2)
       e = -(b2*c3-c2*b3)+(b1*c3-c1*b3)

      alpha = -b/a
      beta = -c/a
      gamma = d/a
      phi  = -e/a
c*******************************************************************
c
      LCOL = 0.0 
      DIAG = 1.0
      RCOL = gamma     
      RSIDE1 = alpha 
      RSIDE2 = beta
      RSIDE3 = phi 
C
c*******************************************************************
c




      RETURN
      END
