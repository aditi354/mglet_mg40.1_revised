










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
      SUBROUTINE RANDWERTDW (KK,K,DX,DDX,STAG,LCOL,DIAG,RCOL,RSIDE1,
     $                        RSIDE2,RSIDE3)
C*MGLET***************************************************************
C     R A N D W E R T D W
C     UNTERPROGRAMM ZU BERECHNUNG DER KOEFFIZIENTEN AM RAND
C     HIER DIE WERTE FUER DEN     A U S S T R O E M R A N D FALLS
C     DORT EINE WAND IST
C     BERECHNUNG DER E R S T E N  A B L E I T U N G ! ! ! ! !
C
C*MGLET***************************************************************
C
C
C DEFINE DIREKTIVEN     : KEINE
C
C UPROG                 : KEINE                                      
C                                                            
C 21.10.1996            : ADNAN MERI (ORIGINAL)
C 21.01.2003            : FLORIAN SCHWERTFIRM
C
C*MGLET***************************************************************
C
      LOGICAL STAG
      REAL LCOL,RCOL,DIAG,RSIDE1, RSIDE2,RSIDE3,
     $     alpha,beta,gamma,phi,h1,h2,h3,h4
       

      REAL a,b,c,d,e,f,g, h,q,l,DX(KK),DDX(KK),
     $     a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4 
     
C

CCCCCC Falls Rand gleich Noslip bestimmung mit cramers. regel
C  Fi+alpha*f(i+1/2)[h1]+beta*f(i-1/2)[h2]+gamma*F(i-1)[h3]+phi*f(i-3/2)[h4] = O**3
C  F unbekannt f bekannt. Taylorreihe -> Koeffizientenvergleich. 
C  Gleichungssystem mit Koeff. als Unbekannte. Loesung mit Cramerscher 
C  Regel


      IF (STAG) THEN
         h1 = DX(K)/2.0
         h2 = DX(K-1)/2.0
         h3 = DX(K-1)
         h4 = DX(K-1)+DX(K-2)/2.0
      ELSE
         h1 =  DX(K-1)/2.0
         h2 =  DX(K-1)/2.0
         h3 = (DX(K-1)+DX(K-2))/2.0
         h4 =  DX(K-1)/2.0 + DX(K-2)
      ENDIF

        a1 = h1 
	a2 = -h2 
	a3 = 1 
	a4 = -h4
        b1 = 0.5*h1**2 
        b2= 0.5*h2**2 
	b3 = -h3 
	b4 = 0.5*h4**2
        c1 = h1**3/6.0 
	c2 = -h2**3/6.0 
	c3 = 0.5*h3**2 
	c4 = -h4**3/6.0


      a = a1*b3*c2-a1*b2*c3+a2*b1*c3-a2*c1*b3-b1*a3*c2+a3*b2*c1 -
     $ a1*b3*c4+a1*b4*c3+b1*a3*c4-b1*a4*c3-a3*c1*b4+c1*a4*b3 +
     $ a2*b3*c4-a2*b4*c3-a3*b2*c4+a3*c2*b4+b2*a4*c3-a4*b3*c2
      b = b2*c3 - b3*c2 + b3*c4 - b4*c3
      c = c1*b3 - b1*c3 - b3*c4 + b4*c3
      d = b1*c2 - b2*c1 - b1*c4 + c1*b4 + b2*c4 - c2*b4
      e = b1*c3 - c1*b3 - b2*c3 + b3*c2
    
      alpha = -b/a
      beta = -c/a
      gamma = d/a
      phi  = -e/a
c
c*******************************************************************
c
      LCOL = gamma        
      DIAG = 1.0
      RCOL = 0.0  
      RSIDE1 = alpha 
      RSIDE2 = beta
      RSIDE3 = phi 
c
c*******************************************************************
c


      RETURN
      END
