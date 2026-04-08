










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
      SUBROUTINE RANDWERTIJD (II,I,DX,DDX,LCOL,DIAG,RCOL,
     $                        RSIDE1,RSIDE2,RSIDE3)
C*MGLET***************************************************************
C     R A N D W E R T I J D
C     UNTERPROGRAMM ZU BERECHNUNG DER KOEFFIZIENTEN AM RAND
C     HIER DIE WERTE FUER DEN     A U S S T R O E M R A N D !
C     BERECHNUNG DER E R S T E N  A B L E I T U N G ! ! ! ! !
C
C      ! ACHTUNG : DIESE FORMULIERUNG IST NUR FUER M G L E T
C      !           (NICHT-GESTAEGERT -------> GESTAEGERT)
C
C*MGLET***************************************************************
C
C DEFINE DIREKTIVEN     : KEINE
C
C UPROG                 : KEINE                                      
C                                                            
C 21.10.1996            : ADNAN MERI (ORIGINAL)
C
C*MGLET***************************************************************
C

      REAL LCOL,RCOL,DIAG,RSIDE1, RSIDE2,RSIDE3,
     $     ALFA,GAMA,BETA,PHI,DX(II),DDX(II)

      REAL a,b,c,d,e,f,g,h,q,l

      l = DX(I) / DX(I-1)

      IF (l.ge.0.9999 .and. l.le.1.0001) THEN
      
      PHI  = 1.0
      GAMA = 0.0
      ALFA = - 1.0
      BETA = - 3.0

      ELSE

      a = 0.5 * DX(I)
      b = DX(I-1) + a
      c = DX(I-2) + b
 
      q = a * (2.0 * DDX(I) - a)/DX(I)
      d = b * (2.0 * DDX(I) - b)/DX(I)             
      e = c * (2.0 * DDX(I) - c)/DX(I) 
 
      f = (a**2.0) * (-3.0 * DDX(I) + 2.0 * a)/DX(I)
      g = (b**2.0) * (-3.0 * DDX(I) + 2.0 * b)/DX(I)
      h = (c**2.0) * (-3.0 * DDX(I) + 2.0 * c)/DX(I)

      PHI = (2.0 * DDX(I)*(g-f))/((h-f)*(d-q) -(e-q)*(g-f))
      GAMA= (-2.0*DDX(I) - (e-q)*PHI)/(d-q)  
      ALFA= - GAMA - PHI
      BETA= -1.0 -((c*PHI +b*GAMA + a*ALFA)/DX(I))       
C
      ENDIF
C
c*******************************************************************
c
      LCOL = BETA
      DIAG = 1.0
      RCOL = 0.0       
      RSIDE1 = ALFA 
      RSIDE2 = GAMA
      RSIDE3 = PHI  
C
c*******************************************************************
c
      RETURN
      END
