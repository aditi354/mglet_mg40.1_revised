










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
      SUBROUTINE GRIDAE  (KK,JJ,II,KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,
     $                    DDX,DDY,DDZ,XTOT,YTOT,ZTOT)
C*STARLET***************************************************************
C        G R I D A E      AUFBAU EINES AEQUIDISTANTEN MASCHENNETZES
C                         OHNE EINBAUTEN
C*STARLET***************************************************************
C
C PARAM: KK,JJ,II       - ARRAYDIMENSIONEN
C        KMX,JMX,IMX    - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        X,Y,Z          + KOORDINATEN DER ZELLMITTELPUNKTE
C        DX,DY,DZ       + ABSTAND DER GITTERPUNKTE
C        DDX,DDY,DDZ    + ZELLABMESSUNGEN
C        XTOT,YTOT,ZTOT - GESAMTLAENGE DES BERECHNUNGSGEBIETES IN X-,
C                         Y- UND Z-RICHTUNG
C
C        05.08.85 (HW)  : ORIGINAL
C
C*STARLET***************************************************************
C
      REAL     DDX(II),     DDY(JJ),     DDZ(KK),
     $          DX(II),      DY(JJ),      DZ(KK),
     $           X(II),       Y(JJ),       Z(KK)
C
      IM1     = IMX-1
      JM1     = JMX-1
      KM1     = KMX-1
C
      DDX(1) = XTOT/FLOAT(IMX-4)
      DDY(1) = YTOT/FLOAT(JMX-4)
      DDZ(1) = ZTOT/FLOAT(KMX-4)
C
      X(1)   = -0.0*XTOT - 1.5*DDX(1)
      Y(1)   = -0.0*YTOT - 1.5*DDY(1)
      Z(1)   = - 1.5*DDZ(1)
C
C                                 BERECHNUNG DER MITTELPUNKTSKOORDINATEN
C                                 AEQUIDISTANTES GITTER !
C
      DO 10 I=1,IMX
   10    X(I) = X(1) + DDX(1)*FLOAT(I-1)
C
      DO 20 J=1,JMX
   20    Y(J) = Y(1) + DDY(1)*FLOAT(J-1)
C
      DO 30 K=1,KMX
   30    Z(K) = Z(1) + DDZ(1)*FLOAT(K-1)
C
C                                 ZELLABMESSUNGEN,
C                                 ABSTAND ZUM NAECHSTEN GITTERPUNKT
C
      DO 40 I=2,IM1
         DDX(I) = 0.5*(X(I+1)-X(I-1))
   40    DX(I)  =      X(I+1)-X(I)
         DX(1)    = DX(2)
         DDX(IMX) = DDX(IM1)
         DX(IMX)  = DX(IM1)
C
      DO 50 J=2,JM1
         DDY(J) = 0.5*(Y(J+1)-Y(J-1))
   50    DY(J)  =      Y(J+1)-Y(J)
         DY(1)    = DY(2)
         DDY(JMX) = DDY(JM1)
         DY(JMX)  = DY(JM1)
C
      DO 60 K=2,KM1
         DDZ(K) = 0.5*(Z(K+1)-Z(K-1))
   60    DZ(K)  =      Z(K+1)-Z(K)
         DZ(1)    = DZ(2)
         DDZ(KMX) = DDZ(KM1)
         DZ(KMX)  = DZ(KM1)
C
      RETURN
      END
