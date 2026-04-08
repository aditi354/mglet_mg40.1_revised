










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
      SUBROUTINE COUSIN (KK,JJ,II,KMX,JMX,IMX,
     $                 ZMIT,YMIT,XMIT,DELTA,ALFA,
     $                  B,HILF,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,BPAR)
C*STARLET***************************************************************
C        C O U S I N    BELEGUNG DES WANDERKENNUNGSFELDES  
C                       FUER EINE COUNINGHAM_FLOSSE
C*STARLET***************************************************************
C
C PARAM: B              + WANDERKENNUNGSFELD FUER KOERPER IM BERECHNUNGS
C                       + GEBIET:
C
C                             B>0.0  ==> IM BERECHNUNGSGEBIET
C                             B=-1.0 ==> NOSLIP-WAND
C                             B=-2.0 ==>   SLIP-WAND
C
C        HILF           - HILFSFELD
C        KK,JJ,II       - ARRAYDIMENSIONEN
C        KMX,JMX,IMX    - GRENZEN DES BERECHNUNGSGEBIETS (MIT BOUND)
C
C        ZMIT           -Z-POSITION DES ELLYPSENMITTELPUNKTES
C        YMIT           -Y-POSITION DES ELLYPSENMITTELPUNKTES
C        XMIT           -X-POSITION DES ELLYPSENMITTELPUNKTES
C        DELTA          - HOEHE DER FLOSSE == GRENZSCHICHTDICKE
C        ALFA           - SPITZENWINKEL
C
C VERS:  12. 5.92 (MM)  - ORIGINAL
C                         _UPWARDS_ EINGEFUEHRT, STROEMUNG NACH OBEN
C
C UPROG                 : ERRR
C
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

      COMMON /CONLES/  CONV2S,CMUE,CAPPA,ECONST
      SAVE   /CONLES/

      REAL        X(II),         Y(JJ),         Z(KK),
     $           DX(II),        DY(JJ),        DZ(KK),
     $          DDX(II),       DDY(JJ),       DDZ(KK),
     $            B(KK,JJ,II),HILF(KK,JJ,II)

      REAL    RRCUB(6)
      INTEGER IICUB(6)

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                     HALBE BREITE DER FLOSSE IN ABHAENGIGKEIT 
C                     DES SPITZENWINKELS UND DER X- Z- POSITION


       WEDGE(WX,WZ) =
     $  TANALF *  
     $  (D2*SQRT(MAX(0.0,1.0-((WZ-Z0)/D1)**2))-X0+WX)

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C             ERMITTELN DER GEOMETRIEBESTIMMENDEN GROESSEN

CSTUFE       Z0   = 1.0
       Z0   = Z(3) - 0.5*DDZ(3)
       X0   = XMIT
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                     HALBACHSEN DER ELLYPSE

       D1 = DELTA
       D2 = 0.5*DELTA
       TANALF = TAN(ALFA*ATAN(1.)/45.)

C      WRITE (6,*) 'COUSIN, X0, Z0, D1, D2:',X0,Z0,D1,D2
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                      BELEGEN DER SCHLEIFENGRENZEN

        IM2  = IMX-2
        JM2  = JMX-2
        KM2  = KMX-2
         
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C

       DO 97 I=3,IM2
       DO 97 J=3,JM2
       DO 97 K=3,KM2

           WIDTH = WEDGE(X(I),Z(K))
C      WRITE (6,*) 'COUSIN, Width:',width

C
C                            ELLYPSE IN X-Z-EBENE
C
           HILF(K,J,I) = 
     $   (Z(K)-Z0)**2/D1**2
     $  +(X(I)-X0)**2/D2**2
     $  -1.0
C
C                            HINTEN ABGESCHNITTEN
C
     $   + SIGN(1.0,(X(I)-X0))+1.0
C
C                            SPITZENWINKEL RECHTS
C
     $   + SIGN(1.0,(Y(J)-(YMIT+WIDTH)))+1.0
C
C                            SPITZENWINKEL LINKS
C
     $   + SIGN(1.0,((YMIT-WIDTH)-Y(J)))+1.0

CCCCC      WRITE (6,*)'COUSIN, B:',i,j,k,hilf(k,j,i)
   97 CONTINUE


CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C                             HINZUFUEGEN DES NEUEN KOERPERS ZU 'B'
C
       DO 145 I=3,IM2
       DO 145 J=3,JM2
       DO 145 K=3,KM2

               B(K,J,I) = -BPAR    *(SIGN(0.5,HILF(K,J,I)) - 0.5)
     $                    +B(K,J,I)*(SIGN(0.5,HILF(K,J,I)) + 0.5)

  145 CONTINUE

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
        RETURN
        END
