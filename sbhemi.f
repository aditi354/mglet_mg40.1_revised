










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
      SUBROUTINE SBHEMI (KK,JJ,II,KMX,JMX,IMX,
     $                  XB1,XB2,YB1,YB2,ZB1,ZB2,
     $                  B,HILF,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,BPAR)
C*STARLET***************************************************************
C        S B H E M I    BELEGUNG DES WANDERKENNUNGSFELDES  
C                       FUER KUGEL
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
C        XB1,XB2        - AUSDEHNUNG IN X-RICHTUNG
C        YB1,YB2        - AUSDEHNUNG IN Y-RICHTUNG
C        ZB1,ZB2        - AUSDEHNUNG IN Z-RICHTUNG
C
C VERS:  20. 7.92 (MM)  - ORIGINAL
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


C
        IM1  = IMX-1
        IM2  = IMX-2
        IM3  = IMX-3
        JM1  = JMX-1
        JM2  = JMX-2
        JM3  = JMX-3
        KM1  = KMX-1
        KM2  = KMX-2
        KM3  = KMX-3
C
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C

       DO 72 I=3,IM2
       DO 72 J=3,JM2
       DO 72 K=3,KM2

           HILF(K,J,I) = 
C
C
     $   (X(I)-XB1)**2/XB2**2
     $  +(Y(J)-YB1)**2/YB2**2
     $  +(Z(K)-ZB1)**2/ZB2**2
     $  -1.0
C

   72   CONTINUE


CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C                             HINZUFUEGEN DES NEUEN KOERPERS ZU 'B'
C
       DO 145 I=3,IM2
       DO 145 J=3,JM2
       DO 145 K=3,KM2

               B(K,J,I) = -BPAR    *(SIGN(0.5,HILF(K,J,I)) - 0.5)
     $                    +B(K,J,I)*(SIGN(0.5,HILF(K,J,I)) + 0.5)

  145   CONTINUE
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
        RETURN
        END
