










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
      SUBROUTINE SBCUB (KK,JJ,II,KMX,JMX,IMX,
     $                  IB1,IB2,JB1,JB2,KB1,KB2,
     $                  XB1,XB2,YB1,YB2,ZB1,ZB2,
     $                  B,HILF,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,BPAR)
C*MGLET*****************************************************************
C        S B C U B      VORBELEGUNG DES WANDERKENNUNGSFELDES       
C                       FUER QUADER
C*MGLET*****************************************************************
C
C PARAM: B              + WANDERKENNUNGSFELD FUER KOERPER IM BERECHNUNGS
C                       + GEBIET:
C
C                             B>0.0  ==> IM BERECHNUNGSGEBIET
C                             B=-1.0 ==> NOSLIP-WAND
C                             B=-2.0 ==>   SLIP-WAND
C
C        JJ,II          - ARRAYDIMENSIONEN
C        JMX,IMX        - GRENZEN DES BERECHNUNGSGEBIETS (MIT BOUND)
C        KB1,KB2        - CUBUS-GRENZE  (INDEX)
C        JB1,JB2        - CUBUS-GRENZE  (INDEX)
C        IB1,IB2        - CUBUS-GRENZE  (INDEX)
C        XB1,XB2        - CUBUS-GRENZEN IN X-RICHTUNG (REAL)
C        YB1,YB2        - CUBUS-GRENZEN IN Y-RICHTUNG (REAL)
C        ZB1,ZB2        - CUBUS-GRENZEN IN Z-RICHTUNG (REAL)
C
C VERS:  12.03.92 (MM)  - STEUERUNG DER BERECHNUG LAEFT ZUNAECHST
C                         MIT IB1,IB2,JB1,JB2,KB
C        23. 4.92 (MM)  - AB JETZT KB1 UND KB2
C        15.12.93 (MM)  - STEUERUNG MIT GRENZEN IN REALER GEOMETRIE 
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
      IF (IB1.GT.0) THEN
C 
C                          STEUERUNG UEBER INDEX-GRENZEN
C
      DO 91 I=IB1,IB2
      DO 91 J=JB1,JB2
      DO 91 K=KB1,KB2

                     B(K,J,I) = BPAR

   91   CONTINUE
C
C
      ELSE
C
C                         STEUERUNG MIT REALER GEOMETRIE
C
C                         ABSTAND VON OBERFLAECHE DES QUADERS
C                         WIRD NEGATIV IM KOERPER
C
        DO I=2,IM1
        DO J=2,JM1
        DO K=2,KM1

             HILF(K,J,I) = MAX ( XB1 - X(I) , X(I) - XB2 ,
     $                           YB1 - Y(J) , Y(J) - YB2 ,
     $                           ZB1 - Z(K) , Z(K) - ZB2 ) - SMALL

         ENDDO
         ENDDO
         ENDDO
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C                             HINZUFUEGEN DES NEUEN KOERPERS ZU 'B'
C
       DO 145 I=2,IM1
       DO 145 J=2,JM1
       DO 145 K=2,KM1

               B(K,J,I) = -BPAR    *(SIGN(0.5,HILF(K,J,I)) - 0.5)
     $                    +B(K,J,I)*(SIGN(0.5,HILF(K,J,I)) + 0.5)

  145   CONTINUE


      ENDIF



        RETURN
        END
