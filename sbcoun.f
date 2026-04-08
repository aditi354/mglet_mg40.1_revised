










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
      SUBROUTINE SBCOUN (KK,JJ,II,KMX,JMX,IMX,
     $                  NCOUN,XMIT,HEIGHT,ALPHA,
     $                  B,HILF,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                  XTOT,YTOT,ZTOT,BPAR)
C*STARLET***************************************************************
C        S B C O U N    BELEGUNG DES WANDERKENNUNGSFELDES  
C                       FUER COUNIHAN-WIRBELERZEUGER
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
C        NCOUN          -ANZAHL AUF KANALBREITE
C        HEIGHT         -HOEHE DER ELEMENTE == GRENZSCHICHTDICKE
C        XMIT           - X-POSITION DER HINTERKANTE
C        ALPHA          - SPITZENWINKEL DER ELEMENTE IN GRAD
C
C VERS:  12. 5.92 (MM)  - ORIGINAL
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


CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C             ERMITTELN DER Y-POSITIONEN DER MITTELPUNKTE

        DELY = YTOT/FLOAT(NCOUN)
        write(6,*) 'Geschatzter lateraler Abstand der Flossen:',DELY
C                                  Anpassung des lateralen Abstands zu 
C                                  den Gitterabstaenden in der Mitte
        DIFFERENCE = ABS(DELY-DDY((JMX+1)/2))
        JMULT = 1
        DO JSTRIDE = 2,(JMX+1)/2

           DELYGRID = DDY((JMX+1)/2)*FLOAT(JSTRIDE)

           IF ( ABS(DELY-DELYGRID) .LT. DIFFERENCE ) THEN
              JMULT = JSTRIDE
              DIFFERENCE = ABS(DELY-DELYGRID)
           ENDIF

        ENDDO
        DELY = DDY((JMX+1)/2)*FLOAT(JMULT)
        write(6,*) 'Gewaehlter  lateraler Abstand der Flossen:',DELY

C                                  ERMITTLUNG DER POSITION DER ERSTEN 
C                                  FLOSSE
        IF (NCOUN .EQ. (NCOUN/2)*2 ) THEN
C                                  GERADE ANZAHL
          YMIT = (Y(JMX-2)+Y(3)) / 2.0 - FLOAT(NCOUN/2)*DELY 
CTEST     $          + DDY((JMX+1)/2)*0.5
     $          - DDY((JMX+1)/2)*0.5 + DELY*0.5

        ELSE
C                                  UNGERADE ANZAHL
          YMIT = (Y(JMX-2)+Y(3)) / 2.0 - FLOAT((NCOUN-1)/2)*DELY 

        ENDIF

       DO 50 I=1,NCOUN

          write(6,*) 'i=',i,'       ymit=',ymit
          CALL COUSIN (KK,JJ,II,KMX,JMX,IMX,
     $                 ZMIT,YMIT,XMIT,HEIGHT,ALPHA,
     $                  B,HILF,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,BPAR)

          YMIT = YMIT + DELY

   50  CONTINUE
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
        RETURN
        END
