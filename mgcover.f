










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
      SUBROUTINE MGCOVER (II1,X1,DX1,II2,X2,DX2,NBND,XS,IPOS1,IPOS2)
C
C--MGLET----------------------------------------------------------------
C
C                    PRUEFT NACH, OB GITTER 2 EIN PARENT SEIN KANN
C                    VON GITTER 1
C
C       XS:          GITTER 1 WIRD UM XSHIFT VERSCHOBEN
C                  
C       IPOS1:       POSITION DER 3. ZELLE DES ERSTEN GITTERS IM ZWEITEN GITTER
C       IPOS2:       POSITION DER IMX-2-TEN ZELLE 
C
C                    SIND 0, FALLS KEINE UEBERLAPPUNG VORHANDEN
C
C
C        17. 3.97 (MM)  : ORIGINAL, AUS MGOVERLAP ABGELEITET
C
C--MGLET----------------------------------------------------------------
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
      REAL     X1(II1),X2(II2),DX1(II1),DX2(II2)
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C

      IF ((X1(1) - X1(II1)) .GE. -SMALL ) CALL ERRR (501,'MGCOVER')
      IF ((X2(1) - X2(II2)) .GE. -SMALL ) CALL ERRR (502,'MGCOVER')

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                     SMALL IST FUER FOLGENDE VERGLEICHE NOCH ZU KLEIN

       SQSM = SQRT(ABS(SMALL))
       SQSM = 0.001
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      IER = 0

      IPOS1 = 0
      IPOS2 = 0
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C                    ANFANGSPUNKT VON GITTER 1 LIEGT INNERHALB 
C                    VON GITTER 2

      IF ((X1(1)+XS .GE. X2(1)) .AND. (X1(1)+XS .LE. X2(II2))) THEN

C                            SUCHE NACH ANFANGSPUNKT IA2
         I1 = 3
         XLEFT = X1(I1)-DX1(I1-1)/2.0 + XS
         DO I2 = 2,II2
            IF ( ABS((X2(I2)-DX2(I2-1)/2.0)-XLEFT) .LE. SQSM ) THEN
               IPOS1 = I2
               GOTO 100
            ENDIF
         ENDDO

         IER = IER + 1
 
  100     CONTINUE

      ENDIF
C                    ENDPUNKT VON GITTER 1 LIEGT INNERHALB
C                    VON GITTER 2

      IF (X1(II1-2)+XS .LE. X2(II2)) THEN

C                            SUCHE NACH ENDPUNKT IPOS2
         I1 = II1-2
         XRIGHT = X1(I1) + DX1(I1)/2.0 + XS
         DO I2 = 1,II2
            IF ( ABS(XRIGHT-(X2(I2)+DX2(I2)/2.0)) .LE. SQSM ) THEN
               IPOS2 = I2
               GOTO 200
            ENDIF
         ENDDO

         IER = IER + 1
         
  200     CONTINUE

      ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C

      RETURN
      END


