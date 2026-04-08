










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
      SUBROUTINE MGOVERLAP (II1,X1,II2,X2,NBND,XS,IA1,IE1,IA2,IE2)
C
C--MGLET----------------------------------------------------------------
C
C                    PRUEFT NACH, OB DIE ZWEI GITTER X1 UND X2
C                    UEBERLAPPEN
C
C       XS:          GITTER 1 WIRD UM XSHIFT VERSCHOBEN
C                  
C       IA1,IA2:     ANFANGSINDIZES DER UEBERLAPPUNG
C       IE1,IE2:     ENDINDIZES DER UEBERLAPPUNG
C
C                    SIND 0, FALLS KEINE UEBERLAPPUNG VORHANDEN
C
C
C        17. 2.94 (MM)  : ORIGINAL
C
C--MGLET----------------------------------------------------------------
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
      REAL     X1(II1),X2(II2)
C
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      IF ((X1(1) - X1(II1)) .GE. -SMALL ) CALL ERRR (501,'MGOVERLAP')
      IF ((X2(1) - X2(II2)) .GE. -SMALL ) CALL ERRR (502,'MGOVERLAP')
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                     SMALL IST FUER FOLGENDE VERGLEICHE NOCH ZU KLEIN

C       SQSM = 10.0*SQRT(ABS(SMALL))
       SQSM = 0.001
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      IER = 0

      IA1 = 0
      IE1 = 0
      IA2 = 0
      IE2 = 0
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C                    ANFANGSPUNKT VON GITTER 1 LIEGT INNERHALB 
C                    VON GITTER 2

      IF ((X1(1)+XS .GE. X2(1)) .AND. (X1(1)+XS .LE. X2(II2))) THEN

C                            SUCHE NACH ANFANGSPUNKT IA2
         DO I1 = 1,2
         DO I2 = 1,II2
            IF ( ABS(X2(I2)-(X1(I1)+XS)) .LE. SQSM ) THEN
               IA2 = I2
               IA1 = I1
               GOTO 100
            ENDIF
         ENDDO
         ENDDO

            IER = IER + 1
 
  100     CONTINUE

C                    ENDPUNKT VON GITTER 1 LIEGT INNERHALB
C                    VON GITTER 2

         IF (X1(II1)+XS .LE. X2(II2)) THEN

C                            SUCHE NACH ENDPUNKT IE2
        DO I1 = II1,II1-1,-1
        DO I2 = 1,II2
            IF ( ABS(X1(I1)+XS-X2(I2)) .LE. SQSM ) THEN
               IE2 = I2
               IE1 = I1
               GOTO 200
            ENDIF
         ENDDO
         ENDDO

            IER = IER + 1
 
  200     CONTINUE

         ENDIF
C
C                    ENDPUNKT VON GITTER 2 LIEGT INNERHALB
C                    VON GITTER 1


         IF (X1(II1)+XS .GT. X2(II2)) THEN

C                            SUCHE NACH ENDPUNKT IE1
         DO I2 = II2,II2-1,-1
         DO I1 = 1,II1
            IF ( ABS(X1(I1)+XS-X2(I2)) .LE. SQSM ) THEN
               IE2 = I2
               IE1 = I1
               GOTO 300
            ENDIF
         ENDDO
         ENDDO

            IER = IER + 1
 
  300     CONTINUE

         ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C                    ANFANGSPUNKT VON GITTER 2 LIEGT INNERHALB 
C                    VON GITTER 1

      ELSEIF ((X1(1)+XS .LE. X2(1)) .AND. (X1(II1)+XS .GE. X2(1))) THEN

C                            SUCHE NACH ANFANGSPUNKT IA1
         DO I2 = 1,2
         DO I1 = 1,II1
            IF ( ABS(X1(I1)+XS-X2(I2)) .LE. SQSM ) THEN
               IA2 = I2
               IA1 = I1
               GOTO 1100
            ENDIF
         ENDDO
         ENDDO

            IER = IER + 1
 
 1100    CONTINUE

C                    ENDPUNKT VON GITTER 1 LIEGT INNERHALB
C                    VON GITTER 2

         IF (X1(II1)+XS .LE. X2(II2)) THEN

C                            SUCHE NACH ENDPUNKT IE2
         DO I1 = II1,II1-1,-1
         DO I2 = 1,II2
            IF ( ABS(X1(I1)+XS-X2(I2)) .LE. SQSM ) THEN
               IE2 = I2
               IE1 = I1
               GOTO 1200
            ENDIF
         ENDDO
         ENDDO

            IER = IER + 1
 
 1200     CONTINUE

         ENDIF
C
C                    ENDPUNKT VON GITTER 2 LIEGT INNERHALB
C                    VON GITTER 1


         IF (X1(II1)+XS .GT. X2(II2)) THEN

C                            SUCHE NACH ENDPUNKT IE1
         DO I2 = II2,II2-1,-1
         DO I1 = 1,II1
            IF ( ABS(X1(I1)+XS-X2(I2)) .LE. SQSM ) THEN
               IE2 = I2
               IE1 = I1
               GOTO 1300
            ENDIF
         ENDDO
         ENDDO

            IER = IER + 1
 
 1300     CONTINUE

         ENDIF

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      ELSE
C                    GITTER UEBERLAPPEN NICHT
         RETURN

      ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C                    CHECK,  OB GITTERPUNKTE UEBEREINSTIMMEN

      DO I = 0,IE1-IA1
         I1 = I+IA1
         I2 = I+IA2

         IF ( ABS(X1(I1)+XS-X2(I2)) .GT. SQSM ) THEN

            IER = IER + 1
            IA1 = 0
            IE1 = 0
            IA2 = 0
            IE2 = 0

         ENDIF

      ENDDO
         
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C

      RETURN
      END


