










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
      SUBROUTINE CALDPDX (II,GRADP,X,NPPHYS,PPHYS,XPPHYS)
C---------------------------------------------------------------------72
C
C     CALCULATES PRESSURE GRADIENTS ON GRID-POINTS (STAGGERED) FROM
C     GIVEN PRESSURE 
C
C      NPPHYS : NUMBER OF GIVEN PRESSURE POINTS
C       PPHYS : GIVEN PRESSURE
C      XPPHYS : POSITIONS OF GIVEN PRESSURE
C
C          II : NUMBER OF GRID POINTS
C      GRADP  : COMPUTED PRESSURE GRADIENT
C           X : POSITIONS OF CELL-CENTERS
C
C      10.02.98 (M.M.) ORIGINAL
C---------------------------------------------------------------------72

      REAL GRADP(II),X(II)

      REAL PPHYS(NPPHYS),XPPHYS(NPPHYS)

C---------------------------------------------------------------------72
C                                FALLS KEINE VERNUENFTIGE EINGABE,
C                                ODER BEI ALTEN STEUERFILES

      IF ( NPPHYS .LT. 2 ) THEN
        DO I = 1,II
          GRADP(I) = 0.0
        ENDDO
      ENDIF

C---------------------------------------------------------------------72

C                                 START- AND STOP-INDICES 
      IX_START = 1
      DO I=2,NPPHYS
         IF ( X(1)   .LE. XPPHYS(I) ) THEN
            IX_START = I-1
            GOTO 30
         ENDIF
      ENDDO
   30 CONTINUE

      IX_STOP = NPPHYS
      DO I=NPPHYS-1,1,-1
         IF ( X(II) .GE. XPPHYS(I) ) THEN
            IX_STOP  = I
            GOTO 31
         ENDIF
      ENDDO
   31 CONTINUE

C---------------------------------------------------------------------72

C                     CORRECTION OF IX_STOP

        IX_STOP = MAX(IX_STOP,IX_START+1)
        IX_STOP = MIN(IX_STOP,NPPHYS)


C---------------------------------------------------------------------72

C                               INTERPOLATION


      DO I=1,II

C                               PRESSURE GRADIENT AT FIRST NODE

         IF ( X(I) .LE. XPPHYS(1     )) THEN
            GRADP(I) = 0.0
         ELSEIF ( X(1) .GE. XPPHYS(NPPHYS)) THEN
            GRADP(I) = 0.0
         ELSE
            
            DO IX=IX_START,IX_STOP
               IF ((X(I)-XPPHYS(IX))*(X(I)-XPPHYS(IX+1)) 
     $              .LE. 0.0) THEN
                  GRADP(I) = 
     $                 ( PPHYS(IX+1) -  PPHYS(IX)) / 
     $                 (XPPHYS(IX+1) - XPPHYS(IX))
               ENDIF
            ENDDO
         ENDIF
      ENDDO

C---------------------------------------------------------------------72

       RETURN
       END
