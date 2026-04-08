










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
      SUBROUTINE BOREC(KK,JJ,II,NBND,WBO,W,HILF,XINNER,XOUTER,
     $     X,XWAND,DELTA,DELTA_ACT,UTAU,UTAU_ACT,UFRCON,UMAX_ACT,
     $     GMOL)
C--MGLET----------------------------------------------------------
C
C     "RECYCLING" BOUNDARY CONDITION (LUND, SPALART) FOR BOTTOM
C     
C
C        UFR(X0)=U_I^INNER + U_I^OUTER
C
C    30.03.98 (MM.):    ORIGINAL AUS SVREC ABGELEITET
C
C-------10--------20--------30--------40--------50--------60--------7072

      REAL
     $     WBO(JJ,II, 2),W(JJ,II),HILF(JJ,II),
     $     X(II),XINNER(II),XOUTER(II)

C---------------------------------------------------------------------72
C                                        BLENDING FACTORS
C                                        WF(BETA) = 0.5
            ALPHA = 4.0
            BETA  = 0.2

C---------------------------------------------------------------------72
      IF (NBND .LT. 1) CALL ERRR(501,"BOREC")
C---------------------------------------------------------------------72
C                               CALCULATION OF INNER COORDINATES
      DO I=1,II
         XINNER(I) = (X(I)-XWAND) * UTAU_ACT / GMOL
      ENDDO
C---------------------------------------------------------------------72
C                               CALCULATION OF OUTER COORDINATES
      DO I=1,II
         XOUTER(I) = (X(I)-XWAND)/DELTA_ACT 
      ENDDO
C---------------------------------------------------------------------72
C                               INTERPOLATION ON INNER COORDINATES

      DO J=NBND,JJ-NBND+1
         DO I=NBND+1,II-NBND
            
            D_INNER = (X(I) - XWAND) * UTAU/GMOL
            D_OUTER = (X(I) - XWAND) / DELTA
 
            IF (DISTANCE .GT. 0.0 ) THEN
               IF (D_INNER .GE. XINNER(II-NBND)) THEN
                  WINNER = UTAU/UTAU_ACT * W(J,II-NBND)
               ELSE
                  DO I2=1,II
                     IF (D_INNER - XINNER(I2+1) .LE. 0.0) THEN
                        
                        WINNER = UTAU/UTAU_ACT * 
     $                       ( W(J,I2) + ((W(J,I2+1) - W(J,I2))*
     $                       (D_INNER - XINNER(I2))/
     $                       (XINNER(I2+1)-XINNER(I2))))
                        
                        GOTO 100
                     ENDIF
                  ENDDO
  100             CONTINUE
               ENDIF
               IF (D_OUTER .GE. XOUTER(II-NBND)) THEN
                  WOUTER = UFRCON/UMAX_ACT * W(J,II-NBND)
               ELSE
                  DO I2=1,II
                     IF (D_OUTER - XOUTER(I2+1) .LE. 0.0) THEN
                        
                        WOUTER = UFRCON/UMAX_ACT * 
     $                       ( W(J,I2) + ((W(J,I2+1) - W(J,I2))*
     $                       (D_OUTER - XOUTER(I2))/
     $                       (XOUTER(I2+1)-XOUTER(I2))))
                        
                        GOTO 200
                     ENDIF
                  ENDDO
  200             CONTINUE
               ENDIF

            ELSE
               WINNER = 0.0
               WOUTER = 0.0
            ENDIF

            WF = 0.5*( 1.0 + 
     $        TANH(ALPHA*(D_OUTER-BETA)/((1.0-2.0*BETA)*D_OUTER + BETA))
     $       /TANH(ALPHA))

            WF = MAX(WF,1.0)

            write (6,*) 'distance, wf',d_outer,wf
 
            WBO(J,I,2) = WINNER * (1.0 - WF) + WOUTER * WF
            
            
         ENDDO
      ENDDO
C---------------------------------------------------------------------72
C---------------------------------------------------------------------72


C---------------------------------------------------------------------72
      RETURN
      END
