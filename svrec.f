










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
      SUBROUTINE SVREC(KK,JJ,II,NBND,UFR,U,HILF,ZINNER,ZOUTER,
     $     Z,ZWAND,DELTA,DELTA_ACT,UTAU,UTAU_ACT,UFRCON,UMAX_ACT,
     $     GMOL)
C--MGLET----------------------------------------------------------
C
C     "RECYCLING" BOUNDARY CONDITION (LUND, SPALART)
C     
C
C        UFR(X0)=U_I^INNER + U_I^OUTER
C
C    24.03.98 (MM.):    ORIGINAL 
C
C-------10--------20--------30--------40--------50--------60--------7072

      REAL
     $     UFR(KK,JJ, 2),U(KK,JJ),HILF(KK,JJ),
     $     Z(KK),ZINNER(KK),ZOUTER(KK)

C---------------------------------------------------------------------72
C                                        BLENDING FACTORS
C                                        WF(BETA) = 0.5
            ALPHA = 4.0
            BETA  = 0.2

C---------------------------------------------------------------------72
      IF (NBND .LT. 1) CALL ERRR(501,"SVREC")
C---------------------------------------------------------------------72
C                               CALCULATION OF INNER COORDINATES
      DO K=1,KK
         ZINNER(K) = (Z(K)-ZWAND) * UTAU_ACT / GMOL
      ENDDO
C---------------------------------------------------------------------72
C                               CALCULATION OF OUTER COORDINATES
      DO K=1,KK
         ZOUTER(K) = (Z(K)-ZWAND)/DELTA_ACT 
      ENDDO
C---------------------------------------------------------------------72
C                               INTERPOLATION ON INNER COORDINATES

      DO J=NBND,JJ-NBND+1
         DO K=NBND+1,KK-NBND
            
            D_INNER = (Z(K) - ZWAND) * UTAU/GMOL
            D_OUTER = (Z(K) - ZWAND) / DELTA
 
            IF (D_OUTER .GT. 0.0 ) THEN
               IF (D_INNER .GE. ZINNER(KK-NBND)) THEN
                  UINNER = UTAU/UTAU_ACT * U(KK-NBND,J)
               ELSE
                  DO K2=1,KK
                     IF (D_INNER - ZINNER(K2+1) .LE. 0.0) THEN
                        
                        UINNER = UTAU/UTAU_ACT * 
     $                       ( U(K2,J) + ((U(K2+1,J) - U(K2,J))*
     $                       (D_INNER - ZINNER(K2))/
     $                       (ZINNER(K2+1)-ZINNER(K2))))
                        
                        GOTO 100
                     ENDIF
                  ENDDO
  100              CONTINUE
               ENDIF
               IF (D_OUTER .GE. ZOUTER(KK-NBND)) THEN
                  UOUTER = UFRCON/UMAX_ACT * U(KK-NBND,J)
               ELSE
                  DO K2=1,KK
                     IF (D_OUTER - ZOUTER(K2+1) .LE. 0.0) THEN
                        
                        UOUTER = UFRCON/UMAX_ACT * 
     $                       ( U(K2,J) + ((U(K2+1,J) - U(K2,J))*
     $                       (D_OUTER - ZOUTER(K2))/
     $                       (ZOUTER(K2+1)-ZOUTER(K2))))
                        
                        GOTO 200
                     ENDIF
                  ENDDO
  200              CONTINUE
               ENDIF

            ELSE
               UINNER = 0.0
               UOUTER = 0.0
            ENDIF

            WF = 0.5*( 1.0 + 
     $        TANH(ALPHA*(D_OUTER-BETA)/((1.0-2.0*BETA)*D_OUTER + BETA))
     $       /TANH(ALPHA))

            WF = MIN(WF,1.0)

            UFR(K,J,2) = UINNER * (1.0 - WF) + UOUTER * WF
            
            
         ENDDO
      ENDDO
C---------------------------------------------------------------------72
C---------------------------------------------------------------------72


C---------------------------------------------------------------------72
      RETURN
      END
