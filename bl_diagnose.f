










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
      SUBROUTINE BL_DIAGNOSE (KK,JJ,II,NBND,U,Z,ZBANF,DELTA,UTAU,UMAX)
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Calculation of boundary layer thickness and utau
C
C     VERSION VOM 23.3.1998 (MM) 
C                          II and JJ are assumed to be 1
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

      REAL U(KK,JJ,II)
      REAL Z(KK)

C---- --------------------  CHECK OF THE DIMENSIONS

      IF (JJ .NE. 1) CALL ERRR (501,"BL_DIAGNOSE")
      IF (II .NE. 1) CALL ERRR (502,"BL_DIAGNOSE")

C---- --------------------  Search of the maximum

      UMAX = -100000

      DO K=NBND+1,KK-NBND
         UMAX = MAX(UMAX,U(K,1,1))
      ENDDO

      U995 = 0.995*UMAX

C---- -------------------- SEARCH OF THE 99.5%-POSITION

      DO K=NBND+1,KK-NBND
         IF ( U(K,1,1) .GT. U995 .AND. U(K-1,1,1) .LE. U995 ) THEN
            DELTA = Z(K-1) +
     $           (Z(K)-Z(K-1))*(U995-U(K-1,1,1))/(U(K,1,1)-U(K-1,1,1))
            GOTO 999
         ENDIF
      ENDDO

  999 CONTINUE

C---- -------------------- CALCULATION OF UTAU

      DO K=NBND+1,KK-NBND
         IF (Z(K) .GT. ZBANF) THEN
C                                                   FUNCTION UTAU ASSUMES
C                                                   THAT DDS IS TWO TIMES
C                                                   THE WALL DISTANCE OF THE
C                                                   DEFINITION POINT
            UTAU = UTAUP (U(K,1,1),2.0*(Z(K)-ZBANF))
            write (6,*) "utaup",utau,U(K,1,1),Z(K)-ZBANF
            GOTO 998
         ENDIF
      ENDDO

  998 CONTINUE

C---- -------------------- CHECK-OUTPUT


      RETURN
      END

