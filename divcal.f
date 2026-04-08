










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
      SUBROUTINE DIVCAL(KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,
     $                  U,V,W,B,DIV,FAK,IREZIP,DIVMAX,
     $                  BU,BV,BW,SDIV)
C*MGLET***************************************************
C  D I V C A L     BERECHNET DIE DIVERGENZ
C
C      FAK:        FALLS NICHT NULL, WIRD DIVERGENZ DAMIT NORMIERT
C   IREZIP:         +1:    DDX,DDY,DDZ SIND ABSTAENDE DER GITTERP.
C                   -1:    DDX,...    SIND SCHON REZIPROKWERTE
C
C    11. 3.94  (MM):   ERWEITERT UM FAK UND IREZIP
C
C*MGLET***************************************************
C
      IMPLICIT NONE

      INTEGER KK,JJ,II,K,J,I,KMX,JMX,IMX,KDM,JDM,IDM,IREZIP
      REAL U(KK,JJ,II), V(KK,JJ,II), W(KK,JJ,II), DIV(KK,JJ,II),
     &     B(KK,JJ,II),BU(KK,JJ,II),BV(KK,JJ,II),BW(KK,JJ,II),
     $     SDIV(KK,JJ,II),
     &     DDX(II),     DDY(JJ),     DDZ(KK)

      REAL FAK,FACTOR,DIVMAX

C
C      DO K=1,KK
C	WRITE(6,*)'U: ',U(K,20,20),K
C      ENDDO
      IF (FAK .EQ. 0.0) FAK = 1.0
C
C---- --------------------------------------------------------
C
      IF (IREZIP .GE. 0) THEN

        DO I=3,IMX-2
        DO J=3,JMX-2
        DO K=3,KMX-2
               FACTOR = B(K,J,I)
          DIV(K,J,I) =  FACTOR*
     &             (( U(K,J,I) - U(K,J,I-1) )/DDX(I) 
     &          +   ( V(K,J,I) - V(K,J-1,I) )/DDY(J) 
     &          +   ( W(K,J,I) - W(K-1,J,I) )/DDZ(K) )*FAK

        ENDDO
        ENDDO
        ENDDO

C
C---- --------------------------------------------------------
C
      ELSEIF (IREZIP .LT. 0) THEN
        IF (ABS(IREZIP).GT.3) THEN
        DO I=3,IMX-2
        DO J=3,JMX-2
        DO K=3,KMX-2
               FACTOR = B(K,J,I)
          DIV(K,J,I) =  FACTOR*
     &          (( U(K,J,I)*BU(K,J,I)
     &           - U(K,J,I-1)*BU(K,J,I-1) )*DDX(I)
     &        +   ( V(K,J,I)*BV(K,J,I)
     &           - V(K,J-1,I)*BV(K,J-1,I) )*DDY(J)
     &       +   ( W(K,J,I)*BW(K,J,I)
     &           - W(K-1,J,I)*BW(K-1,J,I) )*DDZ(K)+
     &             SDIV(K,J,I) )*FAK

        ENDDO
        ENDDO
        ENDDO
        ELSE

        DO I=3,IMX-2
        DO J=3,JMX-2
        DO K=3,KMX-2
               FACTOR = B(K,J,I)
          DIV(K,J,I) =  FACTOR*
     &             (( U(K,J,I) - U(K,J,I-1) )*DDX(I) 
     &          +   ( V(K,J,I) - V(K,J-1,I) )*DDY(J) 
     &          +   ( W(K,J,I) - W(K-1,J,I) )*DDZ(K) )*FAK

        ENDDO
        ENDDO
        ENDDO
        END IF

C
C---- --------------------------------------------------------
C
      ELSE
         CALL ERRR (501,'DIVCAL')
      ENDIF

C
C---- --------------------------------------------------------
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC FOR INFORMATION
C

      DIVMAX = 0.0
      IF (ABS(IREZIP) .EQ. 2 .OR.
     $    ABS(IREZIP) .EQ. 20 ) THEN

      DO I=3,IMX-2
      DO J=3,JMX-2
      DO K=3,KMX-2
         IF (ABS(DIV(K,J,I)) .GT. ABS(DIVMAX)) THEN
                 DIVMAX=ABS(DIV(K,J,I))
         ENDIF
      ENDDO
      ENDDO
      ENDDO

C      WRITE (6,*) 'MAX. DIVERGENZ:',DIVMAX/FAK,IREZIP

      ENDIF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC FOR DETAILED INFORMATION
C
      IF (ABS(IREZIP) .EQ. 3 .OR.
     $    ABS(IREZIP) .EQ. 30 ) THEN

      DO I=3,IMX-2
      DO J=3,JMX-2
      DO K=3,KMX-2

         IF (ABS(DIV(K,J,I)) .GT. ABS(DIVMAX)) THEN

           KDM = K
           JDM = J
           IDM = I
           DIVMAX = ABS(DIV(K,J,I))

         ENDIF
      ENDDO
      ENDDO
      ENDDO

C      WRITE (6,*) 'MAX. DIVERGENZ:',DIVMAX/FAK,KDM,JDM,IDM,IREZIP

      ENDIF
      RETURN
      END
