










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
      SUBROUTINE DFDX (KK,JJ,II,RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                  IVAR,PHI,DPHI)
C------------------------------------------------------------
C
C     GETS DERIVATIVES IN SPACE ON STAGGERED GRIDS
C
C     260600 (SE)  ORIGINAL
C
C------------------------------------------------------------

      REAL PHI(KK,JJ,II), DPHI(KK,JJ,II)

      REAL  RDX(II),  RDY(JJ),  RDZ(KK),
     $     RDDX(II), RDDY(JJ), RDDZ(KK)

      INTEGER I,J,K

      CHARACTER (LEN=3) IVAR

C---------------------------------- COMPUTE DPHIDX

      IF ( IVAR .EQ. 'DXS' ) THEN
         DO I=2,II-1
            DO J=2,JJ-1
               DO K=2,KK-1

                  DPHI(K,J,I) = (PHI(K,J,I+1) - PHI(K,J,I))*RDX(I)

               ENDDO
            ENDDO
         ENDDO

C---------------------------------- COMPUTE DPHIDX ON PRESSURE POINT
      ELSEIF (IVAR .EQ. 'DXT' ) THEN
       DO I=3,II-2
        DO J=3,JJ-2
         DO K=3,KK-2

             DPHI(K,J,I) = (PHI(K,J,I+1) - PHI(K,J,I-1))/
     $                      (1/RDX(I) + 1/RDX(I-1))
    
         ENDDO
        ENDDO
       ENDDO

      ELSEIF ( IVAR .EQ. 'DDX' ) THEN
         DO I=2,II-1
            DO J=2,JJ-1
               DO K=2,KK-1

                  DPHI(K,J,I) = (PHI(K,J,I) - PHI(K,J,I-1))*RDDX(I)

               ENDDO
            ENDDO
         ENDDO

C---------------------------------- COMPUTE DPHIDY

      ELSEIF ( IVAR .EQ. 'DYS' ) THEN
         DO I=2,II-1
            DO J=2,JJ-1
               DO K=2,KK-1

                  DPHI(K,J,I) = (PHI(K,J+1,I) - PHI(K,J,I))*RDY(J)

               ENDDO
            ENDDO
         ENDDO

      ELSEIF (IVAR .EQ. 'DYT' ) THEN
       DO I=3,II-2
        DO J=3,JJ-2
         DO K=3,KK-2

             DPHI(K,J,I) = (PHI(K,J+1,I) - PHI(K,J-1,I))/
     $                      (1/RDY(J) + 1/RDY(J-1))

         ENDDO
        ENDDO
       ENDDO

      ELSEIF ( IVAR .EQ. 'DDY' ) THEN
         DO I=2,II-1
            DO J=2,JJ-1
               DO K=2,KK-1

                  DPHI(K,J,I) = (PHI(K,J,I) - PHI(K,J-1,I))*RDDY(J)

               ENDDO
            ENDDO
         ENDDO

C---------------------------------- COMPUTE DPHIDZ

      ELSEIF ( IVAR .EQ. 'DZS' ) THEN
         DO I=2,II-1
            DO J=2,JJ-1
               DO K=2,KK-1

                  DPHI(K,J,I) = (PHI(K+1,J,I) - PHI(K,J,I))*RDZ(K)

               ENDDO
            ENDDO
         ENDDO

      ELSEIF (IVAR .EQ. 'DZT' ) THEN
       DO I=3,II-2
        DO J=3,JJ-2
             K=3
             DPHI(K,J,I) = (PHI(K+1,J,I) - PHI(K-1,J,I))/
     $                      (1/RDZ(K) + 1/(0.5*RDZ(K-1)))
         DO K=4,KK-3

             DPHI(K,J,I) = (PHI(K+1,J,I) - PHI(K-1,J,I))/
     $                      (1/RDZ(K) + 1/RDZ(K-1))

         ENDDO
            K=KK-2
             DPHI(K,J,I) = (PHI(K+1,J,I) - PHI(K-1,J,I))/
     $                      (1/(0.5*RDZ(K)) + 1/RDZ(K-1))
        ENDDO
       ENDDO


      ELSEIF ( IVAR .EQ. 'DDZ' ) THEN
         DO I=2,II-1
            DO J=2,JJ-1
               DO K=2,KK-1

                  DPHI(K,J,I) = (PHI(K,J,I) - PHI(K-1,J,I))*RDDZ(K)

               ENDDO
            ENDDO
         ENDDO

      ENDIF
   
      RETURN
      END
