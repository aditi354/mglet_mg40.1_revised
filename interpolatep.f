










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
      SUBROUTINE INTERPOLATEP (KK,JJ,II,IVAR,POUT,P)
C------------------------------------------------------------
C
C     INTERPOLATES P ON STAGGERED GRIDS
C
C     13.07.00 (SE) ORIGINAL
C
C------------------------------------------------------------

      REAL POUT(KK,JJ,II), P(KK,JJ,II)

      CHARACTER (LEN=2) IVAR

C------------------------------------------------------------

      IF ( IVAR .EQ. 'XY' ) THEN
         DO I=2,II-1
            DO J=2,JJ-1
               DO K=2,KK-1

                  POUT(K,J,I) = .25
     $                 *( P(K,J,  I) + P(K,J,  I+1)
     $                 +  P(K,J+1,I) + P(K,J+1,I+1))
               ENDDO
            ENDDO
         ENDDO
      ELSEIF ( IVAR .EQ. 'XZ' ) THEN
         DO I=2,II-1
            DO J=2,JJ-1
               DO K=2,KK-1

                  POUT(K,J,I) = .25
     $                 *( P(K,J,I)   + P(K,J,I+1)
     $                 +  P(K+1,J,I) + P(K+1,J,I+1))
               ENDDO
            ENDDO
         ENDDO
      ELSEIF ( IVAR .EQ. 'YZ' ) THEN
         DO I=2,II-1
            DO J=2,JJ-1
               DO K=2,KK-1

                  POUT(K,J,I) = .25
     $                 *( P(K  ,J,I) + P(K  ,J+1,I)
     $                 +  P(K+1,J,I) + P(K+1,J+1,I))
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      RETURN
      END
