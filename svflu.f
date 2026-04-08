










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
      SUBROUTINE SVFLU(KK,JJ,II,UFR,U,HILF,AU,JJA,NBND,
     $                 Z,ZBANF,DELTA)
C--MGLET----------------------------------------------------------
C
C     ADDS FLUCTUATIONS TO FRONT-BUFFERS FOR 
C     "STRUKTUR-PERIODISCHE" BOUNDARY CONDITION
C
C        UFR(X0)=UFR(X0) + (U(X) - <U>(X))
C
C    10.12.96 (MM.):    ORIGINAL
C
C-------10--------20--------30--------40--------50--------60--------7072

      REAL
     $     UFR(KK,JJ, 2),U(KK,JJ),HILF(KK,JJ),Z(KK)
      REAL AU(KK,JJA)

C---------------------------------------------------------------------72
      IF (NBND .LT. 1) CALL ERRR(501,"SVFLU")

C-------------------------------- AVERAGING IN Y ---------------------72

      FAC = 1./FLOAT(JJ - 2*NBND)

C                                        FIRST VALUE
      DO K=NBND+1,KK-NBND
         HILF(K,NBND+1) = FAC*U(K,NBND+1)
      ENDDO

C                                        SUMMING UP


      DO J=NBND+2,JJ-NBND
         DO K= NBND+1,KK-NBND
            HILF(K,J) = HILF(K,J-1) + FAC*U(K,J)
         ENDDO
      ENDDO
C                                      AVERAGE IS IN HILF(K,JJ-NBND)
C------------------------------------- NOW ADDING FLUCTUATIONS -------72

      DO K=NBND+1,KK-NBND
        DISTANCE = (ZBANF+DELTA-Z(K))/DELTA
        FAK = EXP( MIN( 0.0,DISTANCE))**4
CCC		write (6,*) 'svflu:',z(k),distance,fak
        DO J=   NBND+1,JJ-NBND

           UFR(K,J,2)=UFR(K,J,2) + (U(K,J)-HILF(K,JJ-NBND))*FAK

        ENDDO
      ENDDO

C---------------------------------------------------------------------72

      RETURN
      END
