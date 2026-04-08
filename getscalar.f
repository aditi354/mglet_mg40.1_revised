










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
      SUBROUTINE GETSCALAR (KK,JJ,II,X,Y,Z,DX,DY,DZ,RDDX,RDDY,RDDZ,
     $     DDX,DDY,DDZ,
     $     T,H1,HX,HY,HZ,
     $     NPART,XP,YP,ZP,
     $     IINDEX,JINDEX,KINDEX,
     $     FXPART,FYPART,FZPART,PARTSCA,
     $     NFRO,NBAC,NRGT,NLFT,NBOT,NTOP)

C------------------------------------------------------------
C
C     GETS SCALAR CONCENTRATION AT PARTICLE POSITIONS
C
C     08.01.04 (FS) ORIGINAL
C
C------------------------------------------------------------
      IMPLICIT NONE

      INTEGER KK,JJ,II,I,J,K,NPART,NCUB,
     $     NFRO,NBAC,NRGT,NLFT,NBOT,NTOP
  
      REAL T(KK,JJ,II),
     $     H1(KK,JJ,II)

      REAL    X(II),    Y(JJ),    Z(KK),
     $       DX(II),   DY(JJ),   DZ(KK),
     $       DDX(II),  DDY(JJ),  DDZ(KK),
     $     RDDX(II), RDDY(JJ), RDDZ(KK),
     $       HX(II),   HY(JJ),   HZ(KK),
     $     DSUM

      REAL XP(NPART),YP(NPART),ZP(NPART),PARTSCA(NPART)

      REAL FXPART(NPART,4),FYPART(NPART,4),FZPART(NPART,4)
      
      INTEGER IINDEX(NPART), JINDEX(NPART), KINDEX(NPART),
     $     ISUM

C---------------------------------- BOUNDARIES OF PRESSURE CELLS
      DO I=1,II
         HX(I) = X(I) + 0.5*DX(I)
      ENDDO

      DO J=1,JJ
         HY(J) = Y(J) + 0.5*DY(J)
      ENDDO

      DO K=1,KK
         HZ(K) = Z(K) + 0.5*DZ(K)
      ENDDO

C------------------------------------------- FAKTORS FOR INTERPOLATION

      CALL INT3DFAC4(NPART,IINDEX,XP,II,HX,RDDX,
     $              FXPART(1,1),FXPART(1,2),FXPART(1,3),FXPART(1,4))

      CALL INT3DFAC4(NPART,JINDEX,YP,JJ,HY,RDDY,
     $              FYPART(1,1),FYPART(1,2),FYPART(1,3),FYPART(1,4))

      CALL INT3DFAC4(NPART,KINDEX,ZP,KK,HZ,RDDZ,
     $              FZPART(1,1),FZPART(1,2),FZPART(1,3),FZPART(1,4))

C----------- INTERPOLATION OF H1 FIELD TO CORNERS OF PESSURE CELL

      DO I=2,II-2
         DO J=2,JJ-2
            DO K=2,KK-2
               H1(K,J,I) = 0.125*(T(K,J,I) + T(K+1,J,I) +
     $                            T(K,J+1,I) + T(K+1,J+1,I) +
     $                            T(K,J,I+1) + T(K+1,J,I+1) +
     $                            T(K,J+1,I+1) + T(K+1,J+1,I+1))
            ENDDO
         ENDDO
      ENDDO
C----------- ASSUMING: BOUNDARYS CORRECTLY SET IN T-FIELD
      DO J=2,JJ-1
         DO K=2,KK-1
            H1(K,J,1   ) = H1(K,J,II-3)
            H1(K,J,2   ) = H1(K,J,II-2)
            H1(K,J,II-1) = H1(K,J,   3)
            H1(K,J,II  ) = H1(K,J,   4)
C--------- WALL NFRO
C           H1(K,J,1) = 0.0
C           H1(K,J,2) = 0.0
C            H1(K,J,II) = 0.0
C            H1(K,J,II-1) = 0.0
C--------- OP1 NBAC
C            H1(K,J,II-2) = 0.0
         ENDDO
      ENDDO
C---                            BOUNDARY CONDITIONS: PERIODIC/WALL IN Y

      DO I=1,II
         DO K=2,KK-1
C--------- PERIODIC NRGT/NLFT
            H1(K, 1  ,I) = H1(K,  JJ-3, I)
            H1(K, 2  ,I) = H1(K,  JJ-2, I)
            H1(K,JJ-1,I) = H1(K,     3, I)
            H1(K,JJ  ,I) = H1(K,     4, I)
C--------- WALL NRGT
c            H1(K,1,I) = 0.0
c            H1(K,2,I) = 0.0
C--------- WALL NLFT
c            H1(K,JJ,I)= 0.0
c            H1(K,JJ-1,I)= 0.0
c            H1(K,JJ-2,I)= 0.0
         ENDDO
      ENDDO
C---                            BOUNDARY CONDITIONS: ZERO IN K
      DO I=1,II
         DO J=1,JJ
           H1( 1,J,I) = 0.0
           H1( 2,J,I) = 0.0
           H1(KK,J,I) = 0.0
           H1(KK-1,J,I) = 0.0
           H1(KK-2,J,I) = 0.0
         ENDDO
      ENDDO

      CALL INTERPOLATE_3D_CELL3 (KK,JJ,II,H1,RDDX,RDDY,RDDZ,HX,HY,HZ,
     $                          NPART,KINDEX,JINDEX,IINDEX,XP,YP,ZP,
     $                          PARTSCA)

      RETURN
      END
