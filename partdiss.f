










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
      SUBROUTINE PARTDISS (KK,JJ,II,X,Y,Z,DX,DY,DZ,RDDX,RDDY,RDDZ,
     $     DDX,DDY,DDZ,
     $     U,V,W,B,H1,H2,HX,HY,HZ,
     $     NPART,XP,YP,ZP,
     $     IINDEX,JINDEX,KINDEX,DISSI,DISTURB,
     $     FXPART,FYPART,FZPART,
     $     AU,AV,AW,
     $     KKA,JJA,IIA,
     $     FELD1,FELD2,FELD3,
     $     BP,
     $     NFRO,NBAC,NRGT,NLFT,NBOT,NTOP)

C------------------------------------------------------------
C
C     GETS OVERALL DISSIPATION AT PARTICL POSITIONS
C
C     22.05.03 (FS) ORIGINAL
C
C------------------------------------------------------------
      IMPLICIT NONE

      INTEGER KK,JJ,II,I,J,K,NPART,NCUB,
     $     NFRO,NBAC,NRGT,NLFT,NBOT,NTOP
  
      REAL U(KK,JJ,II), V(KK,JJ,II), W(KK,JJ,II),
     $     H1(KK,JJ,II), B(KK,JJ,II), H2(KK,JJ,II)

      REAL    X(II),    Y(JJ),    Z(KK),
     $       DX(II),   DY(JJ),   DZ(KK),
     $       DDX(II),  DDY(JJ),  DDZ(KK),
     $     RDDX(II), RDDY(JJ), RDDZ(KK),
     $       HX(II),   HY(JJ),   HZ(KK),
     $     DSUM

      INTEGER KKA,JJA,IIA
      REAL AU(KKA,JJA,IIA),AV(KKA,JJA,IIA),AW(KKA,JJA,IIA),
     $     FELD1(KK,JJ,II),FELD2(KK,JJ,II),FELD3(KK,JJ,II)

      REAL BP(KK,JJ,II)

      REAL XP(NPART),YP(NPART),ZP(NPART),DISSI(NPART),DISTURB(NPART)

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

C---------------------------------------- CALCULATION OF DISSIPATION
      ISUM = 0
            CALL DISSIPG (KK,JJ,II,KK,JJ,II,DDX,DDY,DDZ,
     $                    DX,DY,DZ,X,Y,Z,
     $                    U,V,W,B,H2,DSUM,
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,ISUM)

C----------- INTERPOLATION OF H1 FIELD TO CORNERS OF PESSURE CELL

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-2
               H1(K,J,I) = 0.125*(H2(K,J,I) + H2(K+1,J,I) +
     $                            H2(K,J+1,I) + H2(K+1,J+1,I) +
     $                            H2(K,J,I+1) + H2(K+1,J,I+1) +
     $                            H2(K,J+1,I+1) + H2(K+1,J+1,I+1))
            ENDDO
         ENDDO
      ENDDO
C----------- BOUNDARY CONDITIONS: X
      DO J=2,JJ-1
         DO K=2,KK-1
C------ WALL NFRO
            H1(K,J,1) = 0.0
            H1(K,J,2) = 0.0
         ENDDO
      ENDDO

C----------- BOUNDARY CONDITIONS Y
      DO I=1,II
         DO K=2,KK-1   
C--------- WALL NRGT
            H1(K,1,I) = 0.0
            H1(K,2,I) = 0.0
C--------- WALL NLFT
            H1(K,JJ,I)= 0.0
            H1(K,JJ-1,I)= 0.0
            H1(K,JJ-2,I)= 0.0
         ENDDO
      ENDDO

C----------- BOUNDARY CONDITIONS Z
      DO I=1,II
         DO J=1,JJ
C----------- PROBLEM: INFLOW REGION
            H1( 1,J,I) = 0.0
            H1( 2,J,I) = 0.0
            H1(KK,J,I) = 0.0
            H1(KK-1,J,I) = 0.0
            H1(KK-2,J,I) = 0.0
         ENDDO
      ENDDO

      CALL INTERPOLATE_3D_CELL3 (KK,JJ,II,H1,RDDX,RDDY,RDDZ,HX,HY,HZ,
     $                          NPART,KINDEX,JINDEX,IINDEX,XP,YP,ZP,
     $                          DISSI)

C-----------------------------------------------------------------------
C---------------------------------- CALCULATION OF TURBULENT DISSIPATION
C-----------------------------------------------------------------------

      IF (KKA .NE. KK .OR. JJA .NE. JJ .OR. IIA .NE. II) THEN
         WRITE(6,*)'ERROR IN PARTDISS'

         RETURN
      ENDIF
C---------------------------------- CALCULATION OF FLUCUATION
      DO K=1,KK
         DO J=1,JJ
            DO I=1,II
               FELD1(K,J,I) = U(K,J,I)-AU(K,J,I)
               FELD2(K,J,I) = V(K,J,I)-AV(K,J,I)
               FELD3(K,J,I) = W(K,J,I)-AW(K,J,I)
            ENDDO
         ENDDO
      ENDDO
      DO K=1,KK
         DO J=1,JJ
            DO I=1,II
               FELD1(K,J,I) = FELD1(K,J,I)*BP(K,J,I)
               FELD2(K,J,I) = FELD2(K,J,I)*BP(K,J,I)
               FELD3(K,J,I) = FELD3(K,J,I)*BP(K,J,I)
            ENDDO
         ENDDO
      ENDDO
      ISUM = 0
      CALL DISSIPG (KK,JJ,II,KK,JJ,II,DDX,DDY,DDZ,
     $     DX,DY,DZ,X,Y,Z,
     $     FELD1,FELD2,FELD3,B,H2,DSUM,
     $     NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,ISUM)
C----------- INTERPOLATION OF H1 FIELD TO CORNERS OF PESSURE CELL

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-2
               H1(K,J,I) = 0.125*(H2(K,J,I) + H2(K+1,J,I) +
     $                            H2(K,J+1,I) + H2(K+1,J+1,I) +
     $                            H2(K,J,I+1) + H2(K+1,J,I+1) +
     $                            H2(K,J+1,I+1) + H2(K+1,J+1,I+1))
            ENDDO
         ENDDO
      ENDDO
C----------- BOUNDARY CONDITIONS: X
      DO J=2,JJ-1
         DO K=2,KK-1
C------ WALL NFRO
            H1(K,J,1) = 0.0
            H1(K,J,2) = 0.0
         ENDDO
      ENDDO

C----------- BOUNDARY CONDITIONS Y
      DO I=1,II
         DO K=2,KK-1   
C--------- WALL NRGT
            H1(K,1,I) = 0.0
            H1(K,2,I) = 0.0
C--------- WALL NLFT
            H1(K,JJ,I)= 0.0
            H1(K,JJ-1,I)= 0.0
            H1(K,JJ-2,I)= 0.0
         ENDDO
      ENDDO

C----------- BOUNDARY CONDITIONS Z
      DO I=1,II
         DO J=1,JJ
C----------- PROBLEM: INFLOW REGION
            H1( 1,J,I) = 0.0
            H1( 2,J,I) = 0.0
            H1(KK,J,I) = 0.0
            H1(KK-1,J,I) = 0.0
            H1(KK-2,J,I) = 0.0
         ENDDO
      ENDDO

      CALL INTERPOLATE_3D_CELL3 (KK,JJ,II,H1,RDDX,RDDY,RDDZ,HX,HY,HZ,
     $                          NPART,KINDEX,JINDEX,IINDEX,XP,YP,ZP,
     $                          DISTURB)
C      WRITE(6,*)'DISS: ',DISTURB,DISSI
      RETURN
      END
