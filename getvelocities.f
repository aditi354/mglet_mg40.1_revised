










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
      SUBROUTINE GETVELOCITIES (KK,JJ,II,X,Y,Z,DX,DY,DZ,RDDX,RDDY,RDDZ,
     $                          U,V,W,H1,H2,HX,HY,HZ,
     $                          NPART,XP,YP,ZP,
     $                          IINDEX,JINDEX,KINDEX,UPART,
     $                          FXPART,FYPART,FZPART)
C------------------------------------------------------------
C
C     GETS NEIGHBOURING VELOCITIES AROUND PARTICLES
C
C     020999 (MM)  ORIGINAL
C
C------------------------------------------------------------

      REAL U(KK,JJ,II), V(KK,JJ,II), W(KK,JJ,II),
     $     H1(KK,JJ,II), H2(KK,JJ,II)

      REAL    X(II),    Y(JJ),    Z(KK),
     $       DX(II),   DY(JJ),   DZ(KK),
     $     RDDX(II), RDDY(JJ), RDDZ(KK),
     $       HX(II),   HY(JJ),   HZ(KK)


      REAL XP(NPART), YP(NPART), ZP(NPART)
      REAL UPART(NPART,3)

      REAL FXPART(NPART,4),FYPART(NPART,4),FZPART(NPART,4)
      
      INTEGER IINDEX(NPART), JINDEX(NPART), KINDEX(NPART)


C---------------------------------- BOUNDARIES OF PRESSURE CELLS
C      WRITE(6,*)'XP: ',XP
C      WRITE(6,*)'YP: ',YP
C      WRITE(6,*)'ZP: ',ZP

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

C---------------------------------- INTERPOLATION TO PARTICLE POSIIONS

C      DO IPART = 1,NPART

C         K=KINDEX(IPART)
C         J=JINDEX(IPART)
C         I=IINDEX(IPART)

C         UPART(IPART,1) = (U(K,J,I-1)*(HX(I    )-XP(IPART)) + 
C     $                     U(K,J,I  )*(XP(IPART)-HX(I-1  )))*RDDX(I)

C         UPART(IPART,2) = (V(K,J-1,I)*(HY(J    )-YP(IPART)) + 
C     $                     V(K,J  ,I)*(YP(IPART)-HY(J-1  )))*RDDY(J)

C         UPART(IPART,3) = (W(K-1,J,I)*(HZ(K    )-ZP(IPART)) + 
C     $                     W(K  ,J,I)*(ZP(IPART)-HZ(K-1  )))*RDDZ(K)
         
C      ENDDO

C----------------------------------------------------------------------
C
C-------------------------------      FIRST STEP: U-VELOCITY
C
C----------------------------------------------------------------------

C------------------------------- INTERPOLATION ON CORNERS OF PRESSURE CELL

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-2
               H1(K,J,I) = 0.25*(U(K  ,J  ,I) + U(K+1,J  ,I) +
     $                           U(K  ,J+1,I) + U(K+1,J+1,I))
            ENDDO
         ENDDO
      ENDDO

C---                            BOUNDARY CONDITIONS: PERIODIC/WALL/OP1 IN X

C      DO J=2,JJ-1
C         DO K=2,KK-1
C            H1(K,J,1   ) = H1(K,J,II-3)
C            H1(K,J,2   ) = H1(K,J,II-2)
C            H1(K,J,II-1) = H1(K,J,   3)
C            H1(K,J,II  ) = H1(K,J,   4)
C--------- WALL NFRO
C           H1(K,J,1) = 0.0
C           H1(K,J,2) = 0.0
C            H1(K,J,II) = 0.0
C            H1(K,J,II-1) = 0.0
C--------- OP1 NBAC
C            H1(K,J,II-2) = 0.0
C         ENDDO
C      ENDDO
C---                            BOUNDARY CONDITIONS: PERIODIC/WALL IN Y

      DO I=1,II
         DO K=2,KK-1
C--------- PERIODIC NRGT/NLFT
C            H1(K, 1  ,I) = H1(K,  JJ-3, I)
C            H1(K, 2  ,I) = H1(K,  JJ-2, I)
C            H1(K,JJ-1,I) = H1(K,     3, I)
C            H1(K,JJ  ,I) = H1(K,     4, I)
C--------- WALL NRGT
            H1(K,1,I) = 0.0
            H1(K,2,I) = 0.0
C--------- WALL NLFT
            H1(K,JJ,I)= 0.0
            H1(K,JJ-1,I)= 0.0
            H1(K,JJ-2,I)= 0.0
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
C      WRITE(6,*)'H1: ',H1(KINDEX(5),JINDEX(5),IINDEX(5))
C      WRITE(6,*)'H1: ',H1(KINDEX(5),JINDEX(5),IINDEX(5)-1)
C      WRITE(6,*)'H1: ',H1(KINDEX(5)-1,JINDEX(5),IINDEX(5))

C------------------------------ INTERPOLATION ON PARTICLE POSITION

      CALL INTERPOLATE_3D_CELL3 (KK,JJ,II,H1,RDDX,RDDY,RDDZ,HX,HY,HZ,
     $                          NPART,KINDEX,JINDEX,IINDEX,XP,YP,ZP,
     $                          UPART(1,1))

C----------------------------------------------------------------------
C
C-------------------------------      SECOND STEP: V-VELOCITY
C
C----------------------------------------------------------------------

C------------------------------- INTERPOLATION ON CORNERS OF PRESSURE CELL

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-3

               H1(K,J,I) = 0.25*(V(K  ,J  ,I  ) + V(K+1,J  ,I  ) +
     $                           V(K  ,J  ,I+1) + V(K+1,J  ,I+1))
            ENDDO
         ENDDO
      ENDDO

C---                            BOUNDARY CONDITIONS: PERIODIC IN X

      DO J=2,JJ-1
         DO K=2,KK-1
C----------- PERIODIC
C           H1(K,J,1   ) = H1(K,J,II-3)
C           H1(K,J,2   ) = H1(K,J,II-2)
C           H1(K,J,II-1) = H1(K,J,   3)
C           H1(K,J,II  ) = H1(K,J,   4)
C----------- WALL
            H1(K,J,1) = 0.0
            H1(K,J,2) = 0.0
            H1(K,J,II) = 0.0
            H1(K,J,II-1) = 0.0
            H1(K,J,II-2) = 0.0
         ENDDO
      ENDDO
C---                            BOUNDARY CONDITIONS: PERIODIC IN Y
C      DO I=1,II
C         DO K=2,KK-1
C            H1(K, 1  ,I) = H1(K,  JJ-3, I)
C            H1(K, 2  ,I) = H1(K,  JJ-2, I)
C            H1(K,JJ-1,I) = H1(K,     3, I)
C            H1(K,JJ  ,I) = H1(K,     4, I)
C            H1(K,1,I) = 0.0
C           H1(K,2,I) = 0.0
C           H1(K,JJ,I)= 0.0
C           H1(K,JJ-1,I)= 0.0
C           H1(K,JJ-2,I)= 0.0
C         ENDDO
C      ENDDO
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

C------------------------------ INTERPOLATION ON PARTICLE POSITION

      CALL INTERPOLATE_3D_CELL3 (KK,JJ,II,H1,RDDX,RDDY,RDDZ,HX,HY,HZ,
     $                          NPART,KINDEX,JINDEX,IINDEX,XP,YP,ZP,
     $                          UPART(1,2))

C----------------------------------------------------------------------
C
C-------------------------------      THIRD STEP: W-VELOCITY
C
C----------------------------------------------------------------------

C------------------------------- INTERPOLATION ON CORNERS OF PRESSURE CELL

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1

               H1(K,J,I) = 0.25*(W(K  ,J  ,I) + W(K  ,J  ,I+1) +
     $                           W(K  ,J+1,I) + W(K  ,J+1,I+1))
            ENDDO
         ENDDO
      ENDDO

C---                            BOUNDARY CONDITIONS: PERIODIC IN X

      DO J=2,JJ-1
         DO K=2,KK-1
C--------- PERIODIC
C            H1(K,J,1   ) = H1(K,J,II-3)
C            H1(K,J,2   ) = H1(K,J,II-2)
C            H1(K,J,II-1) = H1(K,J,   3)
C            H1(K,J,II  ) = H1(K,J,   4)
C--------- WALL
            H1(K,J,1) = 0.0
            H1(K,J,2) = 0.0
            H1(K,J,II) = 0.0
            H1(K,J,II-1) = 0.0
C--------- OP1
C            H1(K,J,II-2) = 0.0
         ENDDO
      ENDDO
C---                            BOUNDARY CONDITIONS: PERIODIC IN Y
      DO I=1,II
         DO K=2,KK-1
C---------- PERIODIC
C            H1(K, 1  ,I) = H1(K,  JJ-3, I)
C            H1(K, 2  ,I) = H1(K,  JJ-2, I)
C            H1(K,JJ-1,I) = H1(K,     3, I)
C            H1(K,JJ  ,I) = H1(K,     4, I)
C---------- WALL
            H1(K,1,I) = 0.0
            H1(K,2,I) = 0.0
            H1(K,JJ,I)= 0.0
            H1(K,JJ-1,I)= 0.0
            H1(K,JJ-2,I)= 0.0
         ENDDO
      ENDDO
C---                            BOUNDARY CONDITIONS: ZERO IN K
C      DO I=1,II
C         DO J=1,JJ
C           H1( 1,J,I) = 0.0
C           H1( 2,J,I) = 0.0
C           H1(KK,J,I) = 0.0
C           H1(KK-1,J,I) = 0.0
C           H1(KK-2,J,I) = 0.0
C           H1( 1,J,I) = H1( 2,J,I)
C           H1( 2,J,I) = 0.0
C           H1(KK,J,I) = H1(KK-2,J,I)
C           H1(KK-1,J,I) = H1(KK-2,J,I)
C           H1(KK-2,J,I) = 0.0
C         ENDDO
C      ENDDO

C------------------------------ INTERPOLATION ON PARTICLE POSITION

      CALL INTERPOLATE_3D_CELL3 (KK,JJ,II,H1,RDDX,RDDY,RDDZ,HX,HY,HZ,
     $                          NPART,KINDEX,JINDEX,IINDEX,XP,YP,ZP,
     $                          UPART(1,3))
C      WRITE(6,*)'UPART: ',UPART(:,3)
C      WRITE(6,*)'H3: ',H1(KINDEX(5),JINDEX(5),IINDEX(5))
C      WRITE(6,*)'H3: ',H1(KINDEX(5),JINDEX(5),IINDEX(5)-1)
C      WRITE(6,*)'H3: ',H1(KINDEX(5)-1,JINDEX(5),IINDEX(5))

C---------------------------------- READY

      RETURN
      END
      SUBROUTINE INTERPOLATE_3D_CELL1 (KK,JJ,II,H2,RDDX,RDDY,RDDZ,
     $                          HX,HY,HZ,
     $                          NPART,KINDEX,JINDEX,IINDEX,XP,YP,ZP,
     $                          DUPART)

C-------------------------------------------------------------
C
C         INTERPOLATION NACH LAGRANGESCHER INTERPOLATION 3. ORDNUNG
C
C-------------------------------------------------------------

      REAL H2(KK,JJ,II)
      REAL RDDX(II),RDDY(JJ),RDDZ(KK),HX(II),HY(JJ),HZ(KK)

      REAL FX(4),FY(4),FZ(4)

      INTEGER KINDEX(NPART),JINDEX(NPART),IINDEX(NPART)

      REAL XP(NPART),YP(NPART),ZP(NPART)
      REAL DUPART(NPART)

      DO IPART = 1,NPART
         
         X = XP(IPART)
         Y = YP(IPART)
         Z = ZP(IPART)

         K=KINDEX(IPART)
         J=JINDEX(IPART)
         I=IINDEX(IPART)

         IF ( K .LT. 3 .OR. K .GT. KK-2 .OR.
     $        J .LT. 3 .OR. J .GT. JJ-2 .OR.
     $        I .LT. 3 .OR. I .GT. II-2      ) THEN

         WRITE (6,*) 'IPART: ',IPART
         WRITE (6,*) 'X,Y,Z: ',X,Y,Z
         WRITE (6,*) 'I,J,K: ',I,J,K
         STOP
         ENDIF
C------------------------------- FAKTOREN IN X-RICHTUNG

         FX(1) = ( (X-HX(I-1)) * (X-HX(I  )) * (X-HX(I+1)) )
     $      /((HX(I-2)-HX(I-1)) * (HX(I-2)-HX(I  )) * (HX(I-2)-HX(I+1)))

         FX(2) = ( (X-HX(I-2)) * (X-HX(I  )) * (X-HX(I+1)) )
     $      /((HX(I-1)-HX(I-2)) * (HX(I-1)-HX(I  )) * (HX(I-1)-HX(I+1)))

         FX(3) = ( (X-HX(I-2)) * (X-HX(I-1)) * (X-HX(I+1)) )
     $      /((HX(I  )-HX(I-2)) * (HX(I  )-HX(I-1)) * (HX(I  )-HX(I+1)))

         FX(4) = ( (X-HX(I-2)) * (X-HX(I-1)) * (X-HX(I  )) )
     $      /((HX(I+1)-HX(I-2)) * (HX(I+1)-HX(I-1)) * (HX(I+1)-HX(I  )))

C------------------------------- FAKTOREN IN Y-RICHTUNG

         FY(1) = ( (Y-HY(J-1)) * (Y-HY(J  )) * (Y-HY(J+1)) )
     $      /((HY(J-2)-HY(J-1)) * (HY(J-2)-HY(J  )) * (HY(J-2)-HY(J+1)))

         FY(2) = ( (Y-HY(J-2)) * (Y-HY(J  )) * (Y-HY(J+1)) )
     $      /((HY(J-1)-HY(J-2)) * (HY(J-1)-HY(J  )) * (HY(J-1)-HY(J+1)))

         FY(3) = ( (Y-HY(J-2)) * (Y-HY(J-1)) * (Y-HY(J+1)) )
     $      /((HY(J  )-HY(J-2)) * (HY(J  )-HY(J-1)) * (HY(J  )-HY(J+1)))

         FY(4) = ( (Y-HY(J-2)) * (Y-HY(J-1)) * (Y-HY(J  )) )
     $      /((HY(J+1)-HY(J-2)) * (HY(J+1)-HY(J-1)) * (HY(J+1)-HY(J  )))


C------------------------------- FAKTOREN IN Z-RICHTUNG

         FZ(1) = ( (Z-HZ(K-1)) * (Z-HZ(K  )) * (Z-HZ(K+1)) )
     $      /((HZ(K-2)-HZ(K-1)) * (HZ(K-2)-HZ(K  )) * (HZ(K-2)-HZ(K+1)))

         FZ(2) = ( (Z-HZ(K-2)) * (Z-HZ(K  )) * (Z-HZ(K+1)) )
     $      /((HZ(K-1)-HZ(K-2)) * (HZ(K-1)-HZ(K  )) * (HZ(K-1)-HZ(K+1)))

         FZ(3) = ( (Z-HZ(K-2)) * (Z-HZ(K-1)) * (Z-HZ(K+1)) )
     $      /((HZ(K  )-HZ(K-2)) * (HZ(K  )-HZ(K-1)) * (HZ(K  )-HZ(K+1)))

         FZ(4) = ( (Z-HZ(K-2)) * (Z-HZ(K-1)) * (Z-HZ(K  )) )
     $      /((HZ(K+1)-HZ(K-2)) * (HZ(K+1)-HZ(K-1)) * (HZ(K+1)-HZ(K  )))


         DUDX = 0.0

         DO IINT=1,4
            DO JINT=1,4
               DO KINT=1,4
                  
                  DUDX = DUDX + 
     $                 H2(K + KINT - 3, J + JINT - 3, I + IINT - 3)*
     $                 FX(IINT)*FY(JINT)*FZ(KINT)

               ENDDO
            ENDDO
         ENDDO
C         DUDX = 
C     $          ((H2(K  ,J-2,I-2)*FX1 +
C     $            H2(K  ,J-2,I-1)*FX2 +
C     $            H2(K  ,J-2,I  )*FX3 +
C     $            H2(K  ,J-2,I+1)*FX4  )* FY1+
C     $           (H2(K  ,J-1,I-2)*FX1 +
C     $            H2(K  ,J-1,I-1)*FX2 +
C     $            H2(K  ,J-1,I  )*FX3 +
C     $            H2(K  ,J-1,I+1)*FX4  )* FY2+
C     $           (H2(K  ,J  ,I-2)*FX1 +
C     $            H2(K  ,J  ,I-1)*FX2 +
C     $            H2(K  ,J  ,I  )*FX3 +
C     $            H2(K  ,J  ,I+1)*FX4  )* FY3+
C     $           (H2(K  ,J-2,I-2)*FX1 +
C     $            H2(K  ,J-2,I-1)*FX2 +
C     $            H2(K  ,J-2,I  )*FX3 +
C     $            H2(K  ,J-2,I+1)*FX4  )* FY4+

         DUPART(IPART) = DUDX
 

       IF ( DUPART(ipart) .GT. 2.0 ) then
             WRITE (6,*) 'GETVEL',ipart,dUPART(ipart)
             WRITE (6,*) 'GETVEL',ipart,k,j,i
             WRITE (6,*) 'GETVEL',ipart,z,y,x
             WRITE (6,*) 'GETVEL',ipart,H2(k,j,i)
         ENDIF


      ENDDO

      RETURN
      END

