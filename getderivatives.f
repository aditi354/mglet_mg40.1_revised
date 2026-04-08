










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
      SUBROUTINE GETDERIVATIVES (KK,JJ,II,X,Y,Z,DX,DY,DZ,
     $                          RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                          U,V,W,H1,H2,HX,HY,HZ,
     $                          NPART,XP,YP,ZP,
     $                          IINDEX,JINDEX,KINDEX,DUPART,
     $                          FXPART,FYPART,FZPART)
C------------------------------------------------------------
C
C     GETS VELOCITY DERIVATIVES AROUND PARTICLES
C
C     020999 (MM)  ORIGINAL
C     170603 (FS)  CHANGED FOR T-MIXER
C
C------------------------------------------------------------

      REAL U(KK,JJ,II), V(KK,JJ,II), W(KK,JJ,II),
     $     H1(KK,JJ,II), H2(KK,JJ,II)

      REAL    X(II),    Y(II),    Z(II),
     $       DX(II),   DY(JJ),   DZ(KK),
     $      RDX(II),  RDY(JJ),  RDZ(KK),
     $     RDDX(II), RDDY(JJ), RDDZ(KK),
     $       HX(II),   HY(JJ),   HZ(KK)


      REAL XP(NPART), YP(NPART), ZP(NPART)
      REAL DUPART(NPART,9)
      
      INTEGER IINDEX(NPART), JINDEX(NPART), KINDEX(NPART)

      REAL FXPART(NPART,4),FYPART(NPART,4),FZPART(NPART,4)

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

C---------------------------------- EVALUATION OF DUDX

C---------------------------------- COMPUTE DUDX ON CELL-CENTERED POSITIONS

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1

               H1(K,J,I) = (U(K,J,I) - U(K,J,I-1))*RDDX(I)

            ENDDO
         ENDDO
      ENDDO
C-----------------  BOUNDARY CONDITION NOTE: ONLY PERIODIC!!!!

C         DO J=2,JJ-1
C            DO K=2,KK-1
C
C               H1(K,J,2) = H1(K,J,II-2)
C
C            ENDDO
C         ENDDO

C--------------------------------- INTERPOLATE DUDX ON CORNERS OF CELLS
C--------------------------------- GRID IS TAKEN TO BE EQUIDISTANT IN ALL
C--------------------------------- THREE DIRECTIONS

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1

               H2(K,J,I) = 0.125 * (H1(K  ,J  ,I  )+
     $                              H1(K+1,J  ,I  )+
     $                              H1(K  ,J+1,I  )+
     $                              H1(K+1,J+1,I  )+
     $                              H1(K  ,J  ,I+1)+
     $                              H1(K+1,J  ,I+1)+
     $                              H1(K  ,J+1,I+1)+
     $                              H1(K+1,J+1,I+1) )

            ENDDO
         ENDDO
      ENDDO

      CALL INTERPOLATE_3D_CELL3 (KK,JJ,II,H2,RDDX,RDDY,RDDZ,HX,HY,HZ,
     $                          NPART,KINDEX,JINDEX,IINDEX,XP,YP,ZP,
     $                          DUPART(1,1))


C------------------------------------------------------------------------
C---------------------------------- EVALUATION OF DUDY ------------------
C---------------------------------- -------------------------------------

C---------------------------------- COMPUTE DUDY ON 
C                                   X- AND Y-STAGGERED POSITIONS

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1

               H1(K,J,I) = (U(K,J+1,I) - U(K,J,I))*RDY(J)

            ENDDO
         ENDDO
      ENDDO

C--------------------------------- INTERPOLATE DUDX ON CORNERS OF CELLS
C--------------------------------- GRID IS TAKEN TO BE EQUIDISTANT IN ALL
C--------------------------------- THREE DIRECTIONS

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1

               H2(K,J,I) = 0.5 * (H1(K  ,J  ,I  )+
     $                            H1(K+1,J  ,I  ))

            ENDDO
         ENDDO
      ENDDO

      CALL INTERPOLATE_3D_CELL3 (KK,JJ,II,H2,RDDX,RDDY,RDDZ,HX,HY,HZ,
     $                          NPART,KINDEX,JINDEX,IINDEX,XP,YP,ZP,
     $                          DUPART(1,2))

C------------------------------------------------------------------------
C---------------------------------- EVALUATION OF DUDZ ------------------
C---------------------------------- -------------------------------------

C---------------------------------- COMPUTE DUDZ ON 
C                                   X- AND Z-STAGGERED POSITIONS

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1

               H1(K,J,I) = (U(K+1,J,I) - U(K,J,I))*RDZ(K)

            ENDDO
         ENDDO
      ENDDO

C--------------------------------- INTERPOLATE DUDZ ON CORNERS OF CELLS
C--------------------------------- GRID IS TAKEN TO BE EQUIDISTANT IN ALL
C--------------------------------- THREE DIRECTIONS

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1

               H2(K,J,I) = 0.5 * (H1(K  ,J  ,I  )+
     $                            H1(K  ,J+1,I  ))

            ENDDO
         ENDDO
      ENDDO

      CALL INTERPOLATE_3D_CELL3 (KK,JJ,II,H2,RDDX,RDDY,RDDZ,HX,HY,HZ,
     $                          NPART,KINDEX,JINDEX,IINDEX,XP,YP,ZP,
     $                          DUPART(1,3))

C------------------------------------------------------------------------
C---------------------------------- EVALUATION OF DVDX ------------------
C---------------------------------- -------------------------------------

C---------------------------------- COMPUTE DUDY ON 
C                                   X- AND Y-STAGGERED POSITIONS

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1

               H1(K,J,I) = (V(K,J,I+1) - V(K,J,I))*RDX(I)

            ENDDO
         ENDDO
      ENDDO

C--------------------------------- INTERPOLATE DVDX ON CORNERS OF CELLS
C--------------------------------- GRID IS TAKEN TO BE EQUIDISTANT IN ALL
C--------------------------------- THREE DIRECTIONS

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1

               H2(K,J,I) = 0.5 * (H1(K  ,J  ,I  )+
     $                            H1(K+1,J  ,I  ))

            ENDDO
         ENDDO
      ENDDO

      CALL INTERPOLATE_3D_CELL3 (KK,JJ,II,H2,RDDX,RDDY,RDDZ,HX,HY,HZ,
     $                          NPART,KINDEX,JINDEX,IINDEX,XP,YP,ZP,
     $                          DUPART(1,4))

C---------------------------------- EVALUATION OF DVDY

C---------------------------------- COMPUTE DUDX ON CELL-CENTERED POSITIONS

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1

               H1(K,J,I) = (V(K,J,I) - V(K,J-1,I))*RDDY(J)

            ENDDO
         ENDDO
      ENDDO

C-----------------  BOUNDARY CONDITION NOTE: ONLY PERIODIC!!!!
C----------- NOT IN T_MIXER CONFIGURATION

C      DO I=2,II-1
C            DO K=2,KK-1
C
C               H1(K,2,I) = H1(K,JJ-2,I)
C
C            ENDDO
C      ENDDO

C--------------------------------- INTERPOLATE DVDY ON CORNERS OF CELLS
C--------------------------------- GRID IS TAKEN TO BE EQUIDISTANT IN ALL
C--------------------------------- THREE DIRECTIONS

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1

               H2(K,J,I) = 0.125 * (H1(K  ,J  ,I  )+
     $                              H1(K+1,J  ,I  )+
     $                              H1(K  ,J+1,I  )+
     $                              H1(K+1,J+1,I  )+
     $                              H1(K  ,J  ,I+1)+
     $                              H1(K+1,J  ,I+1)+
     $                              H1(K  ,J+1,I+1)+
     $                              H1(K+1,J+1,I+1) )

            ENDDO
         ENDDO
      ENDDO

      CALL INTERPOLATE_3D_CELL3 (KK,JJ,II,H2,RDDX,RDDY,RDDZ,HX,HY,HZ,
     $                          NPART,KINDEX,JINDEX,IINDEX,XP,YP,ZP,
     $                          DUPART(1,5))

C------------------------------------------------------------------------
C---------------------------------- EVALUATION OF DVDZ ------------------
C---------------------------------- -------------------------------------

C---------------------------------- COMPUTE DUDY ON 
C                                   Y- AND Z-STAGGERED POSITIONS

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1

               H1(K,J,I) = (V(K+1,J,I) - V(K,J,I))*RDZ(K)

            ENDDO
         ENDDO
      ENDDO

C--------------------------------- INTERPOLATE DVDZ ON CORNERS OF CELLS
C--------------------------------- GRID IS TAKEN TO BE EQUIDISTANT IN ALL
C--------------------------------- THREE DIRECTIONS

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1

               H2(K,J,I) = 0.5 * (H1(K  ,J  ,I  )+
     $                            H1(K  ,J  ,I+1))

            ENDDO
         ENDDO
      ENDDO

      CALL INTERPOLATE_3D_CELL3 (KK,JJ,II,H2,RDDX,RDDY,RDDZ,HX,HY,HZ,
     $                          NPART,KINDEX,JINDEX,IINDEX,XP,YP,ZP,
     $                          DUPART(1,6))

C------------------------------------------------------------------------
C---------------------------------- EVALUATION OF DWDX ------------------
C---------------------------------- -------------------------------------

C---------------------------------- COMPUTE DWDX ON 
C                                   X- AND Z-STAGGERED POSITIONS

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1

               H1(K,J,I) = (W(K,J,I+1) - W(K,J,I))*RDX(I)

            ENDDO
         ENDDO
      ENDDO

C--------------------------------- INTERPOLATE DWDX ON CORNERS OF CELLS
C--------------------------------- GRID IS TAKEN TO BE EQUIDISTANT IN ALL
C--------------------------------- THREE DIRECTIONS

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1

               H2(K,J,I) = 0.5 * (H1(K  ,J  ,I  )+
     $                            H1(K  ,J+1,I  ))

            ENDDO
         ENDDO
      ENDDO

      CALL INTERPOLATE_3D_CELL3 (KK,JJ,II,H2,RDDX,RDDY,RDDZ,HX,HY,HZ,
     $                          NPART,KINDEX,JINDEX,IINDEX,XP,YP,ZP,
     $                          DUPART(1,7))

C------------------------------------------------------------------------
C---------------------------------- EVALUATION OF DWDY ------------------
C---------------------------------- -------------------------------------

C---------------------------------- COMPUTE DWDY ON 
C                                   X- AND Z-STAGGERED POSITIONS

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1

               H1(K,J,I) = (W(K,J+1,I) - W(K,J,I))*RDY(J)

            ENDDO
         ENDDO
      ENDDO

C--------------------------------- INTERPOLATE DWDY ON CORNERS OF CELLS
C--------------------------------- GRID IS TAKEN TO BE EQUIDISTANT IN ALL
C--------------------------------- THREE DIRECTIONS

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1

               H2(K,J,I) = 0.5 * (H1(K  ,J  ,I  )+
     $                            H1(K  ,J  ,I+1))

            ENDDO
         ENDDO
      ENDDO

      CALL INTERPOLATE_3D_CELL3 (KK,JJ,II,H2,RDDX,RDDY,RDDZ,HX,HY,HZ,
     $                          NPART,KINDEX,JINDEX,IINDEX,XP,YP,ZP,
     $                          DUPART(1,8))

C---------------------------------- EVALUATION OF DWDZ

C---------------------------------- COMPUTE DUDX ON CELL-CENTERED POSITIONS

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1

               H1(K,J,I) = (W(K,J,I) - W(K-1,J,I))*RDDZ(K)

            ENDDO
         ENDDO
      ENDDO

C--------------------------------- INTERPOLATE DWDZ ON CORNERS OF CELLS
C--------------------------------- GRID IS TAKEN TO BE EQUIDISTANT IN ALL
C--------------------------------- THREE DIRECTIONS

      DO I=2,II-1
         DO J=2,JJ-1
            DO K=2,KK-1

               H2(K,J,I) = 0.125 * (H1(K  ,J  ,I  )+
     $                              H1(K+1,J  ,I  )+
     $                              H1(K  ,J+1,I  )+
     $                              H1(K+1,J+1,I  )+
     $                              H1(K  ,J  ,I+1)+
     $                              H1(K+1,J  ,I+1)+
     $                              H1(K  ,J+1,I+1)+
     $                              H1(K+1,J+1,I+1) )

            ENDDO
         ENDDO
      ENDDO

      CALL INTERPOLATE_3D_CELL3 (KK,JJ,II,H2,RDDX,RDDY,RDDZ,HX,HY,HZ,
     $                          NPART,KINDEX,JINDEX,IINDEX,XP,YP,ZP,
     $                          DUPART(1,9))

C---------------------------------- READY

C      write (98,*) dupart(100,1),dupart(100,2),dupart(100,3)
C      write (98,*) dupart(100,4),dupart(100,5),dupart(100,6)
C      write (98,*) dupart(100,7),dupart(100,8),dupart(100,9)


      RETURN
      END

      SUBROUTINE INTERPOLATE_3D_CELL (KK,JJ,II,H2,RDDX,RDDY,RDDZ,
     $                          HX,HY,HZ,
     $                          NPART,KINDEX,JINDEX,IINDEX,XP,YP,ZP,
     $                          DUPART)

      REAL H2(KK,JJ,II)
      REAL RDDX(II),RDDY(JJ),RDDZ(KK),HX(II),HY(JJ),HZ(KK)

      INTEGER KINDEX(NPART),JINDEX(NPART),IINDEX(NPART)

      REAL XP(NPART),YP(NPART),ZP(NPART)
      REAL DUPART(NPART)

      DO IPART = 1,NPART

         K=KINDEX(IPART)
         J=JINDEX(IPART)
         I=IINDEX(IPART)

         DELX = (XP(IPART)-HX(I-1  ))*RDDX(I)
         DELY = (YP(IPART)-HY(J-1  ))*RDDY(J)
         DELZ = (ZP(IPART)-HZ(K-1  ))*RDDZ(K)


         DUDX = ((H2(K  ,J  ,I  )*(    DELX)+
     $            H2(K  ,J  ,I-1)*(1.0-DELX) )*(    DELY)+
     $           (H2(K  ,J-1,I  )*(    DELX)+
     $            H2(K  ,J-1,I-1)*(1.0-DELX) )*(1.0-DELY))*(    DELZ)+
     $          ((H2(K-1,J  ,I  )*(    DELX)+
     $            H2(K-1,J  ,I-1)*(1.0-DELX) )*(    DELY)+
     $           (H2(K-1,J-1,I  )*(    DELX)+
     $            H2(K-1,J-1,I-1)*(1.0-DELX) )*(1.0-DELY))*(1.0-DELZ)

         DUPART(IPART) = DUDX

      ENDDO

      RETURN
      END

