










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
      SUBROUTINE PARTICLE_STRESS (KK,JJ,II,X,Y,Z,DX,DY,DZ,
     $                          RDDX,RDDY,RDDZ,
     $                          TAU11,TAU12,TAU13,
     $                          TAU21,TAU22,TAU23,
     $                          TAU31,TAU32,TAU33,
     $                          SAMPLES,HX,HY,HZ,
     $                          NPART,XPART,YPART,ZPART,
     $                          IINDEX,JINDEX,KINDEX,TAU)
C------------------------------------------------------------
C
C     INTERPOLATES PARTICLE STRESSES ON EULERIAN GRID
C
C     020502 (MM)  ORIGINAL
C
C------------------------------------------------------------

      REAL 
     $     TAU11(KK,JJ,II),TAU12(KK,JJ,II),TAU13(KK,JJ,II),
     $     TAU21(KK,JJ,II),TAU22(KK,JJ,II),TAU23(KK,JJ,II),
     $     TAU31(KK,JJ,II),TAU32(KK,JJ,II),TAU33(KK,JJ,II),
     $     SAMPLES(KK,JJ,II)

      REAL    X(II),    Y(II),    Z(II),
     $       DX(II),   DY(JJ),   DZ(KK),
     $     RDDX(II), RDDY(JJ), RDDZ(KK),
     $       HX(II),   HY(JJ),   HZ(KK)


      REAL XPART(NPART), YPART(NPART), ZPART(NPART)
      REAL TAU(NPART,3,3)
      
      INTEGER IINDEX(NPART), JINDEX(NPART), KINDEX(NPART)

      real nverteilung(0:99)
C------------------------------------------------------------
C
C                                      LOOP OVER ALL PARTICLES
C
C------------------------------------------------------------
      
      DO IPART = 1,NPART

         K=KINDEX(IPART)
         J=JINDEX(IPART)
         I=IINDEX(IPART)

C------------------------------------------------------------
C                     DETERMINATION OF WEIGHTING FACOTRS
C                     WE ASSUME LINEAR WEIGHTING 
C                     WE DON T KNOW, IF XPART IS LEFT OR RIGTH FROM
C                     XP -- THEREFORE, WE TEST BOTH DIRECTIONS AND 
C                     SWITCH BY E.G. (FW = 0.5*(ABS(FW) + FW))
C                     WHICH IS ZERO, IF FW IS NEGATIVE


C                                     X-DIRECTION (WEST -- POINT -- EAST)

         XW = X(I - 1)
         XP = X(I    )
         XE = X(I + 1)

         FW = (XP - XPART(IPART)) / (XP - XW)
         FW = 0.5*(ABS(FW) + FW)

         FE = (XPART(IPART) - XP) / (XE - XP)
         FE = 0.5*(ABS(FE) + FE)

         FPX = 1.0 - FW - FE

C                                     Y-DIRECTION (SOUTH -- POINT -- NORTH)

         YS = Y(J - 1)
         YP = Y(J    )
         YN = Y(J + 1)

         FS = (YP - YPART(IPART)) / (YP - YS)
         FS = 0.5*(ABS(FS) + FS)

         FN = (YPART(IPART) - YP) / (YN - YP)
         FN = 0.5*(ABS(FN) + FN)

         FPY = 1.0 - FS - FN

C                                     K-DIRECTION (BOTTOM -- POINT -- TOP)

         ZB = Z(K - 1)
         ZP = Z(K    )
         ZT = Z(K + 1)

         FB = (ZP - ZPART(IPART)) / (ZP - ZB)
         FB = 0.5*(ABS(FB) + FB)

         FT = (ZPART(IPART) - ZP) / (ZT - ZP)
         FT = 0.5*(ABS(FT) + FT)

         FPZ = 1.0 - FB - FT

C------------------------------------------------------------
C              WEIGHTING STRESS TENSOR TO DIFFERENT POINTS IN SPACE
         

C----  point k-1,j-1,i-1

         WEIGHT = FW*FS*FB

         CALL SWEIGHT (KK,JJ,II,NPART,K-1,J-1,I-1,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k  ,j-1,i-1

         WEIGHT = FW*FS*FPZ

         CALL SWEIGHT (KK,JJ,II,NPART,K  ,J-1,I-1,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k+1,j-1,i-1

         WEIGHT = FW*FS*FB

         CALL SWEIGHT (KK,JJ,II,NPART,K+1,J-1,I-1,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k-1,j  ,i-1

         WEIGHT = FW*FPY*FB

         CALL SWEIGHT (KK,JJ,II,NPART,K-1,J  ,I-1,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k  ,j  ,i-1

         WEIGHT = FW*FPY*FPZ

         CALL SWEIGHT (KK,JJ,II,NPART,K  ,J  ,I-1,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k+1,j  ,i-1

         WEIGHT = FW*FPY*FB

         CALL SWEIGHT (KK,JJ,II,NPART,K+1,J  ,I-1,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k-1,j+1,i-1

         WEIGHT = FW*FN*FB

         CALL SWEIGHT (KK,JJ,II,NPART,K-1,J+1,I-1,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k  ,j+1,i-1

         WEIGHT = FW*FN*FPZ

         CALL SWEIGHT (KK,JJ,II,NPART,K  ,J+1,I-1,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k+1,j+1,i-1

         WEIGHT = FW*FN*FB

         CALL SWEIGHT (KK,JJ,II,NPART,K+1,J+1,I-1,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k-1,j-1,i

         WEIGHT = FPX*FS*FB

         CALL SWEIGHT (KK,JJ,II,NPART,K-1,J-1,I  ,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k  ,j-1,i

         WEIGHT = FPX*FS*FPZ

         CALL SWEIGHT (KK,JJ,II,NPART,K  ,J-1,I  ,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k+1,j-1,i

         WEIGHT = FPX*FS*FB

         CALL SWEIGHT (KK,JJ,II,NPART,K+1,J-1,I  ,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k-1,j  ,i

         WEIGHT = FPX*FPY*FB

         CALL SWEIGHT (KK,JJ,II,NPART,K-1,J  ,I  ,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k  ,j  ,i

         WEIGHT = FPX*FPY*FPZ

         CALL SWEIGHT (KK,JJ,II,NPART,K  ,J  ,I  ,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k+1,j  ,i

         WEIGHT = FPX*FPY*FB

         CALL SWEIGHT (KK,JJ,II,NPART,K+1,J  ,I  ,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k-1,j+1,i

         WEIGHT = FPX*FN*FB

         CALL SWEIGHT (KK,JJ,II,NPART,K-1,J+1,I  ,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k  ,j+1,i

         WEIGHT = FPX*FN*FPZ

         CALL SWEIGHT (KK,JJ,II,NPART,K  ,J+1,I  ,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k+1,j+1,i

         WEIGHT = FPX*FN*FB

         CALL SWEIGHT (KK,JJ,II,NPART,K+1,J+1,I  ,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)


C----  point k-1,j-1,i+1

         WEIGHT = FE*FS*FB

         CALL SWEIGHT (KK,JJ,II,NPART,K-1,J-1,I+1,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k  ,j-1,i+1

         WEIGHT = FE*FS*FPZ

         CALL SWEIGHT (KK,JJ,II,NPART,K  ,J-1,I+1,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k+1,j-1,i+1

         WEIGHT = FE*FS*FB

         CALL SWEIGHT (KK,JJ,II,NPART,K+1,J-1,I+1,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k-1,j  ,i+1

         WEIGHT = FE*FPY*FB

         CALL SWEIGHT (KK,JJ,II,NPART,K-1,J  ,I+1,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k  ,j  ,i+1

         WEIGHT = FE*FPY*FPZ

         CALL SWEIGHT (KK,JJ,II,NPART,K  ,J  ,I+1,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k+1,j  ,i+1

         WEIGHT = FE*FPY*FB

         CALL SWEIGHT (KK,JJ,II,NPART,K+1,J  ,I+1,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k-1,j+1,i+1

         WEIGHT = FE*FN*FB

         CALL SWEIGHT (KK,JJ,II,NPART,K-1,J+1,I+1,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k  ,j+1,i+1

         WEIGHT = FE*FN*FPZ

         CALL SWEIGHT (KK,JJ,II,NPART,K  ,J+1,I+1,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

C----  point k+1,j+1,i+1

         WEIGHT = FE*FN*FB

         CALL SWEIGHT (KK,JJ,II,NPART,K+1,J+1,I+1,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)


C------------------------------------------------------------
C                      LOOP OVER PARTICLES FINISHED

      ENDDO



      RETURN

C                   AB HIER NUR DIAGNOSE
C------------------------------------------------------------
      write (6,*) 'particle_stress,',
     $  'tau13,tau31',tau13(5,5,5),tau31(5,5,5),samples(5,5,5)
      write (6,*) 'particle_stress,',
     $  'tau23,tau32',tau23(5,5,5),tau32(5,5,5),samples(5,5,5)
      write (6,*) 'particle_stress,',
     $  'tau12,tau21',tau12(5,5,5),tau21(5,5,5),samples(5,5,5)

      do i=1,ii
      write (6,1000) 'tau22',i,tau22(5,5,i),samples(5,5,i)
      enddo

 1000  FORMAT (A5,I6,2(1X,E12.5E3,1X))
      RN= 0
      N= 0

      do i=0,99
        nverteilung(i) = 0
      enddo

      DO I=2,II-1
       DO J=2,JJ-1
        DO K=2,KK-1

           RN = RN + SAMPLES(K,J,I)
           N = N + int(SAMPLES(K,J,I))

           iv = int(SAMPLES(K,J,I))
           nverteilung(iv) = nverteilung(iv) + 1.0
C           nverteilung(SAMPLES(K,J,I)) = nverteilung(SAMPLES(K,J,I)) + 1

        ENDDO
       ENDDO
      ENDDO

      write (6,*)'gesamtanzahl in particle_stress',rn,n
      do i=0,99
         write (6,*) 'nvert',i,nverteilung(i)
      enddo
C------------------------------------------------------------


      RETURN
      END

         SUBROUTINE SWEIGHT (KK,JJ,II,NPART,K,J,I,IPART,WEIGHT,
     $        TAU11,TAU12,TAU13,TAU22,TAU23,TAU33,TAU,SAMPLES)

      REAL 
     $     TAU11(KK,JJ,II),TAU12(KK,JJ,II),TAU13(KK,JJ,II),
     $                     TAU22(KK,JJ,II),TAU23(KK,JJ,II),
     $                                     TAU33(KK,JJ,II),
     $     SAMPLES(KK,JJ,II)

      REAL TAU(NPART,3,3)

      TAU11(K,J,I) = TAU11(K,J,I) + TAU(IPART, 1, 1)*WEIGHT
      TAU12(K,J,I) = TAU12(K,J,I) + TAU(IPART, 1, 2)*WEIGHT
      TAU13(K,J,I) = TAU13(K,J,I) + TAU(IPART, 1, 3)*WEIGHT
      
      TAU22(K,J,I) = TAU22(K,J,I) + TAU(IPART, 2, 2)*WEIGHT
      TAU23(K,J,I) = TAU23(K,J,I) + TAU(IPART, 2, 3)*WEIGHT
      
      TAU33(K,J,I) = TAU33(K,J,I) + TAU(IPART, 3, 3)*WEIGHT
      
      SAMPLES(K,J,I) = SAMPLES(K,J,I) + WEIGHT

      RETURN
      END
