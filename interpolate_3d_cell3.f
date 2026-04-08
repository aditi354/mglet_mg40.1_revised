










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
      SUBROUTINE INTERPOLATE_3D_CELL3 (KK,JJ,II,H2,RDDX,RDDY,RDDZ,
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

      REAL FX1,FX2,FX3,FX4
      REAL FY1,FY2,FY3,FY4
      REAL FZ1,FZ2,FZ3,FZ4

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

C------------------------------- FAKTOREN IN X-RICHTUNG

         FX1 = ( (X-HX(I-1)) * (X-HX(I  )) * (X-HX(I+1)) )
     $      /((HX(I-2)-HX(I-1)) * (HX(I-2)-HX(I  )) * (HX(I-2)-HX(I+1)))

         FX2 = ( (X-HX(I-2)) * (X-HX(I  )) * (X-HX(I+1)) )
     $      /((HX(I-1)-HX(I-2)) * (HX(I-1)-HX(I  )) * (HX(I-1)-HX(I+1)))

         FX3 = ( (X-HX(I-2)) * (X-HX(I-1)) * (X-HX(I+1)) )
     $      /((HX(I  )-HX(I-2)) * (HX(I  )-HX(I-1)) * (HX(I  )-HX(I+1)))

         FX4 = ( (X-HX(I-2)) * (X-HX(I-1)) * (X-HX(I  )) )
     $      /((HX(I+1)-HX(I-2)) * (HX(I+1)-HX(I-1)) * (HX(I+1)-HX(I  )))

C------------------------------- FAKTOREN IN Y-RICHTUNG

         FY1 = ( (Y-HY(J-1)) * (Y-HY(J  )) * (Y-HY(J+1)) )
     $      /((HY(J-2)-HY(J-1)) * (HY(J-2)-HY(J  )) * (HY(J-2)-HY(J+1)))

         FY2 = ( (Y-HY(J-2)) * (Y-HY(J  )) * (Y-HY(J+1)) )
     $      /((HY(J-1)-HY(J-2)) * (HY(J-1)-HY(J  )) * (HY(J-1)-HY(J+1)))

         FY3 = ( (Y-HY(J-2)) * (Y-HY(J-1)) * (Y-HY(J+1)) )
     $      /((HY(J  )-HY(J-2)) * (HY(J  )-HY(J-1)) * (HY(J  )-HY(J+1)))

         FY4 = ( (Y-HY(J-2)) * (Y-HY(J-1)) * (Y-HY(J  )) )
     $      /((HY(J+1)-HY(J-2)) * (HY(J+1)-HY(J-1)) * (HY(J+1)-HY(J  )))


C------------------------------- FAKTOREN IN Z-RICHTUNG

         FZ1 = ( (Z-HZ(K-1)) * (Z-HZ(K  )) * (Z-HZ(K+1)) )
     $      /((HZ(K-2)-HZ(K-1)) * (HZ(K-2)-HZ(K  )) * (HZ(K-2)-HZ(K+1)))

         FZ2 = ( (Z-HZ(K-2)) * (Z-HZ(K  )) * (Z-HZ(K+1)) )
     $      /((HZ(K-1)-HZ(K-2)) * (HZ(K-1)-HZ(K  )) * (HZ(K-1)-HZ(K+1)))

         FZ3 = ( (Z-HZ(K-2)) * (Z-HZ(K-1)) * (Z-HZ(K+1)) )
     $      /((HZ(K  )-HZ(K-2)) * (HZ(K  )-HZ(K-1)) * (HZ(K  )-HZ(K+1)))

         FZ4 = ( (Z-HZ(K-2)) * (Z-HZ(K-1)) * (Z-HZ(K  )) )
     $      /((HZ(K+1)-HZ(K-2)) * (HZ(K+1)-HZ(K-1)) * (HZ(K+1)-HZ(K  )))


C         DUDX = 0.0

C         DO IINT=1,4
C            DO JINT=1,4
C               DO KINT=1,4
C----------------------------------------------------------- X1:                                 
                  DX1Y1 = 
     $               (H2(K-2, J - 2,  I - 2)*FZ1 +
     $                H2(K-1, J - 2,  I - 2)*FZ2 +
     $                H2(K-0, J - 2,  I - 2)*FZ3 +
     $                H2(K+1, J - 2,  I - 2)*FZ4  )*FY1

                  DX1Y2 = 
     $               (H2(K-2, J - 1,  I - 2)*FZ1 +
     $                H2(K-1, J - 1,  I - 2)*FZ2 +
     $                H2(K-0, J - 1,  I - 2)*FZ3 +
     $                H2(K+1, J - 1,  I - 2)*FZ4  )*FY2


                  DX1Y3 = 
     $               (H2(K-2, J - 0,  I - 2)*FZ1 +
     $                H2(K-1, J - 0,  I - 2)*FZ2 +
     $                H2(K-0, J - 0,  I - 2)*FZ3 +
     $                H2(K+1, J - 0,  I - 2)*FZ4  )*FY3


                  DX1Y4 = 
     $               (H2(K-2, J + 1,  I - 2)*FZ1 +
     $                H2(K-1, J + 1,  I - 2)*FZ2 +
     $                H2(K-0, J + 1,  I - 2)*FZ3 +
     $                H2(K+1, J + 1,  I - 2)*FZ4  )*FY4
C----------------------------------------------------------- X2:                  
                  DX2Y1 = 
     $               (H2(K-2, J - 2,  I - 1)*FZ1 +
     $                H2(K-1, J - 2,  I - 1)*FZ2 +
     $                H2(K-0, J - 2,  I - 1)*FZ3 +
     $                H2(K+1, J - 2,  I - 1)*FZ4  )*FY1

                  DX2Y2 = 
     $               (H2(K-2, J - 1,  I - 1)*FZ1 +
     $                H2(K-1, J - 1,  I - 1)*FZ2 +
     $                H2(K-0, J - 1,  I - 1)*FZ3 +
     $                H2(K+1, J - 1,  I - 1)*FZ4  )*FY2

                  DX2Y3 = 
     $               (H2(K-2, J - 0,  I - 1)*FZ1 +
     $                H2(K-1, J - 0,  I - 1)*FZ2 +
     $                H2(K-0, J - 0,  I - 1)*FZ3 +
     $                H2(K+1, J - 0,  I - 1)*FZ4  )*FY3

                  DX2Y4 = 
     $               (H2(K-2, J + 1,  I - 1)*FZ1 +
     $                H2(K-1, J + 1,  I - 1)*FZ2 +
     $                H2(K-0, J + 1,  I - 1)*FZ3 +
     $                H2(K+1, J + 1,  I - 1)*FZ4  )*FY4
C----------------------------------------------------------- X3:                                 
                  DX3Y1 = 
     $               (H2(K-2, J - 2,  I - 0)*FZ1 +
     $                H2(K-1, J - 2,  I - 0)*FZ2 +
     $                H2(K-0, J - 2,  I - 0)*FZ3 +
     $                H2(K+1, J - 2,  I - 0)*FZ4  )*FY1

                  DX3Y2 = 
     $               (H2(K-2, J - 1,  I - 0)*FZ1 +
     $                H2(K-1, J - 1,  I - 0)*FZ2 +
     $                H2(K-0, J - 1,  I - 0)*FZ3 +
     $                H2(K+1, J - 1,  I - 0)*FZ4  )*FY2

                  DX3Y3 = 
     $               (H2(K-2, J - 0,  I - 0)*FZ1 +
     $                H2(K-1, J - 0,  I - 0)*FZ2 +
     $                H2(K-0, J - 0,  I - 0)*FZ3 +
     $                H2(K+1, J - 0,  I - 0)*FZ4  )*FY3

                  DX3Y4 = 
     $               (H2(K-2, J + 1,  I - 0)*FZ1 +
     $                H2(K-1, J + 1,  I - 0)*FZ2 +
     $                H2(K-0, J + 1,  I - 0)*FZ3 +
     $                H2(K+1, J + 1,  I - 0)*FZ4  )*FY4
C----------------------------------------------------------- X4:                                 
                  DX4Y1 = 
     $               (H2(K-2, J - 2,  I + 1)*FZ1 +
     $                H2(K-1, J - 2,  I + 1)*FZ2 +
     $                H2(K-0, J - 2,  I + 1)*FZ3 +
     $                H2(K+1, J - 2,  I + 1)*FZ4  )*FY1

                  DX4Y2 = 
     $               (H2(K-2, J - 1,  I + 1)*FZ1 +
     $                H2(K-1, J - 1,  I + 1)*FZ2 +
     $                H2(K-0, J - 1,  I + 1)*FZ3 +
     $                H2(K+1, J - 1,  I + 1)*FZ4  )*FY2


                  DX4Y3 = 
     $               (H2(K-2, J - 0,  I + 1)*FZ1 +
     $                H2(K-1, J - 0,  I + 1)*FZ2 +
     $                H2(K-0, J - 0,  I + 1)*FZ3 +
     $                H2(K+1, J - 0,  I + 1)*FZ4  )*FY3


                  DX4Y4 = 
     $               (H2(K-2, J + 1,  I + 1)*FZ1 +
     $                H2(K-1, J + 1,  I + 1)*FZ2 +
     $                H2(K-0, J + 1,  I + 1)*FZ3 +
     $                H2(K+1, J + 1,  I + 1)*FZ4  )*FY4




         DUPART(IPART) = (DX1Y1 + DX1Y2 + DX1Y3 + DX1Y4)*FX1 +
     $                   (DX2Y1 + DX2Y2 + DX2Y3 + DX2Y4)*FX2 +
     $                   (DX3Y1 + DX3Y2 + DX3Y3 + DX3Y4)*FX3 +
     $                   (DX4Y1 + DX4Y2 + DX4Y3 + DX4Y4)*FX4

      ENDDO

      RETURN
      END

