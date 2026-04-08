










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
      SUBROUTINE INTERPOLATE_3D_CELL6 (KK,JJ,II,H2,RDDX,RDDY,RDDZ,
     $                          HX,HY,HZ,
     $                          NPART,KINDEX,JINDEX,IINDEX,XP,YP,ZP,
     $                          DUPART,
     $                          FX,FY,FZ)

C-------------------------------------------------------------
C
C         INTERPOLATION NACH LAGRANGESCHER INTERPOLATION 3. ORDNUNG
C
C-------------------------------------------------------------

      REAL H2(KK,JJ,II)
      REAL RDDX(II),RDDY(JJ),RDDZ(KK),HX(II),HY(JJ),HZ(KK)

      REAL FX(NPART,4),FY(NPART,4),FZ(NPART,4)

      INTEGER KINDEX(NPART),JINDEX(NPART),IINDEX(NPART)

      REAL XP(NPART),YP(NPART),ZP(NPART)
      REAL DUPART(NPART)


      DO IPART = 1,NPART
         DUPART(IPART) = 0.0
      ENDDO

      DO ISHIFT = -2,1
         DO JSHIFT = -2,1

            CALL INTERPOLATE_CELLZ (KK,JJ,II,H2,RDDX,RDDY,RDDZ,
     $           HX,HY,HZ,
     $           NPART,KINDEX,JINDEX,IINDEX,
     $           XP,YP,ZP,
     $           DUPART,
     $           FX(1,ISHIFT+3),FY(1,JSHIFT+3),FZ,ISHIFT,JSHIFT)

         ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE INTERPOLATE_CELLZ (KK,JJ,II,H2,RDDX,RDDY,RDDZ,
     $     HX,HY,HZ,
     $     NPART,KINDEX,JINDEX,IINDEX,
     $     XP,YP,ZP,
     $     DUPART,
     $     FX,FY,FZ,ISHIFT,JSHIFT)

      REAL H2(KK,JJ,II)
      REAL RDDX(II),RDDY(JJ),RDDZ(KK),HX(II),HY(JJ),HZ(KK)

      REAL FX(NPART),FY(NPART),FZ(NPART,4)

      INTEGER KINDEX(NPART),JINDEX(NPART),IINDEX(NPART)

      REAL XP(NPART),YP(NPART),ZP(NPART)
      REAL DUPART(NPART)

      DO IPART = 1,NPART
         
         X = XP(IPART)
         Y = YP(IPART)
         Z = ZP(IPART)

         K=KINDEX(IPART)
         J=JINDEX(IPART) + JSHIFT
         I=IINDEX(IPART) + ISHIFT

C------------------------------- FAKTOREN IN Z-RICHTUNG
                  
                 DUPART(IPART) = DUPART(IPART) +
     $               (H2(K-2, J, I)*FZ(IPART,1) +
     $                H2(K-1, J, I)*FZ(IPART,2) +
     $                H2(K-0, J, I)*FZ(IPART,3) +
     $                H2(K+1, J, I)*FZ(IPART,4)  )*FY(IPART)*FX(IPART)


      ENDDO

      RETURN
      END

      SUBROUTINE INT3DFAC4(NPART,IINDEX,XP,II,HX,RDDX,
     $                    FX1,FX2,FX3,FX4)
      
      REAL RDDX(II),HX(II)
      REAL FX1(NPART),FX2(NPART),FX3(NPART),FX4(NPART)

      INTEGER IINDEX(NPART)

      REAL XP(NPART)

      DO IPART = 1,NPART
         
         X = XP(IPART)

         I=IINDEX(IPART)

C         IF ( I .LT.    3 ) then
C            write (6,*) "INT3DFAC4", i, ii, ipart, x
C         endif

C         IF ( I .GT. II-2 ) then
C            write (6,*) "INT3DFAC4", i, ii, ipart, x
C         endif

         FX1(IPART) = ( (X-HX(I-1)) * (X-HX(I  )) * (X-HX(I+1)) )
     $      /((HX(I-2)-HX(I-1)) * (HX(I-2)-HX(I  )) * (HX(I-2)-HX(I+1)))

         FX2(IPART) = ( (X-HX(I-2)) * (X-HX(I  )) * (X-HX(I+1)) )
     $      /((HX(I-1)-HX(I-2)) * (HX(I-1)-HX(I  )) * (HX(I-1)-HX(I+1)))

         FX3(IPART) = ( (X-HX(I-2)) * (X-HX(I-1)) * (X-HX(I+1)) )
     $      /((HX(I  )-HX(I-2)) * (HX(I  )-HX(I-1)) * (HX(I  )-HX(I+1)))

         FX4(IPART) = ( (X-HX(I-2)) * (X-HX(I-1)) * (X-HX(I  )) )
     $      /((HX(I+1)-HX(I-2)) * (HX(I+1)-HX(I-1)) * (HX(I+1)-HX(I  )))

      ENDDO

      RETURN
      END
