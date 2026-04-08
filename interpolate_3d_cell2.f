










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
      SUBROUTINE INTERPOLATE_3D_CELL2 (KK,JJ,II,H2,RDDX,RDDY,RDDZ,
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
C               DO KINT=1,4
                  
                  DUDX = DUDX  
     $               + H2(K - 2, J + JINT - 3, I + IINT - 3)*
     $                 FX(IINT)*FY(JINT)*FZ(1) 
     $               + H2(K - 1, J + JINT - 3, I + IINT - 3)*
     $                 FX(IINT)*FY(JINT)*FZ(2) 
     $               + H2(K - 0, J + JINT - 3, I + IINT - 3)*
     $                 FX(IINT)*FY(JINT)*FZ(3) 
     $               + H2(K + 1, J + JINT - 3, I + IINT - 3)*
     $                 FX(IINT)*FY(JINT)*FZ(4) 

C               ENDDO
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

