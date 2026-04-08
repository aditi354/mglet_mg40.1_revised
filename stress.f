










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
      SUBROUTINE STRESS (NFAM, RR, RRRREE, DUPART, TAU, IT, DT, 
     $                   RMU0,RMU1,RMU2,RMU3,RMU4,GMOL,VOL_FRAC,LOUT)

      REAL DUPART(NFAM,3,3), RR(NFAM,3,3), RRRREE(NFAM,3,3) 


      REAL DRRRR(3,3), DRR(3,3), RRD(3,3), D(3,3)

      REAL TAU(NFAM,3,3)

      REAL RMU0,RMU1,RMU2,RMU3,RMU4

      LOGICAL LOUT

C      write (6,*)'stress',nfam,it,dt
C      write (6,*)'stress',RMU0,RMU1,RMU2,RMU3,RMU4,gmol


C------------------------------------------- LOOP OVER ALL FAMILIES
      DO IFAM = 1,NFAM
C      DO IFAM = 1,1

C------------------------------------------- SYMMETRIC PART OF GRADIENT
C                                            TENSOR

      DO J = 1,3
         DO I = 1,3
            D(I,J) = 0.5*(DUPART(IFAM,I,J) + DUPART(IFAM,J,I))
         ENDDO
      ENDDO

C---------------------------------------- D:<RR>

      D_RR = 0.0

      DO J = 1,3
         DO I = 1,3
            D_RR = D_RR + D(I,J)*RR(IFAM,I,J)
         ENDDO
      ENDDO


C---------------------------------------- D:<RRRR> = RRRREE

C---------------------------------------- D*<RR> UND <RR>*D

      DO J=1,3
         DO I=1,3
            DRR(I,J) = D(I,1)*RR(IFAM,1,J) + 
     $                 D(I,2)*RR(IFAM,2,J) + 
     $                 D(I,3)*RR(IFAM,3,J)
            RRD(I,J) = RR(IFAM,I,1)*D(1,J) +  
     $                 RR(IFAM,I,2)*D(2,J) +  
     $                 RR(IFAM,I,3)*D(3,J)
         ENDDO
      ENDDO



C---------------------------------------- SPANNUNGSTENSOR

      DO J = 1,3
         DO I = 1,3

            TAU(IFAM,I,J) =   2.0*RMU0 * D(I,J) 
     $                 + RMU2 * RRRREE(IFAM,I,J)
     $                 + 2.0*RMU3 * (DRR(I,J) + RRD(I,J))
     $                 + RMU4 * (3.0*RR(IFAM,I,J))

C            write (6,*) 'i,j,D(i,j)',i,j,d(i,j)
C            write (6,*) 'i,j,RRRREE(i,j)',i,j,RRRREE(i,j)
C            write (6,*) 'i,j,DRR(i,j),RRD',i,j,DRR(i,j),RRD(i,j)
C            write (6,*) 'i,j,RR(i,j)',i,j,RR(i,j)
         ENDDO
      ENDDO
C            write (86,*) 'mu_0,...,D_RR',RMU0,RMU1,RMU2,RMU3,RMU4,D_RR
C      goto 1

      DO I = 1,3
         TAU(IFAM,I,I) = TAU(IFAM,I,I) + RMU1 * D_RR
     $                       - RMU4
      ENDDO

C1     continue
C--------------------------------------- FAMILLIES FINISHED
      ENDDO

C                               SCALAR VERSION FINISHED


      IFAM = 100

C      IF (LOUT) THEN
         WRITE (31,1000) IT*DT,
     $        TAU(IFAM,1,1),TAU(IFAM,1,2),TAU(IFAM,1,3),
     $        TAU(IFAM,2,2),TAU(IFAM,2,3),TAU(IFAM,3,3)
C      ENDIF


      RETURN
 1000 FORMAT (E12.5E3, 6(1X,E12.5E3))
      END
