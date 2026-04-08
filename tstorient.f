










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
      SUBROUTINE TSTORIENT (NPART,NFAM,DUPART,CONF,DT,SHAPE_PART,
     $     IRANDOM,FRANDOM,RANDOM)
C------------------------------------------------------------
C
C     COMPUTES NEW ORIENTATION OF PARTICLES BY 
C
C     OP_NEW = OP_OLD + DT*(OMEGA*N + SHAPE_PART(D*N - (N*D*N)N))
C
C     030999 (MM)  ORIGINAL
C
C------------------------------------------------------------


      REAL CONF(NFAM,NPART,3)

      REAL DUPART(NFAM,9)

      REAL RANDOM(3,NPART),RL

      SQ12DT = SQRT(12.0*DT)

C---------------------------------- TIME STEP


C---------------------------------- all particles
         
      DO IPART = 1,NPART

C---------------------------------- all families

         DO IFAM = 1,NFAM


C------------------------------------ Gradient tensor         
         DUDX = DUPART(IFAM,1)
         DUDY = DUPART(IFAM,2)
         DUDZ = DUPART(IFAM,3)
         
         DVDX = DUPART(IFAM,4)
         DVDY = DUPART(IFAM,5)
         DVDZ = DUPART(IFAM,6)

         DWDX = DUPART(IFAM,7)
         DWDY = DUPART(IFAM,8)
         DWDZ = DUPART(IFAM,9)

C------------------------------------ Vorticity
         OMX = 0.5*(DWDY-DVDZ)
         OMY = 0.5*(DUDZ-DWDX)
         OMZ = 0.5*(DVDX-DUDY)


C------------------------------------ Orientation

            R1 = CONF(IFAM,IPART,1)
            R2 = CONF(IFAM,IPART,2)
            R3 = CONF(IFAM,IPART,3)

C------------------------------------  Vorticity term OMEGA*N

            OM1 =   - OMZ*R2 + OMY*R3
            OM2 =     OMZ*R1 - OMX*R3
            OM3 =   - OMY*R1 + OMX*R2

C-------------------------------------------- D*n

         DR1 =                           DUDX  * R1 +
     $                         0.5*(DUDY+DVDX) * R2 + 
     $                         0.5*(DUDZ+DWDX) * R3  

         DR2 =                 0.5*(DUDY+DVDX) * R1 +
     $                                   DVDY  * R2 + 
     $                         0.5*(DVDZ+DWDY) * R3  

         DR3 =                 0.5*(DUDZ+DWDX) * R1 +
     $                         0.5*(DVDZ+DWDY) * R2 + 
     $                                   DWDZ  * R3  

C-------------------------------------------- n*D*n
         RDR = R1*DR1 + R2*DR2 + R3*DR3

C-------------------------------------------- (n*D*n)n

         DRRR1 = RDR*R1
         DRRR2 = RDR*R2
         DRRR3 = RDR*R3

C-------------------------------------------- time step

         R1 = R1 + DT*(OM1 + SHAPE_PART*(DR1 + DRRR1)) + RANDOM(1,IPART)

         R2 = R2 + DT*(OM2 + SHAPE_PART*(DR2 + DRRR2)) + RANDOM(2,IPART)

         R3 = R3 + DT*(OM3 + SHAPE_PART*(DR3 + DRRR3)) + RANDOM(3,IPART)

C----------------------------- Normierung

C         rL = SQRT(R1*R1 + R2*R2 + R3*R3)
         RL = 1.0/SQRT(R1*R1 + R2*R2 + R3*R3)
C         rl = 1.0

         CONF(IFAM,IPART,1) = R1*RL
         CONF(IFAM,IPART,2) = R2*RL
         CONF(IFAM,IPART,3) = R3*RL

C---------------------------------- PARTICLES READY
      ENDDO
C---------------------------------- FAMILIES READY
      ENDDO

C---------------------------------- READY

      RETURN
      END
