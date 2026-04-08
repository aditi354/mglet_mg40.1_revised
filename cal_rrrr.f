










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
      SUBROUTINE CAL_RRRR(NPART, NFAM, CONF, LOUT, IT, DT, DUPART,
     $           RRRREE)

      REAL CONF(NFAM, NPART, 3), DUPART(NFAM,9), RRRREE(NFAM,9)
      LOGICAL LOUT


      DO IFAM = 1,NFAM

C---------------------------- RATE-OF-STRAIN TENSOR

         DUDX = DUPART(IFAM,1)
         DUDY = DUPART(IFAM,2)
         DUDZ = DUPART(IFAM,3)

         DVDX = DUPART(IFAM,4)
         DVDY = DUPART(IFAM,5)
         DVDZ = DUPART(IFAM,6)

         DWDX = DUPART(IFAM,7)
         DWDY = DUPART(IFAM,8)
         DWDZ = DUPART(IFAM,9)

         E11 = DUDX
         E12 = 0.5*(DUDY+DVDX)
         E13 = 0.5*(DUDZ+DWDX)

         E21 = 0.5*(DUDY+DVDX)
         E22 = DVDY
         E23 = 0.5*(DVDZ+DWDY)

         E31 = 0.5*(DUDZ+DWDX)
         E32 = 0.5*(DVDZ+DWDY)
         E33 = DWDZ

C----------------------------------------

         RRRR_11 = 0.0
         RRRR_12 = 0.0
         RRRR_13 = 0.0
         RRRR_22 = 0.0
         RRRR_23 = 0.0
         RRRR_33 = 0.0

         DO IP = 1,NPART

            R1 = CONF(IFAM,IP,1)
            R2 = CONF(IFAM,IP,2)
            R3 = CONF(IFAM,IP,3)


            RR_EE =  R1 * R1 * E11
     $             + R1 * R2 * E12
     $             + R1 * R3 * E13
     $             + R2 * R1 * E21
     $             + R2 * R2 * E22
     $             + R2 * R3 * E23
     $             + R3 * R1 * E31
     $             + R3 * R2 * E32
     $             + R3 * R3 * E33

            RRRR_11 = RRRR_11 + R1 * R1 * RR_EE
            RRRR_12 = RRRR_12 + R1 * R2 * RR_EE
            RRRR_13 = RRRR_13 + R1 * R3 * RR_EE
            RRRR_22 = RRRR_22 + R2 * R2 * RR_EE
            RRRR_23 = RRRR_23 + R2 * R3 * RR_EE
            RRRR_33 = RRRR_33 + R3 * R3 * RR_EE

         ENDDO

         FNORM = 1.0/FLOAT(NPART)

         RRRR_11 = RRRR_11 * FNORM
         RRRR_12 = RRRR_12 * FNORM
         RRRR_13 = RRRR_13 * FNORM
         RRRR_22 = RRRR_22 * FNORM
         RRRR_23 = RRRR_23 * FNORM
         RRRR_33 = RRRR_33 * FNORM

         RRRREE (IFAM,1) = RRRR_11
         RRRREE (IFAM,2) = RRRR_12
         RRRREE (IFAM,3) = RRRR_13
         RRRREE (IFAM,4) = RRRR_12
         RRRREE (IFAM,5) = RRRR_22
         RRRREE (IFAM,6) = RRRR_23
         RRRREE (IFAM,7) = RRRR_13
         RRRREE (IFAM,8) = RRRR_23
         RRRREE (IFAM,9) = RRRR_33

c---------------------------------------------- families ready

      ENDDO



      RETURN
 1000 FORMAT (E12.5E3, 6(1X,E12.5E3))
      END
