










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
      SUBROUTINE CAL_RR(NPART, NFAM, CONF, LOUT, 
     $                  IT, DT, RR)
      REAL CONF(NFAM,NPART, 3), RR(NFAM,9)

      LOGICAL LDOB, LOUT


      DO IFAM = 1,NFAM

      RR_11 = 0.0
      RR_12 = 0.0
      RR_13 = 0.0
      RR_22 = 0.0
      RR_23 = 0.0
      RR_33 = 0.0

      DO IP = 1,NPART

         RR_11 = RR_11 + CONF(IFAM,IP,1) * CONF(IFAM,IP,1)
         RR_12 = RR_12 + CONF(IFAM,IP,1) * CONF(IFAM,IP,2)
         RR_13 = RR_13 + CONF(IFAM,IP,1) * CONF(IFAM,IP,3)
         RR_22 = RR_22 + CONF(IFAM,IP,2) * CONF(IFAM,IP,2)
         RR_23 = RR_23 + CONF(IFAM,IP,2) * CONF(IFAM,IP,3)
         RR_33 = RR_33 + CONF(IFAM,IP,3) * CONF(IFAM,IP,3)

      ENDDO

      FNORM = 1.0/FLOAT(NPART)

      RR_11 = RR_11 * FNORM
      RR_12 = RR_12 * FNORM
      RR_13 = RR_13 * FNORM
      RR_22 = RR_22 * FNORM
      RR_23 = RR_23 * FNORM
      RR_33 = RR_33 * FNORM

      RR (IFAM,1) = RR_11
      RR (IFAM,2) = RR_12
      RR (IFAM,3) = RR_13
      RR (IFAM,4) = RR_12
      RR (IFAM,5) = RR_22
      RR (IFAM,6) = RR_23
      RR (IFAM,7) = RR_13
      RR (IFAM,8) = RR_23
      RR (IFAM,9) = RR_33


      ENDDO

C---------------------------------- output

C---------------------------------- 

      RETURN
 1000 FORMAT (E12.5E3, 6(1X,E12.5E3))
      END
