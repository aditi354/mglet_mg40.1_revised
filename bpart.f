










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
      SUBROUTINE BPART (KK,JJ,II,X,Y,Z,DX,DY,DZ,
     $                   NPART,XP,YP,ZP,
     $                   NFRO,NBAC,NRGT,NLFT,NBOT,NTOP
     $  )
C------------------------------------------------------------
C
C     TREATS BOUNDARY CONDITIONS FOR PARTICLES
C
C     020999 (MM)  ORIGINAL
C
C------------------------------------------------------------

      REAL    X(II),    Y(II),    Z(II),
     $       DX(II),   DY(JJ),   DZ(KK)

      REAL XP(NPART), YP(NPART), ZP(NPART)
C--------------------------------------- BOUNDARIES OF GRID

      X1 = X(2) + 0.5*DX(2)
      X2 = X(II-2) + 0.5*DX(II-2)

      Y1 = Y(2) + 0.5*DY(2)
      Y2 = Y(JJ-2) + 0.5*DY(JJ-2)

      Z1 = Z(2) + 0.5*DZ(2)
      Z2 = Z(KK-2) + 0.5*DZ(KK-2)

C---------------------------------------- TOP BOUNDARY

      IF (NTOP .EQ. 5) THEN
C                                   NOSLIP: MIRRORING
         DO IPART = 1,NPART
         
            IF (ZP(IPART) .GT. Z2) THEN
               ZP(IPART) = 2.0*Z2 - ZP(IPART)
            ENDIF
            
         ENDDO
      ELSE
         
         CALL ERRR (501,' BPART')

      ENDIF

C---------------------------------------- BOTTOM BOUNDARY

      IF (NBOT .EQ. 5) THEN
C                                   NOSLIP: MIRRORING
         DO IPART = 1,NPART
         
            IF (ZP(IPART) .LT. Z1) THEN
               ZP(IPART) = 2.0*Z1 - ZP(IPART)
            ENDIF
            
         ENDDO
         
      ELSE
         
         CALL ERRR (502,' BPART')

      ENDIF

C---------------------------------------- RIGHT BOUNDARY

      IF (NRGT .EQ. 1) THEN
C                                   PERIODIC
         DO IPART = 1,NPART
         
            IF (YP(IPART) .LT. Y1) THEN
               YP(IPART) = Y2 - (Y1 - YP(IPART))
            ENDIF
            
         ENDDO
         
      ELSEIF (NRGT .EQ. 5) THEN
C                                   NOSLIP: MIRRORING
         DO IPART = 1,NPART
         
            IF (YP(IPART) .LT. Y1) THEN
               YP(IPART) = 2.0*Y1 - YP(IPART)
            ENDIF
            
         ENDDO

      ELSE
         
         CALL ERRR (503,' BPART')

      ENDIF

C---------------------------------------- LEFT BOUNDARY

      IF (NLFT .EQ. 1) THEN
C                                   PERIODIC
         DO IPART = 1,NPART
         
            IF (YP(IPART) .GT. Y2) THEN
               YP(IPART) = Y1 + (YP(IPART) - Y2)
            ENDIF
            
         ENDDO
C                                   MIRRORING: NOSLIP
      ELSEIF(NLFT .EQ. 5) THEN
          
         DO IPART = 1,NPART
         
            IF (YP(IPART) .GT. Y2) THEN
               YP(IPART) = 2.0*Y2 - YP(IPART)
            ENDIF
            
         ENDDO
         
      ELSE
         
         CALL ERRR (504,' BPART')

      ENDIF
C---------------------------------------- FRONT BOUNDARY

      IF (NFRO .EQ. 1) THEN
C                                   PERIODIC
         DO IPART = 1,NPART
         
            IF (XP(IPART) .LT. X1) THEN
               XP(IPART) = X2 - (X1 - XP(IPART))
            ENDIF
            
         ENDDO

      ELSEIF (NFRO .EQ. 5) THEN
C                                   NOSLIP: MIRRORING
         DO IPART = 1,NPART
         
            IF (XP(IPART) .LT. X1) THEN
               XP(IPART) = 2.0*X1 - XP(IPART)
            ENDIF
            
         ENDDO
         
      ELSE
         
         CALL ERRR (505,' BPART')

      ENDIF
C---------------------------------------- BACK BOUNDARY

      IF (NBAC .EQ. 1) THEN
C                                   PERIODIC
         DO IPART = 1,NPART
         
            IF (XP(IPART) .GT. X2) THEN
               XP(IPART) = X1 + (XP(IPART) - X2)
            ENDIF
            
         ENDDO
C                      OP1: LEAVING DOMAIN NEW INITALISED 
C                          NOSLIP: MIRRORING         
      ELSEIF(NBAC .EQ. 5) THEN 
         DO IPART = 1,NPART
            IF (XP(IPART) .GT. X2) THEN
               XP(IPART) = 2*X2 - XP(IPART)
            ENDIF
         ENDDO
      ELSE
         
         CALL ERRR (506,' BPART')

      ENDIF
      RETURN
      END
