










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
      SUBROUTINE COEFSCAL  (II,DX,DDX,COEFTX,NFRO,NBAC)
C*********************************************************************
C                       : 
C
C VERS:  24.11.2003 (FS): ORIGINAL        (KOMPAKT 4. ORDNUNG)
C Es wird davon ausgegangen dass der 2te Ghostcell Punkt (I=2)
C bei allen Ranbedingungen gesetzt wurde und der erste zu berechnende
C Punkt bei 2+1/2 liegt (d.h. direkt auf dem Rand) - Daher aufruf von
C swclesca3d (Korrektur der Randwerte in tstsca4 VOR jeglicher berechnung
C von Interpolationen oder Ableitungen). Dies ist wichtig fuer Fluss-
C Randbedingung  
C
C        13.04.04    FS : Erweiterung Fred Body
C
C*MGLET***************************************************************
      IMPLICIT NONE
      LOGICAL    STAG
      INTEGER    II,NFRO,NBAC,I
      REAL       DX(II), DDX(II), COEFTX(II,12*3)
C----- INTERPOLATION --------------------------------------
      DO I = 2,II-1
         CALL INTERCOEFKOBS0(II,I,DX,DDX,COEFTX(I,1),COEFTX(I,2),
     $        COEFTX(I,3),COEFTX(I,4),COEFTX(I,5),COEFTX(I,6))
      STAG = .FALSE.
c      NFRO, also rechter Koerperrand -->12
         CALL RANDWERTAN(II,I,DX,DDX,STAG,
     $        COEFTX(I,12+1),COEFTX(I,12+2),
     $        COEFTX(I,12+3),COEFTX(I,12+4),
     $        COEFTX(I,12+5),COEFTX(I,12+6),NFRO)
c      NBAC, also linker Koerperrand -->24
         CALL RANDWERTW(II,I,DX,DDX,STAG,
     $        COEFTX(I,24+1),COEFTX(I,24+2),
     $        COEFTX(I,24+3),COEFTX(I,24+4),
     $        COEFTX(I,24+5),COEFTX(I,24+6))
      ENDDO

      IF(NFRO .EQ. 16 .OR. NFRO .EQ. 15) THEN
         I = 3
         WRITE(6,*)'RANDBEDINGUNG INTERPOLATION, SCALAR FIXED VALUE:'
         COEFTX(I,1) = 0.0
	 COEFTX(I,2) = 1.0
	 COEFTX(I,3) = 0.0
         COEFTX(I,4) = 0.0 
	 COEFTX(I,5) = 1.0 
	 COEFTX(I,6) = 0.0
         WRITE(6,6000)I,COEFTX(I,1),COEFTX(I,2),COEFTX(I,3),
     $                  COEFTX(I,4),COEFTX(I,5),COEFTX(I,6)
      ENDIF
      IF(NBAC .EQ. 16 .OR. NBAC .EQ. 15) THEN
         WRITE(6,*)'RANDBEDINGUNG INTERPOLATION, SCALAR FIXED VALUE:'
         I = II-1
         COEFTX(I,1) = 0.0
	 COEFTX(I,2) = 1.0
	 COEFTX(I,3) = 0.0
         COEFTX(I,4) = 0.0
	 COEFTX(I,5) = 1.0
	 COEFTX(I,6) = 0.0
         WRITE(6,6000)I,COEFTX(I,1),COEFTX(I,2),COEFTX(I,3),
     $                  COEFTX(I,4),COEFTX(I,5),COEFTX(I,6)
      ENDIF

C----- ABLEITUNG ------------------------------------------

      DO I = 2,II-1

	 CALL DIFFCOEFKOBS0(II,I,DX,DDX,COEFTX(I,7),COEFTX(I,8),
     $        COEFTX(I,9),COEFTX(I,10),COEFTX(I,11),COEFTX(I,12))

      STAG = .FALSE.
c      NFRO, also rechter Koerperrand -->12
         CALL RANDWERTAND(II,I,DX,STAG,
     $        COEFTX(I,12+7),COEFTX(I,12+8),
     $        COEFTX(I,12+9),COEFTX(I,12+10),
     $        COEFTX(I,12+11),COEFTX(I,12+12))
c      NBAC, also linker Koerperrand -->24
         CALL RANDWERTDW (II,I,DX,DDX,STAG,
     $        COEFTX(I,24+7),COEFTX(I,24+8),
     $        COEFTX(I,24+9),COEFTX(I,24+10),
     $        COEFTX(I,24+11),COEFTX(I,24+12))

      ENDDO

      IF(NFRO .EQ. 16) THEN
      WRITE(6,*)'RANDBEDINGUNG ABLETIUNG, SCALAR FIXED VALUE:'
      I = 3
      CALL WALLKOBDIFFS0A(II,I,DX,DDX,COEFTX(I,7),COEFTX(I,8),
     $     COEFTX(I,9),COEFTX(I,10),COEFTX(I,11),COEFTX(I,12))
      WRITE(6,6000)I,COEFTX(I,7),COEFTX(I,8),COEFTX(I,9),
     $     COEFTX(I,10),COEFTX(I,11),COEFTX(I,12)
      ENDIF

      IF(NBAC .EQ. 16) THEN
      WRITE(6,*)'RANDBEDINGUNG ABLETIUNG, SCALAR FIXED VALUE:'
      I = II-1
      STAG = .FALSE.
      CALL WALLKOBDIFFS0E(II,I,DX,DDX,COEFTX(I,7),COEFTX(I,8),
     $     COEFTX(I,9),COEFTX(I,10),COEFTX(I,11),COEFTX(I,12))
      WRITE(6,6000)I,COEFTX(I,7),COEFTX(I,8),COEFTX(I,9),
     $     COEFTX(I,10),COEFTX(I,11),COEFTX(I,12)
      ENDIF

      IF(NFRO .EQ. 15) THEN
         I = 3
         WRITE(6,*)'RANDBEDINGUNG ABLEITUNG, SCALAR FIXED GRADIENT:'
         COEFTX(I,7) = 0.0
	 COEFTX(I,8) = 1.0
	 COEFTX(I,9) = 0.0
         COEFTX(I,10) = 0.0
	 COEFTX(I,11) = 1.0
	 COEFTX(I,12) = 0.0
         WRITE(6,6000)I,COEFTX(I,7),COEFTX(I,8),COEFTX(I,9),
     $                  COEFTX(I,10),COEFTX(I,11),COEFTX(I,12)
      ENDIF
      IF(NBAC .EQ. 15) THEN
         WRITE(6,*)'RANDBEDINGUNG ABLEITUNG, SCALAR FIXED GRADIENT:'
         I = II-1
         COEFTX(I,7) = 0.0
	 COEFTX(I,8) = 1.0
	 COEFTX(I,9) = 0.0
         COEFTX(I,10) = 0.0
	 COEFTX(I,11) = 1.0
	 COEFTX(I,12) = 0.0
         WRITE(6,6000)I,COEFTX(I,7),COEFTX(I,8),COEFTX(I,9),
     $                  COEFTX(I,10),COEFTX(I,11),COEFTX(I,12)
      ENDIF

 6000 FORMAT(1X,I3,6(F12.6))
      RETURN
      END
