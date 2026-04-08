










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
      SUBROUTINE COEFKON  (II,DX,DDX,NFRO,NBAC,COEFFX)
C*MGLET***************************************************************
C        C A L C O E F F          
C        BERECHNUNG DER KOEFFIZIENTEN FUER DIE KONVEKTIVEN TERME 
C        FUER DAS KOMPAKTVERFAHREN VIERTER ORDNUNG IN EINER
C        KOORDINATENRICHTUNG
C*MGLET***************************************************************
C
C PARAM: COEFFX         - FELD FUER DIE KOEFFIZIENTEN IN X-RICHTUNG 
C 
C      ISTAG = 1        - Interpolation auf Druckpunkt in Koordinaten-
C                         richtung
C      : COEFFX(I,1)    - ENTHAELT KOEFFIZIENTEN DER "LINKEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : COEFFX(I,2)    - ENTHAELT KOEFFIZIENTEN DER HAUPTDIAGONALEN
C      : COEFFX(I,3)    - ENTHAELT KOEFFIZIENTEN DER "RECHTEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : COEFFX(I,4)    - ERSTER KOEFFIZIENT AUF DER RECHTEN SEITE DES 
C                         GLEICHUNGSSYSTEMS
C      : COEFFX(I,5)    - ZWEITER KOEFFIZIENT AUF DER RECHTEN SEITE DES 
C                         GLEICHUNGSSYSTEMS
C      : COEFFX(I,6)    - DRITTER KOEFFIZIENT AUF DER RECHTEN SEITE DES 
C                         GLEICHUNGSSYSTEMS
CCC
C      ISTAG = 0        - Interpolation auf Zellecke in andere Koordinaten
C                         richtung 
C      : COEFFX(I,7)    - ENTHAELT KOEFFIZIENTEN DER "LINKEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : COEFFX(I,8)    - ENTHAELT KOEFFIZIENTEN DER HAUPTDIAGONALEN
C      : COEFFX(I,9)    - ENTHAELT KOEFFIZIENTEN DER "RECHTEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : COEFFX(I,10)    -ERSTER KOEFFIZIENT AUF DER RECHTEN SEITE DES 
C                         GLEICHUNGSSYSTEMS
C      : COEFFX(I,11)   - ZWEITER KOEFFIZIENT AUF DER RECHTEN SEITE DES 
C                         GLEICHUNGSSYSTEMS
C      : COEFFX(I,12)   - DRITTER KOEFFIZIENT AUF DER RECHTEN SEITE DES 
C                         GLEICHUNGSSYSTEMS
C
C      : ISTART         - ERSTER PHYSIKALISCHE PUNKT                   
C      : ISTOP          - LETZTER PHYSIKALISCHE PUNKT                   
C
C      : NFRO           - RANDBEDINGUNG AM ANFANG                       
C      : NBAC           - RANDBEDINGUNG AM ENDE                         
C
C DEFINE DIREKTIVEN     : KEINE 
C
C UNTERPROGRAMME        : INTERCOEF1, INTERCOEF4
C                       : 
C
C VERS:  04.04.2003 (FS): ORIGINAL        (KOMPAKT 4. ORDNUNG)
C
C*MGLET***************************************************************

      IMPLICIT NONE
      LOGICAL    STAG
      INTEGER    II,NFRO,NBAC,I


      REAL       DX(II), DDX(II), COEFFX(II,12*3)


C***********************  ISTAG = 1  ***************************
C****** auf Druckpunkt (in Koordinatenrichtung) ****************

C****** Es werden alle Koeffizienten mit normalem Stenzil  *****
C******         berechnet. Bei Periodizitaet werden diese  *****
C****** nicht mehr veraendert

      WRITE(6,*) ''
      WRITE(6,*)'VORBELEGUNG COEFKON STAG'
      WRITE(6,*)''

      DO I = 2,II-1

         CALL INTERCOEF1(II,I,DX,DDX,COEFFX(I,1),COEFFX(I,2),COEFFX(I,3)
     $        ,COEFFX(I,4),COEFFX(I,5),COEFFX(I,6))
C         WRITE (6,6000) I,COEFFX(I,1),COEFFX(I,2),COEFFX(I,3),
C     $        COEFFX(I,4),COEFFX(I,5),COEFFX(I,6)

c      ENDDO


c      DO I = 2,II-1
      STAG = .TRUE.
c      NFRO, also rechter Koerperrand -->12
         CALL RANDWERTAN(II,I,DX,DDX,STAG,
     $        COEFFX(I,12+1),COEFFX(I,12+2),
     $        COEFFX(I,12+3),COEFFX(I,12+4),
     $        COEFFX(I,12+5),COEFFX(I,12+6),NFRO)
c      NBAC, also linker Koerperrand -->24
         CALL RANDWERTW(II,I,DX,DDX,STAG,
     $        COEFFX(I,24+1),COEFFX(I,24+2),
     $        COEFFX(I,24+3),COEFFX(I,24+4),
     $        COEFFX(I,24+5),COEFFX(I,24+6))

      ENDDO

C             ********** No-SLIP(5) FLU (11)  **********

      IF (NFRO .EQ.5 .OR. NFRO .EQ. 11) THEN

         WRITE(6,*) ''
         WRITE(6,*) ' RANDBEDINGUNG NO-SLIP ODER FLU F. PKT. 3'
         WRITE(6,*) ''

         I = 3
         STAG = .TRUE.
         CALL RANDWERTAN(II,I,DX,DDX,STAG,COEFFX(I,1),COEFFX(I,2),
     $        COEFFX(I,3),COEFFX(I,4),COEFFX(I,5),COEFFX(I,6),NFRO)
         WRITE (6,6000) I,COEFFX(I,1),COEFFX(I,2),COEFFX(I,3),
     $        COEFFX(I,4),COEFFX(I,5),COEFFX(I,6)
         I = 2
         COEFFX(I,1) = 0.0
	 COEFFX(I,2) = 1.0
         COEFFX(I,3) = 0.0
	 COEFFX(I,4) = 0.0
         COEFFX(I,5) = 0.0 
	 COEFFX(I,6) = 0.0

      ENDIF

      IF (NBAC .EQ.5) THEN

         WRITE(6,*) ''
         WRITE(6,*) ' RANDBEDINGUNG NO-SLIP FÜR PKT.',II-2

         I = II-2
         STAG = .TRUE.
         CALL RANDWERTW(II,I,DX,DDX,STAG,COEFFX(I,1),COEFFX(I,2),
     $        COEFFX(I,3),COEFFX(I,4),COEFFX(I,5),COEFFX(I,6))
         WRITE (6,6000) I,COEFFX(I,1),COEFFX(I,2),COEFFX(I,3),
     $        COEFFX(I,4),COEFFX(I,5),COEFFX(I,6)
     
         I = II-1
         COEFFX(I,1) = 0.0
	 COEFFX(I,2) = 1.0
         COEFFX(I,3) = 0.0
	 COEFFX(I,4) = 0.0
         COEFFX(I,5) = 0.0 
	 COEFFX(I,6) = 0.0


      ENDIF

C             ********** OP1 (3)    **********

      IF (NFRO .EQ. 3 .OR. NFRO .EQ. 4) THEN
         WRITE(6,*) ''
         WRITE (6,*) ' RANDBEDIINGUNG  OP1/OP2 FÜR PKT. 3'
         WRITE(6,*) ''

         I = 3
         STAG = .TRUE.
      	 CALL RANDWERTAN(II,I,DX,DDX,STAG,COEFFX(I,1),COEFFX(I,2),
     $         COEFFX(I,3),COEFFX(I,4),COEFFX(I,5),COEFFX(I,6),NFRO)
         WRITE (6,6000) I,COEFFX(I,1),COEFFX(I,2),COEFFX(I,3),
     $        COEFFX(I,4),COEFFX(I,5),COEFFX(I,6)
         WRITE(6,*) 'EXTRAPOLATION NOCH NICHT BERECHNET'
         STOP
      ENDIF

      IF (NBAC .EQ. 3 .OR. NBAC .EQ. 4) THEN
          WRITE(6,*) ''
          WRITE (6,*) ' RANDBEDIINGUNG OP1/OP2 FÜR INTERPOLATION',II-2
          WRITE(6,*) ''

          I = II-2
          STAG = .TRUE.
         CALL RANDWERTW(II,I,DX,DDX,STAG,COEFFX(I,1),COEFFX(I,2),
     $         COEFFX(I,3),COEFFX(I,4),COEFFX(I,5),COEFFX(I,6))
         WRITE (6,6000) I,COEFFX(I,1),COEFFX(I,2),COEFFX(I,3),
     $        COEFFX(I,4),COEFFX(I,5),COEFFX(I,6)

         WRITE(6,*) ''
         WRITE (6,*) ' RANDBEDINGUNG  OP1/OP2 FÜR EXTRAPOL.',II-1
         WRITE(6,*) ''

         I = II-1
         STAG = .TRUE.
         CALL RANDWERT(II,I,DX,STAG,COEFFX(I,1),COEFFX(I,2),
     $        COEFFX(I,3),COEFFX(I,4),COEFFX(I,5),COEFFX(I,6))
         WRITE (6,6000) I,COEFFX(I,1),COEFFX(I,2),COEFFX(I,3),
     $        COEFFX(I,4),COEFFX(I,5),COEFFX(I,6)

      ENDIF


C***********************  ISTAG = 0  ***************************
C****** auf Zelleneckpunkt (andere Koordinatenrichtung) ********

      WRITE(6,*)''
      WRITE(6,*)'VORBELEGUNG COEFKON NONSTAG'
      WRITE(6,*)''

      DO I = 2,II-1

         CALL INTERCOEF4(II,I,DX,DDX,COEFFX(I,7),COEFFX(I,8),
     $        COEFFX(I,9),COEFFX(I,10),COEFFX(I,11),COEFFX(I,12))
C         WRITE (6,6000) I,COEFFX(I,7),COEFFX(I,8),COEFFX(I,9),
C     $        COEFFX(I,10),COEFFX(I,11),COEFFX(I,12)

c      ENDDO 


c      DO I = 2,II-1

      STAG = .FALSE.
c      NFRO, also rechter Koerperrand -->12
         CALL RANDWERTAN(II,I,DX,DDX,STAG,
     $        COEFFX(I,12+7),COEFFX(I,12+8),
     $        COEFFX(I,12+9),COEFFX(I,12+10),
     $        COEFFX(I,12+11),COEFFX(I,12+12),NFRO)
c      NBAC, also linker Koerperrand -->24
         CALL RANDWERTW(II,I,DX,DDX,STAG,
     $        COEFFX(I,24+7),COEFFX(I,24+8),
     $        COEFFX(I,24+9),COEFFX(I,24+10),
     $        COEFFX(I,24+11),COEFFX(I,24+12))

      ENDDO

C             ********** No-SLIP(5) **********

      IF (NFRO .EQ.5) THEN

         I = 3
         COEFFX(I,7) = 0.0 
	COEFFX(I,8) = 1.0
         COEFFX(I,9) = 0.0 
	COEFFX(I,10) = 0.0
         COEFFX(I,11) = 0.0 
	COEFFX(I,12) = 0.0
         WRITE (6,6000) I,COEFFX(I,7),COEFFX(I,8),COEFFX(I,9),
     $        COEFFX(I,10),COEFFX(I,11),COEFFX(I,12)

      ENDIF
      IF (NBAC .EQ.5) THEN

         I = II-1
         COEFFX(I,7) = 0.0 
	COEFFX(I,8) = 1.0
         COEFFX(I,9) = 0.0 
	COEFFX(I,10) = 0.0
         COEFFX(I,11) = 0.0
	 COEFFX(I,12) = 0.0  
         WRITE (6,6000) I,COEFFX(I,7),COEFFX(I,8),COEFFX(I,9),
     $        COEFFX(I,10),COEFFX(I,11),COEFFX(I,12)  

      ENDIF

C             ********** OP1 (3)    **********

      IF (NFRO .EQ. 3 .OR. NFRO .EQ. 4) THEN
         WRITE (6,*) 'ERSTER ZU INTERPOLIERENDER PUNKT', 3
         WRITE (6,*) 'RANDBEDINGUNG  OP1/OP2'
         WRITE (6,*) 'SUBROUTINE RANDWERTIJAN VON ADNAN'
         WRITE (6,*) 'FALSCHER STENZIL IN DOKTORARBEIT!'
         I = 3
         CALL RANDWERTIJAN(II,I,DX,DDX,COEFFX(I,7),COEFFX(I,8),
     $        COEFFX(I,9),COEFFX(I,10),COEFFX(I,11),
     $        COEFFX(I,12),NFRO)
        WRITE (6,6000) I,COEFFX(I,7),COEFFX(I,8),COEFFX(I,9),
     $                COEFFX(I,10),COEFFX(I,11),COEFFX(I,12)

         WRITE (6,*) 'KOEFFIZIENTEN FUER EXTRAPOLATION FUER PUNKT' ,2
         WRITE (6,*) 'NOCH NICHT IMPLEMENTIERT'
      ENDIF

      IF (NBAC .EQ. 3 .OR. NBAC .EQ. 4) THEN
         WRITE (6,*) ' LETZTER ZU INTERPOLIERENDER PUNKT', II-2
         WRITE (6,*) ' RANDBEDINGUNG  OP1/OP2' 
         I = II-2
         STAG = .FALSE.
         CALL RANDWERTW(II,I,DX,DDX,STAG,COEFFX(I,7),COEFFX(I,8),
     $        COEFFX(I,9),COEFFX(I,10),COEFFX(I,11),COEFFX(I,12))
         WRITE (6,6000) I,COEFFX(I,7),COEFFX(I,8),COEFFX(I,9),
     $        COEFFX(I,10),COEFFX(I,11),COEFFX(I,12)

         WRITE (6,*) ' ZU EXTRAPOLIERENDER PUNKT', II-1
         WRITE (6,*) ' RANDBEDINGUNG  OP1/OP2'
         I = II-1
         STAG = .FALSE.
         CALL RANDWERT(II,I,DX,STAG,COEFFX(I,7),COEFFX(I,8),
     $        COEFFX(I,9),COEFFX(I,10),COEFFX(I,11),COEFFX(I,12))
        WRITE (6,6000) I,COEFFX(I,7),COEFFX(I,8),COEFFX(I,9),
     $                COEFFX(I,10),COEFFX(I,11),COEFFX(I,12)
      ENDIF
C ----------------------  FLU (11)   ---------------------------
      IF (NFRO .EQ. 11) THEN
         WRITE (6,*) ' ERSTER ZU INTERPOLIERENDER PUNKT 3'
         WRITE (6,*) ' RANDBEDIINGUNG FLU ZELLE 2 GESETZT'
         I = 3
         STAG = .FALSE.
      	 CALL RANDWERTAN(II,I,DX,DDX,STAG,COEFFX(I,7),COEFFX(I,8),
     $        COEFFX(I,9),COEFFX(I,10),COEFFX(I,11),COEFFX(I,12),NFRO)
         WRITE (6,6000) I,COEFFX(I,7),COEFFX(I,8),COEFFX(I,9),
     $        COEFFX(I,10),COEFFX(I,11),COEFFX(I,12)
         I = 2
         COEFFX(I,7) = 0.0 
	 COEFFX(I,8) = 0.0
         COEFFX(I,9) = 0.0
	  COEFFX(I,10) = 0.0
         COEFFX(I,11) = 0.0 
	COEFFX(I,12) = 0.0
      ENDIF

 6000 FORMAT(1X,I3,6(F12.6))

C***********************  E   N   D  ***************************
      RETURN
      END

