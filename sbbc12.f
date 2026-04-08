










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
       SUBROUTINE SBBC12 (KK,JJ,II,KMX,JMX,IMX,
     $                  B,KC1,KC2,JC1,JC2,IC1,IC2,
     $                  NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

C*STARLET***************************************************************
C        S B B C 1 2      BELEGEN DER GRENZEN DER 
C                         BOUNDING-BOX DES KUBUS
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        B(KK,JJ,II)    + WANDABSTANDSFELD ENTHAELT IM BERECHNUNGSGEBIET
C                         DEN MISCHUNGSWEG
C        IC1,JC1,KC1    + LINKE  RAENDER DER BOUNDING BOX
C        IC2,JC2,KC2    + RECHTE RAENDER DER BOUNDING BOX
C
C UPROG                 : ERRR
C
C VERS:   4. 4.92 (MM)  : ORIGINAL
C         8. 3.95 (MM)  : ES WERDEN AUCH DIE RANDSCHICHTEN BERUECKSICHTIGT
C
C
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

      COMMON /CONLES/  CONV2S,CMUE,CAPPA,ECONST
      SAVE   /CONLES/
C
C
      REAL 
     $      B(KK,JJ,II)
C
C
      KM1    = KMX-1
      JM1    = JMX-1
      IM1    = IMX-1
      KM2    = KMX-2
      JM2    = JMX-2
      IM2    = IMX-2
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

      IC1 = IMX
      JC1 = JMX
      KC1 = KMX

      IC2 = 0
      JC2 = 0
      KC2 = 0

      ISTART = 2
      ISTOP  = IM1
      JSTART = 2
      JSTOP  = JM1
      KSTART = 2
      KSTOP  = KM1

      IF (NFRO .EQ. 5 .OR. NFRO .EQ. 6) ISTART = 3
      IF (NBAC .EQ. 5 .OR. NBAC .EQ. 6) ISTOP  = IM2
      IF (NRGT .EQ. 5 .OR. NRGT .EQ. 6) JSTART = 3
      IF (NLFT .EQ. 5 .OR. NLFT .EQ. 6) JSTOP  = JM2
      IF (NBOT .EQ. 5 .OR. NBOT .EQ. 6) KSTART = 3
      IF (NTOP .EQ. 5 .OR. NTOP .EQ. 6) KSTOP  = KM2

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

      DO 10 I = ISTART,ISTOP
      DO 10 J = JSTART,JSTOP
      DO 10 K = KSTART,KSTOP

            IF (B(K,J,I).LT.1.0) THEN
C                 	write (6,*) 'k,j,i,b:',k,j,i,b(k,j,i)
                  IC1 = MIN (IC1,I-1)
                  JC1 = MIN (JC1,J-1)
                  KC1 = MIN (KC1,K-1)
                  IC2 = MAX (IC2,I+1)
                  JC2 = MAX (JC2,J+1)
                  KC2 = MAX (KC2,K+1)

            ENDIF

   10   CONTINUE

C      stop
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

      IC1 = MAX(IC1,2)
      JC1 = MAX(JC1,2)
      KC1 = MAX(KC1,2)
      IC2 = MIN(IC2,IM1)
      JC2 = MIN(JC2,JM1)
      KC2 = MIN(KC2,KM1)

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

         WRITE(6,*)
         WRITE(6,*)
     $   'MELDUNG AUS SBBC12: GRID,IC1,IC2: ',IGRID,IC1,IC2
         WRITE(6,*)
     $   'MELDUNG AUS SBBC12: GRID,JC1,JC2: ',IGRID,JC1,JC2
         WRITE(6,*)
     $   'MELDUNG AUS SBBC12: GRID,KC1,KC2: ',IGRID,KC1,KC2
         WRITE(6,*)

         IF (IC1.GE.IC2) CALL ERRR ( 435 , ' SBBC12 ' )
         IF (JC1.GE.JC2) CALL ERRR ( 436 , ' SBBC12 ' )
         IF (KC1.GE.KC2) CALL ERRR ( 437 , ' SBBC12 ' )

      RETURN
      END

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

