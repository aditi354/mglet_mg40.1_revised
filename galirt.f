










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
      SUBROUTINE GALIRT  (KK,JJ,II,KMX,JMX,IMX,
     $                    UGRID,PHI,PHIRT,B,IVAR,NFRO,NRGT,NBOT
     $                                          ,NBAC,NLFT,NTOP)
C*STARLET***************************************************************
C        G A L I R T      GALIRT BEWERKSTELLIGT DIE GALILEI-RUECKTRANS-
C                         FORMATION. DAS RUECKTRANSFORMIERTE FELD WIRD
C                         IN DAS FELD 'PHIRT' GESCHRIEBEN; DAS UEBER-
C                         GEBENE FELD 'PHI' BLEIBT UNVERAENDERT.
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        UGRID          - GESCHWINDIGKEIT DES KOORDINATENSYSTEMS
C                         IN X-RICHTUNG (GALILEI-TRANSFORMATION)
C        PHI(KK,JJ,II)  - ALLGEMEINE VARIABLE: 'U', 'V', 'W' ODER
C                         'P' ZULAESSIG
C        PHIRT(KK,JJ,II)+ ENTHAELT DAS RUECKTRANSFORMIERTE FELD
C        IVAR           - CHARACTER-VARIABLE DER FORM '   U', '   V'
C                         ODER '   W'
C
C UPROG                 : ERRR
C
C DEFINE-DIREKTIVEN     : FRFIX, RIFIX
C
C        14.08.86 (HW)  : ORIGINAL
C
C*STARLET***************************************************************
C
      CHARACTER (LEN=4)  IVAR
C
C
      REAL     PHI(KK,JJ,II),       PHIRT(KK,JJ,II),    B(KK,JJ,II)
C
      IF(IVAR .EQ. '   U' .OR. IVAR .EQ. '   V' .OR. IVAR .EQ. '   W'
     $   .OR. IVAR .EQ. '   T' .OR. IVAR .EQ. '   P') GOTO 2000
         CALL ERRR (501,' GALIRT   ')
C
 2000 KSTAG  = 0
      JSTAG  = 0
      ISTAG  = 0
      IFRFIX = 0
      JRIFIX = 0
      KBOFIX = 0
C
      IF(IVAR .EQ. '   P') GOTO 2100
C
      IF(NFRO.EQ.2)   IFRFIX = 1
      IF(NRGT.EQ.2)   JRIFIX = 1
      IF(NBOT.EQ.2)   KBOFIX = 1
C
 2100 UGALIL = 0.0
C
      IF(IVAR .EQ. '   U') THEN
         ISTAG  = 1
         UGALIL = UGRID
      ENDIF
      IF(IVAR .EQ. '   V') THEN
         JSTAG  = 1
      ENDIF
      IF(IVAR .EQ. '   W') THEN
         KSTAG  = 1
      ENDIF
C
C                                 START- UND STOPINDIZES
C
      ISTART = 3-ISTAG-IFRFIX
      ISTOP  = IMX - 2
      JSTART = 3-JSTAG-JRIFIX
      JSTOP  = JMX - 2
      KSTART = 3-KSTAG-KBOFIX
      KSTOP  = KMX - 2
      IF(NFRO .EQ. 5) THEN ISTART = 2
      IF(NBAC .EQ. 5) THEN ISTOP = IMX-1
      IF(NRGT .EQ. 5) THEN JSTART = 2
      IF(NLFT .EQ. 5) THEN JSTOP = JMX-1
      IF(NBOT .EQ. 5) THEN KSTART = 2
      IF(NTOP .EQ. 5) THEN KSTOP = KMX-1 
C
      IF(IVAR .EQ. '   U') THEN

         DO  I = ISTART,ISTOP
         DO  J = JSTART,JSTOP
         DO  K = KSTART,KSTOP
C
               PHIRT(K,J,I) = PHI(K,J,I) 
     $       + UGALIL * 
     $           (SIGN(0.5,B(K,J,I))+0.5) * (SIGN(0.5,B(K,J,I+1))+0.5)
C
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF(IVAR .EQ. '   V') THEN

         DO  I = ISTART,ISTOP
         DO  J = JSTART,JSTOP
         DO  K = KSTART,KSTOP
C
               PHIRT(K,J,I) = PHI(K,J,I)
     $       + UGALIL * 
     $           (SIGN(0.5,B(K,J,I))+0.5) * (SIGN(0.5,B(K,J+1,I))+0.5)
C
         ENDDO
         ENDDO
         ENDDO
      ENDIF
      IF(IVAR .EQ. '   W') THEN

         DO  I = ISTART,ISTOP
         DO  J = JSTART,JSTOP
         DO  K = KSTART,KSTOP
C
               PHIRT(K,J,I) = PHI(K,J,I)
     $       + UGALIL * 
     $           (SIGN(0.5,B(K,J,I))+0.5) * (SIGN(0.5,B(K+1,J,I))+0.5)
C
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF(IVAR .EQ. '   P') THEN

         DO  I = ISTART,ISTOP
         DO  J = JSTART,JSTOP
         DO  K = KSTART,KSTOP
C
               PHIRT(K,J,I) = PHI(K,J,I)
     $       + UGALIL * (SIGN(0.5,B(K,J,I))+0.5) 
C
         ENDDO
         ENDDO
         ENDDO
      ENDIF
C
      RETURN
      END
