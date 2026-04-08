










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
      SUBROUTINE PHIMLT3 (KK,JJ,II,KMX,JMX,IMX,
     $                    PHI123,IVAR,PHI1,PHI2,PHI3,
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
C*STARLET***************************************************************
C        P H I M L T 3    IN  PHIMLT3 WERDEN DREI VARIABLEN MIT EINANDER
C                         MULTIPLIZIERT UM SIE STATISTISCH AUSWERTEN ZU 
C                         KOENNEN ("NEUE" STATISTIK MIT SPEICHERUNG DER
C                         MOMENTANWERTE)
C                         PHI123=PHI1 * PHI2 * PHI3
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        TAUXY(KK,JJ,II)+ PRODUKT DER EINGANGSVARIABLEN
C        PHI1 (KK,JJ,II)- 1. EINGANGSVARIABLE
C        PHI2 (KK,JJ,II)- 2. EINGANGSVARIABLE
C        PHI3 (KK,JJ,II)- 3. EINGANGSVARIABLE
C        IVAR           - CHARACTER-VARIABLE DER FORM 'XXX','UVW' etc.
C                         zur Definition der Interpolationen
C UPROG                 : ERRR
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        11.05.00 (SE)  : Original abgeleitet von TAUFG
C                         UND EINFUEHRUNG DER "NEUEN" STATISTIK
C*STARLET***************************************************************
C
      CHARACTER *3  IVAR
C
C
      REAL     PHI1(KK,JJ,II),PHI2(KK,JJ,II),PHI3(KK,JJ,II),
     &         PHI123(KK,JJ,II)
      INTEGER  I,J,K,I1,I2,I3,J1,J2,J3,K1,K2,K3,IM2,JM2,KM2,IMX,JMX,KMX
      INTEGER  IFRFIX,JRIFIX
      INTEGER  ISTART,JSTART,KSTART,NFRO,NRGT
C
      TOFF   = 0
      IFRFIX = 0
      JRIFIX = 0

      K1     = 0
      J1     = 0
      I1     = 0
      K2     = 0
      J2     = 0
      I2     = 0
      K3     = 0
      J3     = 0
      I3     = 0

C     DEFINITION DER INTERPOLATIONEN FUER UEBRIGE TERME:
C
C     FUER UUU,VVV,WWW KEINE INTERPOLATION NOETIG
C
      IF ( IVAR .EQ. 'XXX' ) THEN

         CONTINUE

      ELSEIF (IVAR .EQ. 'UUV') THEN

         J1     = 1
         J2     = 1
         I3     = 1

      ELSEIF (IVAR .EQ. 'UUW') THEN

         K1     = 1
         K2     = 1
         I3     = 1

      ELSEIF (IVAR .EQ. 'VVU') THEN

         I1     = 1
         I2     = 1
         J3     = 1

      ELSEIF (IVAR .EQ. 'VVW') THEN

         K1     = 1
         K2     = 1
         J3     = 1

      ELSEIF (IVAR .EQ. 'WWU') THEN

         I1     = 1
         I2     = 1
         K3     = 1

      ELSEIF (IVAR .EQ. 'WWV') THEN

         I1     = 1
         I2     = 1
         K3     = 1

      ELSEIF (IVAR .EQ. 'UVW') THEN

         I1     = -1
         J2     = -1
         K3     = -1

      ELSEIF (IVAR .EQ. 'UTT') THEN

          I1    =  1

      ELSEIF (IVAR .EQ. 'VTT') THEN

          J1    =  1

      ELSEIF (IVAR .EQ. 'WTT') THEN

          K1    =  1

      ELSEIF (IVAR .EQ. 'TTT') THEN

        TOFF   = 1

      ELSE

      CALL ERRR (501,' PHIMLT3 ')

      ENDIF
C
C
      IF (NFRO .EQ. 2) THEN
      IF(I1+I2 .EQ. 0) IFRFIX = 1
      ENDIF
      IF (NRGT .EQ. 2) THEN
      IF(J1+J2 .EQ. 0) JRIFIX = 1
      ENDIF
C
C                                 **************************************
C                                 BEHANDLUNG DER RANDFLAECHEN
C                                 **************************************
C

      ISTOP    = IMX - 2 + TOFF
      JSTOP    = JMX - 2 + TOFF
      KSTOP    = KMX - 2 + TOFF

      ISTART = 3 - I1 - I2 - IFRFIX - TOFF
      DO 100 I = ISTART,ISTOP
         JSTART = 3 - J1 - J2 - JRIFIX - TOFF 
         DO 110 J = JSTART,JSTOP
C
            KSTART = 3-K1-K2-TOFF 
C
            DO 120 K = KSTART,KSTOP
            PHI123(K,J,I) = .125
     $                    * (PHI1(K,J,I) + PHI1(K+K1,J+J1,I+I1))
     $                    * (PHI2(K,J,I) + PHI2(K+K2,J+J2,I+I2))
     $                    * (PHI3(K,J,I) + PHI3(K+K3,J+J3,I+I3))

  120       CONTINUE
  110    CONTINUE
  100 CONTINUE

      RETURN
      END
