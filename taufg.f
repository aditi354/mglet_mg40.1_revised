










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
      SUBROUTINE TAUFG   (KK,JJ,II,KMX,JMX,IMX,GMOL,RHO,
     $                    TAUIJ,UIFG,UJFG,IVAR,
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
C*STARLET***************************************************************
C        T A U F G        IN TAUFG WIRD DER ANTEIL DER GROBSTRUKTUR
C                         AN DEN SCHUBSPANNUNGEN
C                           TAUIJ = - RHO * UI''UJ''
C                         BERECHNET.
C                         FOR SCALAR TRANSPORT:
C                           TAUIJ = - UJ''T''
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        GMOL           - MOLEKULARE DYNAMISCHE VISKOSITAET
C        RHO            - DICHTE DES FLUIDS (= CONST)
C        TAUIJ(KK,JJ,II)+ SCHUBSPANNUNG
C        UIFG (KK,JJ,II)- FLUKTUATIONEN DER GESCHWINDIGKEITSKOMP. UI
C        UJFG (KK,JJ,II)- FLUKTUATIONEN DER GESCHWINDIGKEITSKOMP. UJ
C        IVAR           - CHARACTER-VARIABLE DER FORM 'UW', 'VW'
C                         ODER 'UV'. BEISPIELSWEISE BEDEUTET 'UW', DASS
C                         TAUXZ = - RHO * U''W''
C                         AUSGEWERTET WIRD.
C
C UPROG                 : ERRR
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        25.08.86 (HW)  : ORIGINAL
C         1.11.93 (MM)  : RANDBEDINGUNGEN UEBER KOPF
C        12.02.03 (TB)  : MODIFIED FOR USE OF SCALAR HEAT FLUX 
C                         CALCULATION
C
C*STARLET***************************************************************
C
      CHARACTER (LEN=2)  IVAR
C
C
      REAL     UIFG(KK,JJ,II),  UJFG(KK,JJ,II),  TAUIJ(KK,JJ,II)
C
      IM2    = IMX - 2
      JM2    = JMX - 2
      KM2    = KMX - 2
      IFRFIX = 0
      JRIFIX = 0
      RHOFAK = -0.25 * RHO
C
C                                 BESTIMMUNG DER VERSCHIEBUNGSKOEFFI-
C                                 ZIENTEN (DIE SCHUBSPANNUNGEN SIND
C                                 ZWEIFACH 'GESTAGGERT' !)
C
      IF(IVAR .EQ. 'UW') THEN
C                                 BERUECKSICHTIGUNG VON 'W'
         K1     = 1
         J1     = 0
         I1     = 0
C                                 BERUECKSICHTIGUNG VON 'U'
         K2     = 0
         J2     = 0
         I2     = 1
         GOTO 2100
      ENDIF
C
      IF(IVAR .EQ. 'VW') THEN
C                                 BERUECKSICHTIGUNG VON 'W'
         K1     = 1
         J1     = 0
         I1     = 0
C                                 BERUECKSICHTIGUNG VON 'V'
         K2     = 0
         J2     = 1
         I2     = 0
         GOTO 2100
      ENDIF
C
      IF(IVAR .EQ. 'UV') THEN
C                                 BERUECKSICHTIGUNG VON 'V'
         K1     = 0
         J1     = 1
         I1     = 0
C                                 BERUECKSICHTIGUNG VON 'U'
         K2     = 0
         J2     = 0
         I2     = 1
         GOTO 2100
      ENDIF
C     
C      
      CALL ERRR (501,' TAUFG    ')
C
 2100 CONTINUE
C
      IF (NFRO .EQ. 2) THEN
      IF(I1+I2 .EQ. 0) IFRFIX = 1
      ENDIF
      IF (NRGT .EQ. 2) THEN
      IF(J1+J2 .EQ. 0) JRIFIX = 1
      ENDIF
C
C                                 **************************************
C                                 BERECHNUNG DES GROBSTRUKTURANTEILS
C                                 AN DEN SCHUBSPANNUNGEN IM GESAMTEN
C                                 FELD (EINSCHL. DER RANDFLAECHEN)
C                                 **************************************
C
      ISTART = 3 - I1 - I2 - IFRFIX
      DO 100 I = ISTART,IM2
         JSTART = 3 - J1 - J2 - JRIFIX
         DO 110 J = JSTART,JM2
C
            KSTART = 3-K1-K2
C
            DO 120 K = KSTART,KM2
               TAUIJ(K,J,I) = RHOFAK
     $                      * (UIFG(K,J,I) + UIFG(K+K1,J+J1,I+I1))
     $                      * (UJFG(K,J,I) + UJFG(K+K2,J+J2,I+I2))
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
