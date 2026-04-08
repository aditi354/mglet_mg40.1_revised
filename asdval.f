










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
      SUBROUTINE ASDVAL  (FVT,IDIMF,NLARGE,AK,BK,ALARGK,PHIK,ASDK,
     $                    SMAONE,K)
C*STARLET***************************************************************
C        A S D V A L      IN ASDVAL WERDEN DIE FOURIER-KOEFFIZIENTEN
C                         A_K UND B_K GEBILDET (SIEHE BRONSTEIN -
C                         SEMENDJAJEW 17. AUFLAGE (1977) S. 480).
C                         FERNER WIRD DER KOEFFIZIENT GROSS-A_K
C                         (ALARGK) UND DER PHASENWINKEL PHI_K MIT
C                         0 <= PHI_K <= 2*PI BESTIMMT (SIEHE BRONSTEIN
C                         S. 473). AUS DEN KOEFFIZIENTEN A_K UND B_K
C                         WIRD EIN EINZELNER ORDINATENWERT DER
C                         A_UTO-S_PECTRAL-D_ENSITY-FUNCTION GEBILDET.
C                         DIESER K-TE ORDINATENWERT (K = 0,1,2, ... ,
C                         NLARGE/2) GEHOERT ZUR FREQUENZ   F_K =
C                         K / T_PERIODE = K / (NLARGE*DELTA_T)
C
C                         ZUR INTEGRATON WIRD DIE TRAPEZREGEL VER-
C                         WENDET.
C*STARLET***************************************************************
C
C PARAM: FVT (0:IDIMF)  + ORDINATENWERTE DER ZEITREIHE, DIE ZU DEN
C                         ZEITPUNKTEN  T_N = N * DELTA_T  (N = 0,1,2,
C                         ... NLARGE) DEFINIERT SIND. ES MUSS GELTEN:
C                         FVT(0) = FVT(NLARGE). ZUR VORBEUGUNG BEI
C                         FUNKTIONEN MIT UNSTETIGKEITSSTELLEN WIRD STETS
C                         FVT (0)_NEU = FVT (NLARGE)_NEU =
C                         0.5 * (FVT (0)_ALT + FVT (NLARGE)_ALT)
C                         GEBILDET
C        IDIMF          - ARRAYDIMENSION DES FVT-FELDES. ES MUSS
C                         GELTEN: NLARGE <= IDIMF
C        NLARGE         - ANZAHL DER TEILINTERVALLE DER ZEITREIHE.
C                         NLARGE MUSS  G E R A D Z A H L I G  SEIN !
C        AK             + K-TER FOURIERKOEFFIZIENT 'A'
C        BK             + K-TER FOURIERKOEFFIZIENT 'B'
C        ALARGK         + SQRT (AK**2 + BK**2)
C        PHIK           + TAN (PHI_K) = AK / BK
C        ASDK           + K-TER ORDINATENWERT DER AUTO-SPECTRAL-DENSITY-
C                         FUNCTION (NORMIERT MIT DELTA_T !)
C        SMAONE         - SMAONE IST SO GROSS, DASS AUF JEDER MASCHINE
C                         1.0 + SMAONE > 1.0 GERADE NICHT MEHR
C                         ERKANNT WIRD. (MASS FUER DEN RUNDUNGSFEHLER)
C        K              - INDEX
C
C UPROG                 : ERRR
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        08.11.88 (HW)  : ORIGINAL AUS SFTKOE ABGELEITET
C        28.11.88 (HW)  : AUSWERTUNG VON GROSS-A_K UND PHI_K
C
C*STARLET***************************************************************
C
      REAL       FVT (0:IDIMF)
C
C                                 EINIGE KONSTANTEN
C
      PI      = ACOS(-1.0)
      NSMALL  = NLARGE / 2
      AMUL    = 1.0 / FLOAT (NSMALL)
      OMEGA   = 2.0 * PI / FLOAT (NLARGE)
C
C                                 TEST AUF GERADZAHLIGKEIT VON
C                                 NLARGE
C
      AK     = 999.999
      BK     = 999.999
      ALARGK = 999.999
      PHIK   = 999.999
      ASDK   = 999.999
C
      IF(IDIMF     .LT.      NLARGE) THEN
         WRITE (6,6010)
         WRITE (6,6015) NLARGE
         CALL ERRR (501,' ASDVAL   ')
      ENDIF
      IF(NLARGE              .LT. 2) THEN
         WRITE (6,6010)
         WRITE (6,6020) NLARGE
         GOTO 9999
      ENDIF
      IF((NLARGE - 2*NSMALL) .NE. 0) THEN
         WRITE (6,6010)
         WRITE (6,6030) NLARGE
         GOTO 9999
      ENDIF
C
      FVT (0)      = 0.5 * (FVT(0) + FVT(NLARGE))
      FVT (NLARGE) = FVT(0)
C
C                                 BERECHNUNG DES FOURIER-KOEFFIZIENTEN
C                                 A_K
C                                 ------------------------------------
C
      AK    = 0.0
      DO 100 N = 0,NLARGE-1
  100    AK = AK + AMUL * FVT(N) * COS (FLOAT(K*N)*OMEGA)
C
C                                 BERECHNUNG DES FOURIER-KOEFFIZIENTEN
C                                 B_K
C                                 ------------------------------------
C
      IF(K .EQ. 0  .OR.  K .EQ. NSMALL) THEN
         BK    = 0.0
      ELSE
         BK    = 0.0
         DO 200 N = 0,NLARGE-1
  200       BK = BK + AMUL * FVT(N) * SIN (FLOAT(K*N)*OMEGA)
      END IF
C
C                                 BERECHNUNG DES FOURIER-KOEFFIZIENTEN
C                                 GROSS-A_K
C                                 ------------------------------------
C
      ALARGK = SQRT (AK**2 + BK**2)
C
C                                 BERECHNUNG DES PHASENWINKELS PHI_K
C                                 ------------------------------------
C
      IF(ABS (BK) .LT. SMAONE) THEN
         PHIK   = PI * (1.0 - SIGN (0.5, AK))
      ELSE
         PHIK   = ATAN ( ABS(AK) / (AMAX1 (ABS (BK), ABS (SMAONE))))
     $          * SIGN (1.0, AK*BK)
     $          + PI * (1.0 - SIGN(1.0,AK) * (0.5 + SIGN(0.5,BK)))
      END IF
C
      ASDK  = FLOAT (NSMALL) * (AK**2 + BK**2)
C
 9999 RETURN
 6010 FORMAT (/,' ********** FEHLERMELDUNG AUS SUBR. ASDVAL ',
     $        ' **********')
 6015 FORMAT (5X,'IM PARAMETER-STATEMENT MUSS  IDIMF  MINDESTENS',
     $        ' AUF DEN WERT ',I6,' GESETZT WERDEN !')
 6020 FORMAT (5X,'DIE ANZAHL DER TEILINTERVALLE IST KLEINER ALS',
     $        ' 2.  NLARGE = ',I6)
 6030 FORMAT (5X,'DIE ANZAHL DER TEILINTERVALLE IST  N I C H T',
     $        ' GERADZAHLIG (NLARGE = ',I6,') !')
      END
