










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
      SUBROUTINE TAUFS   (KK,JJ,II,KMX,JMX,IMX,IB1,IB2,JB1,JB2,KB,
     $                    G,GMOL,RHO,UGRID,TAUIJ,
     $                    UI,XI,DXI,DDXI,NNI,UJ,XJ,DXJ,DDXJ,NNJ,IVAR
     $                    )
C*STARLET***************************************************************
C        T A U F S        IN TAUFS WIRD DER ANTEIL DES FEINSTRUKTUR-
C                         MODELLS AN DEN SCHUBSPANNUNGEN
C                           TAUIJ = MUE(TUR) * (DUI/DXJ + DUJ/DXI)
C                         BERECHNET.
C                         FOR SCALAR TRANSPORT:
C                           TAUIJ = MUE(TUR)/PRTURB * DT/DXI
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        IB1, IB2       - GRENZE DES KUBUSSES IN X-RI.
C        JB1, JB2       -   ""     "     "     IN Y-RI.
C        KB             -   ""     "     "     IN Z-RI. (TOP-FLAECHE)
C        G(KK,JJ,II)    - EFFEKTIVE DYNAMISCHE VISKOSITAET (= MUE)
C        GMOL           - MOLEKULARE DYNAMISCHE VISKOSITAET
C        RHO            - DICHTE DES FLUIDS (= CONST)
C        UGRID          - GESCHWINDIGKEIT DES KOORDINATENSYSTEMS IN
C                         X-RICHTUNG (GALILEI-TRANSFORMATION)
C        TAUIJ(KK,JJ,II)+ SCHUBSPANNUNG
C        UI   (KK,JJ,II)- GESCHWINDIGKEITSKOMPONENTE UI
C        XI   (NNI)     - KOORDINATEN DER ZELLMITTELPUNKTE IN 'I'-RI.
C        DXI  (NNI)     - ABSTAND DER BASISZELMITTELPUNKTE IN 'I'-RI.
C        DDXI (NNI)     - ABMESSUNGEN DER BASISZELLEN      IN 'I'-RI.
C        NNI            - ARRAYDIMENSION                   IN 'I'-RI.
C        UJ   (KK,JJ,II)- GESCHWINDIGKEITSKOMPONENTE UJ
C        XJ   (NNJ)     - KOORDINATEN DER ZELLMITTELPUNKTE IN 'J'-RI.
C        DXJ  (NNJ)     - ABSTAND DER BASISZELMITTELPUNKTE IN 'J'-RI.
C        DDXJ (NNJ)     - ABMESSUNGEN DER BASISZELLEN      IN 'J'-RI.
C        NNJ            - ARRAYDIMENSION                   IN 'J'-RI.
C        IVAR           - CHARACTER-VARIABLE DER FORM 'UW', 'VW'
C                         ODER 'UV'. BEISPIELSWEISE BEDEUTET 'UW', DASS
C                         TAUXZ = MUE(TUR) * (DU/DZ + DW/DX)
C                         AUSGEWERTET WIRD.
C
C UPROG                 : PRT (FOR _TSCAL_)
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        25.08.86 (HW)  : ORIGINAL
C         1.11.93 (MM)  : RANDBEDINGUNGEN UEBER KOPF
C                         NACHKORREKTUR AN DEN KUBUSKANTEN RAUSGEWORFEN
C        12.02.03 (TB)  : MODIFIED FOR USE OF SCALAR HEAT FLUX 
C                         CALCULATION
C
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
      CHARACTER (LEN=2)  IVAR
C
C
      REAL     UI(KK,JJ,II),  UJ(KK,JJ,II),  G(KK,JJ,II),
     $         TAUIJ(KK,JJ,II)
      REAL     XI(NNI),  DXI(NNI),  DDXI(NNI),
     $         XJ(NNJ),  DXJ(NNJ),  DDXJ(NNJ)
C
      IM2    = IMX - 2
      JM2    = JMX - 2
      KM2    = KMX - 2
      IFRFIX = 0
      JRIFIX = 0
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
      CALL ERRR (501,' TAUFS    ')
C
 2100 CONTINUE
C
C
C                                 **************************************
C                                 BERECHNUNG DES ANTEILS DES FEIN-
C                                 STRUKTURMODELLS AN DEN SCHUBSPANNUNGEN
C                                 IM GESAMTEN FELD (EINSCHL. DER RAND-
C                                 FLAECHEN)
C                                 DER TURBULENTE AUSTAUSCHKOEFFIZIENT
C                                 WIRD AUS KONSISTENZGRUENDEN GEMITTELT
C                                 WIE IN DER IMPULSGLEICHUNG.
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
C

C                                 HIER: MUET * DUI/DXJ
C
               GUI = AMAX1((G(K,J,I)*G(K+K1,J+J1,I+I1)
     $             /       (G(K,J,I)+G(K+K1,J+J1,I+I1))
     $             +        G(K+K2,J+J2,I+I2)*G(K+K1+K2,J+J1+J2,I+I1+I2)
     $             /       (G(K+K2,J+J2,I+I2)+G(K+K1+K2,J+J1+J2,I+I1+I2)
     $                     )) - GMOL,0.0)
               TAUIJ(K,J,I) = GUI* (UI(K+K1,J+J1,I+I1) - UI(K,J,I))
     $                      /       DXJ(K*K1+J*J1+I*I1)
C
C                                 HIER: MUET * DUJ/DXI
C
               GUJ = AMAX1((G(K,J,I)*G(K+K2,J+J2,I+I2)
     $             /       (G(K,J,I)+G(K+K2,J+J2,I+I2))
     $             +        G(K+K1,J+J1,I+I1)*G(K+K1+K2,J+J1+J2,I+I1+I2)
     $             /       (G(K+K1,J+J1,I+I1)+G(K+K1+K2,J+J1+J2,I+I1+I2)
     $                     )) - GMOL,0.0)
               TAUIJ(K,J,I) = TAUIJ(K,J,I)
     $                      + GUJ* (UJ(K+K2,J+J2,I+I2) - UJ(K,J,I))
     $                      /       DXI(K*K2+J*J2+I*I2)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C
C                                 NACHKORREKTUR DER SCHUBSPANNUNGEN AN
C                                 DEN KUBUSKANTEN AUSGESCHALTET (M.M.)
C
C
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
