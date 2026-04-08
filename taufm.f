










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
      SUBROUTINE TAUFM   (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID,TAUIJ,
     $                    UI,XI,DXI,DDXI,NNI,UJ,XJ,DXJ,DDXJ,NNJ,IVAR,
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB
     $                    )     
C*STARLET***************************************************************
C        T A U F M        IN TAUFM WIRD DER MOLEKULARE ANTEIL DER
C                         SCHUBSPANNUNGEN
C                           TAUIJ = MUE(MOL) * (DUI/DXJ + DUJ/DXI)
C                         BERECHNET.
C                         DA DIE WANDSCHUBSPANNUNG MITTELS DES LOG.
C                         WANDGESETZES BESTIMMT WIRD, IST EINE SONDER-
C                         BEHANDLUNG ALLER NOSLIP-SEITENFLAECHEN ER-
C                         FORDERLICH.
C                                                                       
C                                                                       
C                 --------------------------------                      
C                |               |                |                     
C                |    K+K1       |  K+K1+K2       |                     
C                |    J+J1       |  J+J1+J2       |                     
C                |    I+I1       |  I+I1+I2       |   |XJ               
C                |               |                |   |                 
C                |               |                |   |                 
C                 --------------------------------    |                 
C                |               |                |   |                 
C                |      K        |    K+K2        |   |                 
C                |      J        |    J+J2        |    ----------->     
C                |      I        |    I+I2        |            XI       
C                |               |                |                     
C                |               |                |                     
C                 --------------------------------                      
C                           
C                      FOR SCALAR WALL HEAT FLUX
C                      TAUIJ = MUE(MOL)/PRMOL * (DT/DXI)
C                                                                       
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
C                         TAUXZ = MUE(MOL) * (DU/DZ + DW/DX)
C                         AUSGEWERTET WIRD.
C
C UPROG                 : ERRR, TAUWP, QWL
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        25.08.86 (HW)  : ORIGINAL
C        30.12.88 (HW)  : GEPS  WIRD MASCHINEN"UN"ABHAENGIG BELEGT,
C                         UEBERPRUEFUNGEN WURDEN VERBESSERT
C         6.11.92 (MM)  : EINFUEHRUNG VON LENOS UND RINOS
C                         VORSICHT!!!!   NICHT ABGECHECKT!!!!
C         1.11.93 (MM)  : RANDBEDINGUNGEN UEBER KOPF
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
      REAL     UI(KK,JJ,II),  UJ(KK,JJ,II),  G(KK,JJ,II),  B(KK,JJ,II),
     $         TAUIJ(KK,JJ,II)
     
      REAL     XI(NNI),  DXI(NNI),  DDXI(NNI),
     $         XJ(NNJ),  DXJ(NNJ),  DDXJ(NNJ)

C
      IM2    = IMX - 2
      JM2    = JMX - 2
      KM2    = KMX - 2
      IFRFIX = 0
      JRIFIX = 0
      GEPS   = AMAX1 ((GMOL/10000.0) , (10.0 * SMAONE))
      IF(((    GMOL   - GEPS) .LE. 0.0)  .OR.
     $   ((ABS(SMALL) - GEPS) .GE. 0.0)       ) THEN
         CALL ERRR (501,' TAUFM    ')
      END IF
C
C                                 GESCHWINDIGKEIT DES KOORDINATENSYSTEMS IN
C                                 I- UND J-RICHTUNG
        UIGRID = 0.0
        UJGRID = 0.0
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

         UIGRID = UGRID

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

         UIGRID = UGRID

         GOTO 2100
      ENDIF
C
C
      CALL ERRR (510,' TAUFM    ')
C
 2100 CONTINUE
C
C
C
C                                 ***********************************
C                                 BERECHNUNG DES MOLEKULAREN ANTEILS
C                                 DER SCHUBSPANNUNG AN FREIEN GRENZFL.
C                                 DIE SCHUBSPANNUNGEN AN DEN NOSLIP-
C                                 RAENDERN WERDEN DURCH FAKTOREN
C                                 BERUECKSICHTIGT
C                                 ***********************************
C
      ISTART = 3 - I1 - I2 - IFRFIX
      DO 900 I = ISTART,IM2
         JSTART = 3 - J1 - J2 - JRIFIX
         DO 910 J = JSTART,JM2
            KSTART = 3-K1-K2
            DO 920 K = KSTART,KM2
C
               FUI  = SIGN(0.25,(G(K      ,J      ,I      )-GEPS)
     $              *           (G(K+K1   ,J+J1   ,I+I1   )-GEPS))+0.25
C
               FUIP = SIGN(0.25,(G(K   +K2,J   +J2,I   +I2)-GEPS)
     $              *           (G(K+K1+K2,J+J1+J2,I+I1+I2)-GEPS))+0.25
C
               FUJ  = SIGN(0.25,(G(K      ,J      ,I      )-GEPS)
     $              *           (G(K   +K2,J   +J2,I   +I2)-GEPS))+0.25
C
               FUJP = SIGN(0.25,(G(K+K1   ,J+J1   ,I+I1   )-GEPS)
     $              *           (G(K+K1+K2,J+J1+J2,I+I1+I2)-GEPS))+0.25
C     
               DUI = UI(K+K1,J+J1,I+I1) - UI(K,J,I)
               DUJ = UJ(K+K2,J+J2,I+I2) - UJ(K,J,I)
               DUIGRID = DUI + UIGRID*SIGN(1.0,G(K,J,I)-GEPS)
               DUJGRID = DUJ + UJGRID*SIGN(1.0,G(K,J,I)-GEPS)
C
               TAUIJ(K,J,I) = 
C                                        ANTEILE AN FREIEN GRENZFLAECHEN
C                                                             J-RICHTUNG
     $         GMOL*DUI/DXJ(K*K1+J*J1+I*I1) * (FUI+FUIP)
C                                                             I-RICHTUNG
     $       + GMOL*DUJ/DXI(K*K2+J*J2+I*I2) * (FUJ+FUJP)
C
C                                                ANTEILE AN WANDFLAECHEN
C                                                             J-RICHTUNG
     $       + TAUWP(DUIGRID,DXJ(K*K1+J*J1+I*I1))*(1-FUI-FUIP)
C                                                             I-RICHTUNG
     $       + TAUWP(DUJGRID,DXI(K*K2+J*J2+I*I2))*(1-FUJ-FUJP)
  920       CONTINUE
  910    CONTINUE
  900 CONTINUE
C
      RETURN
      END
