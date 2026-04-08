










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
      SUBROUTINE PHIDIDJ (KK,JJ,II,KMX,JMX,IMX,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    B,G,GMOL,RHO,UGRID, DIDJ,
     $                    UI,XJ,DXJ,DDXJ,NNJ,IVAR,
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)
C*MGLET*****************************************************************
C        P H I D I D J    IN PHIDIDJ WIRD DIE  ABLEITUNG DUI/DXJ
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
C*MGLET*****************************************************************
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
C         DIDJ(KK,JJ,II)+ ABLEITUNG DUI/DXJ
C        UI   (KK,JJ,II)- GESCHWINDIGKEITSKOMPONENTE UI
C        XJ   (NNJ)     - KOORDINATEN DER ZELLMITTELPUNKTE IN 'J'-RI.
C        DXJ  (NNJ)     - ABSTAND DER BASISZELMITTELPUNKTE IN 'J'-RI.
C        DDXJ (NNJ)     - ABMESSUNGEN DER BASISZELLEN      IN 'J'-RI.
C        NNJ            - ARRAYDIMENSION                   IN 'J'-RI.
C        IVAR           - CHARACTER-VARIABLE DER FORM 'XX', 'XY'
C                         'XZ',...,'ZY','ZZ'
C                         BEISPIELSWEISE BEDEUTET 'XZ', DASS
C                         DUIDXJ = DU/DZ
C                         AUSGEWERTET WIRD.
C
C UPROG                 : ERRR
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        13.05.94 (MM)  : ORIGINAL AUS TAUFM ABGELEITET
C
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
      CHARACTER (LEN=2)  IVAR
C
C
      REAL     UI(KK,JJ,II),  G(KK,JJ,II),  B(KK,JJ,II),
     $          DIDJ(KK,JJ,II)
      REAL     
     $         XJ(NNJ),  DXJ(NNJ),  DDXJ(NNJ)
C
      I1 = 0
      J1 = 0
      K1 = 0
      I2 = 0
      J2 = 0
      K2 = 0

      IM2    = IMX - 2
      JM2    = JMX - 2
      KM2    = KMX - 2
      IFRFIX = 0
      JRIFIX = 0
      GEPS   = AMAX1 ((GMOL/10000.0) , (10.0 * SMAONE))
      IF(((    GMOL   - GEPS) .LE. 0.0)  .OR.
     $   ((ABS(SMALL) - GEPS) .GE. 0.0)       ) THEN
         CALL ERRR (501,' PHIDIDJ  ')
      END IF
      RGMOL = 1.0/GMOL
C
C                                 GESCHWINDIGKEIT DES KOORDINATENSYSTEMS IN
C                                 I- UND J-RICHTUNG
        UIGRID = 0.0
      IF(IVAR(1:1) .EQ. 'X')  UIGRID = UGRID
C
C                                 BESTIMMUNG DER VERSCHIEBUNGSKOEFFI-
C                                 ZIENTEN (DIE SCHUBSPANNUNGEN SIND
C                                 ZWEIFACH 'GESTAGGERT' !)
C
      IF(IVAR .EQ. 'XX') THEN
C                                 ABLEITUNG VON U IN X-RICHTUNG
         I2     =  1

      ELSEIF(IVAR .EQ. 'XY') THEN
C                                 ABLEITUNG VON U IN Y-RICHTUNG
         J1     =  1
C
      ELSEIF(IVAR .EQ. 'XZ') THEN
C                                 ABLEITUNG VON U IN Z-RICHTUNG
         K1     =  1
C
      ELSEIF(IVAR .EQ. 'YX') THEN
C                                 ABLEITUNG VON V IN X-RICHTUNG
         I1     =  1

      ELSEIF(IVAR .EQ. 'YY') THEN
C                                 ABLEITUNG VON V IN Y-RICHTUNG
         J2     =  1
C
      ELSEIF(IVAR .EQ. 'YZ') THEN
C                                 ABLEITUNG VON V IN Z-RICHTUNG
         K1     =  1
C
      ELSEIF(IVAR .EQ. 'ZX') THEN
C                                 ABLEITUNG VON W IN X-RICHTUNG
         I1     =  1

      ELSEIF(IVAR .EQ. 'ZY') THEN
C                                 ABLEITUNG VON W IN Y-RICHTUNG
         J1     =  1
C
      ELSEIF(IVAR .EQ. 'ZZ') THEN
C                                 ABLEITUNG VON W IN Z-RICHTUNG
         K2     =  1
C
      ELSE 

         CALL ERRR(510,' PHIDIDJ ')

      ENDIF
C
C
       IF (IVAR(1:1) .EQ. IVAR(2:2)) THEN
C                                 ***********************************
C                                 BERECHNUNG DER ABLEITUNGEN, FALLS
C                                 I GLEICH J IN NORMALER RICHTUNG
C                                 ***********************************


      ISTART = 3 - IFRFIX
      DO I = ISTART,IM2
         JSTART = 3 - JRIFIX
         DO J = JSTART,JM2
C
            KSTART = 3
C
            DO K = KSTART,KM2


               DUI = UI(K,J,I) - UI(K-K2,J-J2,I-I2)

                DIDJ(K,J,I) = DUI/DDXJ(K*K2+J*J2+I*I2)

            ENDDO
         ENDDO
       ENDDO
C
       ELSE
C                                 ***********************************
C                                 BERECHNUNG DER ABLEITUNGEN, FALLS
C                                 I UNGLEICH J
C                                 ***********************************
C
      ISTART = 3 - I1  - IFRFIX
      DO 900 I = ISTART,IM2
         JSTART = 3 - J1 - JRIFIX
         DO 910 J = JSTART,JM2
C
            KSTART = 3-K1
C
            DO 920 K = KSTART,KM2
C
               FUI  = SIGN(0.25,(G(K      ,J      ,I      )-GEPS)
     $              *           (G(K+K1   ,J+J1   ,I+I1   )-GEPS))+0.25
C
               FUIP = SIGN(0.25,(G(K   +K2,J   +J2,I   +I2)-GEPS)
     $              *           (G(K+K1+K2,J+J1+J2,I+I1+I2)-GEPS))+0.25
C
C
               DUI = UI(K+K1,J+J1,I+I1) - UI(K+K2,J+J2,I+I2)
               DUIGRID = DUI + UIGRID*SIGN(1.0,G(K,J,I)-GEPS)
               DUIDXJ=DUI/DXJ(K*K1+J*J1+I*I1) * (FUI+FUIP)
               TAU=RGMOL*TAUWP(DUIGRID,DXJ(K*K1+J*J1+I*I1))*(1-FUI-FUIP)
C
                DIDJ(K,J,I) = DUIDXJ + TAU
C                                        ANTEILE AN FREIEN GRENZFLAECHEN
C                                                             J-RICHTUNG
C    $              DUI/DXJ(K*K1+J*J1+I*I1) * (FUI+FUIP)
C
C                                                ANTEILE AN WANDFLAECHEN
C                                                             J-RICHTUNG
C    $       + TAUWP(DUIGRID,DXJ(K*K1+J*J1+I*I1))*(1-FUI-FUIP)

  920       CONTINUE
  910    CONTINUE
  900 CONTINUE
C
      ENDIF
C
      RETURN
      END
