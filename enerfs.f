










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
      SUBROUTINE ENERFS  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,
     $                    DX,DY,DZ,X,Y,Z,B,G,GMOL,RHO,EFS,ESUM,
     $                 IC1,IC2,JC1,JC2,KC1,KC2,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,ISUM,
     $                  CONV1S)
C*STARLET***************************************************************
C        E N E R F S      IN ENERFS WIRD DIE FEINSTRUKTURENERGIE
C                         AUS DEM MOMENTAN VORHANDENEN GESCHWINDIG-
C                         KEITSFELD GEBILDET.
C                         DIE BENOETIGTE CHARAKTERISTISCHE ZELLABMESSUNG
C                         WIRD  N U R  AN DEN  W A N D N A E C H S T E N
C                         ZELLEN DEM PRANDTL''SCHEN MISCHUNGSWEG ANGE-
C                         PASST.
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        DDX,DDY,DDZ    - ABMESSUNGEN DER BASISZELLEN
C        DX,DY,DZ       - ABSTAND DER GITTERPUNKTE
C        X,Y,Z          - KOORDINATEN DER ZELLMITTELPUNKTE
C                         IN WANDNAEHE)
C        B(KK,JJ,II)    - LAENGENMASS DES TURBULENZMODELLS
C        G(KK,JJ,II)    - EFFEKTIVE DYNAMISCHE VISKOSITAET
C                         (MUE-EFF = MUE-TURB + MUE-MOL)
C        GMOL           - MOLEKULARE DYNAMISCHE VISKOSITAET
C        RHO            - DICHTE DES FLUIDS (= CONST.)
C        EFS(KK,JJ,II)  + ENTHAELT DIE FEINSTRUKTURENERGIE
C        ESUM           + ENTHAELT GESAMTENERGIE
C
C        ISUM           - SCHALTER, OB GESAMTENERGIE BERECHNET WERDEN
C                         SOLL:    0 --> NICHT BERECHNET
C                                  1 --> SUMME WIRD GEBILDET
C
C UPROG                 : KEINE
C
C DEFINE-DIREKTIVEN     : DDMIX2, DDMIX3, DDMIX4
C
C        27.08.86 (HW)  : ORIGINAL AUS SUMESG ABGELEITET
C        17.09.86 (HW)  : DDMIX4 KANN DEFINIERT WERDEN. DER MISCHUNGS-
C                         WEG WIRD AUS DEM VERHAELTNIS VON ZELLVOLUMEN
C                         ZU ZELLOBERFLAECHE GEBILDET.
C         1.11.93 (MM)  : RANDBEDINGUNGEN UEBER KOPF
C         7.12.93 (MM)  : ISUM EINGEFUEHRT
C        12.09.96 (A.O.): IN B-FELD STEHT JETZT DIE MISCHUNGSWEGLAENGE SELBER,
C                         SODASS AB JETZT NICHT MEHR MIT DER FEINSTRUKTURKONSTANTE
C                         DIVIDIERT WERDEN MUSS.
C
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

      COMMON /CONLES/  CONV2S,CMUE,CAPPA,ECONST
      SAVE   /CONLES/
C
C
      REAL     G(KK,JJ,II),         EFS(KK,JJ,II),
     $         B(KK,JJ,II)
C
      REAL     DDX(II),   DDY(JJ),  DDZ(KK),
     $          DX(II),    DY(JJ),   DZ(KK),
     $           X(II),     Y(JJ),    Z(KK),
     $     CONV1S(II)
C
C
C                                 KONSTANTEN
C
      IM2    = IMX - 2
      JM2    = JMX - 2
      KM2    = KMX - 2
      IFRFIX = 0
      JRIFIX = 0

      IF ( NFRO .EQ. 2) IFRFIX = 1
      IF ( NRGT .EQ. 2) JRIFIX = 1
C
      SMAL2  = 2.0*SMALL
      THIRD  = 1.0/3.0
C
      ISTART = 3 - IFRFIX
      IEND   = IMX - 2
C
      JSTART = 3 - JRIFIX
      JEND   = JMX - 2
C

      DO 10 I = ISTART,IEND
         DO 20 J = JSTART,JEND
C
            KSTART = 3
            KSTOP  = KM2
C
 2000       DO 30 K = KSTART,KSTOP
C
C                        DELTA  = AMIN1 (DELTA,DELTAG)
                        EFS(K,J,I) = ((G(K,J,I)-GMOL)
     $                        / (RHO*CONV2S*B(K,J,I)))**2
   30       CONTINUE
C
   20    CONTINUE
   10 CONTINUE
C
C
      IF (ISUM.EQ.1) THEN

        DO I = ISTART,IEND
        DO J = JSTART,JEND
        DO K = KSTART,KSTOP

           ESUM = ESUM + EFS(K,J,I)*DDX(I)*DDY(J)*DDZ(K)
  
        ENDDO
        ENDDO
        ENDDO

           XL = X(  IM2   ) + DX(  IM2     )*0.5
     $        -(X(ISTART  ) - DX(ISTART  -1)*0.5)

           YL = Y(  JM2   ) + DY(  JM2     )*0.5
     $        -(Y(JSTART  ) - DY(JSTART  -1)*0.5)

           ZL = Z(  KM2   ) + DZ(  KM2     )*0.5
     $        -(Z(KSTART  ) - DZ(KSTART  -1)*0.5)

           ESUM = ESUM / (XL*YL*ZL)

      ENDIF
      RETURN
      END
