










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
      SUBROUTINE DIC  (CIDENT,IIDENT,RIDENT,KK,JJ,II,KMX,JMX,IMX,NBND,
     $                DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,NVAR,
     $                U,V,W,P,G,B,TK,TE,HILF,IGRID,IDIM3D)
C*STAR******************************************************************
C  D I C          EINGABE VON DATEN (KODIERT,KANAL3)
C*STAR********************************************** ST.FRANK 05.06.85 *
C                                                MIT CYB-STEURUNG FB
C                                                    H.WERNER 19.09.85 *
C                                                    M.MANHART 6. 4.92
C     26. 5.95 (MM) : MPI EINGEFUEHRT, KOPF UM HILF ERWEITERT
C     06.02.03 (TB) : SCALAR TRANSPORT IMPLEMENTED
C
C  PARAMETER
C             KK, JJ, II     - ARRAYGRENZEN
C             KMX, JMX, IMX  - ANZAHL DER GITTERPKTE. IN Z-,Y-,X-DIR.
C             X(II)          - KOORDINATEN IN X-RICHTUNG
C             Y(JJ)          - KOORDINATEN IN Y-RICHTUNG
C             Z(KK)          - KOORDINATEN IN Z-RICHTUNG
C             DDX(II)        - X-KANTENLAENGE DER KONTROLLVOLUMINA
C             DDY(JJ)        - Y-KANTENLAENGE DER KONTROLLVOLUMINA
C             DDZ(KK)        - Z-KANTENLAENGE DER KONTROLLVOLUMINA
C             DX(II)         - X-ABSTAND DER GITTERPUNKTE
C             DY(JJ)         - Y-ABSTAND DER GITTERPUNKTE
C             DZ(KK)         - Z-ABSTAND DER GITTERPUNKTE
C             NVAR           - ANZAHL DER AUSZUGEBENDEN DATENFELDER
C             U(KK,JJ,II)    - GESCHWINDIGKEITSFELD
C             V(KK,JJ,II)    - GESCHWINDIGKEITSFELD
C             W(KK,JJ,II)    - GESCHWINDIGKEITSFELD
C             P(KK,JJ,II)    - DRUCKFELD
C             G(KK,JJ,II)    - GAMMAFELD
C             TK(KK,JJ,II)   - TURBULENZENERGIEFELD
C             TE(KK,JJ,II)   - DISSIPATIONSFELD
C
C  ANMERKUNG  CIDENT(10) ENTHAELT KENNUNG
C
C             CIDENT( 1) = "F.BAETKE"
C             CIDENT( 2) = PROGRAMMNAME
C             CIDENT( 3) = LAUFENDE NUMMER
C             CIDENT( 4) = DATE(DUMMY)
C             CIDENT( 5) = CLOCK(DUMMY)
C
C             IIDENT(100) ENTHAELT INTEGER-KONSTANTEN
C
C             IIDENT( 1) = MTSTEP
C             IIDENT( 2) = MPCORR
C             IIDENT( 5) : VERSIONSNUMMER
C             IIDENT(10) = KB
C             IIDENT(11) = JB1
C             IIDENT(12) = JB2
C             IIDENT(13) = IB1
C             IIDENT(14) = IB2
C             IIDENT(50) = SCALAR FIELD IDENTIFIER
C                           0:= SCALAR FIELD UNSET, FIELD NOT READ
C                           1:= SCALAR FIELD SET, FIELD CAN BE READ & 
C                                                WILL BE WRITTEN)
C
C             RIDENT(100) ENTHAELT REAL-KONSTANTEN
C
C             RIDENT( 1) = RE
C             RIDENT( 2) = DT
C             RIDENT( 3) = EPCORR
C             RIDENT(10) = ZWIDTH
C             RIDENT(11) = YWIDTH
C             RIDENT(12) = XWIDTH
C             RIDENT(13) = BETA
C
C*STAR******************************************************************
      CHARACTER (LEN=8)   CIDENT(10)
      CHARACTER (LEN=80)  CT,  CH1680, TEXT
      INTEGER KK, JJ, II, KMX, JMX, IMX, NVAR, IIDENT(100)
      REAL    RIDENT(100), X(II), Y(JJ), Z(KK),
     $        DX(II), DY(JJ), DZ(KK), DDX(II), DDY(JJ), DDZ(KK),
     $        U( IDIM3D ), V( IDIM3D ), W( IDIM3D ),
     $        P( IDIM3D ), G( IDIM3D ), TK( IDIM3D ), TE( IDIM3D ),
     $        B( IDIM3D )
C
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC FALLS ERGEBNIS AUS MLET
C                                          WIRD FOLGENDES UEBERSPRUNGEN
      IF (IIDENT(5) .GT. 8) THEN

        READ(3,6010) TEXT
        READ(3,6030) IGRIDOLD
        READ(3,6030) KKOLD,JJOLD,IIOLD,KMXOLD,JMXOLD,IMXOLD
C
C                       KKOLD,JJOLD,IIOLD,SIND DUMMYVARIABLE
C                       UND WERDEN NICHT BENOETIGT.
C
 6010   FORMAT (A80)
 6030   FORMAT(6(I9,1X))
C
C                      PLAUSIBILITAETSTEST
C
10003   IF(IGRID.NE.IGRIDOLD) CALL ERRR(700,' DIC      ')
        IF(KMX.NE.KMXOLD) CALL ERRR(701,' DIC      ')
        IF(JMX.NE.JMXOLD) CALL ERRR(702,' DIC      ')
        IF(IMX.NE.IMXOLD) CALL ERRR(703,' DIC      ')

      ENDIF
C
C
      READ(3,6050)(X(I),I=1,IMX),(DX(I),I=1,IMX),(DDX(I),I=1,IMX),
     $             (Y(J),J=1,JMX),(DY(J),J=1,JMX),(DDY(J),J=1,JMX),
     $             (Z(K),K=1,KMX),(DZ(K),K=1,KMX),(DDZ(K),K=1,KMX)
C
C
C                                    **********************************
C
C                                    AB HIER VERSION VOR 6.4.92
C
C                                    **********************************
	IF (IIDENT(5).LT.8) THEN
C
      READ(3,6040) NVAR
      IF (NVAR .LT. 4 .OR. NVAR .GT. 7) CALL ERRR(704,' DIC      ')
 6040 FORMAT(I9,1X)
C
            CALL READ3D(KMX,JMX,IMX,U,'CODIERT ',1)
            CALL READ3D(KMX,JMX,IMX,V,'CODIERT ',1)
            CALL READ3D(KMX,JMX,IMX,W,'CODIERT ',1)
            CALL READ3D(KMX,JMX,IMX,P,'CODIERT ',1)

      IF (NVAR .EQ. 4) RETURN
C
            CALL READ3D(KMX,JMX,IMX,G,'CODIERT ',1)
C
      IF (NVAR .EQ. 5) RETURN
            CALL READ3D(KMX,JMX,IMX,TK,'CODIERT ',1)

      IF (NVAR .EQ. 6) RETURN
            CALL READ3D(KMX,JMX,IMX,TE,'CODIERT ',1)
C
C
 6050 FORMAT(6(E12.5E3,1X))
C
C
C                                    **********************************
C
C                                    AB HIER VERSION NACH 6.4.92
C
C                                    **********************************
      ELSEIF (IIDENT(5).GT.8) THEN

            KANAL = 3

            CT = CH1680 (' U              ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,'CODIERT ',
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   U      ,CT,NBND,HILF,IGRID,IDIM3D)
            CT = CH1680 (' V              ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,'CODIERT ',
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   V      ,CT,NBND,HILF,IGRID,IDIM3D)
            CT = CH1680 (' W              ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,'CODIERT ',
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   W      ,CT,NBND,HILF,IGRID,IDIM3D)
            CT = CH1680 (' P              ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,'CODIERT ',
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   P      ,CT,NBND,HILF,IGRID,IDIM3D)
            CT = CH1680 (' G              ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,'CODIERT ',
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   G      ,CT,NBND,HILF,IGRID,IDIM3D)
            CT = CH1680 (' B              ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,'CODIERT ',
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   B      ,CT,NBND,HILF,IGRID,IDIM3D)

      ELSE
      STOP 'UNBEKANNTE VERSIONSNUMMER IN DIC'
      ENDIF
C
C
      RETURN
C
      END
