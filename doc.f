










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
      SUBROUTINE DOC (CIDENT,IIDENT,RIDENT,KK,JJ,II,KMX,JMX,IMX,NBND,
     $                DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,NVAR,
     $                U,V,W,P,G,B,TK,TE,HILF,IGRID,IDIM3D)
C*STAR******************************************************************
C  D O C          AUSGABE VON DATEN (KODIERT,KANAL4)
C*STAR********************************************** W.STEITZ 20.08.81 *
C                                         GEAENDERT: ST.FRANK 20.10.83 *
C                                                    ST.FRANK 23.05.84 *
C                                                    ST.FRANK 02.07.84 *
C                                                    F.BAETKE 08.10.84 *
C                                                    H.WERNER 19.09.85 *
C                                                    H.WERNER 14.11.85 *
C                                                    M.MANHART 6. 4.92
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
C             TK(KK,JJ,II)   - DISSIPATIONSFELD
C             TE(KK,JJ,II)   - TURBULENZENERGIEFELD
C
C  ANMERKUNG  CIDENT(10) ENTHAELT KENNUNG
C
C             CIDENT( 1) = "F.BAETKE"
C             CIDENT( 2) = PROGRAMMNAME
C             CIDENT( 3) = LAUFENDE NUMMER
C             CIDENT( 4) = DATE()
C             CIDENT( 5) = CLOCK()
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
C  UPROG                 : ERRR
C
C  VERS:  18.12.86 (HW)  : TIDENT IN CIDENT (CHARACTER (LEN=8) !!) UMGE-
C                          WANDELT WEGEN NOS/VE.
C        27.02.92 (MM)  : BELEGUNG der UNDEFINIERTEN WERTE MIT 
C                         GREAT AUSGESCHALTET
C     26. 5.95 (MM) : MPI EINGEFUEHRT, KOPF UM HILF ERWEITERT
C     11.02.03 (TB) : SCALAR TRANSPORT IMPLEMENTED
C
C*STAR******************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
      CHARACTER (LEN=8)   CIDENT(10)
      CHARACTER (LEN=80)  CT,  CH1680, TEXT
      INTEGER KK, JJ, II, KMX, JMX, IMX, NVAR, IIDENT(100)
      REAL    RIDENT(100), X(II), Y(JJ), Z(KK),
     $        DX(II), DY(JJ), DZ(KK), DDX(II), DDY(JJ), DDZ(KK),
     $        U(KK,JJ,II), V(KK,JJ,II), W(KK,JJ,II),
     $        P(KK,JJ,II), G(KK,JJ,II), TK(KK,JJ,II), TE(KK,JJ,II),
     $        B(KK,JJ,II), HILF(KK,JJ,II)
      DATA    TEXT /'GRID            '/
C
C					UPDATE DER VERSIONSNUMMER
C
      IIDENT(5) = IIDENT(6)
C
C
C
C
      WRITE(4,6010) TEXT
      WRITE(4,6030) IGRID
      WRITE(4,6030) KK,JJ,II,KMX,JMX,IMX
 6010 FORMAT (A80)
 6030 FORMAT(6(I9,1X))
C
      WRITE(4,6050)(X(I),I=1,IMX),(DX(I),I=1,IMX),(DDX(I),I=1,IMX),
     $             (Y(J),J=1,JMX),(DY(J),J=1,JMX),(DDY(J),J=1,JMX),
     $             (Z(K),K=1,KMX),(DZ(K),K=1,KMX),(DDZ(K),K=1,KMX)



C
C                              ABPRUEFUNG AUF RINDEF BZW: LEGVAR
C                              AUSGESCHALTET



C
C                                    **********************************
C
C                                    AB HIER VERSION VOR 6.4.92
C
C                                    **********************************
	IF (IIDENT(6).LT.8) THEN
C
      WRITE(4,6040) NVAR
      IF (NVAR .LT. 4 .OR. NVAR .GT. 7) CALL ERRR(704,' DOC      ')
 6040 FORMAT(I9,1X)
C
      DO 10 I=1,IMX
   10 WRITE(4,6050) ((U(K,J,I),K=1,KMX),J=1,JMX)
C
      DO 20 I=1,IMX
   20 WRITE(4,6050) ((V(K,J,I),K=1,KMX),J=1,JMX)
C
      DO 30 I=1,IMX
   30 WRITE(4,6050) ((W(K,J,I),K=1,KMX),J=1,JMX)
C
      DO 40 I=1,IMX
   40 WRITE(4,6050) ((P(K,J,I),K=1,KMX),J=1,JMX)
      IF (NVAR .EQ. 4) RETURN
C
      DO 50 I=1,IMX
   50 WRITE(4,6050) ((G(K,J,I),K=1,KMX),J=1,JMX)
      IF (NVAR .EQ. 5) RETURN
C
      DO 60 I=1,IMX
   60 WRITE(4,6050) ((TK(K,J,I),K=1,KMX),J=1,JMX)
      IF (NVAR .EQ. 6) RETURN
C
      DO 70 I=1,IMX
   70 WRITE(4,6050) ((TE(K,J,I),K=1,KMX),J=1,JMX)
C
C
C
C                                    **********************************
C
C                                    AB HIER VERSION NACH 6.4.92
C
C                                    **********************************
      ELSEIF (IIDENT(6).GE.8) THEN

            KANAL = 4

            CT = CH1680 (' U              ')
            CALL DOBC   (KK ,JJ ,II ,KMX,JMX,IMX,NBND,KANAL,'CODIERT ',
     $                   0    ,0    ,1    ,
     $                   1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0,
     $                   U      ,CT,HILF,IGRID,IDIM3D)
            CT = CH1680 (' V              ')
            CALL DOBC   (KK ,JJ ,II ,KMX,JMX,IMX,NBND,KANAL,'CODIERT ',
     $                   0    ,1    ,0    ,
     $                   1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0,
     $                   V      ,CT,HILF,IGRID,IDIM3D)
            CT = CH1680 (' W              ')
            CALL DOBC   (KK ,JJ ,II ,KMX,JMX,IMX,NBND,KANAL,'CODIERT ',
     $                   1    ,0    ,0    ,
     $                   1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0,
     $                   W      ,CT,HILF,IGRID,IDIM3D)
            CT = CH1680 (' P              ')
            CALL DOBC   (KK ,JJ ,II ,KMX,JMX,IMX,NBND,KANAL,'CODIERT ',
     $                   0    ,0    ,0    ,
     $                   1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0,
     $                   P      ,CT,HILF,IGRID,IDIM3D)
            CT = CH1680 (' G              ')
            CALL DOBC   (KK ,JJ ,II ,KMX,JMX,IMX,NBND,KANAL,'CODIERT ',
     $                   0    ,0    ,0    ,
     $                   1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0,
     $                   G      ,CT,HILF,IGRID,IDIM3D)
            CT = CH1680 (' B              ')
            CALL DOBC   (KK ,JJ ,II ,KMX,JMX,IMX,NBND,KANAL,'CODIERT ',
     $                   0    ,0    ,0    ,
     $                   1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0 ,1.0,
     $                   B      ,CT,HILF,IGRID,IDIM3D)

      ELSE
      STOP 'UNBEKANNTE VERSIONSNUMMER IN DOC'
      ENDIF
C
 6050 FORMAT(6(E12.5E3,1X))
C
      RETURN
C
      END
