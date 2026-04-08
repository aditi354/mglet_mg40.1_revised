










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
      SUBROUTINE BF1CUB  (KK,JJ,II,KMX,JMX,IMX,PHI,
     $                    IB1,IB2,JB1,JB2,KB,IVAR)
C*STARLET***************************************************************
C        B F 1 C U B      BF1CUB SETZT DIE RANDBEDINGUNGEN FUER FLUK-
C                         TUIERENDE UND HIERVON ABGELEITETE GROESSEN
C                         INNERHALB DES KUBUSSES.
C                         FUER DIE U-, V- UND W-KOMPONENTE, SOWIE FUER
C                         DEN DRUCK SIND DIES:
C                         PHIFLUK., PHIRMS, PHISKEW, PHIFLAT.
C                         AUSSERDEM KOENNEN DIE RANDBED. FUER DIE KINET.
C                         ENERGIE DER SCHWANKUNGSGESCHW. AUS GROB- UND
C                         FEINSTRUKTUR GESETZT WERDEN.
C                         FUER VARIABLE, DIE IN MEHR ALS EINER RICHTUNG
C                         IM MASCHENGITTER VERSCHOBEN SIND, KANN DIE
C                         ROUTINE  IM ALLGEMEINEN AUCH  VERWENDET
C                         WERDEN. EINE GENAUE PRUEFUNG IST JEDOCH RATSAM
C
C                         A C H T U N G: SOLLEN RANDBEDINGUNGEN FUER DIE
C                         MOMENTANWERTE DER GESCHW. GESETZT WERDEN, SO
C                         IST ZU PRUEFEN, OB BEI VERWENDUNG DES QUICK-
C                         VERFAHRENS KEINE GROESSEREN FEHLER AUFTRETEN.
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        PHI(KK,JJ,II)  + ALLGEMEINE VARIABLE
C        IB1, IB2       - GRENZE DES KUBUSSES IN X-RI.
C        JB1, JB2       -   "           "     IN Y-RI.
C        KB             -   "           "     IN Z-RI. (TOP-FLAECHE)
C        IVAR           - CHARACTER-VARIABLE, DIE DIE GROESSE PHI NAEHER
C                         SPEZIFIZIERT
C
C UPROG                 : ERRR
C
C DEFINE-DIREKTIVEN     : CUSLI, CUNOS
C
C        25.07.86 (HW)  : ORIGINAL
C        14.08.86 (HW)  : AUCH FUER U, V UND W KOENNEN NOSLIP-RB.
C                         GESETZT WERDEN
C        27.08.86 (HW)  : FUER DIE SCHUBSPANNUNGEN 'UW', 'VW' UND 'UV'
C                         KOENNEN EBENFALLS RANDBEDINGUNGEN GESETZT
C                         WERDEN.
C        11.09.86 (HW)  : KORREKTUREN DER SCHLEIFENGRENZEN
C        19.08.88 (HW)  : BF1CUB AUCH FUER DIE KOMPONENTEN DER VORTICITY
C                         VERWENDBAR
C        19.12.89 (HW)  : BF1CUB SETZT FUER 'P' AN FESTEN WAENDEN
C                         DP/DN = 0.0
C
C         3. 4.93 (MM)  : RANDBEDINGUNGEN WERDEN FUER STATISTISCHE
C                         GROESSEN BIS AUF WEITERES NICHT MEHR GESETZT
C
C
C*STARLET***************************************************************
C
      CHARACTER *6  IVAR
      CHARACTER (LEN=1)  CC1
      CHARACTER (LEN=2)  CC2
C
      REAL     PHI(KK,JJ,II)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                     3. 4.93    RETURN EINGEBAUT!!!!!
      RETURN
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C                                 UNTER BESTIMMTEN (NORMALERWEISE NICHT
C                                 AUFTRETENDEN) BEDINGUNGEN MUESSEN
C                                 DIE GRENZEN DER DO-SCHLEIFEN UEBER-
C                                 PRUEFT WERDEN.
C
      IF(KB      .LE. 3) WRITE (6,6010)
      IF(JB2-JB1 .LT. 1) WRITE (6,6020)
      IF(IB2-IB1 .LT. 1) WRITE (6,6030)
C
C                                 CHECK DER VARIABLEN AUF ZULAESSIGKEIT
C
      CC1    = IVAR(1:1)
      IF(CC1 .EQ. 'U' .OR. CC1 .EQ. 'V' .OR. CC1 .EQ. 'W' .OR.
     $   CC1 .EQ. 'P' .OR. CC1 .EQ. 'E')                       GOTO 2010
         WRITE (6,6040) IVAR
         CALL ERRR (501,' BF1CUB   ')
 2010 CC2    = IVAR(2:3)
      IF(CC2 .EQ. '  ' .OR. CC2 .EQ. 'F ' .OR. CC2 .EQ. 'R ' .OR.
     $   CC2 .EQ. 'FL' .OR. CC2 .EQ. 'SK' .OR. CC2 .EQ. 'U ' .OR.
     $   CC2 .EQ. 'V ' .OR. CC2 .EQ. 'W ')                     GOTO 2020
         WRITE (6,6040) IVAR
         CALL ERRR (502,' BF1CUB   ')
 2020 CC1    = IVAR(4:4)
      IF(CC1 .EQ. ' ' .OR. CC1 .EQ. 'G' .OR. CC1 .EQ. 'S')     GOTO 2030
         WRITE (6,6040) IVAR
         CALL ERRR (503,' BF1CUB   ')
C
 2030 KM2    = KMX - 2
      JM2    = JMX - 2
      IM2    = IMX - 2
C
C                                 BESTIMMUNG DES DEFINITIONSPUNKTES DER
C                                 VARIABLEN INNERHALB DER ZELLE
C                                 (STAGGERED MESH !)
C
      KSTAG  = 0
      JSTAG  = 0
      ISTAG  = 0
C
      IF(IVAR(1:1) .EQ. 'W' ) KSTAG = 1
      IF(IVAR(1:1) .EQ. 'V' ) JSTAG = 1
      IF(IVAR(1:1) .EQ. 'U' ) ISTAG = 1
      IF(IVAR(2:3) .EQ. 'W ') KSTAG = 1
      IF(IVAR(2:3) .EQ. 'V ') JSTAG = 1
      IF(IVAR(2:3) .EQ. 'U ') ISTAG = 1
C
C                                 FALLS SICH DER KUBUS BIS ZUM RAND ER-
C                                 STRECKT MUESSEN DIE GRENZEN DER DO-
C                                 SCHLEIFEN DIESER SITUATION ANGEPASST
C                                 WERDEN
C
      IV1    = 0
      IV2    = 0
      JV1    = 0
      JV2    = 0
      IF(ISTAG .EQ. 1 .AND. IB1 .LE.   3) IV1 = 1
      IF(ISTAG .EQ. 1 .AND. IB2 .GE. IM2) IV2 = 1
      IF(JSTAG .EQ. 1 .AND. JB1 .LE.   3) JV1 = 1
      IF(JSTAG .EQ. 1 .AND. JB2 .GE. JM2) JV2 = 1
C
C                                 RANDBEDINGUNGEN AN DEN SEITENFLAECHEN
C                                 DES KUBUSSES
C                                 -------------------------------------
C
C                                 HIER: FRONT UND BACK
C
      IF(ISTAG .EQ. 0) THEN
C
         DO 10 J = JB1-JV1,JB2-JSTAG+JV2
            IF(IB1 .GT.   3) THEN
               DO 20 K = 3,KB-KSTAG
   20          CONTINUE
            ENDIF
C
            IF(IB2 .LT. IM2) THEN
               DO 30 K = 3,KB-KSTAG
   30          CONTINUE
            ENDIF
   10    CONTINUE
      ENDIF
C
C                                 HIER: LEFT UND RIGHT
C
      IF(JSTAG .EQ. 0) THEN
C
         DO 110 I = IB1-IV1,IB2-ISTAG+IV2
            IF(JB1 .GT.   3) THEN
               DO 120 K = 3,KB-KSTAG
  120          CONTINUE
            ENDIF
C
            IF(JB2 .LT. JM2) THEN
               DO 130 K = 3,KB-KSTAG
  130          CONTINUE
            ENDIF
  110    CONTINUE
      ENDIF
C
C                                 HIER: TOP
C
      IF(KSTAG .EQ. 0) THEN
         IF(KB  .LT. KM2) THEN
            DO 210 I = IB1-IV1,IB2-ISTAG+IV2
               DO 220 J = JB1-JV1,JB2-JSTAG+JV2
  220          CONTINUE
  210       CONTINUE
         ENDIF
      ENDIF
C
C                                 RANDBEDINGUNGEN AN DEN KANTEN DES
C                                 KUBUSSES
C                                 ---------------------------------
C
      IF(ISTAG .EQ. 0 .AND. JSTAG .EQ. 0) THEN
C
C                                 HIER: K-SAEULE MIT I = IB1 UND J = JB1
C
         IF(IB1 .GT.   3 .AND. JB1 .GT.   3) THEN
C
            DO 310 K = 3,KB-KSTAG
  310       CONTINUE
         ENDIF
C
C                                 HIER: K-SAEULE MIT I = IB2 UND J = JB1
C
         IF(IB2 .LT. IM2 .AND. JB1 .GT.   3) THEN
C
            DO 320 K = 3,KB-KSTAG
  320       CONTINUE
         ENDIF
C
C                                 HIER: K-SAEULE MIT I = IB2 UND J = JB2
C
         IF(IB2 .LT. IM2 .AND. JB2 .LT. JM2) THEN
C
            DO 330 K = 3,KB-KSTAG
  330       CONTINUE
         ENDIF
C
C                                 HIER: K-SAEULE MIT I = IB1 UND J = JB2
C
         IF(IB1 .GT.   3 .AND. JB2 .LT. JM2) THEN
C
            DO 340 K = 3,KB-KSTAG
  340       CONTINUE
         ENDIF
      ENDIF
C
      IF(ISTAG .EQ. 0 .AND. KSTAG .EQ. 0) THEN
C
C                                 HIER: J-SAEULE MIT I = IB1 UND K = KB
C
         IF(IB1 .GT.   3 .AND. KB  .LT. KM2) THEN
C
            DO 350 J = JB1-JV1,JB2-JSTAG+JV2
  350       CONTINUE
         ENDIF
C
C                                 HIER: J-SAEULE MIT I = IB2 UND K = KB
C
         IF(IB2 .LT. IM2 .AND. KB  .LT. KM2) THEN
C
            DO 360 J = JB1-JV1,JB2-JSTAG+JV2
  360       CONTINUE
         ENDIF
      ENDIF
C
      IF(JSTAG .EQ. 0 .AND. KSTAG .EQ. 0) THEN
C
C                                 HIER: I-SAEULE MIT J = JB1 UND K = KB
C
         IF(JB1 .GT.   3 .AND. KB  .LT. KM2) THEN
C
            DO 370 I = IB1-IV1,IB2-ISTAG+IV2
  370       CONTINUE
         ENDIF
C
C                                 HIER: I-SAEULE MIT J = JB2 UND K = KB
C
         IF(JB2 .LT. JM2 .AND. KB  .LT. KM2) THEN
C
            DO 380 I = IB1-IV1,IB2-ISTAG+IV2
  380       CONTINUE
         ENDIF
      ENDIF
C
C                                 RANDBEDINGUNGEN AN DEN ECKEN DES
C                                 KUBUSSES
C                                 --------------------------------
C
       IF(ISTAG .EQ. 0 .AND. JSTAG .EQ. 0 .AND. KSTAG .EQ. 0) THEN
C
          THIRD  = 1.0 / 3.0
C
C                                 HIER: ECKE FRONT-RIGHT-TOP
C
          IF(IB1 .GT.   3 .AND. JB1 .GT.   3 .AND. KB  .LT. KM2) THEN
C
          ENDIF
C
C                                 HIER: ECKE BACK-RIGHT-TOP
C
          IF(IB2 .LT. IM2 .AND. JB1 .GT.   3 .AND. KB  .LT. KM2) THEN
C
          ENDIF
C
C                                 HIER: ECKE BACK-LEFT-TOP
C
          IF(IB2 .LT. IM2 .AND. JB2 .LT. JM2 .AND. KB  .LT. KM2) THEN
C
          ENDIF
C
C                                 HIER: ECKE FRONT-LEFT-TOP
C
          IF(IB1 .GT.   3 .AND. JB2 .LT. JM2 .AND. KB  .LT. KM2) THEN
C
          ENDIF
      ENDIF
C
      RETURN
C
 6010 FORMAT (//,2X,' A C H T U N G:  KONTROLLIERE K - SCHLEIFEN',
     $        " IN SUBR. BF1CUB !!!!")
 6020 FORMAT (//,2X,' A C H T U N G:  KONTROLLIERE J - SCHLEIFEN',
     $        " IN SUBR. BF1CUB !!!!")
 6030 FORMAT (//,2X,' A C H T U N G:  KONTROLLIERE I - SCHLEIFEN',
     $        " IN SUBR. BF1CUB !!!!")
 6040 FORMAT (//,2X,' A C H T U N G !!! : FEHLERHAFTE UEBERGABE',
     $        ",  IVAR = ",A6)
      END
