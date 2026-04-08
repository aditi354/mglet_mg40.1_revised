










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
      SUBROUTINE ASDCHE  (KMX, JMX, IMX, KDIR, JDIR, IDIR, NBND,
     $                    LODD, LRETUR, IVAR,NXGRAE,NYGRAE,NZGRAE,
     $                    XHOMOG,YHOMOG,ZHOMOG)
C*STARLET***************************************************************
C        A S D C H E      IN ASDCHE WERDEN EINIGE PARAMETER UEBER-
C                         PRUEFT, DIE IN SUBR. ASDFUN WESENTLICH
C                         SIND.
C*STARLET***************************************************************
C
C PARAM: KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        K-, J-, IDIR   + NDIR = 0 : KEINE AUSWERTUNG DER POWER-SPEKTREN
C                                    IN  N-RICHTUNG
C                         NDIR = 1 : DAS POWER-SPEKTRUM WIRD IN
C                                    N-RICHTUNG AUSGEWERTET
C        NBND           - ANZAHL DER RANDSCHICHTEN DES BERECHNUNGS-
C                         GEBIETES
C        LODD           + LOGISCHE VARIABLE: FALLS DIE ANZAHL DER GITTER
C                         PUNKTE IN DER RICHTUNG, IN DER DAS POWER-
C                         SPEKTRUM GEBILDET WIRD, UNGERADZAHLIG IST,
C                         SO GILT  LODD = .TRUE.
C                         ANDERNFALLS  LODD = .FALSE.
C        LRETUR         + LOGISCHE VARIABLE: FALLS VOM BENUTZER FEHLER
C                         GEMACHT WURDEN, SO WIRD IN SUBR. ASDFUN
C                         EIN RETURN VERANLASST
C        IVAR           - CHARACTER-VARIABLE (DIENT HAUPTSAECHLICH ZUR
C                         IDENTIFIKATION DER RICHTUNG IN DER DAS
C                         POWER-SPEKTRUM GEBILDET WERDEN SOLL.
C
C UPROG                 : ERRR
C
C DEFINE-DIREKTIVEN     : XHOMOG,  YHOMOG,  ZHOMOG
C
C        03.11.88 (HW)  : ORIGINAL
C
C*STARLET***************************************************************
C
C
      CHARACTER (LEN=16)  IVAR
      LOGICAL        LODD,  LRETUR
      LOGICAL  XHOMOG,YHOMOG,ZHOMOG
C
      LRETUR = .FALSE.
C
C                                 CHECK DER CHARACTER-VARIABLEN
C
      IF(IVAR(1:1) .NE. 'S')                CALL ERRR (501,' ASDCHE   ')
      IF(IVAR(2:2) .NE. ' ')                CALL ERRR (502,' ASDCHE   ')
      IF(IVAR(3:3) .NE. 'Z'  .AND.  IVAR(3:3) .NE. 'Y'  .AND.
     $   IVAR(3:3) .NE. 'X')                CALL ERRR (503,' ASDCHE   ')
C
C                                 RICHTUNG FESTSTELLEN
C
      KDIR   = 0
      JDIR   = 0
      IDIR   = 0
      IF(IVAR(3:3) .EQ. 'Z') KDIR = 1
      IF(IVAR(3:3) .EQ. 'Y') JDIR = 1
      IF(IVAR(3:3) .EQ. 'X') IDIR = 1
      IF (XHOMOG) THEN
C
C                                 BEACHTE: DIE X-RICHTUNG IST HOMOGEN !
          NXRIHO = 1
      ELSE
C
C                                 DIE X-RICHTUNG IST  N I C H T  HOMOGEN
          NXRIHO = 0
      ENDIF
C
      IF (YHOMOG) THEN
C
C                                 BEACHTE: DIE Y-RICHTUNG IST HOMOGEN !
      NYRIHO = 1
      ELSE
C
C                                 DIE Y-RICHTUNG IST  N I C H T  HOMOGEN
      NYRIHO = 0
      ENDIF
C
      IF (ZHOMOG) THEN
C
C                                 BEACHTE: DIE Z-RICHTUNG IST HOMOGEN !
      NZRIHO = 1
      ELSE
C
C                                 DIE Z-RICHTUNG IST  N I C H T  HOMOGEN
      NZRIHO = 0
      ENDIF
C
      IF(KDIR .EQ. 1) THEN
C
C                                 HIER: DAS POWER-SPEKTRUM WIRD IN
C                                 Z-RICHTUNG GEBILDET
C
         IF(NZRIHO .NE. 1  .OR.  NZGRAE .NE. 1) THEN
            WRITE (6,6010)
            IF(NZRIHO .NE. 1) THEN
               WRITE (6,6020) IVAR(3:3)
            END IF
            IF(NZGRAE .NE. 1) THEN
               WRITE (6,6030) IVAR(3:3)
            END IF
C
            LRETUR = .TRUE.
         END IF
C
         GOTO 2000
      END IF
C
      IF(JDIR .EQ. 1) THEN
C
C                                 HIER: DAS POWER-SPEKTRUM WIRD IN
C                                 Y-RICHTUNG GEBILDET
C
         IF(NYRIHO .NE. 1  .OR.  NYGRAE .NE. 1) THEN
            WRITE (6,6010)
            IF(NYRIHO .NE. 1) THEN
               WRITE (6,6020) IVAR(3:3)
            END IF
            IF(NYGRAE .NE. 1) THEN
               WRITE (6,6030) IVAR(3:3)
            END IF
C
            LRETUR = .TRUE.
         END IF
C
         GOTO 2000
      END IF
C
      IF(IDIR .EQ. 1) THEN
C
C                                 HIER: DAS POWER-SPEKTRUM WIRD IN
C                                 X-RICHTUNG GEBILDET
C
         IF(NXRIHO .NE. 1  .OR.  NXGRAE .NE. 1) THEN
            WRITE (6,6010)
            IF(NXRIHO .NE. 1) THEN
               WRITE (6,6020) IVAR(3:3)
            END IF
            IF(NXGRAE .NE. 1) THEN
               WRITE (6,6030) IVAR(3:3)
            END IF
C
            LRETUR = .TRUE.
         END IF
C
         GOTO 2000
      END IF
C
      CALL ERRR (510,' ASDCHE   ')
C
 2000 CONTINUE
C
C                                 ANZAHL DER GITTERPUNKTE IN DER
C                                 RICHTUNG, IN DER DAS POWER-SPEKTRUM
C                                 GEBILDET WIRD (OHNE RANDSCHICHTEN)
C
      KJITOT = (KMX - 2*NBND) * KDIR
     $       + (JMX - 2*NBND) * JDIR
     $       + (IMX - 2*NBND) * IDIR
      KJIHAL = KJITOT / 2
C
      IF((KJITOT - 2*KJIHAL) .EQ. 1) THEN
         LODD = .TRUE.
      ELSE
         LODD = .FALSE.
      END IF
C
      RETURN
 6010 FORMAT (/,' ********** FEHLERMELDUNG AUS SUBR. ASDCHE (ASDFUN)',
     $        ' **********')
 6020 FORMAT (' DIE ',A1,'-RICHTUNG IST  K E I N E  HOMOGENE',
     $        ' RICHTUNG !')
 6030 FORMAT (' IN ',A1,'-RICHTUNG IST  DAS MASCHENGITTER',
     $        '  N I C H T  AEQUIDISTANT !')
      END
