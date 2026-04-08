










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
      SUBROUTINE SELACO  (ILIMX, IE1, IE2, IE3, IE4, IE5, IL0, ILX,
     $                    ILY, ILZ, ISLINP, ISLIDI, ISLINI, ILIMXP)
C*STARLET***************************************************************
C        S E L A C O      SELACO ENTHAELT DIE MELDUNGEN AUS SUBR.
C                         SELAUF
C*STARLET***************************************************************
C
C PARAM:                  SIEHE SUBR. SELAUF
C UPROG                 : ERRR,  KJIDCO
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        07.10.88 (HW)  : ORIGINAL
C        24.11.88 (HW)  : KONTROLLAUSDRUCK EINGEBAUT
C        24.02.89 (HW)  : ISLINI(ILIMXP) --> ISLINI(ISLIDI)
C
C*STARLET***************************************************************
C
C
      CHARACTER (LEN=1)  CDIRC
      INTEGER        ISLINP (ISLIDI),  ISLINI (ISLIDI)
      WRITE (6,6010)
      WRITE (6,*)
      WRITE (6,*) '   ******************* INFORMATION AUS SUBR. SELAUF',
     $            ' *******************'
      WRITE (6,*)
      WRITE (6,*) '   IM PARAMETER-STATEMENT KOENNTE  ILIMXP  AUF EINEN'
      WRITE (6,*) '   WERT VON ',ILIMX,' GESETZT WERDEN.'
      WRITE (6,*)
      WRITE (6,*) '   ANZAHL DER AUFPUNKTE MIT RICHTUNGSINFORMATION',
     $            '  0  : ',IL0-1
      WRITE (6,*) '   ANZAHL DER AUFPUNKTE MIT RICHTUNGSINFORMATION',
     $            '  X  : ',ILX-1
      WRITE (6,*) '   ANZAHL DER AUFPUNKTE MIT RICHTUNGSINFORMATION',
     $            '  Y  : ',ILY-1
      WRITE (6,*) '   ANZAHL DER AUFPUNKTE MIT RICHTUNGSINFORMATION',
     $            '  Z  : ',ILZ-1
C
      IF(ILIMX .LE. 1) THEN
         WRITE (6,*) ' PROGRAMMABBRUCH, DA KEINE AUFPUNKTE',
     $               ' IN SUBR. SELAUF DEFINIERT WURDEN !!'
         STOP ' KEINE AUFPUNKTE DEFINIERT (SUBR. SELAUF) !'
      END IF
      IF(IE1 .EQ. 1) THEN
         WRITE (6,*) ' PROGRAMMABBRUCH WEGEN NEGATIVEM K-, J-,',
     $               ' ODER I-INDEXES IN SUBR. SELAUF'
         STOP ' NEGATIVER INDEX IN SUBR. SELAUF !'
      ELSE
         IF(IE1 .NE. 0) CALL ERRR (501,' SELAUF   ')
      END IF
C
      IF(IE2 .EQ. 1) THEN
         WRITE (6,*) ' PROGRAMMABBRUCH WEGEN ZU GROSSEM K-, J-,',
     $               ' ODER I-INDEXES IN SUBR. SELAUF'
         STOP ' ZU GROSSER INDEX IN SUBR. SELAUF !'
      ELSE
         IF(IE2 .NE. 0) CALL ERRR (502,' SELAUF   ')
      END IF
C
      IF(IE3 .EQ. 1) THEN
         WRITE (6,*) ' PROGRAMMABBRUCH, DA DER SPEICHERPLATZ',
     $               ' FUER DIE ANGEGEBENE ZAHL VON AUFPUNKTEN'
         WRITE (6,*) ' (ANZ. DER AUFPUNKTE = ANZ. DER I-LINIEN = ',
     $               ILIMX,') NICHT AUSREICHT !'
         WRITE (6,*) ' IM PARAMETER-STATEMENT MUESSTE  ISLIDI  ',
     $               'MINDESTENS AUF ',ILIMX,' GESETZT WERDEN !'
         WRITE (6,*) ' DIES IST MIT  E R H E B L I C H E M  ',
     $               'UMSTELLUNGSAUFWAND VERBUNDEN.'
         WRITE (6,*) ' DER BENUTZER SOLLTE DAHER DIE ZAHL DER',
     $               ' AUFPUNKTE REDUZIEREN.'
         STOP ' ISLIDI IM PARAMETER-STATEMENT ZU KLEIN !'
      ELSE
         IF(IE3 .NE. 0) CALL ERRR (503,' SELAUF   ')
      END IF
C
      IF(IE4 .EQ. 1) THEN
         WRITE (6,*) ' PROGRAMMABBRUCH, DA DER SPEICHERPLATZ',
     $               ' FUER DIE ANGEGEBENE ZAHL VON AUFPUNKTEN'
         WRITE (6,*) ' (ANZ. DER AUFPUNKTE = ANZ. DER I-LINIEN = ',
     $               ILIMX,') NICHT AUSREICHT !'
         WRITE (6,*) ' IM PARAMETER-STATEMENT MUSS  ILIMXP  MINDESTENS',
     $               ' AUF ',ILIMX,' GESETZT WERDEN !'
         STOP ' ILIMXP IM PARAMETER-STATEMENT ZU KLEIN !'
      ELSE
         IF(IE4 .NE. 0) CALL ERRR (504,' SELAUF   ')
      END IF
C
      IF(IE5 .EQ. 1) THEN
         WRITE (6,*) ' PROGRAMMABBRUCH, DA EINE ODER MEHRERE',
     $               ' RICHTUNGSINFORMATIONEN FALSCH SIND !'
         WRITE (6,*) ' BEI DEN OBEN ANGEGEBENEN AUFPUNKTEN IST DIE',
     $               ' VARIABLE  CD  ZU KORRIGIEREN'
         STOP ' SUBR. SELAUF: RICHTUNGSINFORMATION FALSCH !'
      ELSE
         IF(IE5 .NE. 0) CALL ERRR (505,' SELAUF   ')
      END IF
C
C                                 HIER: KONTROLLAUSDRUCK DER VOM
C                                 PROGRAMM AKZEPTIERTEN AUFPUNKTE
C
      WRITE (6,*)
      WRITE (6,*)
      WRITE (6,6020)
      WRITE (6,6030)
C
C                                 DO 100: DURCHBLAETTERN DES INHALTSVER-
C                                 ZEICHNISSES DER AUFPUNKTE
C
      IF(ISLIDI .GE. 2) THEN
         DO 100 ISLIN = 2,ISLIDI
            CALL KJIDCO (ISLINP(ISLIN),CDIRC,IDIRC,KAUFP,JAUFP,IAUFP)
C
C                                 FALLS DER LETZTE EINTRAG IM INHALTS-
C                                 VERZEICHNIS GELESEN WURDE --> RETURN
C
            IF(CDIRC .EQ. '0'  .AND.  KAUFP .EQ. 0  .AND.
     $         JAUFP .EQ.  0   .AND.  IAUFP .EQ. 0)       GOTO 9999
            IF(CDIRC .EQ. '0') THEN
               WRITE (6,6110) ISLINI (ISLIN), CDIRC, KAUFP, JAUFP,
     $                        IAUFP, ISLIN
            END IF
            IF(CDIRC .EQ. 'X') THEN
               WRITE (6,6120) ISLINI (ISLIN), CDIRC, KAUFP, JAUFP,
     $                        IAUFP, ISLIN
            END IF
            IF(CDIRC .EQ. 'Y') THEN
               WRITE (6,6130) ISLINI (ISLIN), CDIRC, KAUFP, JAUFP,
     $                        IAUFP, ISLIN
            END IF
            IF(CDIRC .EQ. 'Z') THEN
               WRITE (6,6140) ISLINI (ISLIN), CDIRC, KAUFP, JAUFP,
     $                        IAUFP, ISLIN
            END IF
  100    CONTINUE
      END IF
 9999 RETURN
 6010 FORMAT (1H1)
 6020 FORMAT (4X,'----- UEBERSICHT UEBER DIE VOM PROGRAMM',
     $        ' AKZEPTIERTEN AUFPUNKTE -----',/)
 6030 FORMAT (4X,'NR. DES AUFPUNKTES  RICHTUNGS-',8X,'INDIZES',
     $        12X,'I-LINIE',/,5X,'IN SUBR. SELAUF',4X,
     $        'INFORMATION    DES AUFPUNKTES',5X,
     $        'IM ISLINP-FELD',/,39X,'K',5X,'J',5X,'I')
 6110 FORMAT (11X,I4, 8X,A1,12X,3(2X,I4),7X,I4)
 6120 FORMAT (11X,I4,12X,A1, 8X,3(2X,I4),7X,I4)
 6130 FORMAT (11X,I4,16X,A1, 4X,3(2X,I4),7X,I4)
 6140 FORMAT (11X,I4,20X,A1,    3(2X,I4),7X,I4)
      END
