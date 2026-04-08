










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
      SUBROUTINE LINOUT  (DREAD,  DWRITE, IWRB   )
C*STARLET***************************************************************
C        L I N O U T      IN LINOUT WERDEN INFORMATIONEN UEBER DAS EIN-
C                         LESEN UND AUSGEBEN DER DATEN ERZEUGT.
C*STARLET***************************************************************
C
C PARAM: DREAD                    - LOGISCHE VARIABLE:
C                                   .FALSE. : KEIN EINLESEN DER 3D-FEL-
C                                             DER U. D. STATISTIK
C                                   .TRUE.  : EINLESEN DER 3D-FEL-
C                                             DER U. D. STATISTIK
C        DWRITE                   - LOGISCHE VARIABLE:
C                                   .FALSE. : KEINE AUSGABE DER 3D-FEL-
C                                             DER U. D. STATISTIK
C                                   .TRUE.  : AUSGABE DER 3D-FEL-
C                                             DER U. D. STATISTIK
C        IWRB                     - INDEX DER I-SCHEIBE, DEREN VARIABLEN
C                                   AUSGEGEBEN WERDEN UM EINTRITTSPRO-
C                                   FILE ZU ERZEUGEN
C
C UPROG                 : KEINE
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        19.06.89 (HW)  : ORIGINAL
C
C*STARLET***************************************************************
C

      COMMON /CLINOU/

     $                 LDIB,   LDIC,   LDIEIB, LDIEIC,
     $                 LDOB,   LDOC,   LDOEIB, LDOEIC

      SAVE   /CLINOU/
      LOGICAL

     $                 LDIB,   LDIC,   LDIEIB, LDIEIC,
     $                 LDOB,   LDOC,   LDOEIB, LDOEIC
C
      INTEGER       IWRB
      LOGICAL       DREAD,  DWRITE
C
      WRITE (6,6000)
C
      IF(DREAD ) THEN
         IF(LDIB  ) THEN
            WRITE (6,6010) 'UNFORMATIERT (BINAER)',  1
         END IF
C
         IF(LDIC  ) THEN
            WRITE (6,6010) 'FORMATIERT (CODIERT) ',  3
         END IF
C
         IF(LDIB    .AND.  LDIC  ) THEN
            WRITE (6,6020)
         END IF
      END IF
C
C                                 LESEN DER EINTRITTSPROFILE
C
      IF(LDIEIB) THEN
         WRITE (6,6110) 'I',  2,  'UNFORMATIERT (BINAER)',  11
      END IF
C
      IF(LDIEIC) THEN
         WRITE (6,6110) 'I',  2,  'FORMATIERT (CODIERT) ',  13
      END IF
C
      IF(LDIEIB  .AND.  LDIEIC) THEN
         WRITE (6,6020)
      END IF
C
      WRITE (6,6210)
C
C                                 AUSGABE DER EINTRITTSPROFILE
C                                 DIMENSIONSLOS AUSGEGEBEN
C
      IF(LDOEIB) THEN
         WRITE (6,6220) 'I',  IWRB,  IWRB+1,  'UNFORMATIERT (BINAER)',
     $                  12
      END IF
C
      IF(LDOEIC) THEN
         WRITE (6,6220) 'I',  IWRB,  IWRB+1,  'FORMATIERT (CODIERT) ',
     $                  14
      END IF
C
      IF(LDOEIB  .AND.  LDOEIC) THEN
         WRITE (6,6230)
      END IF
C
      IF(DWRITE) THEN
         IF(LDOB  ) THEN
            WRITE (6,6240) 'UNFORMATIERT (BINAER)',  2
         END IF
C
         IF(LDOC  ) THEN
            WRITE (6,6240) 'FORMATIERT (CODIERT) ',  4
         END IF
C
         IF(LDOB    .AND.  LDOC  ) THEN
            WRITE (6,6230)
         END IF
      END IF
C
      WRITE (6,6310)
C
      RETURN
 6000 FORMAT (1H1,3X,10(1H*),'  INFORMATION UEBER EINLESEN UND ',
     $        'AUSGEBEN DER DATEN (SUBR. LINOUT)',//,4X,
     $        'E I N L E S E N :',/)
 6010 FORMAT (/,4X,'DIE DATEN FUER DIE ANFANGSBELEGUNG DER 3D-FELDER ',
     $        'UND DER STATISTIK WERDEN',/,4X,A21,' UEBER KANAL ',I2,
     $        ' EINGELESEN.')
 6020 FORMAT (/,4X,'IST ES UNBEDINGT NOETIG FORMATIERT  U N D  UN',
     $        'FORMATIERT EINZULESEN  ? ')
 6110 FORMAT (/,4X,'DIE DATEN FUER DIE EINTRITTSPROFILE IN DER ',A1,
     $        ' = ',I4,' EBENE WERDEN',/,4X,A21,' UEBER KANAL ',
     $        I2,' EINGELESEN.')
 6210 FORMAT (/,4X,'A U S G A B E  DER DATEN :',/)
 6220 FORMAT (/,4X,'DIE WAEHREND DES LAUFES ERZEUGTEN EINTRITTS',
     $        'PROFILE IN DEN ',A1,'-SCHEIBEN',/,4X,I4,' UND ',
     $        I4,' WERDEN ',A21,' UEBER KANAL ',I2,' AUSGEGEBEN.')
 6230 FORMAT (/,4X,'IST ES UNBEDINGT NOETIG, DIE DATEN FORMATIERT',
     $        '  U N D  UNFORMATIERT AUSZUGEBEN  ? ')
 6240 FORMAT (/,4X,'DIE DATEN DER 3D-FELDER UND DER STATISTIK',
     $        ' WERDEN ',A21,/,4X,'UEBER KANAL ',I2,' AUSGEGEBEN.')
 6310 FORMAT (/,4X,'FALLS DAS PARTIKELPROGRAMM MITLAEUFT, SO WIRD',
     $        ' UEBER',/,4X,'KANAL 24 UND 26 FORMATIERT GELESEN',
     $        ' UND GESCHRIEBEN.',/,4X,'EINZELHEITEN SIEHE ',
     $        'DIPLOMARBEIT VON G. EDER.',//,4X,'STANDARD I N P U T',
     $        '  ERFOLGT UEBER KANAL 5',/,4X,'STANDARD O U T P U T',
     $        '  ERFOLGT UEBER KANAL 6')
      END
