










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
      SUBROUTINE GRDCHE  (XTOT,   XCUB,   SX,     NBX,
     $                    YTOT,   YCUB,   SY,     NBY,
     $                    ZTOT,   ZCUB,   SZ,     NBZ    )
C*STARLET***************************************************************
C        G R D C H E      GRDCHE PRUEFT OB DIE ABMESSUNGEN DES BERECH-
C                         NUNGSGEBIETES WIE SIE DURCH XTOT, YTOT UND
C                         ZTOT GEGEBEN SIND MIT DEN ANGABEN FUER DEN
C                         GITTERGENERATOR GRDFMI UEBEREINSTIMMEN.
C*STARLET***************************************************************
C
C PARAM: STOT           - ABMESSUNG DES BERECHNUNGSGEBIETES IN S-RI.
C                         ENTSPRECHEND DER DATA-ANWEISUNG IM HAUPTPRO-
C                         GRAMM
C        SCUB           - ABMESSUNG DES KUBUSSES IN S-RI.
C                         ENTSPRECHEND DER DATA-ANWEISUNG IM HAUPTPRO-
C                         GRAMM
C        SS   (NBS)     - ENTHAELT FUER JEDEN TEILBEREICH DIE ZU UEBER-
C                         BRUECKENDE STRECKE, GEMESSEN VOM LINKEN RAND
C                         DER ERSTEN ZELLE INNER HALB DES TEILBEREICHES
C                         ZUM RECHTEN RAND DER LETZTEN ZELLE INNERHALB
C                         DES TEILBEREICHES
C        NBS            - ANZAHL DER BLOECKE (TEILBEREICHE) IN DIE DAS
C                         BERECHNUNGSGEBIET AUFGETEILT WURDE
C
C UPROG                 : KEINE
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        20.06.89 (HW)  : ORIGINAL
C
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
      REAL        SX (NBX ),       SY (NBY ),       SZ (NBZ )
C
C                                 UEBERPRUEFUNG DER DREI KOORDINATEN-
C                                 RICHTUNGEN
C
      EPSTOL = 1000.0 * SMAONE
C
      XTOTSO = 0.0
      DO 100 N = 1,NBX
  100    XTOTSO = XTOTSO + SX (N)
C
      IF( ABS (XTOT - XTOTSO) .GT. EPSTOL) THEN
         WRITE (6,6010) 'X', XTOT, 'X', NBX, 'X', XTOTSO
      END IF
C
      YTOTSO = 0.0
      DO 110 N = 1,NBY
  110    YTOTSO = YTOTSO + SY (N)
C
      IF( ABS (YTOT - YTOTSO) .GT. EPSTOL) THEN
         WRITE (6,6010) 'Y', YTOT, 'Y', NBY, 'Y', YTOTSO
      END IF
C
      ZTOTSO = 0.0
      DO 120 N = 1,NBZ
  120    ZTOTSO = ZTOTSO + SZ (N)
C
      IF( ABS (ZTOT - ZTOTSO) .GT. EPSTOL) THEN
         WRITE (6,6010) 'Z', ZTOT, 'Z', NBZ, 'Z', ZTOTSO
      END IF
      RETURN
 6010 FORMAT (//,4X,10(1H!),'  W A R N U N G  (SUBR. GRDCHE)  ',
     $        10(1H!),/,4X,A1,'TOT LAUT DATA-ANWEISUNG = ',1PE12.5,
     $        /,4X,A1,'TOT NACH DEFINITION IN GRDFMI (SUMME 1 ... ',
     $        I2,' S',A1,') = ',1PE12.5,/,4X,'ES  M U S S  UEBER',
     $        'EINSTIMMUNG ERZIELT WERDEN --> LAUF NACH KORREKTUR',
     $        ' WIEDERHOLEN')
      END
