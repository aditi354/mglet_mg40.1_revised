










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
      SUBROUTINE SETKON
C*STAR******************************************************************
C*STAR*  S E T K O N      FESTLEGUNG ALLGEMEINER KONSTANTEN
C*STAR******************************************************************
C
C UPROG                 : KEINE
C
C DEFINE-DIREKTIVEN     : B7800,  CRAY,  CYBER,  C205,  C3820,  IRIS, SP2
C                         LINUX
C
C VERS:  03.01.85 (HW)  : ORIGINAL
C        14.11.85 (HW)  : DEFINE B7800  EINGEFUEHRT
C        30.07.86 (HW)  : DEFINE C205 EINGEFUEHRT
C        21.09.88 (HW)  : SMAONE  EINGEFUEHRT. SMAONE IST GERADE
C                         SO GROSS, DASS AUF JEDER MASCHINE
C                         1.0 + SMAONE > 1.0 GERADE NICHT MEHR
C                         ERKANNT WIRD. (SMAONE IST EIN MASS FUER DEN
C                         RUNDUNGSFEHLER DER MASCHINE)
C        30.12.88 (HW)  : UEBERARBEITUNG, KONTROLLAUSDRUCK
C        27.02.92 (MM)  : EINFUEHRUNG VON PRESET=0.0, AUCH IN KONSTA
C
C*STAR******************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
      DATA RINDE  /1.0E20/
C
      GREAT  = 1.0E40
      SMALL  = 1.0E-10
      RINDEF = RINDE
C
      PRESET = 0.0
C
C                                 MASCHINEN"UN"ABHAENGIGE BESTIMMUNG
C                                 VON SMAONE
C
      SMAONE   = 1.0
 2000 IF(1.0 + SMAONE .GT. 1.0) THEN
         SMAONE = 0.5*SMAONE
         GOTO 2000
      ENDIF
C
      WRITE (6,*)
      WRITE (6,*)
      WRITE (6,*) ' ********** INFORMATION AUS SUBR. SETKON **********'
      WRITE (6,*)
      WRITE (6,*) ' GREAT   (3. WURZEL AUS DER GROESSTEN, AUF DER ',
     $            'MASCHINE DARSTELLBAREN REAL-ZAHL) = ',GREAT
      WRITE (6,*) ' SMALL   (3. WURZEL AUS DER KLEINSTEN, AUF DER ',
     $            'MASCHINE DARSTELLBAREN REAL-ZAHL) = ',SMALL
      WRITE (6,*) ' RINDEF  (REAL-ZAHL MIT DER BEDEUTUNG "UNDEFINI',
     $            'ERT")                             = ',RINDEF
      WRITE (6,*) '         RINDEF IST MASCHINENABHAENGIG !!'
      WRITE (6,*) ' PRESET  (0.0) ZUR VORBELEGUNG ALLER GROESSEN  ',
     $            '                                  = ',PRESET
      WRITE (6,*) 'ALS ERSATZ VON RINDEF AB 27.02.92  '
      WRITE (6,*) ' SMAONE  (DEFINITION: 1.0 + SMAONE = 1.0)      ',
     $            '                                  = ',SMAONE
      WRITE (6,*) '         SMAONE IST MASCHINENABHAENGIG, KANN JE',
     $            'DOCH IM GEGENSATZ ZU "RINDEF"'
      WRITE (6,*) '         AUTOMATISCH BESTIMMT WERDEN'
      WRITE (6,*) ' OBIGE KONSTANTEN STEHEN IM COMMON-BLOCK "KONST',
     $            'A"'
      RETURN
      END
