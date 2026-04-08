










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
      SUBROUTINE VSHIFT  (KK,JJ,II,KMX,JMX,IMX,VAR,ADDKON,ITYP)
C*STARLET***************************************************************
C        V S H I F T      FUER ITYP = 'A' : VAR = VAR + ADDKON
C                         FUER ITYP = 'S' : VAR = VAR - ADDKON
C*STARLET***************************************************************
C
C PARAM: KK,  JJ,  II   - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        VAR(KK,JJ,II)  + ALLGEMEINE VARIABLE Z.B. U, V, W ...
C        ADDKON         - KONSTANTE, DIE ZU VAR ADDIERT BZW.
C                         SUBTRAHIERT WERDEN SOLL
C        ITYP           - HOLLERITH-KONSTANTE
C                         ITYP = 'A' :   ADDKON WIRD ADDIERT
C                         ITYP = 'S' :   ADDKON WIRD SUBTRAHIERT
C
C DEFINE-DIREKTIVEN     : KEINE
C
C VERS:  15.11.85 (HW)  : ORIGINAL
C        21.02.89 (HW)  : NAMENSAENDERUNG VON  SHIFT --> VSHIFT, DA
C                         AUF DER CRAY EINE SYSTEMROUTINE NAMENS SHIFT
C                         EXISTIERT
C
C UPROG                 : ERRR
C
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
      CHARACTER (LEN=1)  ITYP
C
      REAL VAR(KK,JJ,II)
C
      IF(ABS(ADDKON) .LE. SMALL) RETURN
C
      IM2 = IMX-2
      JM2 = JMX-2
      KM2 = KMX-2
C
      IF(ITYP .NE. 'A') GOTO 2000
         SUM =  ADDKON
         GOTO 2100
C
 2000 IF(ITYP .NE. 'S') CALL ERRR(501,' VSHIFT   ')
         SUM = -ADDKON
C
 2100 DO 100 I=3,IM2
         DO 110 J=3,JM2
            DO 120 K=3,KM2
  120          VAR(K,J,I) = VAR(K,J,I) + SUM
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
