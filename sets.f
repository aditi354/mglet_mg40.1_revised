










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
      SUBROUTINE SETS(KK,JJ,II,KMX,JMX,IMX,VAR,CONST)
C*STAR*****************************************************************
C*STAR*  S E T S    VORBELEGUNG SKALARER FELDER MIT CONST
C*STAR*****************************************************************
C
C PARAM :VAR(KK,JJ,II)  + SKALARES FELD (P-,G-,TK-,TE-FELD)
C        KK,JJ,II       - ARRAYDIMENSIONEN
C        KMX,JMX,IMX    - GRENZEN DES BERECHNUNGSGEBIETS (MIT BOUND)
C        CONST          - KONSTANTER ANFANGSWERT
C
C VERS:  21.04.84 (SF)  - ORIGINAL
C        07.12.84 (HW)  : ZUNAECHST VORBELEGUNG DES GESAMTEN FELDES
C        07.01.85 (HW)  : GRENZEN DER DO-SCHLEIFEN GEAENDERT
C                         MIT RINDEF
C        27.02.92 (MM)  : RINDEF AUF PRESET GEAENDERT (0.0)
C
C*STAR*****************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
      REAL VAR(KK,JJ,II)
C

      IM1 = IMX-1
      IM2 = IMX-2
      JM1 = JMX-1
      JM2 = JMX-2
      KM1 = KMX-1
      KM2 = KMX-2
C
      DO 5 I=1,IMX
         DO 6 J=1,JMX
            DO 7 K=1,KMX
    7          VAR(K,J,I) = PRESET
    6    CONTINUE
    5 CONTINUE
C
      DO 40 I=2,IM1
         DO 50 J=2,JM1
            DO 60 K=2,KM1
   60          VAR(K,J,I) = CONST
   50    CONTINUE
   40 CONTINUE
C
      RETURN
      END
