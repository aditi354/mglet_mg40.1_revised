










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
      SUBROUTINE INTLIN  (PHI,KKNL,JJNL,ILIMXP,IBEG,IEND)
C*STARLET***************************************************************
C        I N T L I N      ZUR BESTIMMUNG DES INTEGRALEN LAENGENMASSES
C                         WIRD UEBER EINE "LINIE" MIT HILFE DER TRAPEZ-
C                         REGEL INTEGRIERT.
C*STARLET***************************************************************
C
C PARAM: PHI (KKNL,     + BELIEBIGES "LINIEN"-FELD
C        JJNL,ILIMXP)
C        KKNL,JJNL,     - ARRAYDIMENSIONEN
C        ILIMXP
C        IBEG,IEND      - GRENZEN D. DO-SCHLEIFEN F. I-RICHTUNG
C
C DEFINE DIREKTIVEN     : KEINE
C
C UPROG                 : KEINE
C
C        19.12.88 (HW)  : ORIGINAL
C        27.01.89 (HW)  : KORREKTUR
C
C*STARLET***************************************************************
C
      REAL            PHI (KKNL,JJNL,ILIMXP)
C
      DO 100 IL = IBEG,IEND
         KEND   = IFIX (PHI (KKNL-12, 1, IL)) - 1
         ALINT  = 0.0
         IF(KEND .GE. 1) THEN
            DO 200 K = 1,KEND
               ALINT   = ALINT + 0.5 * (PHI(K  ,1,IL) + PHI(K+1,1,IL))
     $                 *               (PHI(K+1,2,IL) - PHI(K  ,2,IL))
  200       CONTINUE
         END IF
  100 PHI(KKNL-22, 1, IL) = ALINT
C
      RETURN
      END
