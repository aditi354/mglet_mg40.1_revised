










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
      SUBROUTINE AMULT2  (PHI, KKNL, JJNL, ILIMXP,
     $                    CMULT, KB, KE, JB, JE, IB, IE)
C*STARLET***************************************************************
C        A M U L T 2      EIN DEFINIERTER BEREICH INNERHALB DES
C                         PHI-FELDES KANN MIT DER KONSTANTEN CMULT
C                         MULTIPLIZIERT WERDEN.
C*STARLET***************************************************************
C
C PARAM: PHI(KKNL,      + ALLGEMEINE VARIABLE
C        JJNL,ILIMXP)
C        KKNL,JJNL,     - ARRAYDIMENSIONEN
C        ILIMXP
C        CMULT          - MULTIPLIKATOR
C        KB ,KE         - START- UND ENDINDEX FUER DIE "K" DO-SCHLEIFE
C        JB ,JE         - START- UND ENDINDEX FUER DIE "J" DO-SCHLEIFE
C        IB ,IE         - START- UND ENDINDEX FUER DIE "I" DO-SCHLEIFE
C
C VERS:  15.12.88 (HW)  : ORIGINAL
C
C*STARLET***************************************************************
C
      REAL PHI (KKNL,JJNL,ILIMXP)
C
      DO 100 I = IB,IE
         DO 100 J = JB,JE
            DO 100 K = KB,KE
  100          PHI(K,J,I) = PHI(K,J,I) * CMULT
C
      RETURN
      END
