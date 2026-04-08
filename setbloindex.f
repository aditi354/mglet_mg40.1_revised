










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
      SUBROUTINE SETBLOINDEX  (IIJJ,XY,PHYANF,PHYEND,
     $                         IJSTART,IJSTOP)
C*MGLET***************************************************************
C        SETBLOINDEX
C                         BERECHNUNG DER SCHLEIFENINDIZES FUER SCHLEIFE
C                         IN SUBROUTINE BBOBLO ZUR STROEMUNGS-
C                         MANIPULATION. UMSETZUNG DER GEBIETSGRENZEN
C                         IN SCHLEIFENINDIZES.
C
C*MGLET***************************************************************
C
C PARAM:  II     - ARRAYDIMENSIONEN
C         IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        X, Y, Z        - KOORDINATEN DER ZELLDEFINITIONSPUNKTE
C
C VERS:  25.02.94 (RK)  : ORIGINAL
C VERS:  29.09.95 (AO)  : EACH JET IS DEFINED SEPARATELY
C                         EACH DIRECTION IS TREATED WITH THE SAME SUBROUTINE
C   
C*MGLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

      REAL        XY(IIJJ)

C                       	VORBELEGUNG DER VARIABLEN ZUR BILDUNG
C 				VON MINIMALWERTEN

      XMINISTART = 99999999.99
      XMINISTOP  = 99999999.99
C     ADDED FOR CORRECT INITIALISATION
      ISTAMEMO = 1
      ISTOMEMO = IIJJ
C				SUCHE NACH DEN PUNKTEN, DIE XBANF, XBEND
C				YBANF UND YBEND AM MAECHSTEN LIEGEN

      DO I = 1,IIJJ

         DIFF   = XY(I)-PHYANF
         STOREA = ABS(DIFF)
         IF(STOREA.LT.XMINISTART) THEN
            ISTAMEMO    =      I 
            XMINISTART = STOREA
         ENDIF

         DIFF   = XY(I)-PHYEND
         STOREE = ABS(DIFF)
         IF(STOREE.LT.XMINISTOP ) THEN
            ISTOMEMO     =      I 
            XMINISTOP  = STOREE
         ENDIF

      ENDDO


C				FESTLEGUNG DER SCHLEIFENINDIZES
C                               FUER EINEN EINZELNEN JET
       IJSTART = ISTAMEMO
       IJSTOP  = ISTOMEMO


      RETURN
      END
