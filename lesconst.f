










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
       SUBROUTINE LESCONST (KK,JJ,II,KMX,JMX,IMX,
     $                              X,Y,Z,
     $                              DX,DY,DZ,
     $                              DDX,DDY,DDZ,
     $                              CONV1S,CONV1SANF,CONV1SEND,
     $                              TRANSLES1,TRANSLES2)
C*MGLET***************************************************************
C        LESCONST     CALCULATION OF THE SUBGRID SCALE CONSTANT
C                         AND STORING IN CONV1S-FIELD
C*MGLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        X, Y, Z        - KOORDINATEN DER ZELLDEFINITIONSPUNKTE
C        DX,DY,DZ       - ABSTAND DER GITTERPUNKTE
C        DDX,DDY,DDZ    - KANTENLAENGEN DER KONTROLLVOLUMINA
C
C
C VERS:  12.09.96 (A.O.): ORIGINAL
C
C
C*MGLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

      COMMON /CONLES/  CONV2S,CMUE,CAPPA,ECONST
      SAVE   /CONLES/
C
C
      REAL        X(II),         Y(JJ),         Z(KK),
     $           DX(II),        DY(JJ),        DZ(KK),
     $          DDX(II),       DDY(JJ),       DDZ(KK)
C
C       SUBGRIDSCALE CONSTANT
      REAL   CONV1S(II)
C
C
C
      KM1    = KMX-1
      JM1    = JMX-1
      IM1    = IMX-1
      KM2    = KMX-2
      JM2    = JMX-2
      IM2    = IMX-2
C
C
C       GEOMETRY CONSTANT TO CALCULATE THE 
C         TRANSITION INTERPOLATION OF LES

      TRANSCON = 1./(TRANSLES2-TRANSLES1)

      DO I  = 1,IMX
         IF (X(I) .LE. TRANSLES1) THEN
           CONV1S(I) = CONV1SANF
        
         ELSEIF (X(I).LT.TRANSLES2 .AND. X(I).GT.TRANSLES1) THEN
           CONV1S(I) = ((TRANSLES2-X(I))*CONV1SANF+
     &               (X(I)-TRANSLES1)*CONV1SEND)*TRANSCON
         
         ELSEIF (X(I) .GE. TRANSLES2) THEN
           CONV1S(I) = CONV1SEND

         ENDIF
       ENDDO

       WRITE(44,*) '********************'
       WRITE(44,*) 'SUBGRID SCALE CONSTANT CONV1S(I)'
       WRITE(44,*) (CONV1S(I),I=1,IMX,1)
       WRITE(44,*) '********************'

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      RETURN
      END
