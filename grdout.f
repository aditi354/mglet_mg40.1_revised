










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
      SUBROUTINE GRDOUT (IGRID,KMX,JMX,IMX,LEVEL,LCHILD,IPARENT,
     &                   IPOSITION, JPOSITION, KPOSITION,
     &                 XTOT,     YTOT,     ZTOT,
     &                 FRONT,    BACK,     RIGHT,     LEFT,
     &                 BOTTOM,   TOP,      CUBE,
     &                 XHOMOG,   YHOMOG,   ZHOMOG,
     &                 NXGRAE,   NYGRAE,   NZGRAE,
     &                 LPLEVEL)
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C     19. 5.93 (MM) : ORIGINAL
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C

      INTEGER
     &       KMX, JMX, IMX,LEVEL,IPARENT,IPOSITION,JPOSITION,KPOSITION,
     &       NXGRAE, NYGRAE, NZGRAE

      REAL
     &       XTOT,     YTOT,      ZTOT

      CHARACTER (LEN=16)
     &      FRONT,   BACK,
     &      RIGHT,   LEFT,
     &     BOTTOM,    TOP,
     &       CUBE

      LOGICAL
     &      LCHILD,
     &      XHOMOG,   YHOMOG,   ZHOMOG,
     &      LPLEVEL
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      WRITE (6,6000) 'GRD',IGRID,
     $       '      KMX:',KMX,'      JMX:',JMX,'      IMX:',IMX
      WRITE (6,6000) 'GRD',IGRID,
     $       '    LEVEL:',LEVEL

      IF (LCHILD) THEN
          WRITE (6,6002) 'GRD',IGRID,
     $           '    IS CHILD   OF GRID NR:',IPARENT
          WRITE (6,6000) 'GRD',IGRID,
     $           'POSITIONS:'
          WRITE (6,6000) 'GRD',IGRID,
     $           '     KPOS:',KPOSITION,
     $           '     JPOS:',JPOSITION,
     $           '     IPOS:',IPOSITION
      ENDIF

      WRITE (6,*) 
      WRITE (6,*) 
      WRITE (6,*) 

      WRITE (6,6001) 'GRD',IGRID,
     $       '     ZTOT:',ZTOT,'     YTOT:',YTOT,'     XTOT:',XTOT
      WRITE (6,*) 
      WRITE (6,*) 
      WRITE (6,*) 
      IF (XHOMOG) THEN
          WRITE (6,6002) 'GRD',IGRID,
     $       '  IS HOMOGENEOUS IN X-DIRECTION'
      ENDIF
      IF (YHOMOG) THEN
          WRITE (6,6002) 'GRD',IGRID,
     $       '  IS HOMOGENEOUS IN Y-DIRECTION'
      ENDIF
      IF (ZHOMOG) THEN
          WRITE (6,6002) 'GRD',IGRID,
     $       '  IS HOMOGENEOUS IN Z-DIRECTION'
      ENDIF
      WRITE (6,*) 
      WRITE (6,*) 
      WRITE (6,*) 
      IF (NXGRAE.EQ.1) THEN
          WRITE (6,6002) 'GRD',IGRID,
     $       '  IS AEQIDISTANT IN X-DIRECTION'
      ENDIF
      IF (NYGRAE.EQ.1) THEN
          WRITE (6,6002) 'GRD',IGRID,
     $       '  IS AEQIDISTANT IN Y-DIRECTION'
      ENDIF
      IF (NZGRAE.EQ.1) THEN
          WRITE (6,6002) 'GRD',IGRID,
     $       '  IS AEQIDISTANT IN Z-DIRECTION'
      ENDIF

      WRITE (6,*) 
      WRITE (6,*) 
      WRITE (6,*) 

          WRITE (6,6002) 'GRD',IGRID,
     $       '  BOUNDARY CONDITIONS:'
          WRITE (6,6003) 'GRD',IGRID,
     $       '    FRONT:',FRONT,'     BACK:',BACK
          WRITE (6,6003) 'GRD',IGRID,
     $       '    RIGHT:',RIGHT,'     LEFT:',LEFT
          WRITE (6,6003) 'GRD',IGRID,
     $       '   BOTTOM:',BOTTOM,'      TOP:',TOP
          WRITE (6,6003) 'GRD',IGRID,
     $       '     CUBE:',CUBE

      WRITE (6,*) 
      WRITE (6,*) 

C      IF (LPLEVEL) THEN
C          WRITE (6,6002) 'GRD',IGRID,
C     $       '  PLEVEL IS ON !!'
C      ENDIF

      WRITE (6,*) 
      WRITE (6,*) 
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
 6000 FORMAT(1X,A3,I3,A10,I3,A10,I3,A10,I3)
 6001 FORMAT(1X,A3,I3,A10,F8.3,A10,F8.3,A10,F8.3)
 6002 FORMAT(1X,A3,I3,A30,I8)
 6003 FORMAT(1X,A3,I3,A10,A10,A10,A10)


      RETURN
      END
