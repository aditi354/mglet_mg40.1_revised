










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
      SUBROUTINE WRIGEO(KK,JJ,II,KMX,JMX,IMX,X,Y,Z,KANAL,ICOUNT)

      REAL X(II), Y(JJ), Z(KK)

      REWIND (KANAL)

      WRITE (KANAL,'(A)') 'LES        '

      WRITE (KANAL,'(A)') '  X  '
      WRITE (KANAL,6050) (X(I),I=1,IMX)
      WRITE (KANAL,'(A)') '  Y  '
      WRITE (KANAL,6050) (Y(I),I=1,JMX)
      WRITE (KANAL,'(A)') '  Z  '
      WRITE (KANAL,6050) (Z(I),I=1,KMX)
      WRITE (KANAL,'(A)') 'IARR '
      WRITE (KANAL,6060) (I,I=1,IMX)
      WRITE (KANAL,'(A)') 'JARR '
      WRITE (KANAL,6060) (I,I=1,JMX)
      WRITE (KANAL,'(A)') 'KARR '
      WRITE (KANAL,6060) (I,I=1,KMX)

      WRITE (KANAL,'(A)') ' ICOUNT'
      WRITE (KANAL,6060) ICOUNT

 6050  FORMAT (6(E12.5E3,1X))
 6060  FORMAT (4(I9,1X))

      RETURN
      END


