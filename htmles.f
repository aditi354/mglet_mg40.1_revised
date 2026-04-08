










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
      SUBROUTINE HTMLES(IGRDMX,IVOR,INACH,IGROB,IHTMG,EPSHTM,LRESID)

      INTEGER IGRDMX,IVOR,INACH,IGROB,IHTMG
      REAL EPSHTM
      LOGICAL LRESID

C     BEGINN DES EINLESEBLOCKS

C     READ(5,'(A40)') TITLE
C     READ(5,*) NDIM
C     OPEN(35,FILE='HTMG.DAT')
      OPEN(31,FILE='HTMG.OUT')
      OPEN(36,FILE='HTMGDIV.OUT')
      READ(35,*) IGRDMX
      READ(35,*) IVOR,INACH,IGROB,IHTMG
      READ(35,*) EPSHTM
      READ(35,*) LRESID

C    KONTROLLAUSGABE AUF FILE FORT.31

C     OPEN(11,FILE='OPOISOL')
      WRITE(31,*) 'VORGLAETTUNGEN : ',IVOR
      WRITE(31,*) 'NACHGLAETTUNGEN : ',INACH
      WRITE(31,*) 'GROBGITTERGLAETTUNGEN : ',IGROB
      WRITE(31,*) 'ZAHL DER HTMG-ZYKLEN: ',IHTMG
      WRITE(31,*) 'RESIDDUMSBERECHNUNG = ',LRESID
      WRITE(31,*)

      RETURN
      END


