










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
      SUBROUTINE ERRR(NRRR,MRRR)
C***********************************************************************
C  E R R R  FEHLERBEHANDLUNG
C*************************************************** F.BAETKE 05.10.84 *
C
C     NRRR      - FEHLERPARAMETER = 100-499 MELDUNG
C                                 = 500-999 FEHLER ( MIT STOP )
C     MRRR      - MODULNAME ( MAXIMAL 10 ZEICHEN ALS STRING )
C
C***********************************************************************
C
      CHARACTER (LEN=10) MRRR
C                                 TEST OB MELDUNG ODER FEHLER
      IF(NRRR.GT.499)GOTO 2010
C                                 AUSGABE EINER MELDUNG MIT
C                                 PARAMETER UND MODULNAMEN
      WRITE(6,6000)NRRR,MRRR
 6000 FORMAT(/,1X,'  ++++++++++ MELDUNG  ',I3,'  IN MODUL',A8,/)
      GOTO 2900
C                                 AUSGABE EINES FEHLERS MIT
C                                 PARAMETER UND MODULNAMEN
 2010 CONTINUE
      GOTO 2900
 2900 CONTINUE
C                                 TEST, OB PROGRAMMSTOP ODER
C                                 RUECKSPRUNG IN RUFENDE MODUL
       WRITE(6,*)' NRRR= ',NRRR,' MRRR= ',MRRR
      IF(NRRR .GT. 499) THEN
        STOP ' IN ERRR'
      ENDIF
      RETURN
      END
