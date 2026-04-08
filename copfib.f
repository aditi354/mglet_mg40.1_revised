










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
      SUBROUTINE COPFIB  (KK,JJ,II,KMX,JMX,IMX,CIDEND,IIDEND,RIDEND,P1,
     $                    ITSTEP,ITTOT,ITINT,MTSTEP,DT)
C*STARLET***************************************************************
C        C O P F I B      DIE AUF KANAL 11 GELESENEN DATEN WERDEN
C                         AUF KANAL 12 BINAER GESCHRIEBEN (KOPIERT)
C*STARLET***************************************************************
C
C PARAM: KK,JJ,II       - ARRAYGRENZEN
C
C VERS:  06.05.86 (HW)  : ORIGINAL
C
C*STARLET***************************************************************
C
      CHARACTER (LEN=8)     CIDEND(10)
      INTEGER    IIDEND(100)
      REAL       RIDEND(100)
      REAL       P1(KK,JJ, 7)
C
C                                 ABPRUEFUNGEN
C
      IF(ITSTEP .NE. 0) CALL ERRR (501,' COPFIB   ')
      IF(ITTOT  .EQ. 0) CALL ERRR (502,' COPFIB   ')
C
      REWIND 11
      REWIND 12
C
C                                 GELESEN WIRD AUF KANAL 11
C
         READ (11,END=10000) CIDEND
10000    READ (11,END=10001) IIDEND
10001    READ (11,END=10002) RIDEND
10002    READ (11,END=10003) KKOLD,JJOLD,IIOLD,KMXOLD,JMXOLD,IMXOLD
10003    WRITE(12)           CIDEND
         WRITE(12)           IIDEND
C
C                                 ANPASSUNG DER ZYKLUSZEIT
C                                 ZU DER ZAHL DER EINTRITTSPROFILE
C                                 DES VORANGEGANGENEN LAUFES WIRD DIE
C                                 WAEHREND DES AKTUELLEN LAUFES ZU
C                                 ERWARTENDE ANZAHL VON EINTR.-PROF.
C                                 ADDIERT UND DAMIT DIE ZYKLUSZEIT
C                                 GEBILDET.
C
         RIDEND(31) = FLOAT(MIN0((ITTOT+MTSTEP),ITINT)) * DT
C
         WRITE(12)           RIDEND
C
         WRITE(12)           KKOLD,JJOLD,IIOLD,KMXOLD,JMXOLD,IMXOLD
C
      NZEB = 0
C
 2000 READ (11,END=10010) P1
      WRITE(12) (((P1(K,J,I),K=1,KMXOLD),J=1,JMXOLD),I=1,7)
      NZEB = NZEB + 1
      GOTO 2000
C
10010 WRITE (6,6010) ITTOT,NZEB
 6010 FORMAT(//,1X,15(1H*),'  INFORMATION VON COPFIB ',15(1H*),/,
     $       1X,15(1H*),'  ITTOT = ',I7,/,1X,15(1H*),
     $       '  ANZAHL DER GELESENEN ZEITEBENEN : ',I7)
      RETURN
      END
