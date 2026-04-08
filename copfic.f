










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
      SUBROUTINE COPFIC  (KK,JJ,II,KMX,JMX,IMX,CIDEND,IIDEND,RIDEND,P1,
     $                    ITSTEP,ITTOT,ITINT,MTSTEP,DT)
C*STARLET***************************************************************
C        C O P F I C      DIE AUF KANAL 13 GELESENEN DATEN WERDEN
C                         AUF KANAL 14 FORMATIERT GESCHRIEBEN
C                         (KOPIERT)
C*STARLET***************************************************************
C
C PARAM: KK,JJ,II       - ARRAYGRENZEN
C
C VERS:  08.01.87 (HW)  : ORIGINAL
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
      IF(ITSTEP .NE. 0) CALL ERRR (501,' COPFIC   ')
      IF(ITTOT  .EQ. 0) CALL ERRR (502,' COPFIC   ')
C
      REWIND 13
      REWIND 14
C
C                                 GELESEN WIRD AUF KANAL 13
C
         READ (13,6010) (CIDEND(N),N=1,10)
         READ (13,6020) (IIDEND(N),N=1,100)
         READ (13,6030) (RIDEND(N),N=1,100)
C
         READ (13,6040) KKOLD,JJOLD,IIOLD,KMXOLD,JMXOLD,IMXOLD
C
         WRITE(14,6010) (CIDEND(N),N=1,10)
         WRITE(14,6020) (IIDEND(N),N=1,100)
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
         WRITE(14,6030) (RIDEND(N),N=1,100)
C
         WRITE(14,6040) KKOLD,JJOLD,IIOLD,KMXOLD,JMXOLD,IMXOLD
C
      NZEB = 0
C
 2000 READ (13,6030,END=10010) (((P1(K,J,I),K=1,KMXOLD),J=1,JMXOLD),
     $                                      I=1,7)
      WRITE(14,6030)           (((P1(K,J,I),K=1,KMXOLD),J=1,JMXOLD),
     $                                      I=1,7)
      NZEB = NZEB + 1
      GOTO 2000
C
10010 WRITE (6,6000) ITTOT,NZEB
      RETURN
C
 6000 FORMAT (//,1X,15(1H*),'  INFORMATION VON COPFIC ',15(1H*),/,
     $        1X,15(1H*),'  ITTOT = ',I7,/,1X,15(1H*),
     $        '  ANZAHL DER GELESENEN ZEITEBENEN : ',I7)
 6010 FORMAT (5(A8,2X))
 6020 FORMAT (5(I9,1X))
 6030 FORMAT (4(E19.12E3,1X))
 6040 FORMAT (6(I9,1X))
      END
