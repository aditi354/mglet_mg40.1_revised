










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
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C    SUBROUTINE ZUR POSITIONIERUNG AUF DAS FILEENDE
C    FUER RAUSSCHREIBEN VON ZEITRECORDS
C    MM 22.2 1991
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C   NVREC :    ANZAHL DER ZEITRECORDS
C   ITREC :    SPRUNG ZWISCHEN DEN ZEITSCHRITTEN
C   MAXREC:    ANSCHLAG, MAXIMALE ANZAHL DER RECORDS
C   LREC  :    LOGICAL, .TRUE. FALLS ZEITRECORDS RAUSGESCHRIEBEN
C              WERDEN SOLL, .FALSE.  FALLS KANAL 52 LEER
C   NTREC :    ANZAHL DER BIS JETZT GESCHRIEBENEN ZEITRECORDS
C   NPREC:   GESAMTDIMENSIONIERUNG DES FELDES RTREC, AUF DEM DER 
C              ZEITRECORD ZWISCHENGEPUFFERT WIRD
C
C
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
        SUBROUTINE POSREC (NMREC,NVREC,ITREC,MAXREC,LREC,NTREC,
     +                    CIDREC,CIDRE2,
     +                    IVEL,IVOR,KANREC,KANGEO,
     +                    INXREC,INYREC,INZREC,
     +                    XREC,YREC,ZREC,
     +                    IARR,JARR,KARR,
     +                    NPREC,
     +                    KK,JJ,II,
     +                    TIMEPH,IREC)
C
        INTEGER
     +           KANREC(NMREC),KANGEO(NMREC)
     +          ,IVEL(NMREC),IVOR(NMREC)
     +          ,INXREC(NMREC),INYREC(NMREC),INZREC(NMREC)
     +          ,IARR(NMREC,II),JARR(NMREC,JJ),KARR(NMREC,KK)
     +          ,NTREC(NMREC)
C
C
C
        DIMENSION  
     +            XREC(NMREC,II)
     +            ,YREC(NMREC,JJ)
     +            ,ZREC(NMREC,KK)
C
       CHARACTER (LEN=16) CIDREC(NMREC),CIDRE2
C
        LOGICAL LREC
C
C
       OPEN (KANREC(IREC),FORM='UNFORMATTED')
C

C
C                      UEBERSPRINGEN DER GEOMETRIEINFORMATION
C       WRITE(44,*)'AUFRUF VON SKIPGEO'
C       WRITE(44,*)'KANAL:',KANGEO(IREC)
C       WRITE(44,*)'TESTVOLUMEN:',IREC
C
C       CALL SKIPGEO (NMREC,NVREC,ITREC,MAXREC,LREC,NTREC(IREC),
C    +                    IVEL(IREC),IVOR(IREC),KANGEO(IREC),
C    +                    INXREC(IREC),INYREC(IREC),INZREC(IREC),
C    +                    XREC,YREC,ZREC,
C    +                    IARR,JARR,KARR,
C    +                    KK,JJ,II,
C    +                    KMX,JMX,IMX)
C
C       IF(.NOT.LREC) RETURN
C

C
C                      UEBERSPRINGEN SCHON GESCHRIEBENER RECORDS
C
        CALL SKIPREC (NMREC,NVREC,ITREC,MAXREC,LREC,NTREC(IREC),
     +                    CIDREC(IREC),CIDRE2,
     +                    IVEL(IREC),IVOR(IREC),KANREC(IREC),
     +                    INXREC(IREC),INYREC(IREC),INZREC(IREC),
     +                    TIMEPH)
C

C           FEHLERABFRAGE
C
        IF(.NOT.LREC) THEN
        WRITE(44,*)'FEHLER BEIM POSITIONIEREN VON TESTVOLUMEN',IREC
        RETURN
        ENDIF
C
C

C
      WRITE(44,*)'TESTVOLUMEN',IREC,'WURDE POSITIONIERT'
C

C
C
        RETURN
C
        END
C
C
C
C

C
C
C
C
        SUBROUTINE SKIPREC (NMREC,NVREC,ITREC,MAXREC,LREC,NTREC,
     +                    CIDREC,CIDRE2,
     +                    IVEL,IVOR,KANREC,
     +                    INXREC,INYREC,INZREC,
     +                    TIMEPH)
C
        INTEGER
     +           KANREC,IVEL,IVOR
     +          ,INXREC,INYREC,INZREC
     +          ,NTREC
C
C        DIMENSION  RTREC(INZREC,INYREC,INXREC)
C
       CHARACTER (LEN=16) CIDREC,CIDRE2
C
        LOGICAL LREC
C

C                      UEBERPRUEFEN DER IDENTITAET DES RECORDS
       REWIND(KANREC)
       READ(KANREC,END=997) CIDRE2
C

C
       IF (CIDRE2.NE.CIDREC) THEN
       WRITE(44,*)'FALSCHE IDENTITAET VON GESCHWINDIGKEITSRECORD'
       LREC=.FALSE.
       RETURN
       ENDIF
C

C

C
          DO 100  I=1,NTREC
C

C                    ZEITMASSTAB
C
       READ (KANREC,END=998) RECTIME
       WRITE  (44,*) 'TIMEPH',RECTIME
C

C                      DREI GESCHWINDIGKEITSKOMPONENTEN
       READ (KANREC,END=999) 
       READ (KANREC,END=999) 
       READ (KANREC,END=999) 
C

C         FALLS IVEL ODER IVOR=3 WIRD DRUCK MIT RAUSGESCHRIEBEN
C
       IF(IVEL.GE.3.OR.IVOR.GE.3) 
     +READ (KANREC,END=999) 
C

C
  100    CONTINUE
C

C
       WRITE(44,*)'ES KONNTEN',NTREC,'ZEITRECORDS EINGELSESEN WERDEN'
C

C
       RETURN
C

C
  997   REWIND (KANREC)
       WRITE(KANREC) CIDREC
       RETURN
C

C
  998   IF (NTREC.NE.I-1)
     + WRITE(44,*)'WARNUNG ANZAHL DER VORHANDENEN ',
     +'ZEITRECORDS UNVOLLSTAENDIG'
       NTREC=I-1
       RETURN
C

C
  999   WRITE(44,*)'BEI POSITIONIERUNG DER ZEITRECORDS FEHLER ',
     + 'AUFGETRETEN DAHER WERDEN KEINE RECORDS RAUSGESCHRIEBEN'
       LREC = .FALSE.
       RETURN
C

 6050   FORMAT (6(E12.5E3,1X))
 6060   FORMAT (4(I9,1X))
C
        END

        SUBROUTINE SKIPGEO (NMREC,NVREC,ITREC,MAXREC,LREC,NTREC,
     +                    CIDREC,
     +                    IVEL,IVOR,KANGEO,
     +                    INXREC,INYREC,INZREC,
     +                    XREC,YREC,ZREC,
     +                    IARR,JARR,KARR,
     +                    KK,JJ,II,
     +                    KMX,JMX,IMX)
C
        INTEGER
     +           KANGEO,IVEL,IVOR
     +          ,INXREC,INYREC,INZREC
     +          ,IARR(II),JARR(JJ),KARR(KK)
     +          ,NTREC
C
C
C
        DIMENSION  XREC(II)
     +            ,YREC(JJ)
     +            ,ZREC(KK)
C
       CHARACTER (LEN=16) CIDREC(NMREC)
C
        LOGICAL LREC
C

C              SCHLEIFE UEBER TESTVOLUMEN
C
C
C
        WRITE(44,*)'IN SKIPGEO'
        WRITE(44,*)'TESTVOLUMEN:        KANAL:',KANGEO
C

C
C                      MODUL ZUM EINLESEN SCHON ERZEUGTER GEOMETRIEN
C

C
       REWIND KANGEO
C

C

C

C
       READ (KANGEO,6050,END=999)
       READ (KANGEO,6050,END=999)
     +                                 (XREC(I), I=1,INXREC)
       READ (KANGEO,6050,END=999)
       READ (KANGEO,6050,END=999)
     +                                (YREC(I), I=1,INYREC)
       READ (KANGEO,6050,END=999)
       READ (KANGEO,6050,END=999)
     +                                (ZREC(I), I=1,INZREC)
C
       READ (KANGEO,6060,END=999)
       READ (KANGEO,6060,END=999)
     +                                (IARR(I), I=1,INXREC)
       READ (KANGEO,6060,END=999)
       READ (KANGEO,6060,END=999)
     +                                (JARR(I), I=1,INYREC)
       READ (KANGEO,6060,END=999)
       READ (KANGEO,6060,END=999)
     +                                (KARR(I), I=1,INZREC)
C
       READ (KANGEO,6060,END=999)
       READ (KANGEO,6060,END=999)
     +                                  NTREC
C

C               ABFRAGE, OB ANSCHLAG ERREICHT
       IF(NTREC.GE.MAXREC) THEN
       LREC = .FALSE.
       RETURN
       ENDIF
C
C
       RETURN

C
  999   WRITE(44,*)
     + 'BEI UEBERSPRINGEN  DER GEOMETRIEINFORMATION FEHLER ',
     + 'AUFGETRETEN DAHER WERDEN KEINE RECORDS RAUSGESCHRIEBEN'
       LREC = .FALSE.
       RETURN
C

 6050   FORMAT (6(E12.5E3,1X))
 6060   FORMAT (4(I9,1X))
C
        END

