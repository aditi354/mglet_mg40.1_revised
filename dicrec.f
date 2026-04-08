










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

        SUBROUTINE DICREC (LREC)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                           C
C       SUBROUTINE FUER DAS EINLESEN DER STEUERDATEN        C
C       FUER RAUSSCHREIBEN VON GESCHWINDIGKEITSRECORDS
C       MANHART 22.2 1991                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   NVREC :    ANZAHL DER ZEITRECORDS
C   ITREC :    SPRUNG ZWISCHEN DEN ZEITSCHRITTEN
C   MAXREC:    ANSCHLAG, MAXIMALE ANZAHL DER RECORDS
C   LREC  :    LOGICAL, .TRUE. FALLS ZEITRECORDS RAUSGESCHRIEBEN
C              WERDEN SOLL, .FALSE.  FALLS KANAL 52 LEER
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

         PARAMETER   ( NPREC= 16500000 , NMREC= 1 )

      COMMON /CORECDEF/
     &                 NVREC,    ITREC,    MAXREC,
     &                 NGREC,
     &                 CIDREC,   IVEREC,   IVOREC,
     &                 KANREC,   KANGEO,
     &                 INXREC,   INYREC,   INZREC,
     &                 XUGREC,      XOGREC,
     &                 YUGREC,      YOGREC,
     &                 ZUGREC,      ZOGREC

       INTEGER
     +         NGREC(NMREC)
     +        ,KANREC(NMREC),KANGEO(NMREC)
     +        ,IVEREC(NMREC),IVOREC(NMREC)
     +        ,INXREC(NMREC),INYREC(NMREC),INZREC(NMREC)
     +        ,NTREC(NMREC)


       REAL
     +      XUGREC(NMREC),XOGREC(NMREC)
     +     ,YUGREC(NMREC),YOGREC(NMREC)
     +     ,ZUGREC(NMREC),ZOGREC(NMREC)
C
       CHARACTER (LEN=16) CIDREC(NMREC),CIDRE2
C
C
        LOGICAL LREC
C
C
        WRITE (6,*) 'DICREC:  NVREC,ITREC,MAXREC,LREC',
     &                        NVREC,ITREC,MAXREC,LREC
C
        IF(NVREC.GT.NMREC) THEN
                 NVREC = NMREC
                 WRITE(6,*)'ES KOENNEN LEIDER NUR',NMREC
                 WRITE(6,*)'TESTVOLUMEN BERUECKSICHTIGT WERDEN'
                 WRITE(6,*)'ODER NMREC IM HAUPTPROGRAMM ERWEITERN'
        ENDIF

         IF(LREC) THEN
C
C
        DO 20 I=1,NVREC
C
          WRITE (6,*) 'DICREC: CIDREC(I):',CIDREC(I)
          WRITE (6,*) 'DICREC: IVEL,IVOR:',IVEREC(I),IVOREC(I)
C
             IF (IVEREC(I)*IVOREC(I).GE.1) GOTO 998
C
C
          WRITE (6,*) 'DICREC: KANREC,KANGEO:',KANREC(I),KANGEO(I)
          WRITE (6,*) 'DICREC: INXREC,INYREC,INZREC:',
     &                        INXREC(I),INYREC(I),INZREC(I)
          WRITE (6,*) 'DICREC: XUG,XOG:',XUGREC(I),XOGREC(I)
          WRITE (6,*) 'DICREC: YUG,YOG:',YUGREC(I),YOGREC(I)
          WRITE (6,*) 'DICREC: ZUG,ZOG:',ZUGREC(I),ZOGREC(I)
C
        IF(INXREC(I)*INYREC(I)*INZREC(I).GT.NPREC) THEN
            WRITE(6,*)'DIMENSIONIERUNG NPREC IM HAUPTPROGRAMM ZU KLEIN'
            WRITE(6,*)'NOETIG WARE:',INXREC(I)*INYREC(I)*INZREC(I)
            WRITE(6,*)'DESHALB RAUSSCHREIBEN VON RECORDS ABGEBROCHEN'
            LREC = .FALSE.
        ENDIF
C
   20     CONTINUE
C
        RETURN
        ENDIF
C
C                       ENDE DER REGULAEREN EINLESEPHASE
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                       FALLS FEHLER WAEREND EINLESEN AUFGETRETEN,
C                       WERDEN BIS DAHIN FEHLERFREI EINGELESENE
C                       TESTVOLUMEN RAUSGESCHRIEBEN
C
  997   WRITE(6,*)'BEIM EINLESEN VON DATENFILE FEHLER AUFGETRETEN'
C
  998     NVREC = I-1
         IF (NVREC.LT.1) LREC= .FALSE.
        RETURN
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                       FALLS EINGABEFILE KANAL52 LEER ODER NICHT
C                       ANGELINKT WIRD AUCH KEIN TV RAUSGESCHRIEBEN
C                       LREC AUF .FALSE.
C
  999    LREC = .FALSE.
        RETURN
C
C
        END
