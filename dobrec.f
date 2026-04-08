










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
        SUBROUTINE DOBREC (NMREC,NVREC,ITREC,MAXREC,LREC,NTREC,
     +                    CIDREC,CIDRE2,
     +                    IVEL,IVOR,KANREC,KANGEO,
     +                    INXREC,INYREC,INZREC,
     +                    IARR,JARR,KARR,
     +                    XREC,YREC,ZREC,
     +                    NPREC,
     +                    RTREC,TIMEPH,
     +                    U,V,W,P,
     +                    UFG,VFG,WFG,PFG,
     +                    AU,AV,AW,AP,
     +                    OX,OY,OZ,
     +                    OXFG,OYFG,OZFG,
     +                    H1,H2,
     +                    KK,JJ,II,
     +                    DX,DY,DZ,DDX,DDY,DDZ,IREC)
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C
C    SUBROUTINE ZUM RAUSSCHREIBEN DER ZEITRECORDS
C    FUER RAUSSCHREIBEN VON ZEITRECORDS
C    MM 25.2 1991
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C   NVREC :    ANZAHL DER ZEITRECORDS
C   ITREC :    SPRUNG ZWISCHEN DEN ZEITSCHRITTEN
C   MAXREC:    ANSCHLAG, MAXIMALE ANZAHL DER RECORDS
C   LREC  :    LOGICAL, .TRUE. FALLS ZEITRECORDS RAUSGESCHRIEBEN
C              WERDEN SOLL, .FALSE.  FALLS KANAL 52 LEER
C   NTREC :    ANZAHL DER BIS JETZT GESCHRIEBENEN ZEITRECORDS
C   NPREC:   GESAMTDIMENSIONIERUNG DES FELDES RTREC, AUF DEM DER 
C              ZEITRECORD ZWISCHENGEPUFFERT WIRD
C   RTREC  :    PUFFERFELD, AUF DAS GESCHWINDIGK. ODER OMEGA USW.
C              ABGESPEICHERT WERDEN
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C             CONVEX-SPEZIFISCHE FUNKTION "FLUSH" !!!!!!!!!!!
C
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
        INTEGER
     +           KANREC(NMREC),KANGEO(NMREC)
     +          ,IVEL(NMREC),IVOR(NMREC)
     +          ,INXREC(NMREC),INYREC(NMREC),INZREC(NMREC)
     +          ,IARR(NMREC,II),JARR(NMREC,JJ),KARR(NMREC,KK)
     +          ,NTREC(NMREC)

        DIMENSION  XREC(NMREC,II)
     +            ,YREC(NMREC,JJ)
     +            ,ZREC(NMREC,KK)
C
C
        DIMENSION  RTREC(NPREC)
C
      REAL    U(KK,JJ,II),V(KK,JJ,II),W(KK,JJ,II),
     $        AU(KK,JJ,II),AV(KK,JJ,II),AW(KK,JJ,II),
     $        P(KK,JJ,II),AP(KK,JJ,II),PFG(KK,JJ,II),
     $        UFG(KK,JJ,II),VFG(KK,JJ,II),WFG(KK,JJ,II),
     $        OX(KK,JJ,II),OY(KK,JJ,II),OZ(KK,JJ,II),
     $        OXFG(KK,JJ,II),OYFG(KK,JJ,II),OZFG(KK,JJ,II),
     $        DX(II),DY(JJ),DZ(KK),DDX(II),DDY(JJ),DDZ(KK)

      REAL H1(KK,JJ,II),H2(KK,JJ,II)
C

       CHARACTER (LEN=16) CIDREC(NMREC),CIDRE2
C
        LOGICAL LREC
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
       KMX = KK
       JMX = JJ
       IMX = II
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                PHYSIKALISCHE ZEIT
C
       RECTIME = TIMEPH
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                    ABPRUEFEN, OB ANSCLAG SCHON ERREICHT
      IF (NTREC(IREC).GE.MAXREC) GOTO 1000
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                    RAUSCHREIBEN DES ZEITPUNKTES
C
        WRITE (KANREC(IREC)) RECTIME
        WRITE (44,*) 'RAUSSCHREIBEN VON ZEITRECORD BEI TIMEPH=',RECTIME
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C     RAUSSCHREIBEN DER MOMENTANEN GESCHWINDIGKEITEN
C     (AUF DRUCKPUNKT INTERPOLIERT)
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       IF((IVEL(IREC).EQ.1).OR.(IVEL(IREC).EQ.3)
     +     .OR.(IVEL(IREC).EQ.19)) THEN
C
C
        CALL INTVEL (KANREC(IREC),NMREC,IREC,
     +               CIDREC,IVEL(IREC),
     +               INXREC(IREC),INYREC(IREC),INZREC(IREC),
     +               IARR,JARR,KARR,
     +               XREC,YREC,ZREC,
     +               RTREC,
     +               U,V,W,P,
     +               KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ,TIMEPH)
C
C
C
      ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       IF((IVEL(IREC).EQ.2).OR.(IVEL(IREC).EQ.5)) THEN
C
C
        CALL INTVEL (KANREC(IREC),NMREC,IREC,
     +               CIDREC,IVEL(IREC),
     +               INXREC(IREC),INYREC(IREC),INZREC(IREC),
     +               IARR,JARR,KARR,
     +               XREC,YREC,ZREC,
     +               RTREC,
     +               UFG,VFG,WFG,PFG,
     +               KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ,TIMEPH)
C
C
C
      ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       IF(IVEL(IREC).EQ.4) THEN
C
C
        CALL INTVEL (KANREC(IREC),NMREC,IREC,
     +               CIDREC,IVEL(IREC),
     +               INXREC(IREC),INYREC(IREC),INZREC(IREC),
     +               IARR,JARR,KARR,
     +               XREC,YREC,ZREC,
     +               RTREC,
     +               AU,AV,AW,AP,
     +                    KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ,TIMEPH)
C
C
C
      ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       IF(IVEL(IREC).EQ.11) THEN
C
C
        CALL FLTVEL (KANREC(IREC),NMREC,IREC,
     +               CIDREC,IVEL(IREC),
     +               INXREC(IREC),INYREC(IREC),INZREC(IREC),
     +               IARR,JARR,KARR,
     +               RTREC,
     +               U,V,W,P,
     +               KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ)
C
C
C
      ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       IF(IVEL(IREC).EQ.12) THEN
C
C
        CALL FLTVEL (KANREC(IREC),NMREC,IREC,
     +               CIDREC,IVEL(IREC),
     +               INXREC(IREC),INYREC(IREC),INZREC(IREC),
     +               IARR,JARR,KARR,
     +               RTREC,
     +               UFG,VFG,WFG,PFG,
     +               KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ)
C
C
C
      ENDIF
C

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       IF(IVEL(IREC).EQ.14) THEN
C
C
        CALL FLTVEL (KANREC(IREC),NMREC,IREC,
     +               CIDREC,IVEL(IREC),
     +               INXREC(IREC),INYREC(IREC),INZREC(IREC),
     +               IARR,JARR,KARR,
     +               RTREC,
     +               AU,AV,AW,AP,
     +               KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ)
C
C
C
      ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       IF(IVEL(IREC).EQ.21) THEN
C
C
        CALL FLTVEL2 (KANREC(IREC),NMREC,IREC,
     +               CIDREC,IVEL(IREC),
     +               INXREC(IREC),INYREC(IREC),INZREC(IREC),
     +               IARR,JARR,KARR,
     +               RTREC,
     +               U,V,W,P,
     +               KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ)
C
C
C
      ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       IF(IVEL(IREC).EQ.22) THEN
C
C
        CALL FLTVEL2 (KANREC(IREC),NMREC,IREC,
     +               CIDREC,IVEL(IREC),
     +               INXREC(IREC),INYREC(IREC),INZREC(IREC),
     +               IARR,JARR,KARR,
     +               RTREC,
     +               UFG,VFG,WFG,PFG,
     +               KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ)
C
C
C
      ENDIF
C

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       IF(IVEL(IREC).EQ.24) THEN
C
C
        CALL FLTVEL2 (KANREC(IREC),NMREC,IREC,
     +               CIDREC,IVEL(IREC),
     +               INXREC(IREC),INYREC(IREC),INZREC(IREC),
     +               IARR,JARR,KARR,
     +               RTREC,
     +               AU,AV,AW,AP,
     +               KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ)
C
C
C
      ENDIF
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       IF(IVEL(IREC).EQ.31) THEN
C
C
        CALL ACOUSTIC2 (KANREC(IREC),NMREC,IREC,
     +               CIDREC,IVEL(IREC),
     +               INXREC(IREC),INYREC(IREC),INZREC(IREC),
     +               IARR,JARR,KARR,
     +               RTREC,
     +               U,V,W,P,H1,H2,
     +               KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ)
C
C
C
      ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       IF(IVEL(IREC).EQ.7) THEN
C
        IVEL(IREC)=1
        CALL INTVEL (KANREC(IREC),NMREC,IREC,
     +               CIDREC,IVEL(IREC),
     +               INXREC(IREC),INYREC(IREC),INZREC(IREC),
     +               IARR,JARR,KARR,
     +               XREC,YREC,ZREC,
     +               RTREC,
     +               U,V,W,P,
     +               KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ,TIMEPH)
        IVEL(IREC)=7
        CALL Q_KRITERIUM (KANREC(IREC),NMREC,IREC,
     +               CIDREC,IVEL(IREC),
     +               INXREC(IREC),INYREC(IREC),INZREC(IREC),
     +               IARR,JARR,KARR,
     +               RTREC,
     +               U,V,W,P,H1,H2,
     +               KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ)
C
C
C
      ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       IF(IVOR(IREC).EQ.1) THEN
C
C
        CALL INTVOR (KANREC(IREC),NMREC,IREC,
     +                    CIDREC,
     +               INXREC(IREC),INYREC(IREC),INZREC(IREC),
     +               IARR,JARR,KARR,
     +               RTREC,
     +               OX,OY,OZ,
     +                    KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ)
C
C
C
      ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       IF(IVOR(IREC).EQ.2) THEN
C
C
        CALL INTVOR (KANREC(IREC),NMREC,IREC,
     +                    CIDREC,
     +               INXREC(IREC),INYREC(IREC),INZREC(IREC),
     +               IARR,JARR,KARR,
     +               RTREC,
     +               OXFG,OYFG,OZFG,
     +                    KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ)
C
C
C
      ENDIF

      IF(IVEL(IREC).EQ.41) THEN

      ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C             EINTRAGEN DER NEUEN ANZAHL DER ZEITRECORDS AM ENDE
C
        NTREC(IREC) = NTREC(IREC)+1
C
        BACKSPACE (KANGEO(IREC))
        WRITE (KANGEO(IREC),6060) NTREC(IREC)
C
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C
 1000    CONTINUE
C
C
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
        RETURN
C
 6050   FORMAT (6(E12.5E3,1X))
 6060   FORMAT (4(I9,1X))
C
        END
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
        SUBROUTINE INTVEL (KANREC,NMREC,IREC,
     +               CIDREC,IVEL,
     +               INXREC,INYREC,INZREC,
     +               IARR,JARR,KARR,
     +               XREC,YREC,ZREC,
     +               RTREC,
     +               U,V,W,P,
     +                    KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ,TIMEPH)
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
        INTEGER
     +           KANREC,IVEL
     +          ,INXREC,INYREC,INZREC
     +          ,IARR(NMREC,II),JARR(NMREC,JJ),KARR(NMREC,KK)
C
C
        DIMENSION  RTREC(INZREC,INYREC,INXREC)
        DIMENSION  RRTREC(INXREC,INYREC,INZREC)
C
        DIMENSION  XREC(NMREC,II),XXREC(II)
     +            ,YREC(NMREC,JJ),YYREC(JJ)
     +            ,ZREC(NMREC,KK),ZZREC(KK)
C
      REAL    U(KK,JJ,II),V(KK,JJ,II),W(KK,JJ,II),
     $        P(KK,JJ,II),
     $        DX(II),DY(JJ),DZ(KK),DDX(II),DDY(JJ),DDZ(KK)
C
       CHARACTER (LEN=16) CIDREC(NMREC)
       CHARACTER (LEN=32) FNAME
C
        LOGICAL LREC
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
        WRITE(44,*)'RAUSSCHREIBEN EINES RECORDS AUF KANAL',KANREC
C       WRITE(44,*)'DX'
C       WRITE(44,*)(DX(I),I=1,IMX)
C       WRITE(44,*)'DY'
C       WRITE(44,*)(DY(I),I=1,JMX)
C       WRITE(44,*)'DZ'
C       WRITE(44,*)(DZ(I),I=1,KMX)
C       WRITE(44,*)'DDX'
C       WRITE(44,*)(DDX(I),I=1,IMX)
C       WRITE(44,*)'DDY'
C       WRITE(44,*)(DDY(I),I=1,JMX)
C       WRITE(44,*)'DDZ'
C       WRITE(44,*)(DDZ(I),I=1,KMX)
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
        DO I=1,INXREC
          XXREC(I) = XREC(IREC,I)
        ENDDO

        DO J=1,INYREC
	  YYREC(J) = YREC(IREC,J)
        ENDDO

        DO K=1,INZREC
	  ZZREC(K) = ZREC(IREC,K)
        ENDDO
C
      DO 10 I=1,INXREC
      DO 10 J=1,INYREC
      DO 10 K=1,INZREC
C
   10 RTREC(K,J,I) =
     +      (   U   (KARR(IREC,K),JARR(IREC,J),IARR(IREC,I)  )
     +                                 *0.5*DX(IARR(IREC,I)-1)
     +                                    /DDX(IARR(IREC,I)  ) )
     +   +  (   U   (KARR(IREC,K),JARR(IREC,J),IARR(IREC,I)-1)
     +                             *(1.-0.5*DX(IARR(IREC,I)-1)
     +                                    /DDX(IARR(IREC,I)  )))
     +    
C
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C        WRITE (KANREC) RTREC
C

	kanal = KANREC+1
C
	write(FNAME,'("u_vel3_",I3.3,"_",I8.8,F9.8,"_BIN")')
     $		KANREC,int(TIMEPH),TIMEPH-int(TIMEPH)

        call  W6D
     $ (KANAL,0,1,'BIN','STREAM','U','SCALAR',
     $  'ZZREC','YYREC','XXREC','TIMEPH','TIMEPH','TIMEPH',
     $  'RRTREC','RRTREC','RRTREC',
     $  3,INZREC,INYREC,INXREC,1,1,1,
     $  ZZREC,YYREC,XXREC,TIMEPH,TIMEPH,TIMEPH,
     $  RTREC,RTREC,RTREC,FNAME)
C       WRITE (44,*)'U RAUSGESCHRIEBEN'
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      DO 20 I=1,INXREC
      DO 20 J=1,INYREC
      DO 20 K=1,INZREC
C
   20 RTREC(K,J,I) =
     +      (   V   (KARR(IREC,K),JARR(IREC,J)  ,IARR(IREC,I))
     +                    *0.5*DY(JARR(IREC,J)-1)
     +                       /DDY(JARR(IREC,J)  )              )
     +   +  (   V   (KARR(IREC,K),JARR(IREC,J)-1,IARR(IREC,I))
     +                *(1.-0.5*DY(JARR(IREC,J)-1)
     +                       /DDY(JARR(IREC,J)  )             ))
     +    
C
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C        WRITE (KANREC) RTREC
C       WRITE (44,*)'V RAUSGESCHRIEBEN'
C
	kanal = KANREC + 1 
C
	write(FNAME,'("v_vel3_",I3.3,"_",I8.8,F9.8,"_BIN")')
     $		KANREC,int(TIMEPH),TIMEPH-int(TIMEPH)

        call  W6D
     $ (KANAL,0,1,'BIN','STREAM','V','SCALAR',
     $  'ZZREC','YYREC','XXREC','TIMEPH','TIMEPH','TIMEPH',
     $  'RRTREC','RRTREC','RRTREC',
     $  3,INZREC,INYREC,INXREC,1,1,1,
     $  ZZREC,YYREC,XXREC,TIMEPH,TIMEPH,TIMEPH,
     $  RTREC,RTREC,RTREC,FNAME)
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      DO 30 I=1,INXREC
      DO 30 J=1,INYREC
      DO 30 K=1,INZREC
C
   30 RTREC(K,J,I) =
     +      (   W   (KARR(IREC,K)  ,JARR(IREC,J),IARR(IREC,I))
     +       *0.5*DZ(KARR(IREC,K)-1)
     +          /DDZ(KARR(IREC,K)  )                           )
     +   +  (   W   (KARR(IREC,K)-1,JARR(IREC,J),IARR(IREC,I))
     +   *(1.-0.5*DZ(KARR(IREC,K)-1)
     +          /DDZ(KARR(IREC,K)  )                          ))
C
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C        WRITE (KANREC) RTREC
C       WRITE (44,*)'W RAUSGESCHRIEBEN'
C
	kanal = KANREC + 1
C
	write(FNAME,'("w_vel3_",I3.3,"_",I8.8,F9.8,"_BIN")')
     $		KANREC,int(TIMEPH),TIMEPH-int(TIMEPH)

        call  W6D
     $ (KANAL,0,1,'BIN','STREAM','W','SCALAR',
     $  'ZZREC','YYREC','XXREC','TIMEPH','TIMEPH','TIMEPH',
     $  'RRTREC','RRTREC','RRTREC',
     $  3,INZREC,INYREC,INXREC,1,1,1,
     $  ZZREC,YYREC,XXREC,TIMEPH,TIMEPH,TIMEPH,
     $  RTREC,RTREC,RTREC,FNAME)
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       IF(IVEL.GE.3) THEN
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      DO 40 I=1,INXREC
      DO 40 J=1,INYREC
      DO 40 K=1,INZREC
C
   40 RTREC(K,J,I) =
     +          P   (KARR(IREC,K)  ,JARR(IREC,J),IARR(IREC,I))
C
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C        WRITE (KANREC) RTREC
C
	kanal = KANREC + 1
C
	write(FNAME,'("p_vel3_",I3.3,"_",I8.8,F9.8,"_BIN")')
     $		KANREC,int(TIMEPH),TIMEPH-int(TIMEPH)

        call  W6D
     $ (KANAL,0,1,'BIN','STREAM','P','SCALAR',
     $  'ZZREC','YYREC','XXREC','TIMEPH','TIMEPH','TIMEPH',
     $  'RRTREC','RRTREC','RRTREC',
     $  3,INZREC,INYREC,INXREC,1,1,1,
     $  ZZREC,YYREC,XXREC,TIMEPH,TIMEPH,TIMEPH,
     $  RTREC,RTREC,RTREC,FNAME)
C
C       WRITE (44,*)'P RAUSGESCHRIEBEN'
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
        ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      DO  ILAUF=1,INXREC
      DO  JLAUF=1,INYREC
      DO  KLAUF=1,INZREC

          K = KARR(IREC,KLAUF)
          J = JARR(IREC,JLAUF)
          I = IARR(IREC,ILAUF)

	  IF (I.EQ.1) THEN
	     D2X =  
     $           ( 2.0*P( K,J,1) 
     $            -5.0*P( K,J,2)
     $            +4.0*P( K,J,3)
     $            -1.0*P( K,J,4) )
     $          /(DX(1)**2)
	  ENDIF

          IF (I.EQ.II) THEN
	     D2X =  
     $           ( 2.0*P(K,J,II  ) 
     $            -5.0*P(K,J,II-1)
     $            +4.0*P(K,J,II-2)
     $            -1.0*P(K,J,II-3) )
     $          /(DX(II-1)**2)
	  ENDIF

	  IF ((I.GE.2).AND.(I.LE.II-1)) THEN
	      D2X = 
     $            (P(K,J,I+1)+P(K,J,I-1)-2.0*P(K,J,I))/(DX(I)*DX(I-1))
          ENDIF

CCCCCCCCCCCCCCCCCCCCCCC
	  
	  IF (J.EQ.1) THEN
	     D2Y =  
     $           ( 2.0*P( K,1,I) 
     $            -5.0*P( K,2,I)
     $            +4.0*P( K,3,I)
     $            -1.0*P( K,4,I) )
     $          /(DY(1)**2)
	  ENDIF

          IF (J.EQ.JJ) THEN
	     D2Y =  
     $           ( 2.0*P(K,J ,I) 
     $            -5.0*P(K,J-1,I)
     $            +4.0*P(K,J-2,I)
     $            -1.0*P(K,J-3,I) )
     $          /(DY(JJ-1)**2)
	  ENDIF

	  IF ((J.GE.2).AND.(J.LE.JJ-1)) THEN
	      D2Y = 
     $            (P(K,J+1,I)+P(K,J-1,I)-2.0*P(K,J,I))/(DY(J)*DY(J-1))
          ENDIF

CCCCCCCCCCCCCCCCCCCCCCC
	  
	  IF (K.EQ.1) THEN
	     D2Z =  
     $           ( 2.0*P( 1,J,I) 
     $            -5.0*P( 2,J,I)
     $            +4.0*P( 3,J,I)
     $            -1.0*P( 4,J,I) )
     $          /(DZ(1)**2)
	  ENDIF

          IF (K.EQ.KK) THEN
	     D2Z =  
     $           ( 2.0*P(K  ,J,I) 
     $            -5.0*P(K-1,J,I)
     $            +4.0*P(K-2,J,I)
     $            -1.0*P(K-3,J,I) )
     $          /(DZ(KK-1)**2)
	  ENDIF

	  IF ((K.GE.2).AND.(K.LE.KK-1)) THEN
	     D2Z = 
     $            (P(K+1,J,I)+P(K-1,J,I)-2.0*P(K,J,I))/(DZ(K)*DZ(K-1))
          ENDIF
	  
CCCCCCCCCCCCCCCCCCCCCCC
C On veut MOINS le lapacien !!!

          RTREC(KLAUF,JLAUF,ILAUF) =  -(D2X +  D2Y +  D2Z)

CCCCCCCCCCCCCCCCCCCCCCC
 
      ENDDO
      ENDDO
      ENDDO

	kanal = KANREC + 1
C
	write(FNAME,'("dd_Tij_",I3.3,"_",I8.8,F9.8,"_BIN")')
     $		KANREC,int(TIMEPH),TIMEPH-int(TIMEPH)

        call  W6D
     $ (KANAL,0,1,'BIN','STREAM','dd_Tij','SCALAR',
     $  'ZZREC','YYREC','XXREC','TIMEPH','TIMEPH','TIMEPH',
     $  'RRTREC','RRTREC','RRTREC',
     $  3,INZREC,INYREC,INXREC,1,1,1,
     $  ZZREC,YYREC,XXREC,TIMEPH,TIMEPH,TIMEPH,
     $  RTREC,RTREC,RTREC,FNAME)

       RETURN
       END
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
        SUBROUTINE INTVOR (KANREC,NMREC,IREC,
     +                    CIDREC,
     +               INXREC,INYREC,INZREC,
     +               IARR,JARR,KARR,
     +               RTREC,
     +               OX,OY,OZ,
     +                    KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ)
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
        INTEGER
     +           KANREC
     +          ,INXREC,INYREC,INZREC
     +          ,IARR(NMREC,II),JARR(NMREC,JJ),KARR(NMREC,KK)
C
C
        DIMENSION  RTREC(INZREC,INYREC,INXREC)
C
      REAL    OX(KK,JJ,II),OY(KK,JJ,II),OZ(KK,JJ,II),
     $        DX(II),DY(JJ),DZ(KK),DDX(II),DDY(JJ),DDZ(KK)
C
       CHARACTER (LEN=16) CIDREC(NMREC)
C
        LOGICAL LREC
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
        WRITE(44,*)'RAUSSCHREIBEN EINES RECORDS AUF KANAL',KANREC
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      DO 10 I=1,INXREC
      DO 10 J=1,INYREC
      DO 10 K=1,INZREC
C
   10 RTREC(K,J,I) =
     +     (   OX   (KARR(IREC,K),JARR(IREC,J),IARR(IREC,I)  )
     +   +     OX   (KARR(IREC,K),JARR(IREC,J),IARR(IREC,I)+1))*0.5
     +    
C
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
        WRITE (KANREC) RTREC
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      DO 20 I=1,INXREC
      DO 20 J=1,INYREC
      DO 20 K=1,INZREC
C
   20 RTREC(K,J,I) =
     +     (   OY   (KARR(IREC,K),JARR(IREC,J)  ,IARR(IREC,I))
     +   +     OY   (KARR(IREC,K),JARR(IREC,J)+1,IARR(IREC,I)))*0.5
     +    
C
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
        WRITE (KANREC) RTREC
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      DO 30 I=1,INXREC
      DO 30 J=1,INYREC
      DO 30 K=1,INZREC
C
   30 RTREC(K,J,I) =
     +     (   OZ   (KARR(IREC,K)  ,JARR(IREC,J),IARR(IREC,I))
     +   +     OZ   (KARR(IREC,K)+1,JARR(IREC,J),IARR(IREC,I)))*0.5
C
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
        WRITE (KANREC) RTREC
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       RETURN
       END
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
        SUBROUTINE FLTVEL (KANREC,NMREC,IREC,
     +               CIDREC,IVEL,
     +               INXREC,INYREC,INZREC,
     +               IARR,JARR,KARR,
     +               RTREC,
     +               U,V,W,P,
     +               KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ)
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
        INTEGER
     +           KANREC,IVEL
     +          ,INXREC,INYREC,INZREC
     +          ,IARR(NMREC,II),JARR(NMREC,JJ),KARR(NMREC,KK)
C
C
        DIMENSION  RTREC(INZREC,INYREC,INXREC)
C
      REAL    U(KK,JJ,II),V(KK,JJ,II),W(KK,JJ,II),
     $        P(KK,JJ,II),
     $        DX(II),DY(JJ),DZ(KK),DDX(II),DDY(JJ),DDZ(KK)
C
       CHARACTER (LEN=16) CIDREC(NMREC)
C
        LOGICAL LREC
C
C                     U U U U U U U U U U U U U U U U U U U U U U U U
C
C
      DO  ILAUF=1,INXREC
      DO  JLAUF=1,INYREC
      DO  KLAUF=1,INZREC

          K = KARR(IREC,KLAUF)
          J = JARR(IREC,JLAUF)
          I = IARR(IREC,ILAUF)

        VOL = (DX(I-1)+DX(I)) * (DY(J-1)+DY(J)) * (DZ(K-1)+DZ(K))

      RTREC(KLAUF,JLAUF,ILAUF) =

C            K - SCHEIBE 

     +    (
     +     (DDZ(K)*
     +     (DDY(J  )*(DX(I)*U(K  ,J  ,I  ) + DX(I-1)*U(K  ,J  ,I-1))
     +  +0.5*DY(J-1)*(DX(I)*U(K  ,J-1,I  ) + DX(I-1)*U(K  ,J-1,I-1))
     +  +0.5*DY(J  )*(DX(I)*U(K  ,J+1,I  ) + DX(I-1)*U(K  ,J+1,I-1))))

C            (K-1) - SCHEIBE 

     +  +(0.5*DZ(K-1)*
     +     (DDY(J  )*(DX(I)*U(K-1,J  ,I  ) + DX(I-1)*U(K-1,J  ,I-1))
     +  +0.5*DY(J-1)*(DX(I)*U(K-1,J-1,I  ) + DX(I-1)*U(K-1,J-1,I-1))
     +  +0.5*DY(J  )*(DX(I)*U(K-1,J+1,I  ) + DX(I-1)*U(K-1,J+1,I-1))))

C            (K+1) - SCHEIBE 

     +  +(0.5*DZ(K )*
     +     (DDY(J  )*(DX(I)*U(K+1,J  ,I  ) + DX(I-1)*U(K+1,J  ,I-1))
     +  +0.5*DY(J-1)*(DX(I)*U(K+1,J-1,I  ) + DX(I-1)*U(K+1,J-1,I-1))
     +  +0.5*DY(J  )*(DX(I)*U(K+1,J+1,I  ) + DX(I-1)*U(K+1,J+1,I-1))))

C           ... GETEILT DURCH DAS GESAMTVOLUMEN

     +   )/VOL
C
      ENDDO
      ENDDO
      ENDDO
C
C
        WRITE (KANREC) RTREC
C
C       WRITE (44,*)'U RAUSGESCHRIEBEN'
C
C                     V V V V V V V V V V V V V V V V V V V V V V V V
C
C
      DO  ILAUF=1,INXREC
      DO  JLAUF=1,INYREC
      DO  KLAUF=1,INZREC

          K = KARR(IREC,KLAUF)
          J = JARR(IREC,JLAUF)
          I = IARR(IREC,ILAUF)

        VOL = (DX(I-1)+DX(I)) * (DY(J-1)+DY(J)) * (DZ(K-1)+DZ(K))

      RTREC(KLAUF,JLAUF,ILAUF) =

C            K - SCHEIBE 

     +    (
     +     (DDZ(K)*
     +     (DDX(I  )*(DY(J)*V(K  ,J  ,I  ) + DY(J-1)*V(K  ,J-1,I  ))
     +  +0.5*DX(I-1)*(DY(J)*V(K  ,J  ,I-1) + DY(J-1)*V(K  ,J-1,I-1))
     +  +0.5*DX(I  )*(DY(J)*V(K  ,J  ,I+1) + DY(J-1)*V(K  ,J-1,I+1))))

C            (K-1) - SCHEIBE 

     +  +(0.5*DZ(K-1)*
     +     (DDX(I  )*(DY(J)*V(K-1,J  ,I  ) + DY(J-1)*V(K-1,J-1,I  ))
     +  +0.5*DX(I-1)*(DY(J)*V(K-1,J  ,I-1) + DY(J-1)*V(K-1,J-1,I-1))
     +  +0.5*DX(I  )*(DY(J)*V(K-1,J  ,I+1) + DY(J-1)*V(K-1,J-1,I+1))))

C            (K+1) - SCHEIBE 

     +  +(0.5*DZ(K )*
     +     (DDX(I  )*(DY(J)*V(K+1,J  ,I  ) + DY(J-1)*V(K+1,J-1,I  ))
     +  +0.5*DX(I-1)*(DY(J)*V(K+1,J  ,I-1) + DY(J-1)*V(K+1,J-1,I-1))
     +  +0.5*DX(I  )*(DY(J)*V(K+1,J  ,I+1) + DY(J-1)*V(K+1,J-1,I+1))))

C           ... GETEILT DURCH DAS GESAMTVOLUMEN

     +   )/VOL
C
      ENDDO
      ENDDO
      ENDDO
C
C
        WRITE (KANREC) RTREC
C
C       WRITE (44,*)'V RAUSGESCHRIEBEN'
C
C
C
C                     W W W W W W W W W W W W W W W W W W W W W W W W
C
C
C
      DO  ILAUF=1,INXREC
      DO  JLAUF=1,INYREC
      DO  KLAUF=1,INZREC

          K = KARR(IREC,KLAUF)
          J = JARR(IREC,JLAUF)
          I = IARR(IREC,ILAUF)

        VOL = (DX(I-1)+DX(I)) * (DY(J-1)+DY(J)) * (DZ(K-1)+DZ(K))

      RTREC(KLAUF,JLAUF,ILAUF) =

C            I - SCHEIBE 

     +    (
     +     (DDX(I)*
     +     (DDY(J  )*(DZ(K)*W(K  ,J  ,I  ) + DZ(K-1)*W(K-1,J  ,I  ))
     +  +0.5*DY(J-1)*(DZ(K)*W(K  ,J-1,I  ) + DZ(K-1)*W(K-1,J-1,I  ))
     +  +0.5*DY(J  )*(DZ(K)*W(K  ,J+1,I  ) + DZ(K-1)*W(K-1,J+1,I  ))))

C            (I-1) - SCHEIBE 

     +  +(0.5*DX(I-1)*
     +     (DDY(J  )*(DZ(K)*W(K  ,J  ,I-1) + DZ(K-1)*W(K-1,J  ,I-1))
     +  +0.5*DY(J-1)*(DZ(K)*W(K  ,J-1,I-1) + DZ(K-1)*W(K-1,J-1,I-1))
     +  +0.5*DY(J  )*(DZ(K)*W(K  ,J+1,I-1) + DZ(K-1)*W(K-1,J+1,I-1))))

C            (I+1) - SCHEIBE 

     +  +(0.5*DX(I )*
     +     (DDY(J  )*(DZ(K)*W(K  ,J  ,I+1) + DZ(K-1)*W(K-1,J  ,I+1))
     +  +0.5*DY(J-1)*(DZ(K)*W(K  ,J-1,I+1) + DZ(K-1)*W(K-1,J-1,I+1))
     +  +0.5*DY(J  )*(DZ(K)*W(K  ,J+1,I+1) + DZ(K-1)*W(K-1,J+1,I+1))))

C           ... GETEILT DURCH DAS GESAMTVOLUMEN

     +   )/VOL
C
      ENDDO
      ENDDO
      ENDDO
C
        WRITE (KANREC) RTREC
C
C       WRITE (44,*)'W RAUSGESCHRIEBEN'
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       IF(IVEL.GE.3) THEN
C
C
C
C                     P P P P P P P P P P P P P P P P P P P P P P P P
C
C
      DO  ILAUF=1,INXREC
      DO  JLAUF=1,INYREC
      DO  KLAUF=1,INZREC

          K = KARR(IREC,KLAUF)
          J = JARR(IREC,JLAUF)
          I = IARR(IREC,ILAUF)

        VOL = (DX(I-1)+DX(I)) * (DY(J-1)+DY(J)) * (DZ(K-1)+DZ(K))

      RTREC(KLAUF,JLAUF,ILAUF) =

C                             I - SCHEIBE 

     +    (
     +     (DDX(I)*

C                             J 
     +     (DDY(J  )*
     +               (DDZ(K  )*P(K  ,J  ,I  ) 
     +               + DZ(K-1)*P(K-1,J  ,I  ) 
     +               + DZ(K+1)*P(K+1,J  ,I  ))

C                             J - 1
     +     +0.5*DY(J-1)*
     +               (DDZ(K  )*P(K  ,J-1,I  ) 
     +               + DZ(K-1)*P(K-1,J-1,I  ) 
     +               + DZ(K+1)*P(K+1,J-1,I  ))

C                             J + 1
     +     +0.5*DY(J  )*
     +               (DDZ(K  )*P(K  ,J+1,I  ) 
     +               + DZ(K-1)*P(K-1,J+1,I  ) 
     +               + DZ(K+1)*P(K+1,J+1,I  ))))

C            (I-1) - SCHEIBE 

     +  +(0.5*DX(I-1)*

C                             J 
     +     (DDY(J  )*
     +               (DDZ(K  )*P(K  ,J  ,I-1) 
     +               + DZ(K-1)*P(K-1,J  ,I-1) 
     +               + DZ(K+1)*P(K+1,J  ,I-1))

C                             J - 1
     +     +0.5*DY(J-1)*
     +               (DDZ(K  )*P(K  ,J-1,I-1) 
     +               + DZ(K-1)*P(K-1,J-1,I-1) 
     +               + DZ(K+1)*P(K+1,J-1,I-1))

C                             J + 1
     +     +0.5*DY(J  )*
     +               (DDZ(K  )*P(K  ,J+1,I-1) 
     +               + DZ(K-1)*P(K-1,J+1,I-1) 
     +               + DZ(K+1)*P(K+1,J+1,I-1))))



C            (I+1) - SCHEIBE 

     +  +(0.5*DX(I )*

C                             J 
     +     (DDY(J  )*
     +               (DDZ(K  )*P(K  ,J  ,I+1) 
     +               + DZ(K-1)*P(K-1,J  ,I+1) 
     +               + DZ(K+1)*P(K+1,J  ,I+1))

C                             J - 1
     +     +0.5*DY(J-1)*
     +               (DDZ(K  )*P(K  ,J-1,I+1) 
     +               + DZ(K-1)*P(K-1,J-1,I+1) 
     +               + DZ(K+1)*P(K+1,J-1,I+1))

C                             J + 1
     +     +0.5*DY(J  )*
     +               (DDZ(K  )*P(K  ,J+1,I+1)
     +               + DZ(K-1)*P(K-1,J+1,I+1)
     +               + DZ(K+1)*P(K+1,J+1,I+1))))



C           ... GETEILT DURCH DAS GESAMTVOLUMEN

     +   )/VOL
C
      ENDDO
      ENDDO
      ENDDO
C
C
        WRITE (KANREC) RTREC
C
C       WRITE (44,*)'P RAUSGESCHRIEBEN'
C
C
        ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       RETURN
       END
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
        SUBROUTINE FLTVEL2 (KANREC,NMREC,IREC,
     +               CIDREC,IVEL,
     +               INXREC,INYREC,INZREC,
     +               IARR,JARR,KARR,
     +               RTREC,
     +               U,V,W,P,
     +               KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ)
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
        INTEGER
     +           KANREC,IVEL
     +          ,INXREC,INYREC,INZREC
     +          ,IARR(NMREC,II),JARR(NMREC,JJ),KARR(NMREC,KK)
C
C
        DIMENSION  RTREC(INZREC,INYREC,INXREC)
C
      REAL    U(KK,JJ,II),V(KK,JJ,II),W(KK,JJ,II),
     $        P(KK,JJ,II),
     $        DX(II),DY(JJ),DZ(KK),DDX(II),DDY(JJ),DDZ(KK)
C
       CHARACTER (LEN=16) CIDREC(NMREC)
C
        LOGICAL LREC
C
C                     U U U U U U U U U U U U U U U U U U U U U U U U
C
C
      DO  ILAUF=1,INXREC
      DO  JLAUF=1,INYREC
      DO  KLAUF=1,INZREC

          K = KARR(IREC,KLAUF)
          J = JARR(IREC,JLAUF)
          I = IARR(IREC,ILAUF)

          VOL = (DDX(I)+DDX(I+1))*(DDY(J)+DDY(J+1))*(DDZ(K)+DDZ(K+1))

      RTREC(KLAUF,JLAUF,ILAUF) =

C            (I  ) - SCHEIBE 

     $      (0.5*(DDX(I)+DDX(I+1))*
     $              (DDY(J  )*
     $                        (DDZ(K  )*U(K  ,J  ,I  ) + 
     $                         DDZ(K+1)*U(K+1,J  ,I  )  ) +
     $               DDY(J+1)*
     $                        (DDZ(K  )*U(K  ,J+1,I  ) + 
     $                         DDZ(K+1)*U(K+1,J+1,I  )  )   )

C            (I-1) - SCHEIBE 

     $  +  0.5*DDX(I)*
     $              (DDY(J  )*
     $                        (DDZ(K  )*U(K  ,J  ,I-1) + 
     $                         DDZ(K+1)*U(K+1,J  ,I-1)  ) +
     $               DDY(J+1)*
     $                        (DDZ(K  )*U(K  ,J+1,I-1) + 
     $                         DDZ(K+1)*U(K+1,J+1,I-1)  )   )


C            (I+1) - SCHEIBE 

     $  +  0.5*DDX(I+1)*
     $              (DDY(J  )*
     $                        (DDZ(K  )*U(K  ,J  ,I+1) + 
     $                         DDZ(K+1)*U(K+1,J  ,I+1)  ) +
     $               DDY(J+1)*
     $                        (DDZ(K  )*U(K  ,J+1,I+1) + 
     $                         DDZ(K+1)*U(K+1,J+1,I+1)  )   )  

C           ... GETEILT DURCH DAS GESAMTVOLUMEN

     $   )/VOL
C
      ENDDO
      ENDDO
      ENDDO
C
C
        WRITE (KANREC) RTREC
C
C       WRITE (44,*)'U RAUSGESCHRIEBEN'
C
C                     V V V V V V V V V V V V V V V V V V V V V V V V
C
C
      DO  ILAUF=1,INXREC
      DO  JLAUF=1,INYREC
      DO  KLAUF=1,INZREC

          K = KARR(IREC,KLAUF)
          J = JARR(IREC,JLAUF)
          I = IARR(IREC,ILAUF)

          VOL = (DDX(I)+DDX(I+1))*(DDY(J)+DDY(J+1))*(DDZ(K)+DDZ(K+1))

      RTREC(KLAUF,JLAUF,ILAUF) =

C            (J  ) - SCHEIBE 

     $      (0.5*(DDY(J)+DDY(J+1))*
     $              (DDX(I  )*
     $                        (DDZ(K  )*V(K  ,J  ,I  ) + 
     $                         DDZ(K+1)*V(K+1,J  ,I  )  ) +
     $               DDX(I+1)*
     $                        (DDZ(K  )*V(K  ,J  ,I+1) + 
     $                         DDZ(K+1)*V(K+1,J  ,I+1)  )   )  

C            (J-1) - SCHEIBE 

     $  +0.5*DDY(J)*
     $              (DDX(I  )*
     $                        (DDZ(K  )*V(K  ,J-1,I  ) + 
     $                         DDZ(K+1)*V(K+1,J-1,I  )  ) +
     $               DDX(I+1)*
     $                        (DDZ(K  )*V(K  ,J-1,I+1) + 
     $                         DDZ(K+1)*V(K+1,J-1,I+1)  )   )  


C            (J+1) - SCHEIBE 

     $  +0.5*DDY(J+1)*
     $              (DDX(I  )*
     $                        (DDZ(K  )*V(K  ,J+1,I  ) + 
     $                         DDZ(K+1)*V(K+1,J+1,I  )  ) +
     $               DDX(I+1)*
     $                        (DDZ(K  )*V(K  ,J+1,I+1) + 
     $                         DDZ(K+1)*V(K+1,J+1,I+1)  )   )  

C           ... GETEILT DURCH DAS GESAMTVOLUMEN

     +   )/VOL
C
      ENDDO
      ENDDO
      ENDDO
C
C
        WRITE (KANREC) RTREC
C
C       WRITE (44,*)'V RAUSGESCHRIEBEN'
C
C
C
C                     W W W W W W W W W W W W W W W W W W W W W W W W
C
C
C
      DO  ILAUF=1,INXREC
      DO  JLAUF=1,INYREC
      DO  KLAUF=1,INZREC

          K = KARR(IREC,KLAUF)
          J = JARR(IREC,JLAUF)
          I = IARR(IREC,ILAUF)

          VOL = (DDX(I)+DDX(I+1))*(DDY(J)+DDY(J+1))*(DDZ(K)+DDZ(K+1))

      RTREC(KLAUF,JLAUF,ILAUF) =

C            (K  ) - SCHEIBE 

     $      (0.5*(DDZ(K)+DDZ(K+1))*
     $              (DDX(I  )*
     $                        (DDY(J  )*W(K  ,J  ,I  ) + 
     $                         DDY(J+1)*W(K  ,J+1,I  )  ) +
     $               DDX(I+1)*
     $                        (DDY(J  )*W(K  ,J  ,I+1) + 
     $                         DDY(J+1)*W(K  ,J+1,I+1)  )   )  

C            (K-1) - SCHEIBE 

     $  +0.5*DDZ(K)*
     $              (DDX(I  )*
     $                        (DDY(J  )*W(K-1,J  ,I  ) + 
     $                         DDY(J+1)*W(K-1,J+1,I  )  ) +
     $               DDX(I+1)*
     $                        (DDY(J  )*W(K-1,J  ,I+1) + 
     $                         DDY(J+1)*W(K-1,J+1,I+1)  )   )  


C            (K+1) - SCHEIBE 

     $  +0.5*DDZ(K+1)*
     $              (DDX(I  )*
     $                        (DDY(J  )*W(K+1,J  ,I  ) + 
     $                         DDY(J+1)*W(K+1,J+1,I  )  ) +
     $               DDX(I+1)*
     $                        (DDY(J  )*W(K+1,J  ,I+1) + 
     $                         DDY(J+1)*W(K+1,J+1,I+1)  )   )  

C           ... GETEILT DURCH DAS GESAMTVOLUMEN

     $   )/VOL
C
      ENDDO
      ENDDO
      ENDDO
C
        WRITE (KANREC) RTREC
C
C       WRITE (44,*)'W RAUSGESCHRIEBEN'
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       IF(IVEL.GE.3) THEN
C
C
C
C                     P P P P P P P P P P P P P P P P P P P P P P P P
C
C
      DO  ILAUF=1,INXREC
      DO  JLAUF=1,INYREC
      DO  KLAUF=1,INZREC

          K = KARR(IREC,KLAUF)
          J = JARR(IREC,JLAUF)
          I = IARR(IREC,ILAUF)

C                       BEIM DRUCK EINFACHE MITTELUNG, 
C                       DA PUNKT GENAU IN DER MITTE LIEGT

      RTREC(KLAUF,JLAUF,ILAUF) = 0.125*

     $         P(K  ,J  ,I  )+
     $         P(K+1,J  ,I  )+
     $         P(K  ,J+1,I  )+
     $         P(K+1,J+1,I  )+
     $         P(K  ,J  ,I+1)+
     $         P(K+1,J  ,I+1)+
     $         P(K  ,J+1,I+1)+
     $         P(K+1,J+1,I+1)
C
      ENDDO
      ENDDO
      ENDDO
C
C
        WRITE (KANREC) RTREC
C
C       WRITE (44,*)'P RAUSGESCHRIEBEN'
C
C
        ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       RETURN
       END
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
        SUBROUTINE ACOUSTICTERM (KANREC,NMREC,IREC,
     +               CIDREC,IVEL,
     +               INXREC,INYREC,INZREC,
     +               IARR,JARR,KARR,
     +               RTREC,
     +               U,V,W,P,H1,H2,
     +               KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ)
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
        INTEGER
     +           KANREC,IVEL
     +          ,INXREC,INYREC,INZREC
     +          ,IARR(NMREC,II),JARR(NMREC,JJ),KARR(NMREC,KK)
C
C
        DIMENSION  RTREC(INZREC,INYREC,INXREC)
C
      REAL    U(KK,JJ,II),V(KK,JJ,II),W(KK,JJ,II),
     $        P(KK,JJ,II),
     $        DX(II),DY(JJ),DZ(KK),DDX(II),DDY(JJ),DDZ(KK)
C
       CHARACTER (LEN=16) CIDREC(NMREC)
C
        LOGICAL LREC
C
C---- -------------------------------------------- 
C
C
      DO  ILAUF=1,INXREC
      DO  JLAUF=1,INYREC
      DO  KLAUF=1,INZREC

          K = KARR(IREC,KLAUF)
          J = JARR(IREC,JLAUF)
          I = IARR(IREC,ILAUF)

          UUE = 0.25*(U(K  ,J  ,I+1) + U(K  ,J+1,I+1)+
     $                U(K+1,J  ,I+1) + U(K+1,J+1,I+1) )**2
          UUP = 0.25*(U(K  ,J  ,I  ) + U(K  ,J+1,I  )+
     $                U(K+1,J  ,I  ) + U(K+1,J+1,I  ) )**2
          UUW = 0.25*(U(K  ,J  ,I-1) + U(K  ,J+1,I-1)+
     $                U(K+1,J  ,I-1) + U(K+1,J+1,I-1) )**2

          D2UUDX2  = (UUE - 2.0*UUP + UUW) / ((DDX(I) * DDX(I+1)))

C---------------------------------------------------------------------

          VVN = 0.25*(V(K  ,J+1,I  ) + V(K  ,J+1,I+1)+
     $                V(K+1,J+1,I  ) + V(K+1,J+1,I+1) )**2
          VVP = 0.25*(V(K  ,J  ,I  ) + V(K  ,J  ,I+1)+
     $                V(K+1,J  ,I  ) + V(K+1,J  ,I+1) )**2
          VVS = 0.25*(V(K  ,J-1,I  ) + V(K  ,J-1,I+1)+
     $                V(K+1,J-1,I  ) + V(K+1,J-1,I+1) )**2
          
          D2VVDY2  = (VVN - 2.0*VVP + VVS) / ((DDY(J) * DDY(J+1)))

C---------------------------------------------------------------------

          WWT = 0.25*(W(K+1,J  ,I  ) + W(K+1,J  ,I+1)+
     $                W(K+1,J+1,I  ) + W(K+1,J+1,I+1) )**2
          WWP = 0.25*(W(K  ,J  ,I  ) + W(K  ,J  ,I+1)+
     $                W(K  ,J+1,I  ) + W(K  ,J+1,I+1) )**2
          WWB = 0.25*(W(K-1,J  ,I  ) + W(K-1,J  ,I+1)+
     $                W(K-1,J+1,I  ) + W(K-1,J+1,I+1) )**2

          D2WWDZ2  = (WWT - 2.0*WWP + WWB) / ((DDZ(K) * DDZ(K+1)))

C---------------------------------------------------------------------

          UVNE = 0.5*(( 0.5*(U(K  ,J+1,I+1) + U(K  ,J+2,I+1))*
     $                  0.5*(V(K  ,J+1,I+1) + V(K  ,J+1,I+2)) )+
     $                ( 0.5*(U(K+1,J+1,I+1) + U(K+1,J+2,I+1))*
     $                  0.5*(V(K+1,J+1,I+1) + V(K+1,J+1,I+2)) ))
          UVSE = 0.5*(( 0.5*(U(K  ,J-1,I+1) + U(K  ,J  ,I+1))*
     $                  0.5*(V(K  ,J-1,I+1) + V(K  ,J-1,I+2)) )+
     $                ( 0.5*(U(K+1,J-1,I+1) + U(K+1,J  ,I+1))*
     $                  0.5*(V(K+1,J-1,I+1) + V(K+1,J-1,I+2)) ))

          DUVDYE = (UVNE - UVSE)/(DDY(J) + DDY(J+1))


          UVNW = 0.5*(( 0.5*(U(K  ,J+1,I-1) + U(K  ,J+2,I-1))*
     $                  0.5*(V(K  ,J+1,I-1) + V(K  ,J+1,I  )) )+
     $                ( 0.5*(U(K+1,J+1,I-1) + U(K+1,J+2,I-1))*
     $                  0.5*(V(K+1,J+1,I-1) + V(K+1,J+1,I  )) ))
          UVSW = 0.5*(( 0.5*(U(K  ,J-1,I-1) + U(K  ,J  ,I-1))*
     $                  0.5*(V(K  ,J-1,I-1) + V(K  ,J-1,I  )) )+
     $                ( 0.5*(U(K+1,J-1,I-1) + U(K+1,J  ,I-1))*
     $                  0.5*(V(K+1,J-1,I-1) + V(K+1,J-1,I  )) ))

          DUVDYW = (UVNW - UVSW)/(DDY(J) + DDY(J+1))

          D2UVDXDY = (DUVDYE - DUVDYW) / (DDX(I) + DDX(I+1))

C---------------------------------------------------------------------

          UWTE = 0.5*(( 0.5*(U(K+1,J  ,I+1) + U(K+2,J  ,I+1))*
     $                  0.5*(W(K+1,J  ,I+1) + W(K+1,J  ,I+2)) )+
     $                ( 0.5*(U(K+1,J+1,I+1) + U(K+2,J+1,I+1))*
     $                  0.5*(W(K+1,J+1,I+1) + W(K+1,J+1,I+2)) ))
          UWBE = 0.5*(( 0.5*(U(K-1,J  ,I+1) + U(K  ,J  ,I+1))*
     $                  0.5*(W(K-1,J  ,I+1) + W(K-1,J  ,I+2)) )+
     $                ( 0.5*(U(K-1,J+1,I+1) + U(K  ,J+1,I+1))*
     $                  0.5*(W(K-1,J+1,I+1) + W(K-1,J+1,I+2)) ))

          DUWDZE = (UWTE - UWBE)/(DDZ(K) + DDZ(K+1))


          UWTW = 0.5*(( 0.5*(U(K+1,J  ,I-1) + U(K+2,J  ,I-1))*
     $                  0.5*(W(K+1,J  ,I-1) + W(K+1,J  ,I  )) )+
     $                ( 0.5*(U(K+1,J+1,I-1) + U(K+2,J+1,I-1))*
     $                  0.5*(W(K+1,J+1,I-1) + W(K+1,J+1,I  )) ))
          UWBW = 0.5*(( 0.5*(U(K-1,J  ,I-1) + U(K  ,J  ,I-1))*
     $                  0.5*(W(K-1,J  ,I-1) + W(K-1,J  ,I  )) )+
     $                ( 0.5*(U(K-1,J+1,I-1) + U(K  ,J+1,I-1))*
     $                  0.5*(W(K-1,J+1,I-1) + W(K-1,J+1,I  )) ))

          DUWDZW = (UWTW - UWBW)/(DDZ(K) + DDZ(K+1))


          D2UWDXDZ = (DUWDZE - DUWDZW) / (DDX(I) + DDX(I+1))

C---------------------------------------------------------------------

          D2VWDYDZ = (DVWDZN - DVWDZS) / (DDY(J) + DDY(J+1))

          VWTN = 0.5*(( 0.5*(V(K+1,J+1,I  ) + V(K+2,J+1,I  ))*
     $                  0.5*(W(K+1,J+1,I  ) + W(K+1,J+2,I  )) )+
     $                ( 0.5*(V(K+1,J+1,I+1) + V(K+2,J+1,I+1))*
     $                  0.5*(W(K+1,J+1,I+1) + W(K+1,J+2,I+1)) ))
          VWBN = 0.5*(( 0.5*(V(K-1,J+1,I  ) + V(K  ,J+1,I  ))*
     $                  0.5*(W(K-1,J+1,I  ) + W(K-1,J+2,I  )) )+
     $                ( 0.5*(V(K-1,J+1,I+1) + V(K  ,J+1,I+1))*
     $                  0.5*(W(K-1,J+1,I+1) + W(K-1,J+2,I+1)) ))

          DVWDZN = (VWTN - VWBN)/(DDZ(K) + DDZ(K+1))


          UWTW = 0.5*(( 0.5*(V(K+1,J-1,I  ) + V(K+2,J-1,I  ))*
     $                  0.5*(W(K+1,J-1,I  ) + W(K+1,J  ,I  )) )+
     $                ( 0.5*(V(K+1,J-1,I+1) + V(K+2,J-1,I+1))*
     $                  0.5*(W(K+1,J-1,I+1) + W(K+1,J  ,I+1)) ))
          UWBW = 0.5*(( 0.5*(V(K-1,J-1,I  ) + V(K  ,J-1,I  ))*
     $                  0.5*(W(K-1,J-1,I  ) + W(K-1,J  ,I  )) )+
     $                ( 0.5*(V(K-1,J-1,I+1) + V(K  ,J-1,I+1))*
     $                  0.5*(W(K-1,J-1,I+1) + W(K-1,J  ,I+1)) ))

          DVWDZS = (VWTS - VWBS)/(DDZ(K) + DDZ(K+1))


          D2VWDYDZ = (DVWDZN - DVWDZS) / (DDY(J) + DDY(J+1))

C---------------------------------------------------------------------

      RTREC(KLAUF,JLAUF,ILAUF) = D2UUDX2  + D2UVDXDY + D2UWDXDZ +
     $                           D2UVDXDY + D2VVDY2  + D2VWDYDZ +
     $                           D2UWDXDZ + D2VWDYDZ + D2WWDZ2 

      ENDDO
      ENDDO
      ENDDO
C
C
        WRITE (KANREC) RTREC
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       RETURN
       END
C
        SUBROUTINE Q_KRITERIUM (KANREC,NMREC,IREC,
     +               CIDREC,IVEL,
     +               INXREC,INYREC,INZREC,
     +               IARR,JARR,KARR,
     +               RTREC,
     +               U,V,W,P,H1,H2,
     +               KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ)
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
        INTEGER
     +           KANREC,IVEL
     +          ,INXREC,INYREC,INZREC
     +          ,IARR(NMREC,II),JARR(NMREC,JJ),KARR(NMREC,KK)
C
C
        DIMENSION  RTREC(INZREC,INYREC,INXREC)
C
      REAL    U(KK,JJ,II),V(KK,JJ,II),W(KK,JJ,II),
     $        P(KK,JJ,II),
     $        DX(II),DY(JJ),DZ(KK),DDX(II),DDY(JJ),DDZ(KK)
C
        CHARACTER*16 CIDREC(NMREC)
C
      REAL    D2X,D2Y,D2Z
C
        LOGICAL LREC
C
C---- -------------------------------------------- 
C
C

      DO  ILAUF=1,INXREC
      DO  JLAUF=1,INYREC
      DO  KLAUF=1,INZREC

          K = KARR(IREC,KLAUF)
          J = JARR(IREC,JLAUF)
          I = IARR(IREC,ILAUF)

	  IF (I.EQ.1) THEN
	     D2X =  
     $           ( 2.0*P( K,J,1) 
     $            -5.0*P( K,J,2)
     $            +4.0*P( K,J,3)
     $            -1.0*P( K,J,4) )
     $          /(DX(1)**2)
	  ENDIF

          IF (I.EQ.II) THEN
	     D2X =  
     $           ( 2.0*P(K,J,II  ) 
     $            -5.0*P(K,J,II-1)
     $            +4.0*P(K,J,II-2)
     $            -1.0*P(K,J,II-3) )
     $          /(DX(II-1)**2)
	  ENDIF

	  IF ((I.GE.2).AND.(I.LE.II-1)) THEN
	      D2X = 
     $            (P(K,J,I+1)+P(K,J,I-1)-2.0*P(K,J,I))/(DX(I)*DX(I-1))
          ENDIF

CCCCCCCCCCCCCCCCCCCCCCC
	  
	  IF (J.EQ.1) THEN
	     D2Y =  
     $           ( 2.0*P( K,1,I) 
     $            -5.0*P( K,2,I)
     $            +4.0*P( K,3,I)
     $            -1.0*P( K,4,I) )
     $          /(DY(1)**2)
	  ENDIF

          IF (J.EQ.JJ) THEN
	     D2Y =  
     $           ( 2.0*P(K,J ,I) 
     $            -5.0*P(K,J-1,I)
     $            +4.0*P(K,J-2,I)
     $            -1.0*P(K,J-3,I) )
     $          /(DY(JJ-1)**2)
	  ENDIF

	  IF ((J.GE.2).AND.(J.LE.JJ-1)) THEN
	      D2Y = 
     $            (P(K,J+1,I)+P(K,J-1,I)-2.0*P(K,J,I))/(DY(J)*DY(J-1))
          ENDIF

CCCCCCCCCCCCCCCCCCCCCCC
	  
	  IF (K.EQ.1) THEN
	     D2Z =  
     $           ( 2.0*P( 1,J,I) 
     $            -5.0*P( 2,J,I)
     $            +4.0*P( 3,J,I)
     $            -1.0*P( 4,J,I) )
     $          /(DZ(1)**2)
	  ENDIF

          IF (K.EQ.KK) THEN
	     D2Z =  
     $           ( 2.0*P(K  ,J,I) 
     $            -5.0*P(K-1,J,I)
     $            +4.0*P(K-2,J,I)
     $            -1.0*P(K-3,J,I) )
     $          /(DZ(KK-1)**2)
	  ENDIF

	  IF ((K.GE.2).AND.(K.LE.KK-1)) THEN
	     D2Z = 
     $            (P(K+1,J,I)+P(K-1,J,I)-2.0*P(K,J,I))/(DZ(K)*DZ(K-1))
          ENDIF
	  
CCCCCCCCCCCCCCCCCCCCCCC

          RTREC(KLAUF,JLAUF,ILAUF) =  D2X +  D2Y +  D2Z

CCCCCCCCCCCCCCCCCCCCCCC
 
      ENDDO
      ENDDO
      ENDDO
C
C
        WRITE (KANREC) RTREC
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       RETURN
       END
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
        SUBROUTINE ACOUSTIC2 (KANREC,NMREC,IREC,
     +               CIDREC,IVEL,
     +               INXREC,INYREC,INZREC,
     +               IARR,JARR,KARR,
     +               RTREC,
     +               U,V,W,P,H1,H2,
     +               KK,JJ,II,
     +               KMX,JMX,IMX,
     +               DX,DY,DZ,DDX,DDY,DDZ)
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
        INTEGER
     +           KANREC,IVEL
     +          ,INXREC,INYREC,INZREC
     +          ,IARR(NMREC,II),JARR(NMREC,JJ),KARR(NMREC,KK)
C
C
        DIMENSION  RTREC(INZREC,INYREC,INXREC)
C
      REAL    U(KK,JJ,II),V(KK,JJ,II),W(KK,JJ,II),
     $        P(KK,JJ,II),H1(KK,JJ,II),H2(KK,JJ,II),
     $        DX(II),DY(JJ),DZ(KK),DDX(II),DDY(JJ),DDZ(KK)
C
       CHARACTER (LEN=16) CIDREC(NMREC)
C
        LOGICAL LREC
C
C---- -------------------------------------------- 
C        RÜCKSETZEN VON RTREC

        CALL SETS(INZREC,INYREC,INXREC,INZREC,INYREC,INXREC,RTREC,0.0)
C
C---- -------------------------------------------- 
C                    BERECHNUNG VON D2UUDX2
C
	CALL NONSTAGX (KK,JJ,II,U,H1,DX,DDX)
 	CALL MULTIPLY (KK,JJ,II,H1,H1,H2)
 	CALL D2DX2 (KK,JJ,II,DX,DDX,H2,H1)
        CALL ADDFLT (KK,JJ,II,INXREC,INYREC,INZREC,NMREC,IREC,
     $              IARR,JARR,KARR,
     $              RTREC,H1)
C
C---- -------------------------------------------- 
C                    BERECHNUNG VON D2VVDY2
C
	CALL NONSTAGY (KK,JJ,II,V,H1,DY,DDY)
 	CALL MULTIPLY (KK,JJ,II,H1,H1,H2)
 	CALL D2DY2 (KK,JJ,II,DY,DDY,H2,H1)
        CALL ADDFLT (KK,JJ,II,INXREC,INYREC,INZREC,NMREC,IREC,
     $              IARR,JARR,KARR,
     $              RTREC,H1)
C
C---- -------------------------------------------- 
C                    BERECHNUNG VON D2WWDZ2
C
	CALL NONSTAGY (KK,JJ,II,W,H1,DZ,DDZ)
 	CALL MULTIPLY (KK,JJ,II,H1,H1,H2)
 	CALL D2DZ2 (KK,JJ,II,DZ,DDZ,H2,H1)
        CALL ADDFLT (KK,JJ,II,INXREC,INYREC,INZREC,NMREC,IREC,
     $              IARR,JARR,KARR,
     $              RTREC,H1)
C
C---- -------------------------------------------- 
C                    BERECHNUNG VON D2UVDXDY

	CALL NONSTAGX (KK,JJ,II,U,H1,DX,DDX)
	CALL NONSTAGY (KK,JJ,II,V,H2,DY,DDY)
 	CALL MULTIPLY (KK,JJ,II,H1,H2,H2)
 	CALL DPHIDX   (KK,JJ,II,DX,DDX,H2,H1)
 	CALL DPHIDY   (KK,JJ,II,DY,DDY,H1,H2)
        CALL ADDFLT (KK,JJ,II,INXREC,INYREC,INZREC,NMREC,IREC,
     $              IARR,JARR,KARR,
     $              RTREC,H2)
        CALL ADDFLT (KK,JJ,II,INXREC,INYREC,INZREC,NMREC,IREC,
     $              IARR,JARR,KARR,
     $              RTREC,H2)
C
C---- -------------------------------------------- 
C                    BERECHNUNG VON D2UWDXDZ

	CALL NONSTAGX (KK,JJ,II,U,H1,DX,DDX)
	CALL NONSTAGZ (KK,JJ,II,W,H2,DZ,DDZ)
 	CALL MULTIPLY (KK,JJ,II,H1,H2,H2)
 	CALL DPHIDX   (KK,JJ,II,DX,DDX,H2,H1)
 	CALL DPHIDZ   (KK,JJ,II,DZ,DDZ,H1,H2)
        CALL ADDFLT (KK,JJ,II,INXREC,INYREC,INZREC,NMREC,IREC,
     $              IARR,JARR,KARR,
     $              RTREC,H2)
        CALL ADDFLT (KK,JJ,II,INXREC,INYREC,INZREC,NMREC,IREC,
     $              IARR,JARR,KARR,
     $              RTREC,H2)
C
C---- -------------------------------------------- 
C                    BERECHNUNG VON D2VWDYDZ

	CALL NONSTAGY (KK,JJ,II,V,H2,DY,DDY)
	CALL NONSTAGZ (KK,JJ,II,W,H1,DZ,DDZ)
 	CALL MULTIPLY (KK,JJ,II,H1,H2,H2)
 	CALL DPHIDY   (KK,JJ,II,DY,DDY,H2,H1)
 	CALL DPHIDZ   (KK,JJ,II,DZ,DDZ,H1,H2)
        CALL ADDFLT (KK,JJ,II,INXREC,INYREC,INZREC,NMREC,IREC,
     $              IARR,JARR,KARR,
     $              RTREC,H2)
        CALL ADDFLT (KK,JJ,II,INXREC,INYREC,INZREC,NMREC,IREC,
     $              IARR,JARR,KARR,
     $              RTREC,H2)
C
C---- -------------------------------------------- 
C
        WRITE (KANREC) RTREC
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
       RETURN
       END

C---------------------------------------------------------------------------
      SUBROUTINE NONSTAGX(KK,JJ,II,U,H1,DX,DDX)
      REAL U(KK,JJ,II),H1(KK,JJ,II)
      REAL DX(II),DDX(II)
      
      DO I=1,2
         DO J=1,JJ
            DO K=1,KK

               H1(K,J,I) = U(K,J,2)

            ENDDO
         ENDDO
      ENDDO

      DO I=3,II-1
         DO J=1,JJ
            DO K=1,KK

               H1(K,J,I) = (DX(I)*U(K,J,I-1) + DX(I-1)*U(K,J,I))/
     $                      (DX(I)+DX(I-1))

            ENDDO
         ENDDO
      ENDDO
      DO I=II,II
         DO J=1,JJ
            DO K=1,KK

               H1(K,J,I) = U(K,J,II-1)

            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
C---------------------------------------------------------------------------
      SUBROUTINE NONSTAGY(KK,JJ,II,V,H1,DY,DDY)
      REAL V(KK,JJ,II),H1(KK,JJ,II)
      REAL DY(II),DDY(II)
      
      DO J=1,2
         DO I=1,II
            DO K=1,KK

               H1(K,J,I) = V(K,2,I)

            ENDDO
         ENDDO
      ENDDO

      DO J=3,JJ-1
         DO I=1,II
            DO K=1,KK

               H1(K,J,I) = (DY(J)*V(K,J-1,I) + DY(J-1)*V(K,J,I))/
     $                      (DY(J)+DY(J-1))

            ENDDO
         ENDDO
      ENDDO
      DO J=JJ,JJ
         DO I=1,II
            DO K=1,KK

               H1(K,J,I) = V(K,JJ-1,I)

            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
C---------------------------------------------------------------------------
      SUBROUTINE NONSTAGZ(KK,JJ,II,W,H1,DZ,DDZ)
      REAL W(KK,JJ,II),H1(KK,JJ,II)
      REAL DZ(KK),DDZ(KK)
      
      DO K=1,2
         DO J=1,JJ
            DO I=1,II

               H1(K,J,I) = W(2,J,I)

            ENDDO
         ENDDO
      ENDDO

      DO I=1,II
         DO J=1,JJ
            DO K=3,KK-1

               H1(K,J,I) = (DZ(K)*W(K-1,J,I) + DZ(K-1)*W(K,J,I))/
     $                      (DZ(K)+DZ(K-1))

            ENDDO
         ENDDO
      ENDDO
      DO K=KK,KK
         DO J=1,JJ
            DO I=1,II

               H1(K,J,I) = W(KK-1,J,I)

            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
C---------------------------------------------------------------------------
      SUBROUTINE MULTIPLY (KK,JJ,II,H1,H2,H3)
      REAL H1(KK*JJ*II),H2(KK*JJ*II),H3(KK*JJ*II)
      
      DO I=1,KK*JJ*II
         H3(I) = H1(I)*H2(I)
      ENDDO
      RETURN
      END
C---------------------------------------------------------------------------
      SUBROUTINE D2DX2 (KK,JJ,II,DX,DDX,H1,H2)
      REAL H1(KK,JJ,II),H2(KK,JJ,II)
      REAL DX(II),DDX(II)

      DO I=2,II-1
         DO J=1,JJ
            DO K=1,KK
               
               H2(K,J,I) = (H1(K,J,I-1) + H1(K,J,I+1) - 2.0*H1(K,J,I))/
     $                     (DX(I) * DX(I-1))
               
            ENDDO
         ENDDO
      ENDDO

      DO I=1,1
         DO J=1,JJ
            DO K=1,KK
               
               H2(K,J,I) = H2(K,J,I+1)
               
            ENDDO
         ENDDO
      ENDDO
      
      DO I=II,II
         DO J=1,JJ
            DO K=1,KK
               
               H2(K,J,I) = H2(K,J,I-1)
               
            ENDDO
         ENDDO
      ENDDO
      
      RETURN
      END
C---------------------------------------------------------------------------
      SUBROUTINE D2DY2 (KK,JJ,II,DY,DDY,H1,H2)
      REAL H1(KK,JJ,II),H2(KK,JJ,II)
      REAL DY(JJ),DDY(JJ)

      DO J=2,JJ-1
         DO I=1,II
            DO K=1,KK
               
               H2(K,J,I) = (H1(K,J-1,I) + H1(K,J+1,I) - 2.0*H1(K,J,I))/
     $                     (DY(J) * DY(J-1))
               
            ENDDO
         ENDDO
      ENDDO

      DO J=1,1
         DO I=1,II
            DO K=1,KK
               
               H2(K,J,I) = H2(K,J+1,I)
               
            ENDDO
         ENDDO
      ENDDO
      
      DO J=JJ,JJ
         DO I=1,II
            DO K=1,KK
               
               H2(K,J,I) = H2(K,J-1,I)
               
            ENDDO
         ENDDO
      ENDDO
      
      RETURN
      END
C---------------------------------------------------------------------------
      SUBROUTINE D2DZ2 (KK,JJ,II,DZ,DDZ,H1,H2)
      REAL H1(KK,JJ,II),H2(KK,JJ,II)
      REAL DZ(KK),DDZ(KK)

      DO I=1,II
         DO J=1,JJ
            DO K=2,KK-1
               
               H2(K,J,I) = (H1(K-1,J,I) + H1(K+1,J,I) - 2.0*H1(K,J,I))/
     $                     (DZ(K) * DZ(K-1))
               
            ENDDO
         ENDDO
      ENDDO

      DO K=1,1
         DO J=1,JJ
            DO I=1,II
               
               H2(K,J,I) = H2(K+1,J,I)
               
            ENDDO
         ENDDO
      ENDDO
      
      DO K=KK,KK
         DO J=1,JJ
            DO I=1,II
               
               H2(K,J,I) = H2(K-1,J,I)
               
            ENDDO
         ENDDO
      ENDDO
      
      RETURN
      END
      SUBROUTINE DPHIDX (KK,JJ,II,DX,DDX,H1,H2)
      REAL H1(KK,JJ,II),H2(KK,JJ,II)
      REAL DX(II),DDX(II)

      DO I=2,II-1
         DO J=1,JJ
            DO K=1,KK
               
               H2(K,J,I) = (H1(K,J,I+1) - H1(K,J,I-1))/(DX(I) + DX(I-1))
               
            ENDDO
         ENDDO
      ENDDO

      DO I=1,1
         DO J=1,JJ
            DO K=1,KK
               
               H2(K,J,I) = H2(K,J,I+1)
               
            ENDDO
         ENDDO
      ENDDO
      
      DO I=II,II
         DO J=1,JJ
            DO K=1,KK
               
               H2(K,J,I) = H2(K,J,I-1)
               
            ENDDO
         ENDDO
      ENDDO
      
      RETURN
      END
C---------------------------------------------------------------------------
      SUBROUTINE DPHIDY (KK,JJ,II,DY,DDY,H1,H2)
      REAL H1(KK,JJ,II),H2(KK,JJ,II)
      REAL DY(JJ),DDY(JJ)

      DO J=2,JJ-1
         DO I=1,II
            DO K=1,KK
               
               H2(K,J,I) = (H1(K,J+1,I) - H1(K,J-1,I))/(DY(J) + DY(J-1))
               
            ENDDO
         ENDDO
      ENDDO

      DO J=1,1
         DO I=1,II
            DO K=1,KK
               
               H2(K,J,I) = H2(K,J+1,I)
               
            ENDDO
         ENDDO
      ENDDO
      
      DO J=JJ,JJ
         DO I=1,II
            DO K=1,KK
               
               H2(K,J,I) = H2(K,J-1,I)
               
            ENDDO
         ENDDO
      ENDDO
      
      RETURN
      END
C---------------------------------------------------------------------------
      SUBROUTINE DPHIDZ (KK,JJ,II,DZ,DDZ,H1,H2)
      REAL H1(KK,JJ,II),H2(KK,JJ,II)
      REAL DZ(KK),DDZ(KK)

      DO I=1,II
         DO J=1,JJ
            DO K=2,KK-1
               
               H2(K,J,I) = (H1(K+1,J,I) - H1(K-1,J,I))/(DZ(K) + DZ(K-1))
               
            ENDDO
         ENDDO
      ENDDO

      DO K=1,1
         DO J=1,JJ
            DO I=1,II
               
               H2(K,J,I) = H2(K+1,J,I)
               
            ENDDO
         ENDDO
      ENDDO
      
      DO K=KK,KK
         DO J=1,JJ
            DO I=1,II
               
               H2(K,J,I) = H2(K-1,J,I)
               
            ENDDO
         ENDDO
      ENDDO
      
      RETURN
      END
      SUBROUTINE ADDFLT (KK,JJ,II,INXREC,INYREC,INZREC,NMREC,IREC,
     $     IARR,JARR,KARR,
     $     RTREC,H1)
      REAL RTREC(INZREC,INYREC,INXREC)
      INTEGER IARR(NMREC,II),JARR(NMREC,JJ),KARR(NMREC,KK)
      REAL H1(KK,JJ,II)

      DO  ILAUF=1,INXREC
      DO  JLAUF=1,INYREC
      DO  KLAUF=1,INZREC

          K = KARR(IREC,KLAUF)
          J = JARR(IREC,JLAUF)
          I = IARR(IREC,ILAUF)

          RTREC(KLAUF,JLAUF,ILAUF) = RTREC(KLAUF,JLAUF,ILAUF) + 
     $                      0.125*
     $                              ( H1(K  ,J  ,I  ) + H1(K+1,J  ,I  )
     $                             +  H1(K  ,J+1,I  ) + H1(K+1,J+1,I  )
     $                             +  H1(K  ,J  ,I+1) + H1(K+1,J  ,I+1)
     $                             +  H1(K  ,J+1,I+1) + H1(K+1,J+1,I+1))

      ENDDO
      ENDDO
      ENDDO


      RETURN
      END

      SUBROUTINE W6D
     * (KANAL,IAPPEND,ICLOSE,MODIO,CLOOP,CNAME,CTYPE,
     *  COORD1,COORD2,COORD3,COORD4,COORD5,COORD6,C1,C2,C3,
     *  NDIMS,N1,N2,N3,N4,N5,N6,
     *  X1,X2,X3,X4,X5,X6,D1,D2,D3,FNAME)
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                 UNTERPROGRAMM ZUM SCHREIBEN VON
C                 ALLGEMEINEN, BIS ZU 6-DIMENSIONALEN, FELDERN
C     
C           BINAER (MODIO='BIN')          KODIERT     (MODIO='COD')
C           STREAM (CLOOP='STREAM')       LOOP        (CLOOP='LOOP')
C
C           DATENTYP (CNAME='EIGENVALUE' 
C                        OR 'SOLUTION' 
C                        OR 'SNAP-SHOT-MATRIX' ETC.)
C
C           VEKTOR   (CTYPE='VECTOR')     SKALAR      (CTYPE='SCALAR')
C           COORD1...6:   NAME OF COORDINATES, ONLY FIRST NDIM WILL BE
C                         WRITTEN
C           C1,C2,C3 :    NAME OF DATA, IF CTYPE='SCALAR', ONLY
C                         C1 WILL BE WRITTEN
C
C           NDIMS    :    NUMBER OF DIMENSIONS
C           N1,...,N6:    DIMENSIONS
C           X1,...,X6:    COORDINATE-VECTORS, ONLY FIRST NDIM WILL BE
C                         WRITTEN
C           D1,...,D3:    DATA-FIELDS
C                         IF CTYPE='VECTOR', ALL DATA-FIELDS MUST BE
C                         PROVIDED
C                         IF CTYPE='SCALAR', ONLY ONE DATA-FIELD 
C                         REQUIRED
C
C           IAPPEND  :    0 --> NEW FILE IS WRITTEN
C                               WITH HEADER AND GEOMETRY
C                         1 --> DATA WILL BE APPENDED TO EXISTING FILE
C                               FILE MUST BE OPENED AND POSITIONED
C           ICLOSE   :    0 --> FILE REMAINS OPENED,
C                               DATA CAN BE APPENDED
C                         1 --> FILE WILL BE CLOSED,
C                               DATA CANNOT BE APPENDED
C
C        2.6.94 : ORIGINAL
C        4.11.94: 
C       
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

      DIMENSION 
     *         D1(N1,N2,N3,N4,N5,N6)
     *        ,D2(N1,N2,N3,N4,N5,N6)
     *        ,D3(N1,N2,N3,N4,N5,N6)
     *        ,X1(N1)
     *        ,X2(N2)
     *        ,X3(N3)
     *        ,X4(N4)
     *        ,X5(N5)
     *        ,X6(N6)

      CHARACTER (LEN=80) MODIO,CLOOP,CNAME,CTYPE,CTEST,
     *  COORD1,COORD2,COORD3,COORD4,COORD5,COORD6
      CHARACTER (LEN=16) C1,C2,C3
      CHARACTER (LEN=32) FNAME

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                       CHECK OF THE DIMENSIONS

      IF ((NDIMS .EQ. 0) .OR. (NDIMS .GT. 6)) STOP 'NDIMS ERROR IN W6D'
      IF ((N1*N2*N3*N4*N5*N6) .EQ. 0) STOP 'DIMENSIONS ERROR IN W6D'

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                       OPENING OF THE FILE

      IF (IAPPEND .EQ. 0) THEN
         IF (MODIO(1:3) .EQ. 'BIN' ) THEN
            OPEN (KANAL,FILE=FNAME,FORM='UNFORMATTED')
         ELSEIF ((MODIO(1:3) .EQ. 'COD') 
     *      .OR. (MODIO(1:3) .EQ. 'KOD')) THEN
            OPEN (KANAL,FILE=FNAME,FORM='FORMATTED')
         ELSE
            STOP 'MODIO NICHT SPEZIFIZIERT, ERROR IN W6D'
         ENDIF
         REWIND (KANAL)
      ENDIF

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                       WRITING OF THE HEADER OF FILE

C              HEADER WIRD JEDESMAL GESCHRIEBEN (A.O. 909.08.1995)
C     IF (IAPPEND .EQ. 0) THEN
         IF (MODIO(1:3) .EQ. 'BIN' ) THEN
            CTEST = '6D-FILE-BIN'
            WRITE (KANAL) CTEST
            WRITE (KANAL) MODIO
            WRITE (KANAL) CLOOP
            WRITE (KANAL) CNAME
         ELSEIF ((MODIO(1:3) .EQ. 'COD')
     *      .OR. (MODIO(1:3) .EQ. 'KOD')) THEN
            WRITE (KANAL,'(A)') '6D-FILE-COD'
            WRITE (KANAL,'(A)') MODIO
            WRITE (KANAL,'(A)') CLOOP
            WRITE (KANAL,'(A)') CNAME
         ENDIF
C     ENDIF

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                       WRITING OF THE HEADER OF DATA

         IF (MODIO(1:3) .EQ. 'BIN' ) THEN
            CTEST = 'NUMBER OF DIMENSIONS'
            WRITE (KANAL) CTEST
            WRITE (KANAL) NDIMS
            CTEST = 'DIMENSIONS'
            WRITE (KANAL) CTEST
            WRITE (KANAL) N1,N2,N3,N4,N5,N6
         ELSEIF ((MODIO(1:3) .EQ. 'COD')
     *      .OR. (MODIO(1:3) .EQ. 'KOD')) THEN
            WRITE (KANAL,'(A)') 'NUMBER OF DIMENSIONS'
            WRITE (KANAL,'(I6)') NDIMS
            WRITE (KANAL,'(A)') 'DIMENSIONS'
            WRITE (KANAL,'(6(1X,I6))') N1,N2,N3,N4,N5,N6
         ENDIF


CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                             WRITING OF COORDINATES

         IF (NDIMS .GE. 1) THEN
            CALL WCOORD (KANAL,MODIO,COORD1,N1,X1)
         ENDIF
         IF (NDIMS .GE. 2) THEN
            CALL WCOORD (KANAL,MODIO,COORD2,N2,X2)
         ENDIF
         IF (NDIMS .GE. 3) THEN
            CALL WCOORD (KANAL,MODIO,COORD3,N3,X3)
         ENDIF
         IF (NDIMS .GE. 4) THEN
            CALL WCOORD (KANAL,MODIO,COORD4,N4,X4)
         ENDIF
         IF (NDIMS .GE. 5) THEN
            CALL WCOORD (KANAL,MODIO,COORD5,N5,X5)
         ENDIF
         IF (NDIMS .GE. 6) THEN
            CALL WCOORD (KANAL,MODIO,COORD6,N6,X6)
         ENDIF

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                             WRITING OF DATA

      IF (CTYPE(1:3) .EQ. 'SKA' .OR. CTYPE(1:3) .EQ. 'SCA') THEN
         CALL WSCALAR (KANAL,MODIO,CLOOP,C1,N1,N2,N3,N4,N5,N6,D1)
      ELSEIF (CTYPE(1:3) .EQ. 'VEK' .OR. CTYPE(1:3) .EQ. 'VEC') THEN
         CALL WVECTOR (KANAL,MODIO,CLOOP,
     *                 C1,C2,C3,N1,N2,N3,N4,N5,N6,D1,D2,D3)
      ELSE
         STOP 'CTYPE ERROR IN W6D'
      ENDIF

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

      RETURN
      END

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

      SUBROUTINE WCOORD(KANAL,MODIO,COORD,N,X)
      DIMENSION X(N)
      CHARACTER (LEN=80) MODIO,COORD,CTEST

      IF (MODIO(1:3) .EQ. 'BIN' ) THEN
         CTEST = 'COORDINATE'
         WRITE (KANAL) CTEST
         WRITE (KANAL) COORD
         WRITE(KANAL) (X(I),I=1,N)
      ELSE
         WRITE (KANAL,'(A)') 'COORDINATE'
         WRITE (KANAL,'(A)') COORD
         WRITE(KANAL,*) (X(I),I=1,N)
      ENDIF

      RETURN
      END
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

      SUBROUTINE WSCALAR (KANAL,MODIO,CLOOP,C1,N1,N2,N3,N4,N5,N6,D1)

      DIMENSION D1(N1,N2,N3,N4,N5,N6)
      CHARACTER (LEN=80) MODIO,CLOOP,CTEST
      CHARACTER (LEN=16) C1, C2, C3

      IF (MODIO(1:3) .EQ. 'BIN' ) THEN
         CTEST = 'SCALAR'
         WRITE (KANAL) CTEST
         WRITE (KANAL) C1

         IF (CLOOP(1:6) .EQ. 'STREAM') THEN

            WRITE (KANAL) D1

         ELSEIF (CLOOP(1:4) .EQ. 'LOOP') THEN

            DO I6=1,N6
            DO I5=1,N5
            DO I4=1,N4
            DO I3=1,N3
            DO I2=1,N2
            DO I1=1,N1

               WRITE (KANAL) D1(I1,I2,I3,I4,I5,I6)

            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO

         ELSE

            STOP 'CLOOP ERROR IN WSCALAR'

         ENDIF

      ELSEIF ((MODIO(1:3) .EQ. 'COD')
     *   .OR. (MODIO(1:3) .EQ. 'KOD')) THEN

         WRITE (KANAL,'(A)') 'SCALAR'
         WRITE (KANAL,'(A)') C1

         IF (CLOOP(1:6) .EQ. 'STREAM') THEN

            WRITE (KANAL,*) D1

         ELSEIF (CLOOP(1:4) .EQ. 'LOOP') THEN

            DO I6=1,N6
            DO I5=1,N5
            DO I4=1,N4
            DO I3=1,N3
            DO I2=1,N2
            DO I1=1,N1

               WRITE (KANAL,*) D1(I1,I2,I3,I4,I5,I6)

            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO

         ELSE
  
            STOP 'CLOOP ERROR IN WSCALAR'

         ENDIF

      ELSE

         STOP 'MODIO ERROR IN WSCALAR'

      ENDIF

      RETURN
      END

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

      SUBROUTINE WVECTOR (KANAL,MODIO,CLOOP,
     *                    C1,C2,C3,N1,N2,N3,N4,N5,N6,D1,D2,D3)

      DIMENSION D1(N1,N2,N3,N4,N5,N6)
      DIMENSION D2(N1,N2,N3,N4,N5,N6)
      DIMENSION D3(N1,N2,N3,N4,N5,N6)

      CHARACTER (LEN=80) MODIO,CLOOP,CTEST
      CHARACTER (LEN=16) C1,C2,C3

      IF (MODIO(1:3) .EQ. 'BIN' ) THEN
         CTEST = 'VECTOR'
         WRITE (KANAL) CTEST

         IF (CLOOP(1:6) .EQ. 'STREAM') THEN

            WRITE (KANAL) C1
            WRITE (KANAL) D1

            WRITE (KANAL) C2
            WRITE (KANAL) D2

            WRITE (KANAL) C3
            WRITE (KANAL) D3

         ELSEIF (CLOOP(1:4) .EQ. 'LOOP') THEN

               WRITE (KANAL) C1,C2,C3

            DO I6=1,N6
            DO I5=1,N5
            DO I4=1,N4
            DO I3=1,N3
            DO I2=1,N2
            DO I1=1,N1

               WRITE (KANAL) D1(I1,I2,I3,I4,I5,I6), 
     &                       D2(I1,I2,I3,I4,I5,I6),
     &                       D3(I1,I2,I3,I4,I5,I6)

            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO

         ELSE

            STOP 'CLOOP ERROR IN WVECTOR'

         ENDIF

      ELSEIF ((MODIO(1:3) .EQ. 'COD')
     *   .OR. (MODIO(1:3) .EQ. 'KOD')) THEN

         WRITE (KANAL,'(A)') 'VECTOR'

         IF (CLOOP(1:6) .EQ. 'STREAM') THEN

            WRITE (KANAL,'(A)') C1
            WRITE (KANAL,*) D1

            WRITE (KANAL,'(A)') C2
            WRITE (KANAL,*) D2

            WRITE (KANAL,'(A)') C3
            WRITE (KANAL,*) D3

         ELSEIF (CLOOP(1:4) .EQ. 'LOOP') THEN

               WRITE (KANAL,'(3A)') C1,C2,C3

            DO I6=1,N6
            DO I5=1,N5
            DO I4=1,N4
            DO I3=1,N3
            DO I2=1,N2
            DO I1=1,N1

               WRITE (KANAL,*) D1(I1,I2,I3,I4,I5,I6), 
     &                         D2(I1,I2,I3,I4,I5,I6),
     &                         D3(I1,I2,I3,I4,I5,I6)

            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            ENDDO

         ELSE
  
            STOP 'CLOOP ERROR IN WVECTOR'

         ENDIF

      ELSE

         STOP 'MODIO ERROR IN WVECTOR'

      ENDIF

      RETURN
      END

C------ Subroutine for writing scalar in w6d-format without anything else
      SUBROUTINE W6DSCALAR(KANREC,NMREC,IREC,INXREC,INYREC,INZREC,
     $ IARR,JARR,KARR,XREC,YREC,ZREC,RTREC,
     $ PHI,
     $ KK,JJ,II,KMX,JMX,IMX,
     $ DX,DY,DZ,DDX,DDY,DDZ,TIMEPH)

C IMPLICIT NONE

      INTEGER KANREC,IREC,INXREC,INYREC,INZREC,
     $ II,JJ,KK,IMX,JMX,KMX,NMREC,K,J,I,
     $ ILAUF,JLAUF,KLAUF

      INTEGER IARR(NMREC,II),JARR(NMREC,JJ),KARR(NMREC,KK)

        REAL RTREC(INZREC,INYREC,INXREC)
        REAL XREC(NMREC,II),XXREC(II)
     + ,YREC(NMREC,JJ),YYREC(JJ)
     + ,ZREC(NMREC,KK),ZZREC(KK)

      REAL TIMEPH
      REAL PHI(KK,JJ,II)
      REAL DX(II),DY(JJ),DZ(KK),DDX(II),DDY(JJ),DDZ(KK)

      CHARACTER (LEN=16) CIDREC(NMREC)
      CHARACTER (LEN=18) FNAME

      WRITE(44,*)'RAUSSCHREIBEN EINES RECORDS AUF KANAL',KANREC

      DO I=1,INXREC
        XXREC(I) = XREC(IREC,I)
      ENDDO

      DO J=1,INYREC
        YYREC(J) = YREC(IREC,J)
      ENDDO

      DO K=1,INZREC
        ZZREC(K) = ZREC(IREC,K)
      ENDDO
C      write(6,*)'HIER bin ich',INXREC,INYREC,INZREC,KARR(IREC,INZREC),
C     $ JARR(IREC,INYREC),IARR(IREC,INXREC),KARR(IREC,1),
C     $ JARR(IREC,1),IARR(IREC,1),KK,JJ,II,
C     $ SIZE(PHI),SIZE(RTREC)

      DO ILAUF=1,INXREC
       DO JLAUF=1,INYREC
        DO KLAUF=1,INZREC
C
       I = IARR(IREC,ILAUF)
       J = JARR(IREC,JLAUF)
       K = KARR(IREC,KLAUF)

       RTREC(KLAUF,JLAUF,ILAUF) =
     + PHI(K,J,I)
C          WRITE(6,*)KLAUF,JLAUF,ILAUF,K,J,I
        ENDDO
       ENDDO
      ENDDO

        write(FNAME,'("t_",I3.3,"_",I4.4,F4.3,"_BIN")')
     $ KANREC,int(TIMEPH),TIMEPH-int(TIMEPH)

C       write(6,*)'HIER: ',FNAME


        call W6D
     $ (KANAL,0,1,'BIN','STREAM','T','SCALAR',
     $ 'ZZREC','YYREC','XXREC','TIMEPH','TIMEPH','TIMEPH',
     $ 'RRTREC','RRTREC','RRTREC',
     $ 3,INZREC,INYREC,INXREC,1,1,1,
     $ ZZREC,YYREC,XXREC,TIMEPH,TIMEPH,TIMEPH,
     $ RTREC,RTREC,RTREC,FNAME)

        RETURN
        END

