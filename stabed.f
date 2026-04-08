










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
      SUBROUTINE STABED  (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,U,V,W,G,B,     
     $                    RHO,DT,XTOT,YTOT,ZTOT,IGRID)
C*STARLET***************************************************************
C        S T A B E D      ERMITTLUNG DES MAXIMALEN ZEITSCHRITTES NACH
C                         EINER VON NEUMANN''SCHEN STABILITAETSANALYSE
C                         (LITERATUR SIEHE SUBR. 'STABHI'), DER MAX.
C                         COURANT-ZAHL IM FELD, DER MAX. UND MIN. ZELL-
C                         PECLET-ZAHL, SOWIE DER MAX. DIVERGENZ.
C                         DER ORT DES JEWEILIGEN EXTREMUMS WIRD DURCH
C                         DIE ZUGEHOERIGEN INDIZES I,J,K ANGEGEBEN.
C        BEACHTE:         DIE VON NEUMANN''SCHE STABILITAETSANALYSE IST
C                         NUR KORREKT, WENN DIE ZU UNTERSUCHENDE TRANS-
C                         PORTGLEICHUNG PERIODISCHE RANDBEDINGUNGEN HAT!
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        DDX,DDY,DDZ    - KANTENLAENGEN DER KONTROLLVOLUMINA IN X-, Y-
C                         UND Z-RICHTUNG
C        U,V,W          - GESCHWINDIGKEITSFELDER
C        G(KK,JJ,II)    - EFFEKTIVE DYNAMISCHE VISKOSITAET
C        RHO            - DICHTE DES FLUIDS (= CONST)
C        DT             - ZEITSCHRITT
C        XTOT,YTOT,ZTOT - GESAMTABMESSUNG DES BERECHNUNGSGEBIETES IN
C                         X-, Y- UND Z-RICHTUNG
C
C UPROG                 : STABHI, PAGE6, ERRR
C
C DEFINE-DIREKTIVEN     : UPW, ZEN, QUD, LEAPF
C
C        02.08.85 (HW)  : ORIGINAL
C        10.11.87 (HW)  : STABED LAEUFT AUCH MIT EINER SCHEIBE
C                         IM BERECHNUNGSGEBIET (Z.B. JMX = 5)
C
C        13.05.02 (TB)  : DIV OUTPUT FOR INDICATED IJK DOMAIN INTEGRATED
C
C*STARLET***************************************************************
C
      REAL       U(KK,JJ,II),   V(KK,JJ,II),   W(KK,JJ,II),
     $           G(KK,JJ,II),
     $           B(KK,JJ,II)
      REAL       DDX(II),       DDY(JJ),       DDZ(KK)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      WRITE (6,*) '  '
      WRITE (6,*) '  '
      WRITE (6,*) ' ************************************************'
      WRITE (6,*) '  '
      WRITE (6,*) '       STABILITAETSANALYSE FUER GITTER:',IGRID
      WRITE (6,*) '  '

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IM2 = IMX-2
      JM2 = JMX-2
      KM2 = KMX-2
C
C                                 STABILITAETSUNTERSUCHUNG FUER DIE X-,
C                                 Y- UND Z-RICHTUNG ERFOLGT NACHEINANDER
C
C      IF( (IMX-5) .LE. 0) GOTO 2001
      CALL STABHI  (KK,JJ,II,KMX,JMX,IMX,'U',DDX,II,U,G,B,RHO,
     $              DT,XTOT,YTOT,ZTOT,COUMAU,ICMAU,JCMAU,KCMAU,DTMAXU,
     $              IDTU,JDTU,KDTU,PECMAU,IPMAU,JPMAU,KPMAU,PECMIU,
     $              IPMIU,JPMIU,KPMIU)
C 2001 IF( (JMX-5) .LE. 0) GOTO 2002
      CALL STABHI  (KK,JJ,II,KMX,JMX,IMX,'V',DDY,JJ,V,G,B,RHO,
     $              DT,XTOT,YTOT,ZTOT,COUMAV,ICMAV,JCMAV,KCMAV,DTMAXV,
     $              IDTV,JDTV,KDTV,PECMAV,IPMAV,JPMAV,KPMAV,PECMIV,
     $              IPMIV,JPMIV,KPMIV)
C 2002 IF( (KMX-5) .LE. 0) GOTO 2003
      CALL STABHI  (KK,JJ,II,KMX,JMX,IMX,'W',DDZ,KK,W,G,B,RHO,
     $              DT,XTOT,YTOT,ZTOT,COUMAW,ICMAW,JCMAW,KCMAW,DTMAXW,
     $              IDTW,JDTW,KDTW,PECMAW,IPMAW,JPMAW,KPMAW,PECMIW,
     $              IPMIW,JPMIW,KPMIW)
C
 2003 CONTINUE
C
C                                 ERMITTLUNG DER MAX. DIVERGENZ
C
CA13  TB130502 CCCCC WRITE DIVG FOR ALL INDICATED IJK CCCCCCCCCCCCCCCCCC
      GRID = 0
      IANF = 0
      JANF = 0
      KANF = 0
      IEND = 0
      JEND = 0
      KEND = 0
      OPEN(80,ERR=2005,FILE='DIVIJK.dat',STATUS='OLD')
      READ (80,*,ERR=2005) GRID,IANF,IEND,JANF,JEND,KANF,KEND
      CLOSE(80,ERR=2005)
      IF (IGRID .EQ. GRID) THEN
       IF (IEND .GT. 0 .OR. JEND .GT. 0 .OR. KEND .GT. 0) THEN
       OPEN(81,ERR=2005,FILE='DIVIJK.div',STATUS='NEW',FORM='FORMATTED')
       WRITE (81,*) 'GRID: ',IGRID    	
       ENDIF
      ENDIF
 2005 CONTINUE     
CENDE TB130502 CCCCC WRITE DIVG FOR ALL INDICATED IJK CCCCCCCCCCCCCCCCCC
      DIVGMX = 0.0
C
      DO 10 I=3,IM2
         DO 20 J=3,JM2
            DO 30 K=3,KM2

            IF (B(K,J,I).GT.0.0) THEN

               DIVG = ABS((U(K,J,I) - U(K,J,I-1))/DDX(I)
     $              +     (V(K,J,I) - V(K,J-1,I))/DDY(J)
     $              +     (W(K,J,I) - W(K-1,J,I))/DDZ(K))
CA6   TB130502 CCCCC WRITE DIVG FOR ALL INDICATED IJK CCCCCCCCCCCCCCCCCC
	      IF (IGRID .EQ. GRID .AND. I .GE. IANF .AND. J .GE. JANF 
     $           .AND. K .GE. KANF .AND. I .LE. IEND
     $           .AND. J .LE. JEND .AND. K .LE. KEND) THEN
      	        WRITE (81,6001) DIVG,I,J,K    	
     	      ENDIF
CENDE TB130502 CCCCC WRITE DIVG FOR ALL INDICATED IJK CCCCCCCCCCCCCCCCCC
    
               IF(DIVG .LT. DIVGMX) GOTO 2010
                  DIVGMX = DIVG
                  IDIV   = I
                  JDIV   = J
                  KDIV   = K
 2010          CONTINUE

            ENDIF

   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
CA6   TB130502 CCCCC WRITE DIVG FOR ALL INDICATED IJK CCCCCCCCCCCCCCCCCC
      IF (IGRID .EQ. GRID) THEN
       IF (IEND .GT. 0 .OR. JEND .GT. 0 .OR. KEND .GT. 0) THEN
        CLOSE(81,ERR=2011)
       ENDIF
      ENDIF
 2011 CONTINUE     
CENDE TB130502 CCCCC WRITE DIVG FOR ALL INDICATED IJK CCCCCCCCCCCCCCCCCC
C
 2200 CONTINUE
C
C                                 AUSGABE DER ERGEBNISSE
C
      CALL PAGE6
C
CA2   TB130502 CCCCC WRITE DIVG FOR ALL INDICATED IJK CCCCCCCCCCCCCCCCCC
 6001 FORMAT(1X,'DIV :  ',1PE11.4,'  I=',I3,'  J=',I3,
     $       '  K=',I3) 
CENDE TB130502 CCCCC WRITE DIVG FOR ALL INDICATED IJK CCCCCCCCCCCCCCCCCC     
      WRITE (6,6010)
 6010 FORMAT(/,1X,52(1H-),'  STABILITAETSUNTERSUCHUNG  ',52(1H-),/,
     $       28X,'X-RICHTUNG  (U-KOMP.)',16X,'Y-RICHTUNG  (V-KOMP.)',
     $       16X,'Z-RICHTUNG  (W-KOMP.)',/,1X,'MAX. ZEITSCHRITT :')
      WRITE (6,6020) DTMAXU,IDTU,JDTU,KDTU,DTMAXV,IDTV,JDTV,KDTV,
     $               DTMAXW,IDTW,JDTW,KDTW
 6020 FORMAT(1H+,20X,3(1PE11.4,'  I=',I3,'  J=',I3,'  K=',I3,5X))
      WRITE (6,6030)
 6030 FORMAT(/,1X,'MAX. COURANT-ZAHL:')
      WRITE (6,6020) COUMAU,ICMAU,JCMAU,KCMAU,COUMAV,ICMAV,JCMAV,KCMAV,
     $               COUMAW,ICMAW,JCMAW,KCMAW
      WRITE (6,6040)
 6040 FORMAT(/,1X,'MAX. ZELL-PECLET :')
      WRITE (6,6020) PECMAU,IPMAU,JPMAU,KPMAU,PECMAV,IPMAV,JPMAV,KPMAV,
     $               PECMAW,IPMAW,JPMAW,KPMAW
      WRITE (6,6050)
 6050 FORMAT(/,1X,'MIN. ZELL-PECLET :')
      WRITE (6,6020) PECMIU,IPMIU,JPMIU,KPMIU,PECMIV,IPMIV,JPMIV,KPMIV,
     $               PECMIW,IPMIW,JPMIW,KPMIW
      WRITE (6,6060) IGRID,DIVGMX,IDIV,JDIV,KDIV
 6060 FORMAT(/,I3,1X,'MAX. DIVERGENZ   :  ',1PE11.4,'  I=',I3,'  J=',I3,
     $       '  K=',I3,/,1X,130(1H-))
C
      RETURN
      END
      SUBROUTINE STABHI  (KK,JJ,II,KMX,JMX,IMX,IC,DDS,LL,A,G,B,RHO,
     $                    DT,XTOT,YTOT,ZTOT,COURMA,ICMA,JCMA,KCMA,DTMAX,
     $                    IDT,JDT,KDT,PECMA,IPMA,JPMA,KPMA,PECMI,IPMI,
     $                    JPMI,KPMI)
C*STARLET***************************************************************
C        S T A B H I      HILFSROUTINE ZUR ERMITTLUNG DES MAXIMALEN
C                         ZEITSCHRITTES, SOWIE WEITERER CHARAKTERIST.
C                         GROESSEN DES STROEMUNGSFELDES (MAX. COURANT-
C                         ZAHL, MAX. UND MIN. ZELL-PECLETZAHL)
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        IC             - HOLLERITH-KONSTANTE ('U', 'V', 'W')
C        DDS(LL)        - KANTENLAENGE DER KONTROLLVOLUMINA IN DER
C                         KOORD.-RICHTUNG, DIE DER VARIABLEN "IC" ZUGE-
C                         ORDNET IST
C        A(KK,JJ,II)    - GESCHWINDIGKEITSFELD (U-, V-, ODER W-FELD)
C                         A  MUSS MIT "IC" UEBEREINSTIMMEN !
C        G(KK,JJ,II)    - EFFEKTIVE DYNAMISCHE VISKOSITAET
C        RHO            - DICHTE DES FLUIDS (= CONST)
C        DT             - ZEITSCHRITT
C        XTOT,YTOT,ZTOT - GESAMTABMESSUNG DES BERECHNUNGSGEBIETES IN
C                         X-, Y- UND Z-RICHTUNG
C        COURMA         + MAXIMALE COURANT-ZAHL IM FELD (FUER "IC"!)
C        ICMA,JCMA,KCMA + INDIZES DES ORTES DER MAX. COURANT-ZAHL
C        DTMAX          + MAXIMAL MOEGLICHER ZEITSCHRITT, ERMITTELT
C                         NACH EINER VON NEUMANN''SCHEN STABILITAETS-
C                         ANALYSE (GUELTIG FUER UPWIND-, ZENTRALE- UND
C                         QUICK-DIFFERENZEN)
C                         LITERATUR: PAOLUCCI, S.; CHENOWETH, D.R.:
C                                   ""STABILITY OF THE EXPLICIT FINITE
C                                     DIFFERENCED TRANSPORT EQUATION"",
C                                    J. OF COMP. PHYSICS 47, 489-496
C                                    (1982)
C        IDT,JDT,KDT    + INDIZES DER ZELLE, DIE DIE GROESSTE RESTRIK-
C                         TION AUF DEN ZEITSCHRITT AUSUEBT
C        PECMA          + MAXIMALE ZELL-PECLET-ZAHL IM FELD
C                         PECMA = MAX(RHO*U*DDX/GEFF)
C        IPMA,JPMA,KPMA + INDIZES DER ZELLE MIT DER GROESSTEN PECLET-Z.
C        PECMI          + MINIMALE ZELL-PECLET-ZAHL IM FELD
C        IPMI,JPMI,KPMI + INDIZES DER ZELLE MIT DER KLEINSTEN PECLET-Z.
C
C        01.08.85 (HW)  : ORIGINAL
C
C*STARLET***************************************************************
C
      CHARACTER (LEN=1)  IC
C
      REAL       A(KK,JJ,II),   
     $           G(KK,JJ,II),
     $           B(KK,JJ,II),
     $           DDS(LL)
C
      CST1(PHI1) = (1.0-COS(PHI1)) - CONJ/12.0*PECLS*(4.0*COS(PHI1)
     $                             - COS(2.0*PHI1) - 3.0)
C
      CSTAB(PHI) = (CST1(PHI) / (CST1(PHI)**2 + (0.5*PECLS*(CONJ/6.0
     $           * (2.0*SIN(PHI) - SIN(2.0*PHI)) + SIN(PHI)))**2))
     $           / ALPHS*DIST**2
C
      IM2    = IMX-2
      JM2    = JMX-2
      KM2    = KMX-2
C
C                                 KONSTANTEN FUER DIE STABILITAETS-
C                                 ANALYSE
C
      CONI   = -1.1E+10
C
      CONI   = 1.0
      CONJ   = 0.0
C
C                                 VORBELEGUNG
C
      COURMA = -1.0E+10
      DTMAX  =  1.0E+10
      PECMA  = -1.0E+10
      PECMI  =  1.0E+10
C
      PI     =  ACOS(-1.0)
C
      IDT  = -999
      JDT  = IDT
      KDT  = IDT
      IPMA = IDT
      JPMA = IDT
      KPMA = IDT
      IPMI = IDT
      JPMI = IDT
      KPMI = IDT
      ICMA = IDT
      JCMA = IDT
      KCMA = IDT
C
C                                 MULTIPLIKATOREN ZUR KORREKTEN BE-
C                                 STIMMUNG DER MASCHENWEITE
C
      IF(IC .NE. 'U') GOTO 2010
         IMULT = 1
         JMULT = 0
         KMULT = 0
         GOTO 2030
 2010 IF(IC .NE. 'V') GOTO 2020
         IMULT = 0
         JMULT = 1
         KMULT = 0
         GOTO 2030
 2020 IF(IC .NE. 'W') CALL ERRR(501,' STABHI   ')
         IMULT = 0
         JMULT = 0
         KMULT = 1
C
 2030 CONTINUE
C
      DO 10 I=3,IM2
         DO 20 J=3,JM2
            DO 30 K=3,KM2

            IF (B(K,J,I).GT.0.0) THEN

               DIST   = DDS(I*IMULT + J*JMULT + K*KMULT)
               UQUER  = ABS(0.5*(A(K,J,I)+A(K-KMULT,J-JMULT,I-IMULT)))
               COUR   = UQUER*DT/DIST
               PECL   = RHO*UQUER*DIST/G(K,J,I)
               IF(CONI .LT. -1.0E+10) GOTO 2110
                  FAKTS  = 1.0 + 0.5*(1.0-CONI)*PECL
                  ALPHS  = G(K,J,I)/RHO*FAKTS
                  PECLS  = PECL/FAKTS
                  PHIST  = 2.0*DIST*PI/(XTOT*FLOAT(IMULT) + YTOT
     $                   * FLOAT(JMULT) + ZTOT*FLOAT(KMULT))
                  DT1    = CSTAB(PHIST)
                  DT2    = CSTAB(PI)
                  IF((DT1 .GT. DTMAX) .AND. (DT2 .GT. DTMAX)) GOTO 2110
                     DTMAX  = AMIN1(DT1,DT2)
                     IDT    = I
                     JDT    = J
                     KDT    = K
 2110          IF(COUR .LT. COURMA)                           GOTO 2120
                  COURMA = COUR
                  ICMA   = I
                  JCMA   = J
                  KCMA   = K
 2120          IF(PECL .LT. PECMA )                           GOTO 2130
                  PECMA  = PECL
                  IPMA   = I
                  JPMA   = J
                  KPMA   = K
 2130          IF(PECL .GE. PECMI )                           GOTO 2140
                  PECMI  = PECL
                  IPMI   = I
                  JPMI   = J
                  KPMI   = K
 2140          CONTINUE

            ENDIF

   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
      RETURN
      END


