










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
      SUBROUTINE TSTLE2  (KK,JJ,II,KMX,JMX,IMX,
     $                    X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                    RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                    U,V,W,UO,VO,WO,
     $                    UP,VP,WP,WCU,WCV,WCW,P,G,
     $                    B,GRADP,DT,
     $                    GRADPX,ZTOT,NBUF,ITSTEP,WPHI,WKON,WDIF,
     $                    WSOR,IDUZ,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,
     $                   RSGS3,RSGS2,FAKTOR,
     $                   UINI,VINI,WINI,FDUI,FDVI,FDWI,
     $                   UINJ,VINJ,WINJ,FDUJ,FDVJ,FDWJ,
     $                   UINK,VINK,WINK,FDUK,FDVK,FDWK,
     $                   COEFFX,COEFFY,COEFFZ,COEFDX,COEFDY,COEFDZ,
     $                   LCL,DIL,RCL,UZ,RSP,LCR,DIR,RCR,LCOL,DIAG,RCOL,
     $                   FAKTOR2
     $                    )
C*STARLET***************************************************************
C        T S T L E 2      TIMESTEP - TURBULENT - ALLE VERFAHREN - LES
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        DX,DY,DZ       - ABSTAND DER GITTERPUNKTE
C        DDX,DDY,DDZ    - KANTENLAENGEN DER KONTROLLVOLUMINA
C        U(KK,JJ,II)    - NEUES GESCHWINDIGKEITSFELD
C        V(KK,JJ,II)    - NEUES GESCHWINDIGKEITSFELD
C        W(KK,JJ,II)    - NEUES GESCHWINDIGKEITSFELD
C        UO(KK,JJ,II)   + ALTES GESCHWINDIGKEITSFELD (FELD WIRD F.
C                         IDUZ = 1 UND IDUZ = 2 VERAENDERT !!)
C        VO(KK,JJ,II)   + ALTES GESCHWINDIGKEITSFELD (  "    "  "    "
C        WO(KK,JJ,II)   + ALTES GESCHWINDIGKEITSFELD (  "    "  "    "
C        UP(KK,JJ,NBUF) + PUFFERARRAYS
C        VP(KK,JJ,NBUF) + PUFFERARRAYS
C        WP(KK,JJ,NBUF) + PUFFERARRAYS
C        P(KK,JJ,II)    - ALTES DRUCKFELD
C        G(KK,JJ,II)    - EFFEKTIVE DYNAMISCHE VISKOSITAET (= MUE)
C        GMOL           - MOLEKULARE DYNAMISCHE VISKOSITAET
C        DT             - ZEITSCHRITT
C        RHO            - DICHTE (= CONST)
C        UGRID          - GESCHWINDIGKEIT DES KOORDINATENSYSTEMS IN
C                         X-RICHTUNG (GALILEI-TRANSFORMATION)
C        WCU(KK,JJ)     + WANDKORREKTUR F. U-KOMP. (DIFFUSIVER IMPULSVER
C                                                   LUST INFOLGE REIBUNG
C        WCV(KK,JJ)     + WANDKORREKTUR F. V-KOMP.           ""
C        WCW(KK,JJ)     + WANDKORREKTUR F. W-KOMP.           ""
C        GRADPX         - DRUCKGRADIENT IN X-RICHTUNG
C        ZTOT           - SENKRECHTER ABSTAND DER KANALWAENDE
C        MTURB          - SCHALTER: 0 = LAMINAR,  1 = TURBULENT
C        NBUF           - ANZAHL DER PUFFERSCHICHTEN (F. QUICK: NBUF=3)
C        ITSTEP         - ZEITSCHRITTZAEHLER
C        WPHI           - WICHTUNGSFAKTOR FUER DEN PUNKTWERT
C        WKON           -         "        "   DIE KONVEKTIVEN TERME
C        WDIF           -         "        "   DIE DIFFUSIVEN TERME
C        WSOR           -         "        "   DEN QUELLTERM
C        IDUZ           - FUER DEN ZEITSCHRITT 2. ORDNUNG MUSS DER ZEIT-
C                         SCHRITT ZWEIMAL DURCHLAUFEN WERDEN. FUER IDUZ
C                         = 1 WERDEN GESCHW.KOMP VOM AELTESTEN ZEIT-
C                         NIVEAU UEBERGEBEN; FUER IDUZ = 2 DIE DES
C                         AKTUELLEREN ZEITNIVEAUS
C        IC1,JC1,KC1    - LINKE  RAENDER DER BOUNDING BOX
C        IC2,JC2,KC2    - RECHTE RAENDER DER BOUNDING BOX
C
C     UINI,VINI,WINI    - INTERPOLIERTE GESCHW. IN I-RICHTUNG
C     UINJ,VINJ,WINJ    - INTERPOLIERTE GESCHW. IN J-RICHTUNG
C     UINK,VINK,WINK    - INTERPOLIERTE GESCHW. IN K-RICHTUNG
C
C     FDUI,FDVI,FDWI    - ERSTE ABLEITUNGEN     IN I-RICHTUNG           
C     FDUJ,FDVJ,FDWJ    - ERSTE ABLEITUNGEN     IN J-RICHTUNG           
C     FDUK,FDVK,FDWK    - ERSTE ABLEITUNGEN     IN K-RICHTUNG           
C
C                  X    - X-KOORDINATE DES ZELLMITTELPUNKTES            
C
C DEFINE DIREKTIVEN     : ZEN, UPW, QUD, EULER, LEAPF, ADBA,
C                         FRPER, RIPER
C UPROG                 : SWCLE1
C
C        23.08.85 (HW)  : TSTLE2 : ZEITSCHRITT WAHLWEISE MITTELS EULER-,
C                                  LEAPFROG- ODER ADAMS/BASHFORTH-VER-
C                                  FAHREN
C                         U-SCHLEIFE LAEUFT BIS I=IM2 EINSCHLIESSLICH !
C                         V-   "       "     "  J=JM2        "        !
C                         DER DRUCK  P  IST ALS EFFEKTIVER DRUCK ZU BE-
C                         TRACHTEN.  P = P + 1/3UI*UI
C                         IN X-RI. WIRD DER DRUCKGRADIENT EXPLIZIT VOR-
C                         GEGEBEN.
C        26.03.86 (HW)  : LAUFLAENGE DER SCHLEIFEN DURCH DEFINE-
C                         DIREKTIVEN GESTEUERT
C        02.04.92 (MM)  : SCHLEIFEN LAUFEN UEBER GANZES GEBIET
C         4. 4.92 (MM)  : SCHLEIFEN ZUR WANDKORREKTUR AM KOERPER
C                         LAUFEN NUR INNERHALB DER
C                         BOUNDING BOX (IC1,IC2,JC1,JC2,KC1,KC2)
C
C        10.10.96 (AM)  : KOMPAKTVERFAHREN VIERTER ORDNUNG    
C                         BERECHNUNG DER INTERPOLATIONEN UND ABLEITUNGEN
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

      PARAMETER (NPPHYS_MAX=1)
      COMMON /COPHYSPAR/
     $                  MTURB,  TU_LEVEL,   RHO,    GMOL,   UGRID,
     $                  VREF,   EXPON,  UTAUX,
     $                  CIDUFR, UFRCON, UFRFREQ, DELTA,  XREF,
     $                  IPP,    JPP,    KPP,
     $                  XPER,   YPER,   ZPER, 
     $                  NXPER,  NYPER,  NZPER,
     $                  NPPHYS, PPHYS, XPPHYS

      INTEGER           IPP,  JPP,  KPP,   NPPHYS
     
      REAL              TU_LEVEL,  RHO,   GMOL,   UGRID,  VREF,  
     $                  EXPON, UFRCON, 
     $                  PPHYS(NPPHYS_MAX), XPPHYS(NPPHYS_MAX)

      CHARACTER (LEN=16)      CIDUFR

C
      REAL        X(II),         Y(JJ),         Z(KK),
     &           DX(II),        DY(JJ),        DZ(KK),
     &          DDX(II),       DDY(JJ),       DDZ(KK),
     $        RDX(II),RDY(JJ),RDZ(KK),RDDX(II),RDDY(JJ),RDDZ(KK)
C
      REAL        U(KK,JJ,II),   V(KK,JJ,II),   W(KK,JJ,II),
     $           UO(KK,JJ,II),  VO(KK,JJ,II),  WO(KK,JJ,II),
     &           UP(KK,JJ,NBUF),VP(KK,JJ,NBUF),WP(KK,JJ,NBUF)
C
      REAL        P(KK,JJ,II),   G(KK,JJ,II)
      REAL        B(KK,JJ,II) 

      REAL        GRADP(II)


C
      REAL      WCU(KK,JJ),    WCV(KK,JJ),    WCW(KK,JJ)
C 
      REAL      RSGS3(KK,JJ,II), RSGS2(KK,JJ)
C
      REAL      COEFFX(II,12),COEFFY(JJ,12),COEFFZ(KK,12),
     $          COEFDX(II,12),COEFDY(JJ,12),COEFDZ(KK,12)
C
      REAL      LCOL(II),DIAG(II),RCOL(II),UZ(II),RSP(II),
     $          FAKTOR(KK,JJ,II),FAKTOR2(II)
C
      REAL      LCL(KK,JJ,II),DIL(KK,JJ,II),RCL(KK,JJ,II)
      REAL      LCR(KK,JJ,II),DIR(KK,JJ,II),RCR(KK,JJ,II)

C
      REAL      UINI(KK,JJ,II),VINI(KK,JJ,II),WINI(KK,JJ,II)
      REAL      UINJ(KK,JJ),VINJ(KK,JJ),WINJ(KK,JJ)
      REAL      UINK(KK,JJ),VINK(KK,JJ),WINK(KK,JJ)

      REAL      FDUI(KK,JJ,II),FDVI(KK,JJ,II),FDWI(KK,JJ,II)
      REAL      FDUJ(KK,JJ),FDVJ(KK,JJ),FDWJ(KK,JJ)
      REAL      FDUK(KK,JJ),FDVK(KK,JJ),FDWK(KK,JJ)
C
C                                  *************************************
C                                  ZEN ZEN ZEN ZEN ZEN ZEN ZEN ZEN ZEN
C                                  *************************************
C
              WF(FFF)  = 0.5*ABS(FFF)+AMIN1(FFF,0.0)
C
C
C
              KSTART = 3
              JSTART = 2
              ISTART = 2
              KSTOP = KMX - 2
              JSTOP = JMX - 2
              ISTOP = IMX - 2
C
      HP    =  0.5
      HN    = -0.5
C
C
C                                 BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 UND BEI GITTERKOPPLUNG
C                                 SOWIE BEI AUSFLUSSBEDINGUNG
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT

      NBU = 0
      IF ((NBAC.EQ.1).OR.(NBAC.EQ.7))  NBU = 1
      IF ((NBAC.EQ.3).OR.(NBAC.EQ.4))  NBU = 1
      IF (NBAC.EQ.8)  NBU = 1
      NLV = 0
      IF ((NLFT.EQ.1).OR.(NLFT.EQ.7).OR.(NLFT.EQ.8))  NLV = 1
      IF ((NLFT.EQ.3).OR.(NLFT.EQ.4))  NLV = 1
      NTW = 0
      IF ((NTOP.EQ.1).OR.(NTOP.EQ.7).OR.(NTOP.EQ.3))  NTW = 1
      IF (NTOP.EQ.8)  NTW = 1
C
C                                 BEI LOKAL VERFEINERTEM GITTER
C                                 WIRD DIE NORMALKOMPONENTE DER ERSTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT

      NFU = 0
      IF ( NFRO .EQ. 8 .OR. NFRO .EQ. 3 ) NFU = 1
      NRV = 0
      IF ( NRGT .EQ. 8 .OR. NRGT .EQ. 3 ) NRV = 1
      NBW = 0
      IF ( NBOT .EQ. 8 .OR. NBOT .EQ. 3 ) NBW = 1

C
C***********************************************************************
C           HIER KOMPAKTVERFAHREN VIERTER ORDNUNG IN X-RICHTUNG        *
C***********************************************************************
CC                                PERIODISCHE RANDBEDINGUNGEN 
CC                                NUR AEQUIDISTANTE GITTER !! 
      IF (NFRO .EQ.1) THEN
C
C
C
C
C
CC                                NICHT-PERIODISCHE RANDBED.
      ELSE
C
C                                 BERECHNUNG DER INTERPOLIERTEN GESCH.
C
C
C                                 BERECHNUNG DER ERSTEN ABLEITUNG       
C
      ENDIF
C
C
C***********************************************************************
C                                 10 10 10 10 10 10 10 10 10 10 10 10 10
      DO  10  I = ISTART, ISTOP
C
C***********************************************************************
C           HIER KOMPAKTVERFAHREN VIERTER ORDNUNG IN Y-RICHTUNG        *
C***********************************************************************
C
CC                                PERIODISCHE RANDBEDINGUNGEN 
CC                                NUR AEQUIDISTANTE GITTER !! 
C     IF (NRGT .EQ.1) THEN
C
C
C
C
CC                                NICHT-PERIODISCHE RANDBED.
C     ELSE
C                                 BERECHNUNG DER INTERPOLIERTEN GESCH.
C#ifdef _KOMYK_
C   
C      CALL INTERPOLATEUWY(KK,JJ,II,KSTART,JSTART,ISTART,KSTOP,JSTOP,
C     $                    ISTOP,DY,DDY,COEFFY,RSGS2,FAKTOR,U(1,1,I),
C     $                    UINJ,NRGT,NLFT,LCOL,DIAG,RCOL)
CC
C      CALL INTERPOLATEVY (KK,JJ,II,KSTART,JSTART,ISTART,KSTOP,JSTOP,
C    $                    ISTOP,DY,DDY,COEFFY,RSGS2,FAKTOR,V(1,1,I),
C    $                    VINJ,NRGT,NLFT,LCOL,DIAG,RCOL)
C    
C     CALL INTERPOLATEUWY(KK,JJ,II,KSTART,JSTART,ISTART,KSTOP,JSTOP,
C    $                    ISTOP,DY,DDY,COEFFY,RSGS2,FAKTOR,W(1,1,I),
C    $                    WINJ,NRGT,NLFT,LCOL,DIAG,RCOL)
C
C#endif
C
C                                 BERECHNUNG DER ERSTEN ABLEITUNG       
C
C#ifdef _KOMYD_
C
C     CALL FDERFOUWY(KK,JJ,II,KSTART,JSTART,ISTART,KSTOP,JSTOP,
C    $               ISTOP,RDY,RDDY,COEFDY,RSGS2,FAKTOR,
C    $               U(1,1,I),UINJ,FDUJ,NRGT,NLFT,LCOL,DIAG,RCOL)
C
C     CALL FDERFOVY (KK,JJ,II,KSTART,JSTART,ISTART,KSTOP,JSTOP,
C    $               ISTOP,RDY,RDDY,COEFDY,RSGS2,FAKTOR,
C    $               V(1,1,I),VINJ,FDVJ,NRGT,NLFT,LCOL,DIAG,RCOL)
C
C     CALL FDERFOUWY(KK,JJ,II,KSTART,JSTART,ISTART,KSTOP,JSTOP,
C    $               ISTOP,RDY,RDDY,COEFDY,RSGS2,FAKTOR,
C    $               W(1,1,I),WINJ,FDWJ,NRGT,NLFT,LCOL,DIAG,RCOL)
C
C#endif
C
C     ENDIF
C
C
C***********************************************************************
C           HIER KOMPAKTVERFAHREN VIERTER ORDNUNG IN Z-RICHTUNG        *
C***********************************************************************
C
C                                 BERECHNUNG DER ERSTEN ABLEITUNG       
C
C
C***********************************************************************
C 
      ICOM = 1 + MOD(I,NBUF)
C
C                                  WALL KORRECTION FUER U,V,W
C
      IF(ABS(WDIF) .GT. SMALL) THEN
C     $   CALL SWCLE1  (KK,JJ,II,KMX,JMX,IMX,DX,DY,DZ,DDX,DDY,DDZ,
C     $                 U,V,W,P,G,TK,
C     $                 B,
C     $                 WCU,WCV,WCW,GMOL,RHO,UGRID,MTURB,ZTOT,GRADPX,I,
C     $                              IC1,IC2,JC1,JC2,KC1,KC2,
C     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,
C     $                 RDDX,RDDY,RDDZ)
      ENDIF
C
C                            BERECHNUNG DER DRUCKGRADIENTEN IN Z-RICHTUNG
C
C
C                                  GEOMETRISCHE KONSTANTEN (I)
                DXI   =  DX(I)
                DXIM  =  DX(I-1)
                DXIP  =  DX(I+1)
               DDXI   = DDX(I)
               DDXIM  = DDX(I-1)
               DDXIP  = DDX(I+1)

               RDXI   =  RDX(I)
               RDXIM  =  RDX(I-1)
               RDXIP  =  RDX(I+1)
              RDDXI   = RDDX(I)
              RDDXIM  = RDDX(I-1)
              RDDXIP  = RDDX(I+1)
C
C                                 20 20 20 20 20 20 20 20 20 20 20 20 20
      DO  20  J = JSTART, JSTOP
C                                  GEOMETRISCHE KONTANTEN (J)
                DYJ   =  DY(J)
                DYJM  =  DY(J-1)
                DYJP  =  DY(J+1)
               DDYJ   = DDY(J)
               DDYJM  = DDY(J-1)
               DDYJP  = DDY(J+1)
C
               RDYJ   =  RDY(J)
               RDYJM  =  RDY(J-1)
               RDYJP  =  RDY(J+1)
              RDDYJ   = RDDY(J)
              RDDYJM  = RDDY(J-1)
              RDDYJP  = RDDY(J+1)
C
                AUZ   =  DXI*DDYJ
                AVZ   = DDXI* DYJ
                AWZ   = DDXI*DDYJ

C
C                                 PUFFERFELDER WERDEN MIT NULL BELEGT
C
      DO  30 K = 1, KMX
          UP(K,J,ICOM) = 0.0
          VP(K,J,ICOM) = 0.0
   30     WP(K,J,ICOM) = 0.0
C
      IF(ABS(WKON) .LE. SMALL) GOTO 2200
C
C                                 **************************************
C                                 KONVEKTIVE TERME   KONVEKTIVE TERME
C                                 **************************************
C
C                             KOEFFIZIENTEN BEIM KONVEKTIVEN ZEITSCHRITT
C
		FKDTU = -1.0*DT*RDDY(J)* RDX(I)*WKON
		FKDTV = -1.0*DT*RDDX(I)* RDY(J)*WKON
		FKDTW = -1.0*DT*RDDX(I)*RDDY(J)*WKON
C
C                                 UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C                                 --------------------------------------
C
C                                 NUR BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 UND BEI GITTERKOPPLUNG
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT

      IF ((NBU.EQ.0).AND.(I .EQ. ISTOP))  GOTO 2110
      IF ((NFU.EQ.0).AND.(I .EQ. ISTART)) GOTO 2110
      IF ( J .EQ. JSTART ) GOTO 2110

C
C
      DO 100  K = KSTART, KSTOP
C                                  GEOMETRISCHE KONSTANTEN (K)
                DZK   =  DZ(K)
                DZKM  =  DZ(K-1)
                DZKMM =  DZ(K-2)
                DZKP  =  DZ(K+1)
               DDZK   =  DDZ(K)
              RDDZK   = RDDZ(K)
C
                AUX   = DDYJ*DDZK
                AUY   =  DXI*DDZK
C                                 KONVEKTIVE TERME
C
C#if defined _KOMXK_ | defined _KOMXK6O_ | defined _KOMXK8O_ | defined _KOU3K_
CC
CC
C      FUE  =     AUX* UINI(K,J,I+1)
C      FUW  =     AUX* UINI(K,J,I)
C      FUN  =     AUY* VINI(K,J,I+1)
C      FUS  =     AUY* VINI(K,J-1,I+1)
C      FUT  =     AUZ* WINI(K,J,I+1)
C      FUB  =     AUZ* WINI(K-1,J,I+1)
CC                                 
C#else
C
      FUE  =     AUX*(U(K,J,I)  +(U(K,J,I+1)-U(K,J,I)  )*0.5*DXI /DDXIP)
      FUW  =     AUX*(U(K,J,I-1)+(U(K,J,I)  -U(K,J,I-1))*0.5*DXIM/DDXI)
      FUN  =     AUY*(V(K,J,I)  + V(K,J,I+1)  )         *0.5
      FUS  =     AUY*(V(K,J-1,I)+ V(K,J-1,I+1))         *0.5
      FUT  =     AUZ*(W(K,J,I)  + W(K,J,I+1)  )         *0.5
      FUB  =     AUZ*(W(K-1,J,I)+ W(K-1,J,I+1))         *0.5
C
C#endif
C
C
C#ifdef 
C                                 JETZT ZENTRALE INTERPOLATION
C
C     QKUE = 0.5*FUE*(U(K,J,I)   + U(K,J,I+1))
C     QKUW = 0.5*FUW*(U(K,J,I-1) + U(K,J,I)  )
C     QKUN = 0.5*FUN*(U(K,J,I)   + U(K,J+1,I))
C     QKUS = 0.5*FUS*(U(K,J-1,I) + U(K,J,I)  )
C     QKUT = 0.5*FUT*(U(K,J,I)   + U(K+1,J,I))
C     QKUB = 0.5*FUB*(U(K-1,J,I) + U(K,J,I)  )
C
C#endif
C                                 JETZT KOMPAKTE INTERPOLATION
C                                 U IN I,J UND K-RICHTUNG
      QKUE = 0.5*FUE*(U(K,J,I)   + U(K,J,I+1))
      QKUW = 0.5*FUW*(U(K,J,I-1) + U(K,J,I)  )
      QKUN = 0.5*FUN*(U(K,J,I)   + U(K,J+1,I))
      QKUS = 0.5*FUS*(U(K,J-1,I) + U(K,J,I)  )
      QKUT = 0.5*FUT*(U(K,J,I)   + U(K+1,J,I)) 
      QKUB = 0.5*FUB*(U(K-1,J,I) + U(K,J,I)  )
C
C
C                                 ZEITSCHRITT
C
C     UP(K,J,ICOM) = UP(K,J,ICOM) - DT/(    AUX*DXI)
C    $             * WKON         * (QKUE-QKUW+QKUN-QKUS+QKUT-QKUB)
      UP(K,J,ICOM) =
     $                FKDTU * RDDZK * (QKUE-QKUW+QKUN-QKUS+QKUT-QKUB)   
C
  100 CONTINUE
C
 2110 CONTINUE
C
C                                 VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
C                                 --------------------------------------
C
C                                 NUR BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT

      IF ((NLV.EQ.0).AND.(J .EQ. JSTOP))  GOTO 2120
      IF ((NRV.EQ.0).AND.(J .EQ. JSTART)) GOTO 2120
      IF ( I .EQ. ISTART ) GOTO 2120

C
C
      DO 110  K = KSTART, KSTOP
C                                  GEOMETRISCHE KONSTANTEN (K)
                DZK   =  DZ(K)
                DZKM  =  DZ(K-1)
                DZKMM =  DZ(K-2)
                DZKP  =  DZ(K+1)
               DDZK   =  DDZ(K)
              RDDZK   = RDDZ(K)
C
                AVX   =  DYJ*DDZK
                AVY   = DDXI*DDZK
C
C                                 KONVEKTIVE F-TERME
C
C
      FVE  =     AVX*(U(K,J,I)  + U(K,J+1,I)  )         *0.5
      FVW  =     AVX*(U(K,J,I-1)+ U(K,J+1,I-1))         *0.5
      FVN  =     AVY*(V(K,J,I)  +(V(K,J+1,I)-V(K,J,I) ) *0.5*DYJ /DDYJP)
      FVS  =     AVY*(V(K,J-1,I)+(V(K,J,I)  -V(K,J-1,I))*0.5*DYJM/DDYJ )
      FVT  =     AVZ*(W(K,J,I)  + W(K,J+1,I)  )         *0.5
      FVB  =     AVZ*(W(K-1,J,I)+ W(K-1,J+1,I))         *0.5
C
C
C
C#ifdef 
C                                 JETZT ZENTRALE INTERPOLATION
C
C     QKVE = 0.5*FVE*(V(K,J,I)   + V(K,J,I+1))
C     QKVW = 0.5*FVW*(V(K,J,I-1) + V(K,J,I)  )
C     QKVN = 0.5*FVN*(V(K,J,I)   + V(K,J+1,I))
C     QKVS = 0.5*FVS*(V(K,J-1,I) + V(K,J,I)  )
C     QKVT = 0.5*FVT*(V(K,J,I)   + V(K+1,J,I))
C     QKVB = 0.5*FVB*(V(K-1,J,I) + V(K,J,I)  )
C
C#endif
C                                 JETZT KOMPAKTE INTERPOLATION
C                                 V IN I,J UND K-RICHTUNG
C
      QKVE = 0.5*FVE*(V(K,J,I)   + V(K,J,I+1))
      QKVW = 0.5*FVW*(V(K,J,I-1) + V(K,J,I)  )
      QKVN = 0.5*FVN*(V(K,J,I)   + V(K,J+1,I))
      QKVS = 0.5*FVS*(V(K,J-1,I) + V(K,J,I)  )
      QKVT = 0.5*FVT*(V(K,J,I)   + V(K+1,J,I))
      QKVB = 0.5*FVB*(V(K-1,J,I) + V(K,J,I)  )
C
C
C                                 ZEITSCHRITT
C
C     VP(K,J,ICOM) = VP(K,J,ICOM) - DT/(    AVY*DYJ)
C    $             * WKON         * (QKVE-QKVW+QKVN-QKVS+QKVT-QKVB)
      VP(K,J,ICOM) =
     $               FKDTV * RDDZK * (QKVE-QKVW+QKVN-QKVS+QKVT-QKVB)
C
  110 CONTINUE
C
 2120 CONTINUE
C
C                                 WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
C                                 --------------------------------------
C
      IF ( I .EQ. ISTART ) GOTO 2130
      IF ( J .EQ. JSTART ) GOTO 2130

      KSTOPM = KSTOP - 1 + NTW
      KSTARTM= KSTART - NBW
C
      DO 120  K = KSTARTM, KSTOPM
C                                  GEOMETRISCHE KONSTANTEN (K)
                DZK   =  DZ(K)
               RDZK   = RDZ(K)
                DZKM  =  DZ(K-1)
                DZKMM =  DZ(K-2)
                DZKP  =  DZ(K+1)
               DDZK   = DDZ(K)
               DDZKM  = DDZ(K-1)
               DDZKP  = DDZ(K+1)
               DDZKPP = DDZ(K+2)
C
                AWX   = DDYJ*DZK
                AWY   = DDXI*DZK
C
C                                 KONVEKTIVE F-TERME
C
C
      FWE  =     AWX*(U(K,J,I)  + U(K+1,J,I)  )         *0.5
      FWW  =     AWX*(U(K,J,I-1)+ U(K+1,J,I-1))         *0.5
      FWN  =     AWY*(V(K,J,I)  + V(K+1,J,I)  )         *0.5
      FWS  =     AWY*(V(K,J-1,I)+ V(K+1,J-1,I))         *0.5
      FWT  =     AWZ*(W(K,J,I)  +(W(K+1,J,I)-W(K,J,I)  )*0.5*DZK /DDZKP)
      FWB  =     AWZ*(W(K-1,J,I)+(W(K,J,I)  -W(K-1,J,I))*0.5*DZKM/DDZK )
C
C
C
C#ifdef 
C                                 JETZT ZENTRALE INTERPOLATION
C 
C     QKWE = 0.5*FWE*(W(K,J,I)   + W(K,J,I+1))
C     QKWW = 0.5*FWW*(W(K,J,I-1) + W(K,J,I)  )
C     QKWN = 0.5*FWN*(W(K,J,I)   + W(K,J+1,I))
C     QKWS = 0.5*FWS*(W(K,J-1,I) + W(K,J,I)  )
C     QKWT = 0.5*FWT*(W(K,J,I)   + W(K+1,J,I))
C     QKWB = 0.5*FWB*(W(K-1,J,I) + W(K,J,I)  )
C
C#endif
C                                 JETZT KOMPAKTE INTERPOLATION
C                                 W IN I,J UND K-RICHTUNG
      QKWE = 0.5*FWE*(W(K,J,I)   + W(K,J,I+1))
      QKWW = 0.5*FWW*(W(K,J,I-1) + W(K,J,I)  )
      QKWN = 0.5*FWN*(W(K,J,I)   + W(K,J+1,I))
      QKWS = 0.5*FWS*(W(K,J-1,I) + W(K,J,I)  )
      QKWT = 0.5*FWT*(W(K,J,I)   + W(K+1,J,I))
      QKWB = 0.5*FWB*(W(K-1,J,I) + W(K,J,I)  )
C
C
C                                 ZEITSCHRITT
C     WP(K,J,ICOM) = WP(K,J,ICOM) - DT/(    AWZ*DZK)
C    &             * WKON         * (QKWE-QKWW+QKWN-QKWS+QKWT-QKWB)
      WP(K,J,ICOM) =
     $               FKDTW * RDZK  * (QKWE-QKWW+QKWN-QKWS+QKWT-QKWB)
C
  120 CONTINUE
C
 2130 CONTINUE
 2200 IF(ABS(WDIF) .LE. SMALL) GOTO 2300
C
C                                 **************************************
C                                 DIFFUSIVE TERME   DIFFUSIVE TERME   DI
C                                 **************************************
C
C                             KOEFFIZIENTEN BEIM DIFFUSIVEN ZEITSCHRITT
C
		FDDTU = -1.0*DT/RHO*RDDY(J)* RDX(I)*WDIF
		FDDTV = -1.0*DT/RHO*RDDX(I)* RDY(J)*WDIF
		FDDTW = -1.0*DT/RHO*RDDX(I)*RDDY(J)*WDIF
C
C                                 UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C                                 --------------------------------------
C
C
C                                 NUR BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT

      IF ((NBU.EQ.0).AND.(I .EQ. ISTOP))  GOTO 2210
      IF ((NFU.EQ.0).AND.(I .EQ. ISTART)) GOTO 2210
      IF ( J .EQ. JSTART ) GOTO 2210
C
C
      DO 200  K = KSTART, KSTOP
C                                  GEOMETRISCHE KONSTANTEN (K)
               RDZK   =  RDZ(K)
               RDZKM  =  RDZ(K-1)
               DDZK   =  DDZ(K)
              RDDZK   = RDDZ(K)
C
                AUX   = DDYJ*DDZK
                AUY   =  DXI*DDZK
C
C                                 DIFFUSIVE D-TERME MIT
C                                 TRANSPORTKOEFFIZIENTEN
      GUE  = G(K,J,I+1)
      GUW  = G(K,J,I)
      GUN  = G(K,J,I)  *G(K,J+1,I)   / (G(K,J,I)  +G(K,J+1,I)  )
     $     + G(K,J,I+1)*G(K,J+1,I+1) / (G(K,J,I+1)+G(K,J+1,I+1))
      GUS  = G(K,J,I)  *G(K,J-1,I)   / (G(K,J,I)  +G(K,J-1,I)  )
     $     + G(K,J,I+1)*G(K,J-1,I+1) / (G(K,J,I+1)+G(K,J-1,I+1))
      GUT  = G(K,J,I)  *G(K+1,J,I)   / (G(K,J,I)  +G(K+1,J,I)  )
     $     + G(K,J,I+1)*G(K+1,J,I+1) / (G(K,J,I+1)+G(K+1,J,I+1))
      GUB  = G(K,J,I)  *G(K-1,J,I)   / (G(K,J,I)  +G(K-1,J,I)  )
     $     + G(K,J,I+1)*G(K-1,J,I+1) / (G(K,J,I+1)+G(K-1,J,I+1))
C
C                                     HIER KOMPAKTE BERECHNUNG
C                                     DER ERSTEN ABLEITUNG U
      QDUE = -GUE * AUX*RDDXIP * (U(K,J,I+1) - U(K,J,I))
      QDUW = -GUW * AUX*RDDXI  * (U(K,J,I)   - U(K,J,I-1))
      QDUN = -GUN * AUY*RDYJ   * (U(K,J+1,I) - U(K,J,I))
      QDUS = -GUS * AUY*RDYJM  * (U(K,J,I)   - U(K,J-1,I))
      QDUT = -GUT * AUZ*RDZK   * (U(K+1,J,I) - U(K,J,I))
      QDUB = -GUB * AUZ*RDZKM  * (U(K,J,I)   - U(K-1,J,I))
C
C     QDUE = -GUE * AUX*RDDXIP * (U(K,J,I+1) - U(K,J,I))
C     QDUW = -GUW * AUX*RDDXI  * (U(K,J,I)   - U(K,J,I-1))
C     QDUN = -GUN * AUY*RDYJ   * (U(K,J+1,I) - U(K,J,I))
C     QDUS = -GUS * AUY*RDYJM  * (U(K,J,I)   - U(K,J-1,I))
C     QDUT = -GUT * AUZ*RDZK   * (U(K+1,J,I) - U(K,J,I))
C     QDUB = -GUB * AUZ*RDZKM  * (U(K,J,I)   - U(K-1,J,I))
C
C
      QDUC = WCU(K,J)
C
C
C                                 ZEITSCHRITT
C
C     UP(K,J,ICOM) = UP(K,J,ICOM) - DT/(RHO*AUX*DXI)
C    $             * WDIF         * (QDUE-QDUW+QDUN-QDUS+QDUT-QDUB-QDUC)
      UP(K,J,ICOM) = UP(K,J,ICOM) 
     $            + FDDTU * RDDZK * (QDUE-QDUW+QDUN-QDUS+QDUT-QDUB-QDUC) 
C
  200 CONTINUE
C
 2210 CONTINUE
C
C                                 VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
C                                 --------------------------------------
C
C                                 NUR BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT

      IF ((NLV.EQ.0).AND.(J .EQ. JSTOP))  GOTO 2220
      IF ((NRV.EQ.0).AND.(J .EQ. JSTART)) GOTO 2220
      IF ( I .EQ. ISTART ) GOTO 2220
C
C
      DO 210  K = KSTART, KSTOP
C                                  GEOMETRISCHE KONSTANTEN (K)
                RDZK   =  RDZ(K)
                RDZKM  =  RDZ(K-1)
               DDZK   = DDZ(K)
              RDDZK   =RDDZ(K)
C
                AVX   =  DYJ*DDZK
                AVY   = DDXI*DDZK
C
C                                 DIFFUSIVE D-TERME MIT
C                                 TRANSPORTKOEFFIZIENTEN
C
      GVE  = G(K,J,I)  *G(K,J,I+1)   / (G(K,J,I)+  G(K,J,I+1)  )
     $     + G(K,J+1,I)*G(K,J+1,I+1) / (G(K,J+1,I)+G(K,J+1,I+1))
      GVW  = G(K,J,I)  *G(K,J,I-1)   / (G(K,J,I)  +G(K,J,I-1)  )
     $     + G(K,J+1,I)*G(K,J+1,I-1) / (G(K,J+1,I)+G(K,J+1,I-1))
      GVN  = G(K,J+1,I)
      GVS  = G(K,J,I)
      GVT  = G(K,J,I)  *G(K+1,J,I)   / (G(K,J,I)  +G(K+1,J,I)  )
     $     + G(K,J+1,I)*G(K+1,J+1,I) / (G(K,J+1,I)+G(K+1,J+1,I))
      GVB  = G(K,J,I)  *G(K-1,J,I)   / (G(K,J,I)  +G(K-1,J,I)  )
     $     + G(K,J+1,I)*G(K-1,J+1,I) / (G(K,J+1,I)+G(K-1,J+1,I))
C
C                                     HIER KOMPAKTE BERECHNUNG
C                                     DER ERSTEN ABLEITUNG V
      QDVE = -GVE * AVX*RDXI   * (V(K,J,I+1) - V(K,J,I))
      QDVW = -GVW * AVX*RDXIM  * (V(K,J,I)   - V(K,J,I-1))
      QDVN = -GVN * AVY*RDDYJP * (V(K,J+1,I) - V(K,J,I))
      QDVS = -GVS * AVY*RDDYJ  * (V(K,J,I)   - V(K,J-1,I))
      QDVT = -GVT * AVZ*RDZK   * (V(K+1,J,I) - V(K,J,I))
      QDVB = -GVB * AVZ*RDZKM  * (V(K,J,I)   - V(K-1,J,I))

C     QDVE = -GVE * AVX*RDXI   * (V(K,J,I+1) - V(K,J,I))
C     QDVW = -GVW * AVX*RDXIM  * (V(K,J,I)   - V(K,J,I-1))
C     QDVN = -GVN * AVY*RDDYJP * (V(K,J+1,I) - V(K,J,I))
C     QDVS = -GVS * AVY*RDDYJ  * (V(K,J,I)   - V(K,J-1,I))
C     QDVT = -GVT * AVZ*RDZK   * (V(K+1,J,I) - V(K,J,I))
C     QDVB = -GVB * AVZ*RDZKM  * (V(K,J,I)   - V(K-1,J,I))
C
C
      QDVC =  WCV(K,J)
C
C
C                                  ZEITSCHRITT
C
C     VP(K,J,ICOM) = VP(K,J,ICOM) - DT/(RHO*AVY*DYJ)
C    $             * WDIF         * (QDVE-QDVW+QDVN-QDVS+QDVT-QDVB-QDVC)
      VP(K,J,ICOM) = VP(K,J,ICOM) 
     $            + FDDTV * RDDZK * (QDVE-QDVW+QDVN-QDVS+QDVT-QDVB-QDVC)
C
  210 CONTINUE
C
 2220 CONTINUE
C
C                                 WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
C                                 --------------------------------------
C
      IF ( I .EQ. ISTART ) GOTO 2230
      IF ( J .EQ. JSTART ) GOTO 2230
C
      KSTOPM = KSTOP - 1 + NTW
      KSTARTM= KSTART - NBW
C
      DO 220  K = KSTARTM, KSTOPM
C                                  GEOMETRISCHE KONSTANTEN (K)
                DZK   =   DZ(K)
               RDZK   =  RDZ(K)
              RDDZK   = RDDZ(K)
              RDDZKP  = RDDZ(K+1)
C
                AWX   = DDYJ*DZK
                AWY   = DDXI*DZK
C
C                                 DIFFUSIVE D-TERME MIT
C                                 TRANSPORTKOEFFIZIENTEN
C
      GWE  = G(K,J,I)  *G(K,J,I+1)   / (G(K,J,I)  +G(K,J,I+1)  )
     $     + G(K+1,J,I)*G(K+1,J,I+1) / (G(K+1,J,I)+G(K+1,J,I+1))
      GWW  = G(K,J,I)  *G(K,J,I-1)   / (G(K,J,I)  +G(K,J,I-1)  )
     $     + G(K+1,J,I)*G(K+1,J,I-1) / (G(K+1,J,I)+G(K+1,J,I-1))
      GWN  = G(K,J,I)  *G(K,J+1,I)   / (G(K,J,I)  +G(K,J+1,I)  )
     $     + G(K+1,J,I)*G(K+1,J+1,I) / (G(K+1,J,I)+G(K+1,J+1,I))
      GWS  = G(K,J,I)  *G(K,J-1,I)   / (G(K,J,I)  +G(K,J-1,I)  )
     $     + G(K+1,J,I)*G(K+1,J-1,I) / (G(K+1,J,I)+G(K+1,J-1,I))
      GWT  = G(K+1,J,I)
      GWB  = G(K,J,I)
C
C                                     HIER KOMPAKTE BERECHNUNG
C                                     DER ERSTEN ABLEITUNG W
      QDWE = -GWE * AWX*RDXI   * (W(K,J,I+1) - W(K,J,I))
      QDWW = -GWW * AWX*RDXIM  * (W(K,J,I)   - W(K,J,I-1))
      QDWN = -GWN * AWY*RDYJ   * (W(K,J+1,I) - W(K,J,I))
      QDWS = -GWS * AWY*RDYJM  * (W(K,J,I)   - W(K,J-1,I))
      QDWT = -GWT * AWZ*RDDZKP * (W(K+1,J,I) - W(K,J,I))
      QDWB = -GWB * AWZ*RDDZK  * (W(K,J,I)   - W(K-1,J,I))
C
C     QDWE = -GWE * AWX*RDXI   * (W(K,J,I+1) - W(K,J,I))
C     QDWW = -GWW * AWX*RDXIM  * (W(K,J,I)   - W(K,J,I-1))
C     QDWN = -GWN * AWY*RDYJ   * (W(K,J+1,I) - W(K,J,I))
C     QDWS = -GWS * AWY*RDYJM  * (W(K,J,I)   - W(K,J-1,I))
C     QDWT = -GWT * AWZ*RDDZKP * (W(K+1,J,I) - W(K,J,I))
C     QDWB = -GWB * AWZ*RDDZK  * (W(K,J,I)   - W(K-1,J,I))
C
C
      QDWC = WCW(K,J)
C
C
C                                 ZEITSCHRITT
C     WP(K,J,ICOM) = WP(K,J,ICOM) - DT/(RHO*AWZ*DZK)
C    $             * WDIF         * (QDWE-QDWW+QDWN-QDWS+QDWT-QDWB-QDWC)
      WP(K,J,ICOM) = WP(K,J,ICOM) 
     $            + FDDTW * RDZK  * (QDWE-QDWW+QDWN-QDWS+QDWT-QDWB-QDWC)
C
  220 CONTINUE
C
 2230 CONTINUE
C
 2300 IF(IDUZ .EQ. 1) GOTO 2400
C
C                                 **************************************
C                                 QUELL-TERME   QUELL-TERME   QUELL-TERM
C                                 **************************************
C
C                                 UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C                                 --------------------------------------
C
C
C                                 NUR BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT

      IF ((NBU.EQ.0).AND.(I .EQ. ISTOP))  GOTO 2310
      IF ((NFU.EQ.0).AND.(I .EQ. ISTART)) GOTO 2310
      IF ( J .EQ. JSTART ) GOTO 2310

C
C
      DO 300  K = KSTART, KSTOP
C
C                                 ACHTUNG: IN X-RICHTUNG KANN EIN
C                                 MITTLERER DRUCKGRADIENT VORGEGEBEN
C                                 WERDEN
C
C                                 ZEITSCHRITT
C
      UP(K,J,ICOM) = UP(K,J,ICOM) - DT/(RHO*DXI) * WSOR
     $             * (P(K,J,I+1)-P(K,J,I) + GRADPX*DXI)  
C
  300 CONTINUE
C
 2310 CONTINUE
C
C                                 VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
C                                 --------------------------------------
C                                 NUR BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT

      IF ((NLV.EQ.0).AND.(J .EQ. JSTOP))  GOTO 2320
      IF ((NRV.EQ.0).AND.(J .EQ. JSTART)) GOTO 2320
      IF (I .EQ. ISTART) GOTO 2320

C
C
      DO 310  K = KSTART, KSTOP
C
C                                  ZEITSCHRITT
C
      VP(K,J,ICOM) = VP(K,J,ICOM) - DT/(RHO*DYJ) * WSOR
     $             * (P(K,J+1,I)-P(K,J,I))
C
  310 CONTINUE
C
 2320 CONTINUE
C
C                                 WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
C                                 --------------------------------------

      IF (I .EQ. ISTART) GOTO 2400
      IF (J .EQ. JSTART) GOTO 2400

      KSTOPM = KSTOP - 1 + NTW
      KSTARTM= KSTART - NBW
C
      DO 320  K = KSTARTM, KSTOPM
C
C                                 ZEITSCHRITT
C
      WP(K,J,ICOM) = WP(K,J,ICOM) - DT/(RHO*DZ(K)) * WSOR
     $             * (P(K+1,J,I)-P(K,J,I))
C
  320 CONTINUE
C
 2400 IF(ABS(WPHI) .LE. SMALL) GOTO 2500
C
C                                 **************************************
C                                 ANTEIL DER PUNKTWERTE   ANTEIL DER PUN
C                                 **************************************
C
C                                 UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C                                 --------------------------------------
C
C                                 NUR BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT

      IF ((NBU.EQ.0).AND.(I .EQ. ISTOP))  GOTO 2410
      IF ((NFU.EQ.0).AND.(I .EQ. ISTART)) GOTO 2410
      IF ( J .EQ. JSTART ) GOTO 2410
C
C
      DO 400  K = KSTART, KSTOP
C
C                                 ZEITSCHRITT
C
      UP(K,J,ICOM) = UP(K,J,ICOM) + WPHI * U(K,J,I) 
  400 CONTINUE
C
 2410 CONTINUE
C
C                                 VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
C                                 --------------------------------------
C                                 NUR BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT

      IF ((NLV.EQ.0).AND.(J .EQ. JSTOP))  GOTO 2420
      IF ((NRV.EQ.0).AND.(J .EQ. JSTART)) GOTO 2420
      IF (I .EQ. ISTART) GOTO 2420

C
C
      DO 410  K = KSTART, KSTOP
C
C                                  ZEITSCHRITT
C
  410 VP(K,J,ICOM) = VP(K,J,ICOM) + WPHI * V(K,J,I)
C
 2420 CONTINUE
C
C                                 WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
C                                 --------------------------------------
      IF (I .EQ. ISTART) GOTO 2500
      IF (J .EQ. JSTART) GOTO 2500
      KSTOPM = KSTOP - 1 + NTW
      KSTARTM= KSTART - NBW
C
      DO 420  K = KSTARTM, KSTOPM
C
C                                 ZEITSCHRITT
C
  420 WP(K,J,ICOM) = WP(K,J,ICOM) + WPHI * W(K,J,I)
C
 2500 CONTINUE
C
   20 CONTINUE
C                                 20 20 20 20 20 20 20 20 20 20 20 20 20
C
C
C                                 JETZT ZURUECKSPEICHERN DER TEMPORAEREN
C                                 PUFFER. NUR TATSAECHLICH BERECHNETE
C                                 PUFFERWERTE WERDEN ZURUECKGESPEICHERT
      IF(I .LT. (ISTART+NBUF-1)) GOTO 3010

      INEW   = I
      IREPMX = 1
      IF(I .EQ. ISTOP) IREPMX = NBUF
C
      DO 810 IREP=1,IREPMX
         IF(I .LT. ISTOP) GOTO 3020
            INEW = ISTOP + IREP - 1
 3020    IBUFF = 1 + MOD(INEW+1,NBUF)
         IBACK = INEW - NBUF + 1
C

         DO 820 J=JSTART,JSTOP
C
C                                 NUR BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT

      IF ((NBU.EQ.0).AND.(INEW .EQ. (ISTOP+NBUF-1)))  GOTO 3030
      IF ((NFU.EQ.0).AND.(I .EQ. ISTART)) GOTO 3030
      IF (J .EQ. JSTART) GOTO 3030

            DO 830 K=KSTART,KSTOP              
               UO(K,J,IBACK) = UO(K,J,IBACK)*ZERONE + UP(K,J,IBUFF)
  830       CONTINUE
 3030       CONTINUE

C                                 NUR BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT

      IF ((NLV.EQ.0).AND.(J .EQ. JSTOP))  GOTO 3040
      IF ((NRV.EQ.0).AND.(J .EQ. JSTART)) GOTO 3040
      IF (I .EQ. ISTART) GOTO 3040

            DO 840 K=KSTART,KSTOP
  840          VO(K,J,IBACK) = VO(K,J,IBACK)*ZERONE + VP(K,J,IBUFF)
 3040               CONTINUE
      IF (I .EQ. ISTART) GOTO 3050
      IF (J .EQ. JSTART) GOTO 3050
            DO 850 K=KSTARTM,KSTOPM
  850          WO(K,J,IBACK) = WO(K,J,IBACK)*ZERONE + WP(K,J,IBUFF)
C
 3050               CONTINUE
  820    CONTINUE
  810 CONTINUE
 3010 CONTINUE
C
   10 CONTINUE
C                                 10 10 10 10 10 10 10 10 10 10 10 10 10
C
      RETURN
      END
