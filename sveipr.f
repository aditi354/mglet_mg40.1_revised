










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

      SUBROUTINE SVEIPR   (KK,JJ,II,KMX,JMX,IMX,
     $                    UFR,VFR,WFR,PFR,GFR,Y,Z,DZ,VCON,WCON,
     $                    YBANF,YBEND,ZBANF,ZBEND,FREQB,ANIVEAU,
     $                    AUB,AVB,AWB,TIMEPH,DT,
     $                    UBO,VBO,WBO)
C*STARLET***************************************************************
C        S V E I P R      VORBELEGUNG DER GESCHWINDIGKEITSFELDER 
C                         UFR
C                         FUER DAS U-FELD WIRD EINE VOLL ENTWICKELTE
C                         KANALSTROEMUNG ANGENOMMEN. NACH DER ERZEUGUNG
C                         DER ZEITLICH GEMITTELTEN PROFILE WERDEN DIESEN
C                         WERTEN SCHWANKUNGSGESCHWINDIGKEITEN UEBERLA-
C                         GERT.
C*STARLET***************************************************************
C
C PARAM: KK,  JJ,  II   - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        TK(KK,JJ,II)   + TURBULENZENERGIE (= K)
C        UFR(KK,JJ,2)   + EINSTROEMPROFIL FUER _FRFIX_
C        MTURB          - SCHALTER (1 : TURBULENT,  0 : LAMINAR)
C        VCON           - V-GESCHWINDIGKEITSKOMP. AM EINTRITTSRAND
C        WCON           - W-GESCHWINDIGKEITSKOMP. AM EINTRITTSRAND
C        VREF           - UEBER DIE KANALHOEHE GEMITTELTE GESCHWINDIGKEI
C        EXPON          - EXPONENT DES GESCHW.-PROFILES
C        GMOL           - MOLEKULARE DYN. VISKOSITAET
C        RHO            - DICHTE DES FLUIDS (= CONST)
C        Z (KK)         - Z-KOORDINATEN DER ZELLMITTELPUNKTE
C        DZ(KK)         - ABSTAND DER GITTERPUNKTE IN Z-RICHTUNG
C        TU             - TURBULENZGRAD DER ANSTROEMUNG
C        UFRFREQ        - FREQUENZ DER ZUSTROEMBEDINGUNG
C        UFRCON         - AMPLITUDE DER ZUSTROEMBEDINGUNG
C        UBO(2,JJ,II)  + U-GESCHWINDIGKEITSKOMP.
C        VBO(2,JJ,II)  + V-GESCHWINDIGKEITSKOMP.
C
C VERS:  14. 4.92 (MM)  : AUS SVLE1.SR11 ABGELEITET
C        21. 7.95 (MM)  : 'DUCT' EINGEFUEHRT
C         8. 3.96 (MM)  : YBANF...AWB EINGEFUEHRT
C
C UPROG: FUNCTION VORZ
C        FUNCTION RANEQU
C
C*STARLET***************************************************************
C

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

      COMMON /CONLES/  CONV2S,CMUE,CAPPA,ECONST
      SAVE   /CONLES/
C
      REAL 
     $           UFR(KK,JJ, 2),  VFR(KK,JJ, 2),  WFR(KK,JJ, 2),
     $           UBO(2,JJ,II),  VBO(2,JJ,II),  WBO(2,JJ,II),
     $           PFR(KK,JJ, 2),  GFR(KK,JJ, 2),
     $            Y(JJ),
     $           DZ(KK),         Z(KK),
     $           ETA(41),        F(41),          DF(41),
     $           ETAN,           FN,             DFN,
     $           XIDELTA,        KDELTAI,    WFREND,
     $           FAKTOR,         REX,         XIKDELTA
     $           R
C

      UEIN  (ZWAND)  = (1.0+EXPON)/(2.0-RTURB)*VREF*((1.0-RTURB)
     $               + SIGN(1.0,(RTURB-0.5))*(ABS(HK2*(1.0-RTURB)
     $               - ABS(ZWAND))/HK2)**EXPON) * SIGN(1.0,ZWAND)
      TKEIN (ZMITTE) = CK1 + CK2 * ABS(ZMITTE)**CK3
C
       DATA           REV         /20.0/
       DATA           TU          / 0.10/
       DATA           PI          /3.1415927/
C         DATA-STATEMENT FOR BLASIUS-PROFILE
       DATA           ETA(1)      /0.0/
       DATA           ETA(2)      /0.1/
       DATA           ETA(3)      /0.2/
       DATA           ETA(4)      /0.3/
       DATA           ETA(5)      /0.4/
       DATA           ETA(6)      /0.5/
       DATA           ETA(7)      /0.6/
       DATA           ETA(8)      /0.7/
       DATA           ETA(9)      /0.8/
       DATA           ETA(10)     /0.9/
       DATA           ETA(11)     /1.0/
       DATA           ETA(12)     /1.1/
       DATA           ETA(13)     /1.2/
       DATA           ETA(14)     /1.3/
       DATA           ETA(15)     /1.4/
       DATA           ETA(16)     /1.5/
       DATA           ETA(17)     /1.6/
       DATA           ETA(18)     /1.7/
       DATA           ETA(19)     /1.8/
       DATA           ETA(20)     /1.9/
       DATA           ETA(21)     /2.0/
       DATA           ETA(22)     /2.2/
       DATA           ETA(23)     /2.4/
       DATA           ETA(24)     /2.6/
       DATA           ETA(25)     /2.8/
       DATA           ETA(26)     /3.0/
       DATA           ETA(27)     /3.2/
       DATA           ETA(28)     /3.4/
       DATA           ETA(29)     /3.6/
       DATA           ETA(30)     /3.8/
       DATA           ETA(31)     /4.0/
       DATA           ETA(32)     /4.2/
       DATA           ETA(33)     /4.4/
       DATA           ETA(34)     /4.6/
       DATA           ETA(35)     /4.8/
       DATA           ETA(36)     /5.0/
       DATA           ETA(37)     /5.2/
       DATA           ETA(38)     /5.4/
       DATA           ETA(39)     /5.6/
       DATA           ETA(40)     /5.8/
       DATA           ETA(41)     /6.0/
       DATA           F(1)        /0.0/
       DATA           F(2)        /0.00235/
       DATA           F(3)        /0.00939/
       DATA           F(4)        /0.02113/
       DATA           F(5)        /0.03755/
       DATA           F(6)        /0.05864/
       DATA           F(7)        /0.08439/
       DATA           F(8)        /0.11474/
       DATA           F(9)        /0.14967/
       DATA           F(10)       /0.18911/
       DATA           F(11)       /0.23299/
       DATA           F(12)       /0.28121/
       DATA           F(13)       /0.33366/
       DATA           F(14)       /0.39021/
       DATA           F(15)       /0.45072/
       DATA           F(16)       /0.51503/
       DATA           F(17)       /0.58296/
       DATA           F(18)       /0.65430/
       DATA           F(19)       /0.72887/
       DATA           F(20)       /0.80644/
       DATA           F(21)       /0.88680/
       DATA           F(22)       /1.05495/
       DATA           F(23)       /1.23153/
       DATA           F(24)       /1.41482/
       DATA           F(25)       /1.60328/
       DATA           F(26)       /1.79557/
       DATA           F(27)       /1.99058/
       DATA           F(28)       /2.18747/
       DATA           F(29)       /2.38559/
       DATA           F(30)       /2.58450/
       DATA           F(31)       /2.78388/
       DATA           F(32)       /2.98355/
       DATA           F(33)       /3.18338/
       DATA           F(34)       /3.38329/
       DATA           F(35)       /3.58325/
       DATA           F(36)       /3.78323/
       DATA           F(37)       /3.98322/
       DATA           F(38)       /4.18322/
       DATA           F(39)       /4.38322/
       DATA           F(40)       /4.58322/
       DATA           F(41)       /4.78322/
       DATA           DF(1)        /0.0/
       DATA           DF(2)        /0.04696/
       DATA           DF(3)        /0.09391/
       DATA           DF(4)        /0.14081/
       DATA           DF(5)        /0.18761/
       DATA           DF(6)        /0.23423/
       DATA           DF(7)        /0.28058/
       DATA           DF(8)        /0.32653/
       DATA           DF(9)        /0.37196/
       DATA           DF(10)        /0.41672/
       DATA           DF(11)       /0.46063/
       DATA           DF(12)       /0.50354/
       DATA           DF(13)       /0.54525/
       DATA           DF(14)       /0.58559/
       DATA           DF(15)       /0.62439/
       DATA           DF(16)       /0.66147/
       DATA           DF(17)       /0.69670/
       DATA           DF(18)       /0.72993/
       DATA           DF(19)       /0.76106/
       DATA           DF(20)       /0.79000/
       DATA           DF(21)       /0.81669/
       DATA           DF(22)       /0.86330/
       DATA           DF(23)       /0.90107/
       DATA           DF(24)       /0.93060/
       DATA           DF(25)       /0.95288/
       DATA           DF(26)       /0.96905/
       DATA           DF(27)       /0.98037/
       DATA           DF(28)       /0.98797/
       DATA           DF(29)       /0.99289/
       DATA           DF(30)       /0.99594/
       DATA           DF(31)       /0.99777/
       DATA           DF(32)       /0.99882/
       DATA           DF(33)       /0.99940/
       DATA           DF(34)       /0.99970/
       DATA           DF(35)       /0.99986/
       DATA           DF(36)       /0.99994/
       DATA           DF(37)       /0.999971/
       DATA           DF(38)       /0.999988/
       DATA           DF(39)       /0.999995/
       DATA           DF(40)       /0.999998/
       DATA           DF(41)       /0.999999/
C
CTEST      WRITE (0,*) 
CTEST     $                    YBANF,YBEND,ZBANF,ZBEND,FREQB,ANIVEAU,
CTEST     $                    AUB,AVB,AWB
C
      IF (UFRFREQ .LT. SMALL) THEN

          AMPLITUDE = 1.0

      ELSE

          AMPLITUDE = SIN( 2.0 * PI * TIMEPH * UFRFREQ )

      ENDIF
C
C
      IM2 = IMX-2
      JM1 = JMX-1
      JM2 = JMX-2
      JM3 = JMX-3
      KM1 = KMX-1
      KM2 = KMX-2
      KM3 = KMX-3
C
C                                 EINIGE KONSTANTEN
C
      CM25   = CMUE**0.25
      SQCMUE = CM25**2
      CM25RE = CM25* REV
      RTURB  = FLOAT(MTURB)
C
C                                 VORBELEGUNG DER GESCHWINDIGKEITSFELDER
C                                 (MIT PRESET !!)
C
      DO 10 I=1,2
         DO 10 J=1,JMX
            DO 10 K=1,KMX
               UFR(K,J,I) = PRESET
               VFR(K,J,I) = PRESET
               WFR(K,J,I) = PRESET
               PFR(K,J,I) = PRESET
               GFR(K,J,I) = PRESET
   10 CONTINUE
C
C                                 **************************************
C                                 ZEITLICH GEMITTELTE GESCHW.PROFILE
C                                 **************************************
C
      IF(CIDUFR(1:7).EQ.'UNIFORM') THEN
C
C    					UNIFORME GESCHWINDIGKEITSVERT.
C					IM GANZEN FELD
C
      DO 70 I=1,2
      DO 70 J=2,JM1
      DO 70 K=2,KM1
            UFR(K,J,I) = UFRCON
   70 CONTINUE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   BEI 'UNIFORM' KEIN RANDOM!
      RETURN
C      IF(MTURB .EQ. 0) RETURN

      DO  71 I=1,2
      DO  72 J=3,JM2
         DO  73 K=3,KM1
            UFR(K,J,I) = UFR(K,J,I)+SQRT(2.0*RANEQU(0.123)*TU)*VORZ(1)
   73    CONTINUE
   72 CONTINUE
   71 CONTINUE

      RETURN
C
      ELSEIF(CIDUFR(1:7).EQ.'CHANNEL') THEN
C
C                                 BEI LAMINARER STROEMUNG WIRD EIN PARA-
C                                 BOLISCHES GESCHWINDIGKEITSPROFIL ER-
C                                 ZEUGT. DIESES PROFIL WIRD ZWISCHEN
C                                 ZBANF UND ZBEND GELEGT.
C
      IF(MTURB .EQ. 0) THEN
          EXPON = 2.0

          HEIGHT = ZBEND-ZBANF
          HEIGHT2 = HEIGHT**2
          HK2   = 0.5*(ZBEND-ZBANF)

            DO K=3,KM2
               IF ( Z(K) .GE. ZBANF .AND. Z(K) .LE. ZBEND ) THEN
                 DO I=1,2
                 DO J=3,JM2
                    DISTANCE = MIN(Z(K)-ZBANF,ZBEND-Z(K))
                    UFR(K,J,I) = 
     $                   UFRCON*DISTANCE*(HEIGHT - DISTANCE)/(HEIGHT2)
C TEST                    write (6,*)
C TEST     $                  'zbanf,...',ZBANF,ZBEND,Z(K),DISTANCE,UFR(K,J,I)
                 ENDDO
                 ENDDO
               ELSE
                 DO I=1,2
                 DO J=3,JM2
                    UFR(K,J,I) = 0.0
                 ENDDO
                 ENDDO
               ENDIF
            ENDDO
          RETURN
      ELSE
C                                 BEI TURBULENTER STROEMUNG, WIE VOR (8.3.96)
C
C                                 GEOMETRISCHE GROESSEN, INDIZES DES
C                                 ERSTEN UND LETZTEN PUNKTES DES FREI
C                                 DURCHSTROEMTEN KANALS
C
      NSTART = 3
      NSTOP  = KMX-2
      NSTUE  = (NSTOP-NSTART)/2
      NMITTE = (NSTOP+NSTART)/2
      NZUS   =  NSTOP+NSTART-2*NMITTE
C
      ZBOT   = 0.5*(Z(NSTART) + Z(NSTART-1))
      ZTOP   = 0.5*(Z(NSTOP+1)+ Z(NSTOP)   )
      ZMIT   = 0.5*(ZBOT+ZTOP)
C
C                                 HALBE KANALHOEHE
C
      HK2    = 0.5*(ZTOP-ZBOT)
C
      NSTP   = NSTUE + 1
C
      DO 100 I=1,2
      DO 110 J=3,JM2
         DO 120 N=1,NSTP
            NMIN        = NMITTE + 1 - N
            NPLU        = NMITTE - 1 + NZUS + N
C
            UFR(NMIN,J,I) = UEIN(Z(NMIN)-ZBOT)
  120       UFR(NPLU,J,I) = UEIN(ZTOP-Z(NPLU))
C
C
  110 CONTINUE
  100 CONTINUE
C
C
C
C                                  *************************************
C                                  ERMITTLUNG DER SCHWANKUNGSGESCHWIN-
C                                  DIGKEITEN MITTELS ZUFALLSZAHLENGENE-
C                                  RATOR. DIE SCHWANKUNGEN WERDEN SO BE-
C                                  STIMMT, DASS EINE VORGEGEBENE TURBU-
C                                  LENZENERGIEVERTEILUNG EINGEHALTEN
C                                  WIRD.
C                                  *************************************
C
C                                  KONSTANTEN DER K-VERTEILUNG (AN-
C                                  NAEHERUNG EINER MIT DEM K-EPS-MODELL
C                                  BERECHNETEN K-VERTEILUNG)
C
      CK1    = 0.0015
      CK2    = 0.017
      CK3    = 1.443
C
C                                  NAEHERUNGSWEISE BEST. DER DICKE DER
C                                  VISKOSEN UNTERSCHICHT
C
      DYV    = (HK2**EXPON * CM25RE*GMOL*LOG(ECONST*CM25RE)
     $       /  (CAPPA*RHO*(1.0+EXPON)))**(1.0/(1.0+EXPON))
      TKV    = TKEIN(HK2-DYV)
      GRADKV = TKV/DYV
C
      NSTP = NSTUE + 1
      DO 170 I=1,2
      DO 180 J=3,JM2
         DO 190 N=1,NSTP
            NMIN         = NMITTE + 1 - N
            NPLU         = NMITTE - 1 + NZUS + N
            ZMMIN        = ZMIT - Z(NMIN)
            ZMPLU        = Z(NPLU) - ZMIT
            TKMIN        = TKEIN(ZMMIN)
            TKPLU        = TKEIN(ZMPLU)
            IF((HK2-ZMMIN) .LE. DYV) TKMIN = GRADKV*(HK2-ZMMIN)
            IF((HK2-ZMPLU) .LE. DYV) TKPLU = GRADKV*(HK2-ZMPLU)
C
C                                 ERSTE SCHAETZUNG FUER DIE QUADRATE DER
C                                 SCHWANKUNGSGESCHWINDIGKEITEN
C
            US2M  = RANEQU(0.123)
            US2P  = RANEQU(0.123)
            VS2M  = RANEQU(0.123)
            VS2P  = RANEQU(0.123)
            WS2M  = RANEQU(0.123)
            WS2P  = RANEQU(0.123)
C                                 ERMITTLUNG EINES FAKTORS, DER DIE
C                                 SCHAETZWERTE D. SCHWANKUNGSGESCHW. AN
C                                 DIE GEFORDERTE TURBULENZENERGIE AN-
C                                 PASST
C
            FAKM  = 2.0*TKMIN/(US2M+VS2M+WS2M)
            FAKP  = 2.0*TKPLU/(US2P+VS2P+WS2P)
C
            UFR(NMIN,J,I) = UFR(NMIN,J,I)+SQRT(FAKM*US2M)*VORZ(N)
            UFR(NPLU,J,I) = UFR(NPLU,J,I)+SQRT(FAKP*US2P)*VORZ(N+1)
  190    CONTINUE
  180 CONTINUE
  170 CONTINUE

      ENDIF
C
      ELSEIF(CIDUFR(1:7).EQ.'BLASIUS') THEN
C
      DO I=1,2
	 DO K=3,KM1
         ETAN = 3.5 * Z(K)/DELTA
C
C                       UMRECHNUG VON Z(K) ZUR ETAN 
C
C
      DO N=1,40
      IF (ETAN.GE.ETA(N).AND.ETAN.LT.ETA(N+1)) THEN
      FN = (F(N+1)-F(N)) * (ETAN - ETA(N))/(ETA(N+1)-ETA(N)) + F(N)
      DFN = (DF(N+1)-DF(N)) * (ETAN - ETA(N))/(ETA(N+1)-ETA(N)) + DF(N)
C                      INTERPOLATION FN UND DFN	  
      ENDIF
      ENDDO
        DO J=3,JM2
          UFR(K,J,I) = DFN * UFRCON
	  VFR(K,J,I) = 0.0
	  WFR(K,J,I) = 0.0
C	  WFR(K,J,I) = (ETAN*DFN-FN)*3.5*GMOL/DELTA
        ENDDO
	 ENDDO
      ENDDO


C                         EXPERIMENTAL PROFILE IS SET AS INLET

      ELSEIF(CIDUFR(1:8).EQ.'EXPERIME') THEN
      WRITE(18,*) 'EXPERIMENT'
	 ZSCA=0.025
	 USCA=0.9
	 ZEXP=0.0
	 UEXP=0.0
	 REWIND(16)
	 READ(16,*) DUMMY,ZEXPNEU,DUMMY,UEXPNEU
	 ZEXPNEU=ZEXPNEU*ZSCA
	 UEXPNEU=UEXPNEU*USCA
	 WRITE(18,*) ZEXPNEU, UEXPNEU
      DO  K=3,KM1

  333    IF (ZEXPNEU.LE.Z(K)) THEN
	 ZEXP=ZEXPNEU
	 UEXP=UEXPNEU
	 READ(16,*) DUMMY,ZEXPNEU,DUMMY,UEXPNEU
	 ZEXPNEU=ZEXPNEU*ZSCA
	 UEXPNEU=UEXPNEU*USCA
	 WRITE(18,*) ZEXPNEU, UEXPNEU
	 GOTO 333
	 ENDIF

      DO  I=1,2
      DO  J=3,JM2
	 UFR(K,J,I) = UEXP+(Z(K)-ZEXP)*(UEXPNEU-UEXP)/(ZEXPNEU-ZEXP)
      VFR(K,J,I) = 0.0
      WFR(K,J,I) = 0.0

      ENDDO
      ENDDO
      WRITE(18,*) 'Gitterpunktwerte:','Z(',K,')=',Z(K),
     &            'U=',UFR(K,3,2)
      ENDDO

      RETURN



       ELSEIF(CIDUFR(1:10).EQ.'SPALART300') THEN
C
C                                   MITTLERES PROFIL VON SPALART (1987)
C
          CALL SPALART300 (KK,JJ,2,UFR,Z,ZBANF,DELTA,UFRCON)
C
C

      RETURN
C
       ELSEIF(CIDUFR(1:10).EQ.'SPALART670') THEN
C
C                                   MITTLERES PROFIL VON SPALART (1987)
C
          CALL SPALART670 (KK,JJ,2,UFR,Z,ZBANF,DELTA,UFRCON)
C
C

      RETURN
C
       ELSEIF(CIDUFR(1:11).EQ.'SPALART1410') THEN
C
C                                   MITTLERES PROFIL VON SPALART (1987)
C
          CALL SPALART1410 (KK,JJ,2,UFR,Z,ZBANF,DELTA,UFRCON)
C
C

      RETURN
C
      ELSEIF(CIDUFR(1:8).EQ.'BOUNDARY') THEN

C
C                                 BEI LAMINARER STROEMUNG WIRD EIN PARA-
C                                 BOLISCHES GESCHWINDIGKEITSPROFIL ER-
C                                 ZEUGT.
C
      IF(MTURB .EQ. 0) THEN
      EXPON = 2.0
C
C
      ZBOT   = 0.5*(Z(3) + Z(2))
C
      DO 200 I=1,2
      DO 210 J=3,JM2
         DO 220 K=3,KM1
C
          ZDELTA=(DELTA+ZBOT-Z(K))/DELTA
          Verd=MAX(0.0,ZDELTA)**EXPON
          UFR(K,J,I)= VREF*(1.0 - Verd)
C
  220    CONTINUE
  210 CONTINUE
  200 CONTINUE
C
C
C                             UM BEI GEBIETSZERLEGUNG IN Y-RICHTUNG KEINE
C                             SCHWIERIGKEITEN ZU BEKOMMEN:
      RETURN
      ELSE
C                             HIER TURBULENTE GRENZSCHICHT
C
      ZBOT = 0.5*(Z(3) + Z(2))
C     ZBOT = ZBANF
C     
      DO I=1,2
      DO K=3,KM1
        IF(Z(K).LT.ZBOT) THEN
         fak =0.0
        ELSE
         fak=(MIN(1.0,(Z(k)-ZBOT)/DELTA))**EXPON*VREF
        ENDIF
         DO J=2,JM1
C    
               UFR(K,J,I) = fak
C	 UFR(K,J,I) = (MIN(1.0,(Z(k)-ZBOT)/DELTA))**EXPON*VREF
	 ENDDO
      ENDDO
      ENDDO

      RETURN

C
C
C                                  *************************************
C                                  ERMITTLUNG DER SCHWANKUNGSGESCHWIN-
C                                  DIGKEITEN MITTELS ZUFALLSZAHLENGENE-
C                                  RATOR. DIE SCHWANKUNGEN WERDEN SO BE-
C                                  STIMMT, DASS EINE VORGEGEBENE TURBU-
C                                  LENZENERGIEVERTEILUNG EINGEHALTEN
C                                  WIRD.
C                                  *************************************
C
C                                  KONSTANTEN DER K-VERTEILUNG (AN-
C                                  NAEHERUNG EINER MIT DEM K-EPS-MODELL
C                                  BERECHNETEN K-VERTEILUNG)
C
      CK1    = 0.0015
      CK2    = 0.017
      CK3    = 1.443
C
C                                  NAEHERUNGSWEISE BEST. DER DICKE DER
C                                  VISKOSEN UNTERSCHICHT
C
      HK2 = DELTA
      DYV    = (HK2**EXPON * CM25RE*GMOL*LOG(ECONST*CM25RE)
     $       /  (CAPPA*RHO*(1.0+EXPON)))**(1.0/(1.0+EXPON))
      TKV    = TKEIN(HK2-DYV)
      GRADKV = TKV/DYV
C
      DO 270 I=1,2
      DO 280 J=3,JM2
         DO 290 N=3,KM1
	    IF (Z(N).GT.DELTA) GOTO 280
            ZMMIN        = DELTA - Z(N)
            TKMIN        = TKEIN(ZMMIN)
            IF((HK2-ZMMIN) .LE. DYV) TKMIN = GRADKV*(HK2-ZMMIN)
C
C                                 ERSTE SCHAETZUNG FUER DIE QUADRATE DER
C                                 SCHWANKUNGSGESCHWINDIGKEITEN
C
            US2M  = RANEQU(0.123)
            VS2M  = RANEQU(0.123)
            WS2M  = RANEQU(0.123)
C                                 ERMITTLUNG EINES FAKTORS, DER DIE
C                                 SCHAETZWERTE D. SCHWANKUNGSGESCHW. AN
C                                 DIE GEFORDERTE TURBULENZENERGIE AN-
C                                 PASST
C
            FAKM  = 2.0*TKMIN/(US2M+VS2M+WS2M)
C
            UFR(N,J,I) = UFR(N,J,I)+SQRT(FAKM*US2M)*VORZ(N)
  290    CONTINUE
  280 CONTINUE
  270 CONTINUE
      ENDIF
C
      ELSEIF(CIDUFR(1:4).EQ.'DUCT') THEN
C
      ZBOT   = 0.5*(Z(3) + Z(2))
      ZTOP   = 0.5*(Z(KM2) + Z(KM1))
      YBOT   = 0.5*(Y(3) + Y(2))
      YTOP   = 0.5*(Y(JM2) + Y(JM1))
C
C
C                                 HALBE KANALHOEHE
C
      HKZ    = 0.5*(ZTOP-ZBOT)
      HKY    = 0.5*(YTOP-YBOT)
      HK2    = 0.5*(HKY + HKZ)
      RHK    = 1./(HKZ**2 * HKY**2)
C
C                                 BEI LAMINARER STROEMUNG WIRD EIN PARA-
C                                 BOLISCHES GESCHWINDIGKEITSPROFIL ER-
C                                 ZEUGT.
C
      IF(MTURB .EQ. 0) THEN

          DO I=1,2
           DO J=3,JM2
            DO K=3,KM1

               UFR(K,J,I) = UFRCON*
     $             (Z(K)-ZBOT)*(Y(J)-YBOT)*(ZTOP-Z(K))*(YTOP-Y(J))*
     $             AMPLITUDE * RHK

            ENDDO
           ENDDO
          ENDDO
          RETURN
C
      ELSE
C
      DO 300 I=1,2
      DO 310 J=3,JM2
         DO 320 K=3,KM1
C
            UFR(K,J,I) = (MIN(1.0              ,
     +                     (Z(K)-ZBOT)/HKZ,
     +                     (ZTOP-Z(K))/HKZ,
     +                     (Y(J)-YBOT)/HKY,
     +                     (YTOP-Y(J))/HKY))**EXPON*VREF
                           
C
  320    CONTINUE
  310 CONTINUE
  300 CONTINUE
C
      ENDIF
C
      ELSEIF(CIDUFR(1:5).EQ.'PDE1D') THEN
C
          DO I=1,2
           DO J=3,JM2
            DO K=3,KM1

               VFR(K,J,I) = UFRCON*AMPLITUDE
               UFR(K,J,I) = 1.0

            ENDDO
           ENDDO
          ENDDO
          RETURN
C
       ELSEIF(CIDUFR(1:3).EQ.'BLA') THEN
      RETURN
C
       ELSEIF(CIDUFR(1:8).EQ.'ROUNDJET') THEN
C
C                                   RUNDER JET UM DIE X-ACHSE
C                                          RADIUS 0.625
      DO I=1,2
      DO J=2,JM1
      DO K=2,KM1
	 R = SQRT( Y(J)**2 + Z(K)**2)
	 IF ( R .LE. 0.625 ) THEN
            UFR(K,J,I) = UFRCON
	 ELSE 
	    UFR(K,J,I) = 0.0
	 ENDIF
	 VFR(K,J,I) = 0.0
	 WFR(K,J,I) = 0.0
      END DO
      END DO
      END DO
C
      RETURN	
CCC
      ELSEIF(CIDUFR(1:8).EQ.'INBOTTOM') THEN
C
C      DO I=2,IM1
C      DO J=2,JM1
C      DO K=1,2
C         UFR(K,J,I) = 0.0
C         VFR(K,J,I) = 0.0
C         WFR(K,J,I) = UFRCON
C      END DO
C      END DO
C      END DO
C
      RETURN
C
CCC
C
      ELSE 
          CALL ERRR (745,'SVEIPR ')
C
      ENDIF
C
      RETURN
      END
      SUBROUTINE SPALART300(KK,JJ,II,UFR,Z,ZBANF,DELTA,UFRCON)
C---------------------------------------------------------------------72
C
C                      BELEGT DAS FELD UFRCON MIT DEM MITTLEREM
C                      PROFIL AUS SPALART (1987), RE_DELTA_2=300
C
C    22.11.96     (MM)
C---------------------------------------------------------------------72

      REAL UFR(KK,JJ,II),Z(KK)

      REAL ZWAND(48),UMEAN(48)


      DATA    ZWAND(  1)  /1.5089E-03/
      DATA    ZWAND(  2)  /3.7109E-03/
      DATA    ZWAND(  3)  /6.8970E-03/
      DATA    ZWAND(  4)  /1.1073E-02/
      DATA    ZWAND(  5)  /1.6248E-02/
      DATA    ZWAND(  6)  /2.2431E-02/
      DATA    ZWAND(  7)  /2.9635E-02/
      DATA    ZWAND(  8)  /3.7874E-02/
      DATA    ZWAND(  9)  /4.7164E-02/
      DATA    ZWAND( 10)  /5.7526E-02/
      DATA    ZWAND( 11)  /6.8980E-02/
      DATA    ZWAND( 12)  /8.1551E-02/
      DATA    ZWAND( 13)  /9.5267E-02/
      DATA    ZWAND( 14)  /1.1016E-01/
      DATA    ZWAND( 15)  /1.2626E-01/
      DATA    ZWAND( 16)  /1.4361E-01/
      DATA    ZWAND( 17)  /1.6226E-01/
      DATA    ZWAND( 18)  /1.8224E-01/
      DATA    ZWAND( 19)  /2.0361E-01/
      DATA    ZWAND( 20)  /2.2644E-01/
      DATA    ZWAND( 21)  /2.5078E-01/
      DATA    ZWAND( 22)  /2.7672E-01/
      DATA    ZWAND( 23)  /3.0432E-01/
      DATA    ZWAND( 24)  /3.3370E-01/
      DATA    ZWAND( 25)  /3.6494E-01/
      DATA    ZWAND( 26)  /3.9817E-01/
      DATA    ZWAND( 27)  /4.3352E-01/
      DATA    ZWAND( 28)  /4.7113E-01/
      DATA    ZWAND( 29)  /5.1119E-01/
      DATA    ZWAND( 30)  /5.5389E-01/
      DATA    ZWAND( 31)  /5.9946E-01/
      DATA    ZWAND( 32)  /6.4816E-01/
      DATA    ZWAND( 33)  /7.0030E-01/
      DATA    ZWAND( 34)  /7.5627E-01/
      DATA    ZWAND( 35)  /8.1648E-01/
      DATA    ZWAND( 36)  /8.8148E-01/
      DATA    ZWAND( 37)  /9.5190E-01/
      DATA    ZWAND( 38)  /1.0286E+00/
      DATA    ZWAND( 39)  /1.1124E+00/
      DATA    ZWAND( 40)  /1.2048E+00/
      DATA    ZWAND( 41)  /1.3074E+00/
      DATA    ZWAND( 42)  /1.4225E+00/
      DATA    ZWAND( 43)  /1.5532E+00/
      DATA    ZWAND( 44)  /1.7041E+00/
      DATA    ZWAND( 45)  /1.8822E+00/
      DATA    ZWAND( 46)  /2.0990E+00/
      DATA    ZWAND( 47)  /2.3756E+00/
      DATA    ZWAND( 48)  /2.7574E+00/

C---------------------------------------------------------------------72


      DATA    UMEAN(  1)  /0.0124532/
      DATA    UMEAN(  2)  /0.0306272/
      DATA    UMEAN(  3)  /0.0569061/
      DATA    UMEAN(  4)  /0.0912614/
      DATA    UMEAN(  5)  /0.133475/
      DATA    UMEAN(  6)  /0.18289/
      DATA    UMEAN(  7)  /0.238096/
      DATA    UMEAN(  8)  /0.296842/
      DATA    UMEAN(  9)  /0.356312/
      DATA    UMEAN( 10)  /0.413774/
      DATA    UMEAN( 11)  /0.467129/
      DATA    UMEAN( 12)  /0.515216/
      DATA    UMEAN( 13)  /0.557727/
      DATA    UMEAN( 14)  /0.59497/
      DATA    UMEAN( 15)  /0.627517/
      DATA    UMEAN( 16)  /0.656016/
      DATA    UMEAN( 17)  /0.681222/
      DATA    UMEAN( 18)  /0.703676/
      DATA    UMEAN( 19)  /0.724024/
      DATA    UMEAN( 20)  /0.742754/
      DATA    UMEAN( 21)  /0.760296/
      DATA    UMEAN( 22)  /0.777082/
      DATA    UMEAN( 23)  /0.793491/
      DATA    UMEAN( 24)  /0.809629/
      DATA    UMEAN( 25)  /0.825768/
      DATA    UMEAN( 26)  /0.842014/
      DATA    UMEAN( 27)  /0.858369/
      DATA    UMEAN( 28)  /0.874939/
      DATA    UMEAN( 29)  /0.891672/
      DATA    UMEAN( 30)  /0.908512/
      DATA    UMEAN( 31)  /0.92519/
      DATA    UMEAN( 32)  /0.941329/
      DATA    UMEAN( 33)  /0.956226/
      DATA    UMEAN( 34)  /0.969558/
      DATA    UMEAN( 35)  /0.980677/
      DATA    UMEAN( 36)  /0.989097/
      DATA    UMEAN( 37)  /0.994548/
      DATA    UMEAN( 38)  /0.997625/
      DATA    UMEAN( 39)  /0.99919/
      DATA    UMEAN( 40)  /0.999784/
      DATA    UMEAN( 41)  /0.999946/
      DATA    UMEAN( 42)  /1/
      DATA    UMEAN( 43)  /1/
      DATA    UMEAN( 44)  /1/
      DATA    UMEAN( 45)  /1/
      DATA    UMEAN( 46)  /1/
      DATA    UMEAN( 47)  /1/
      DATA    UMEAN( 48)  /1/

C---------------------------------------------------------------------72

      DO K=1,KK

         DISTANCE = (Z(K)-ZBANF)/DELTA

C----------------------------------------- INTERPOLATION TO GRID

         IF (DISTANCE .GT. 0.0 ) THEN

            DO K2=1,47
               
               IF ((DISTANCE-ZWAND(K2+1)) .LE. 0.0) THEN

                  INDEX = K2

                  UINTER = UFRCON*(UMEAN(K2) + 
     $                 ((UMEAN(K2+1)-UMEAN(K2))* 
     $                  (DISTANCE   -ZWAND(K2))/ 
     $                  (ZWAND(K2+1)-ZWAND(K2))))

                  GOTO 100

               ENDIF
               
            ENDDO

  100       CONTINUE

C--------------------------------------- BELEGEN DER WERTE

            DO I=1,II
               DO J=1,JJ
                  
                  UFR(K,J,I) = UINTER
                  
               ENDDO
            ENDDO

         endif
      enddo

C-------------------------------------- READY

      RETURN
      END
      SUBROUTINE SPALART670(KK,JJ,II,UFR,Z,ZBANF,DELTA,UFRCON)
C---------------------------------------------------------------------72
C
C                      BELEGT DAS FELD UFRCON MIT DEM MITTLEREM
C                      PROFIL AUS SPALART (1987), RE_DELTA_2=670
C
C    22.11.96     (MM)
C---------------------------------------------------------------------72

      REAL UFR(KK,JJ,II),Z(KK)

      REAL ZWAND(62),UMEAN(62)



      DATA    ZWAND(1)  /7.9361E-04/
      DATA    ZWAND(2)  /1.9512E-03/
      DATA    ZWAND(3)  /3.6251E-03/
      DATA    ZWAND(4)  /5.8172E-03/
      DATA    ZWAND(5)  /8.5301E-03/
      DATA    ZWAND(6)  /1.1767E-02/
      DATA    ZWAND(7)  /1.5532E-02/
      DATA    ZWAND(8)  /1.9829E-02/
      DATA    ZWAND(9)  /2.4664E-02/
      DATA    ZWAND(10)  /3.0043E-02/
      DATA    ZWAND(11)  /3.5972E-02/
      DATA    ZWAND(12)  /4.2459E-02/
      DATA    ZWAND(13)  /4.9512E-02/
      DATA    ZWAND(14)  /5.7140E-02/
      DATA    ZWAND(15)  /6.5353E-02/
      DATA    ZWAND(16)  /7.4163E-02/
      DATA    ZWAND(17)  /8.3581E-02/
      DATA    ZWAND(18)  /9.3621E-02/
      DATA    ZWAND(19)  /1.0430E-01/
      DATA    ZWAND(20)  /1.1562E-01/
      DATA    ZWAND(21)  /1.2762E-01/
      DATA    ZWAND(22)  /1.4030E-01/
      DATA    ZWAND(23)  /1.5369E-01/
      DATA    ZWAND(24)  /1.6781E-01/
      DATA    ZWAND(25)  /1.8268E-01/
      DATA    ZWAND(26)  /1.9833E-01/
      DATA    ZWAND(27)  /2.1479E-01/
      DATA    ZWAND(28)  /2.3208E-01/
      DATA    ZWAND(29)  /2.5025E-01/
      DATA    ZWAND(30)  /2.6932E-01/
      DATA    ZWAND(31)  /2.8935E-01/
      DATA    ZWAND(32)  /3.1036E-01/
      DATA    ZWAND(33)  /3.3242E-01/
      DATA    ZWAND(34)  /3.5557E-01/
      DATA    ZWAND(35)  /3.7987E-01/
      DATA    ZWAND(36)  /4.0539E-01/
      DATA    ZWAND(37)  /4.3221E-01/
      DATA    ZWAND(38)  /4.6040E-01/
      DATA    ZWAND(39)  /4.9005E-01/
      DATA    ZWAND(40)  /5.2127E-01/
      DATA    ZWAND(41)  /5.5418E-01/
      DATA    ZWAND(42)  /5.8889E-01/
      DATA    ZWAND(43)  /6.2557E-01/
      DATA    ZWAND(44)  /6.6438E-01/
      DATA    ZWAND(45)  /7.0552E-01/
      DATA    ZWAND(46)  /7.4921E-01/
      DATA    ZWAND(47)  /7.9574E-01/
      DATA    ZWAND(48)  /8.4540E-01/
      DATA    ZWAND(49)  /8.9859E-01/
      DATA    ZWAND(50)  /9.5575E-01/
      DATA    ZWAND(51)  /1.0174E+00/
      DATA    ZWAND(52)  /1.0843E+00/
      DATA    ZWAND(53)  /1.1573E+00/
      DATA    ZWAND(54)  /1.2374E+00/
      DATA    ZWAND(55)  /1.3261E+00/
      DATA    ZWAND(56)  /1.4254E+00/
      DATA    ZWAND(57)  /1.5379E+00/
      DATA    ZWAND(58)  /1.6676E+00/
      DATA    ZWAND(59)  /1.8204E+00/
      DATA    ZWAND(60)  /2.0062E+00/
      DATA    ZWAND(61)  /2.2429E+00/
      DATA    ZWAND(62)  /2.5695E+00/
C---------------------------------------------------------------------72
      DATA    UMEAN(1)  /0.0127183/
      DATA    UMEAN(2)  /0.0312694/
      DATA    UMEAN(3)  /0.0580678/
      DATA    UMEAN(4)  /0.0930285/
      DATA    UMEAN(5)  /0.135782/
      DATA    UMEAN(6)  /0.185334/
      DATA    UMEAN(7)  /0.239828/
      DATA    UMEAN(8)  /0.296554/
      DATA    UMEAN(9)  /0.352522/
      DATA    UMEAN(10)  /0.405172/
      DATA    UMEAN(11)  /0.452842/
      DATA    UMEAN(12)  /0.494838/
      DATA    UMEAN(13)  /0.531268/
      DATA    UMEAN(14)  /0.562635/
      DATA    UMEAN(15)  /0.589626/
      DATA    UMEAN(16)  /0.612979/
      DATA    UMEAN(17)  /0.633284/
      DATA    UMEAN(18)  /0.651131/
      DATA    UMEAN(19)  /0.66706/
      DATA    UMEAN(20)  /0.681416/
      DATA    UMEAN(21)  /0.694592/
      DATA    UMEAN(22)  /0.706883/
      DATA    UMEAN(23)  /0.718486/
      DATA    UMEAN(24)  /0.729597/
      DATA    UMEAN(25)  /0.740315/
      DATA    UMEAN(26)  /0.750836/
      DATA    UMEAN(27)  /0.761308/
      DATA    UMEAN(28)  /0.771731/
      DATA    UMEAN(29)  /0.782203/
      DATA    UMEAN(30)  /0.792625/
      DATA    UMEAN(31)  /0.803196/
      DATA    UMEAN(32)  /0.813864/
      DATA    UMEAN(33)  /0.824582/
      DATA    UMEAN(34)  /0.835447/
      DATA    UMEAN(35)  /0.846559/
      DATA    UMEAN(36)  /0.857866/
      DATA    UMEAN(37)  /0.86942/
      DATA    UMEAN(38)  /0.881318/
      DATA    UMEAN(39)  /0.893363/
      DATA    UMEAN(40)  /0.905556/
      DATA    UMEAN(41)  /0.917748/
      DATA    UMEAN(42)  /0.929892/
      DATA    UMEAN(43)  /0.94174/
      DATA    UMEAN(44)  /0.953245/
      DATA    UMEAN(45)  /0.96411/
      DATA    UMEAN(46)  /0.973845/
      DATA    UMEAN(47)  /0.982055/
      DATA    UMEAN(48)  /0.988545/
      DATA    UMEAN(49)  /0.993314/
      DATA    UMEAN(50)  /0.99646/
      DATA    UMEAN(51)  /0.998328/
      DATA    UMEAN(52)  /0.999263/
      DATA    UMEAN(53)  /0.999705/
      DATA    UMEAN(54)  /0.999902/
      DATA    UMEAN(55)  /0.999951/
      DATA    UMEAN(56)  /1/
      DATA    UMEAN(57)  /1/
      DATA    UMEAN(58)  /1/
      DATA    UMEAN(59)  /1/
      DATA    UMEAN(60)  /1/
      DATA    UMEAN(61)  /1/
      DATA    UMEAN(62)  /1/

C---------------------------------------------------------------------72

      DO K=1,KK

         DISTANCE = (Z(K)-ZBANF)/DELTA

C----------------------------------------- INTERPOLATION TO GRID

         IF (DISTANCE .GT. 0.0 ) THEN

            DO K2=1,61
               
               IF ((DISTANCE-ZWAND(K2+1)) .LE. 0.0) THEN

                  INDEX = K2

                  UINTER = UFRCON*(UMEAN(K2) + 
     $                 ((UMEAN(K2+1)-UMEAN(K2))* 
     $                  (DISTANCE   -ZWAND(K2))/ 
     $                  (ZWAND(K2+1)-ZWAND(K2))))

                  GOTO 100

               ENDIF
               
            ENDDO

  100       CONTINUE

C--------------------------------------- BELEGEN DER WERTE

            DO I=1,II
               DO J=1,JJ
                  
                  UFR(K,J,I) = UINTER
                  
               ENDDO
            ENDDO

         endif
      enddo

C-------------------------------------- READY

      RETURN
      END
      SUBROUTINE SPALART1410(KK,JJ,II,UFR,Z,ZBANF,DELTA,UFRCON)
C---------------------------------------------------------------------72
C
C                      BELEGT DAS FELD UFRCON MIT DEM MITTLEREM
C                      PROFIL AUS SPALART (1987), RE_DELTA_2=300
C
C    22.11.96     (MM)
C---------------------------------------------------------------------72

      REAL UFR(KK,JJ,II),Z(KK)

      REAL ZWAND(78),UMEAN(78)

      DATA    ZWAND( 1 )  /4.7397E-04/
      DATA    ZWAND( 2 )  /1.1652E-03/
      DATA    ZWAND( 3 )  /2.1642E-03/
      DATA    ZWAND( 4 )  /3.4719E-03/
      DATA    ZWAND( 5 )  /5.0892E-03/
      DATA    ZWAND( 6 )  /7.0174E-03/
      DATA    ZWAND( 7 )  /9.2578E-03/
      DATA    ZWAND( 8 )  /1.1812E-02/
      DATA    ZWAND( 9 )  /1.4683E-02/
      DATA    ZWAND( 10 )  /1.7872E-02/
      DATA    ZWAND( 11 )  /2.1381E-02/
      DATA    ZWAND( 12 )  /2.5214E-02/
      DATA    ZWAND( 13 )  /2.9374E-02/
      DATA    ZWAND( 14 )  /3.3863E-02/
      DATA    ZWAND( 15 )  /3.8686E-02/
      DATA    ZWAND( 16 )  /4.3847E-02/
      DATA    ZWAND( 17 )  /4.9349E-02/
      DATA    ZWAND( 18 )  /5.5198E-02/
      DATA    ZWAND( 19 )  /6.1399E-02/
      DATA    ZWAND( 20 )  /6.7956E-02/
      DATA    ZWAND( 21 )  /7.4876E-02/
      DATA    ZWAND( 22 )  /8.2164E-02/
      DATA    ZWAND( 23 )  /8.9829E-02/
      DATA    ZWAND( 24 )  /9.7875E-02/
      DATA    ZWAND( 25 )  /1.0631E-01/
      DATA    ZWAND( 26 )  /1.1515E-01/
      DATA    ZWAND( 27 )  /1.2439E-01/
      DATA    ZWAND( 28 )  /1.3405E-01/
      DATA    ZWAND( 29 )  /1.4413E-01/
      DATA    ZWAND( 30 )  /1.5466E-01/
      DATA    ZWAND( 31 )  /1.6563E-01/
      DATA    ZWAND( 32 )  /1.7706E-01/
      DATA    ZWAND( 33 )  /1.8896E-01/
      DATA    ZWAND( 34 )  /2.0136E-01/
      DATA    ZWAND( 35 )  /2.1425E-01/
      DATA    ZWAND( 36 )  /2.2767E-01/
      DATA    ZWAND( 37 )  /2.4162E-01/
      DATA    ZWAND( 38 )  /2.5613E-01/
      DATA    ZWAND( 39 )  /2.7121E-01/
      DATA    ZWAND( 40 )  /2.8689E-01/
      DATA    ZWAND( 41 )  /3.0319E-01/
      DATA    ZWAND( 42 )  /3.2013E-01/
      DATA    ZWAND( 43 )  /3.3775E-01/
      DATA    ZWAND( 44 )  /3.5607E-01/
      DATA    ZWAND( 45 )  /3.7512E-01/
      DATA    ZWAND( 46 )  /3.9495E-01/
      DATA    ZWAND( 47 )  /4.1558E-01/
      DATA    ZWAND( 48 )  /4.3706E-01/
      DATA    ZWAND( 49 )  /4.5944E-01/
      DATA    ZWAND( 50 )  /4.8277E-01/
      DATA    ZWAND( 51 )  /5.0709E-01/
      DATA    ZWAND( 52 )  /5.3249E-01/
      DATA    ZWAND( 53 )  /5.5901E-01/
      DATA    ZWAND( 54 )  /5.8675E-01/
      DATA    ZWAND( 55 )  /6.1577E-01/
      DATA    ZWAND( 56 )  /6.4619E-01/
      DATA    ZWAND( 57 )  /6.7810E-01/
      DATA    ZWAND( 58 )  /7.1163E-01/
      DATA    ZWAND( 59 )  /7.4692E-01/
      DATA    ZWAND( 60 )  /7.8413E-01/
      DATA    ZWAND( 61 )  /8.2344E-01/
      DATA    ZWAND( 62 )  /8.6505E-01/
      DATA    ZWAND( 63 )  /9.0923E-01/
      DATA    ZWAND( 64 )  /9.5627E-01/
      DATA    ZWAND( 65 )  /1.0065E+00/
      DATA    ZWAND( 66 )  /1.0604E+00/
      DATA    ZWAND( 67 )  /1.1184E+00/
      DATA    ZWAND( 68 )  /1.1811E+00/
      DATA    ZWAND( 69 )  /1.2495E+00/
      DATA    ZWAND( 70 )  /1.3244E+00/
      DATA    ZWAND( 71 )  /1.4073E+00/
      DATA    ZWAND( 72 )  /1.4999E+00/
      DATA    ZWAND( 73 )  /1.6047E+00/
      DATA    ZWAND( 74 )  /1.7253E+00/
      DATA    ZWAND( 75 )  /1.8674E+00/
      DATA    ZWAND( 76 )  /2.0400E+00/
      DATA    ZWAND( 77 )  /2.2598E+00/
      DATA    ZWAND( 78 )  /2.5629E+00/

C---------------------------------------------------------------------72

      DATA    UMEAN( 1 )  /0.0139919/
      DATA    UMEAN( 2 )  /0.034394/
      DATA    UMEAN( 3 )  /0.0638296/
      DATA    UMEAN( 4 )  /0.102076/
      DATA    UMEAN( 5 )  /0.148357/
      DATA    UMEAN( 6 )  /0.20084/
      DATA    UMEAN( 7 )  /0.256572/
      DATA    UMEAN( 8 )  /0.312059/
      DATA    UMEAN( 9 )  /0.364233/
      DATA    UMEAN( 10 )  /0.411168/
      DATA    UMEAN( 11 )  /0.452119/
      DATA    UMEAN( 12 )  /0.487253/
      DATA    UMEAN( 13 )  /0.517153/
      DATA    UMEAN( 14 )  /0.54269/
      DATA    UMEAN( 15 )  /0.564592/
      DATA    UMEAN( 16 )  /0.58354/
      DATA    UMEAN( 17 )  /0.600172/
      DATA    UMEAN( 18 )  /0.614894/
      DATA    UMEAN( 19 )  /0.628253/
      DATA    UMEAN( 20 )  /0.640431/
      DATA    UMEAN( 21 )  /0.651791/
      DATA    UMEAN( 22 )  /0.662424/
      DATA    UMEAN( 23 )  /0.672557/
      DATA    UMEAN( 24 )  /0.682191/
      DATA    UMEAN( 25 )  /0.69146/
      DATA    UMEAN( 26 )  /0.700503/
      DATA    UMEAN( 27 )  /0.709228/
      DATA    UMEAN( 28 )  /0.717725/
      DATA    UMEAN( 29 )  /0.725995/
      DATA    UMEAN( 30 )  /0.734083/
      DATA    UMEAN( 31 )  /0.74199/
      DATA    UMEAN( 32 )  /0.749851/
      DATA    UMEAN( 33 )  /0.757621/
      DATA    UMEAN( 34 )  /0.765346/
      DATA    UMEAN( 35 )  /0.773025/
      DATA    UMEAN( 36 )  /0.780705/
      DATA    UMEAN( 37 )  /0.788293/
      DATA    UMEAN( 38 )  /0.795927/
      DATA    UMEAN( 39 )  /0.803561/
      DATA    UMEAN( 40 )  /0.811286/
      DATA    UMEAN( 41 )  /0.819101/
      DATA    UMEAN( 42 )  /0.826917/
      DATA    UMEAN( 43 )  /0.834824/
      DATA    UMEAN( 44 )  /0.842776/
      DATA    UMEAN( 45 )  /0.850819/
      DATA    UMEAN( 46 )  /0.858907/
      DATA    UMEAN( 47 )  /0.867086/
      DATA    UMEAN( 48 )  /0.875447/
      DATA    UMEAN( 49 )  /0.884126/
      DATA    UMEAN( 50 )  /0.892987/
      DATA    UMEAN( 51 )  /0.901984/
      DATA    UMEAN( 52 )  /0.911208/
      DATA    UMEAN( 53 )  /0.920569/
      DATA    UMEAN( 54 )  /0.929975/
      DATA    UMEAN( 55 )  /0.939336/
      DATA    UMEAN( 56 )  /0.948469/
      DATA    UMEAN( 57 )  /0.957239/
      DATA    UMEAN( 58 )  /0.9656/
      DATA    UMEAN( 59 )  /0.973188/
      DATA    UMEAN( 60 )  /0.979868/
      DATA    UMEAN( 61 )  /0.985594/
      DATA    UMEAN( 62 )  /0.990228/
      DATA    UMEAN( 63 )  /0.993773/
      DATA    UMEAN( 64 )  /0.996272/
      DATA    UMEAN( 65 )  /0.997908/
      DATA    UMEAN( 66 )  /0.998908/
      DATA    UMEAN( 67 )  /0.999544/
      DATA    UMEAN( 68 )  /0.999862/
      DATA    UMEAN( 69 )  /0.999953/
      DATA    UMEAN( 70 )  /0.999953/
      DATA    UMEAN( 71 )  /0.999953/
      DATA    UMEAN( 72 )  /0.999953/
      DATA    UMEAN( 73 )  /0.999953/
      DATA    UMEAN( 74 )  /0.999953/
      DATA    UMEAN( 75 )  /0.999998/
      DATA    UMEAN( 76 )  /0.999998/
      DATA    UMEAN( 77 )  /1.00004/
      DATA    UMEAN( 78 )  /0.999998/


C---------------------------------------------------------------------72

      DO K=1,KK

         DISTANCE = (Z(K)-ZBANF)/DELTA

C----------------------------------------- INTERPOLATION TO GRID

         IF (DISTANCE .GT. 0.0 ) THEN

            DO K2=1,77
               
               IF ((DISTANCE-ZWAND(K2+1)) .LE. 0.0) THEN

                  INDEX = K2

                  UINTER = UFRCON*(UMEAN(K2) + 
     $                 ((UMEAN(K2+1)-UMEAN(K2))* 
     $                  (DISTANCE   -ZWAND(K2))/ 
     $                  (ZWAND(K2+1)-ZWAND(K2))))

                  GOTO 100

               ENDIF
               
            ENDDO

  100       CONTINUE

C--------------------------------------- BELEGEN DER WERTE

            DO I=1,II
               DO J=1,JJ
                  
                  UFR(K,J,I) = UINTER
                  
               ENDDO
            ENDDO

         endif
      enddo

C-------------------------------------- READY

      RETURN
      END
