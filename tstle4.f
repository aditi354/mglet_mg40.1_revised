










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
      SUBROUTINE TSTLE4  (KK,JJ,II,KMX,JMX,IMX,
     $                    X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                    RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                    U,V,W,UO,VO,WO,
     $                    WCU,WCV,WCW,P,G,B,
     $                    BP,BU,BV,BW,
     $                    DT,GRADPX,ZTOT,NBUF,ITSTEP,
     $                    WPHI,WKON,WDIF,WSOR,IDUZ,
     $                    IC1,IC2,JC1,JC2,KC1,KC2,
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,
     $                    RSGS3,FAKTOR,
     $                    FUI,FVI,FWI,
     $                    FUJ,FVJ,FWJ,
     $                    FUK,FVK,FWK,
     $                    COEFFX,COEFFY,COEFFZ,COEFDX,COEFDY,COEFDZ,
     $                    LCOL,DIAG,RCOL,UZ,RSP,
     $                    GI,GJ,GK
CC#if defined _PREPROC_
C     $        ,HPX
C#endif
     $                    )

C*STARLET***************************************************************
C        T S T L E 1      TIMESTEP - TURBULENT - ALLE VERFAHREN - LES
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
C        FUI,FVI,FWI    - ERGEBNIS (INT/ABL) IN I-RICHTUNG           
C        FUJ,FVJ,FWJ    - ERGEBNIS (INT/ABL) IN J-RICHTUNG           
C        FUK,FVK,FWK    - ERGEBNIS (INT/ABL) IN K-RICHTUNG           
C
C                  X    - X-KOORDINATE DES ZELLMITTELPUNKTES            
C
C DEFINE DIREKTIVEN     : ZEN, UPW, QUD, EULER, LEAPF, ADBA,
C                         FRPER, RIPER
C UPROG                 : SWCLE1
C
C      10.03.03 (SE,FS) : Original von TSTLE2 abgeleitet
C      30.06.04 (FS)    : von TSTLE3 abgeleitet, fuer RK
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


      INTEGER    KK,JJ,II,KMX,JMX,IMX,K,J,I,NBUF,ITSTEP,
     $           IC1,IC2,JC1,JC2,KC1,KC2,
     $           KSTART,KSTOP,JSTART,JSTOP,ISTART,ISTOP,
     $           KSTARTM,KSTOPM,
     $           NBU,NFU,NRV,NBW,NTW,NLV,
     $           NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,
     $           NXPER,NYPER,NZPER,
     $           IDUZ,ICOM,IST8,IEND8,MTURB

      REAL        X(II),         Y(JJ),         Z(KK),
     &           DX(II),        DY(JJ),        DZ(KK),
     &          DDX(II),       DDY(JJ),       DDZ(KK),
     &        RDX(II),RDY(JJ),RDZ(KK),RDDX(II),RDDY(JJ),RDDZ(KK)

      REAL       U  (KK,JJ,II),V  (KK,JJ,II),W  (KK,JJ,II),
     $           P  (KK,JJ,II),G  (KK,JJ,II),B  (KK,JJ,II),
     $           UO (KK,JJ,II),VO (KK,JJ,II),WO (KK,JJ,II),
     $           GI (KK,JJ,II),GJ (KK,JJ,II),GK (KK,JJ,II),
     $           WCU(KK,JJ,II),WCV(KK,JJ,II),WCW(KK,JJ,II)

      REAL  QKUW,QKUE,QKUT,QKUB,QKUN,QKUS,FKDTU,FKDTV,FKDTW,
     $      QKVW,QKVE,QKVT,QKVB,QKVN,QKVS,FDDTU,FDDTV,FDDTW,
     $      QKWW,QKWE,QKWT,QKWB,QKWN,QKWS,
     $      QDUW,QDUE,QDUT,QDUB,QDUN,QDUS,QDUC,
     $      QDVW,QDVE,QDVT,QDVB,QDVN,QDVS,QDVC,
     $      QDWW,QDWE,QDWT,QDWB,QDWN,QDWS,QDWC,
     $      GUW, GUE, GUT, GUB, GUN, GUS,
     $      GVW, GVE, GVT, GVB, GVN, GVS,
     $      GWW, GWE, GWT, GWB, GWN, GWS,
     $      FUW, FUE, FUT, FUB, FUN, FUS, AUX,  AUY,  AUZ,
     $      FVW, FVE, FVT, FVB, FVN, FVS, AVX,  AVY,  AVZ,
     $      FWW, FWE, FWT, FWB, FWN, FWS, AWX,  AWY,  AWZ

      REAL   DXI, DXIM, DXIP, DDXI, DDXIM, DDXIP,
     $       DYJ, DYJM, DYJP, DDYJ, DDYJM, DDYJP,
     $       DZK, DZKM, DZKP, DDZK, DDZKM, DDZKP,DDZKPP,DZKMM,
     $      RDXI,RDXIM,RDXIP,RDDXI,RDDXIM,RDDXIP,
     $      RDYJ,RDYJM,RDYJP,RDDYJ,RDDYJM,RDDYJP,
     $      RDZK,RDZKM,RDZKP,RDDZK,RDDZKM,RDDZKP

      REAL  WSOR,WKON,WPHI,WDIF,HP,HN,WF,ZERONE,DT,ZTOT,DELTA,
     $      GRADPX,PRESET,UFRFREQ,SMAONE,GREAT,SMALL,RINDEF,
     $      XPER,YPER,ZPER,XREF,FFF
      REAL        GRADP(II)


C                             FELDER FUER DIE INTERPOLIERTEN WERTE
C                             FUER DAS KOMPAKTVERFAHREN


      REAL  BU(KK,JJ,II),BV(KK,JJ,II),BW(KK,JJ,II),BP(KK,JJ,II)

      REAL  FUI(1),FVI(1),FWI(1)
      REAL FUJ(1),FVJ(1),FWJ(1)
      REAL FUK(1),FVK(1),FWK(1)
C#if defined _PREPROC_
C      REAL HPX(KK,JJ,II)
C#endif



C                                  *************************************
C                                  ZEN ZEN ZEN ZEN ZEN ZEN ZEN ZEN ZEN
C                                  *************************************

              WF(FFF)  = 0.5*ABS(FFF)+MIN(FFF,0.0)




              KSTART = 3
              JSTART = 2
              ISTART = 2
              KSTOP = KMX - 2
              JSTOP = JMX - 2
              ISTOP = IMX - 2

      HP    =  0.5
      HN    = -0.5


C
C                                 BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 UND BEI GITTERKOPPLUNG
C                                 SOWIE BEI AUSFLUSSBEDINGUNG
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT

      NBU = 0
      IF ((NBAC.EQ.1).OR.(NBAC.EQ.7))  NBU = 1
      IF ((NBAC.EQ.3).OR.(NBAC.EQ.4).OR. (NBAC.EQ.13))  NBU = 1
      IF (NBAC.EQ.8)  NBU = 1
      NLV = 0
      IF ((NLFT.EQ.1).OR.(NLFT.EQ.7).OR.(NLFT.EQ.8))  NLV = 1
      IF ((NLFT.EQ.3).OR.(NLFT.EQ.4))  NLV = 1
      NTW = 0
      IF ((NTOP.EQ.1).OR.(NTOP.EQ.7).OR.(NTOP.EQ.3))  NTW = 1
      IF (NTOP.EQ.8)  NTW = 1

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
C***********************************************************************
C*     DECONVOLUTION OF FIELD FOR IMPROVED RESOLUTION KONVEKTIVE TERM  *
C***********************************************************************
C***********************************************************************
C           HIER KOMPAKTVERFAHREN VIERTER ORDNUNG IN X-RICHTUNG        *
C***********************************************************************





      IF(ABS(WKON) .LE. SMALL) GOTO 2200

C***********************************************************************
C                                 10 10 10 10 10 10 10 10 10 10 10 10 10

      DO  10  I = ISTART, ISTOP


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

C                                 20 20 20 20 20 20 20 20 20 20 20 20 20
      DO  20  J = JSTART, JSTOP
C                                  GEOMETRISCHE KONTANTEN (J)
                DYJ   =  DY(J)
                DYJM  =  DY(J-1)
                DYJP  =  DY(J+1)
               DDYJ   = DDY(J)
               DDYJM  = DDY(J-1)
               DDYJP  = DDY(J+1)

               RDYJ   =  RDY(J)
               RDYJM  =  RDY(J-1)
               RDYJP  =  RDY(J+1)
              RDDYJ   = RDDY(J)
              RDDYJM  = RDDY(J-1)
              RDDYJP  = RDDY(J+1)

              AUZ   =  DXI*DDYJ
                AVZ   = DDXI* DYJ
                AWZ   = DDXI*DDYJ


C                                 **************************************
C                                 KONVEKTIVE TERME   KONVEKTIVE TERME
C                                 **************************************

C                             KOEFFIZIENTEN BEIM KONVEKTIVEN ZEITSCHRITT


                FKDTU = -1.0*RDDY(J)* RDX(I)*WKON
                FKDTV = -1.0*RDDX(I)* RDY(J)*WKON
                FKDTW = -1.0*RDDX(I)*RDDY(J)*WKON




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



      DO 100  K = KSTART, KSTOP
C                                  GEOMETRISCHE KONSTANTEN (K)
                DZK   =  DZ(K)
                DZKM  =  DZ(K-1)
                DZKMM =  DZ(K-2)
                DZKP  =  DZ(K+1)
               DDZK   =  DDZ(K)
              RDDZK   = RDDZ(K)

                AUX   = DDYJ*DDZK
                AUY   =  DXI*DDZK

      FUE  =  AUX*(U(K,J,I)  +(U(K,J,I+1)-U(K,J,I)  )*0.5*DXI /DDXIP)
      FUW  =  AUX*(U(K,J,I-1)+(U(K,J,I)  -U(K,J,I-1))*0.5*DXIM/DDXI)
      FUN  =  AUY*(V(K,J,I)  + V(K,J,I+1)  )         *0.5
      FUS  =  AUY*(V(K,J-1,I)+ V(K,J-1,I+1))         *0.5
      FUT  =  AUZ*(W(K,J,I)  + W(K,J,I+1)  )         *0.5
      FUB  =  AUZ*(W(K-1,J,I)+ W(K-1,J,I+1))         *0.5





C                                 JETZT KOMPAKTE INTERPOLATION
C                                 U IN I,J UND K-RICHTUNG
       QKUE =FUE*(U(K,J,I)+(U(K,J,I+1)-U(K,J,I)  )*0.5*DXI /DDXIP)
       QKUW =FUW*(U(K,J,I-1)+(U(K,J,I)-U(K,J,I-1))*0.5*DXIM/DDXI)

      QKUN = 0.5*FUN*(U(K,J,I)   + U(K,J+1,I))
      QKUS = 0.5*FUS*(U(K,J-1,I) + U(K,J,I)  )

      QKUT = 0.5*FUT*(U(K,J,I)   + U(K+1,J,I)) 
      QKUB = 0.5*FUB*(U(K-1,J,I) + U(K,J,I)  )



C                                 ZEITSCHRITT


      UO(K,J,I) = FKDTU * RDDZK * (QKUE-QKUW+QKUN-QKUS+QKUT-QKUB)



  100 CONTINUE

 2110 CONTINUE

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



      DO 110  K = KSTART, KSTOP
C                                  GEOMETRISCHE KONSTANTEN (K)
                DZK   =  DZ(K)
                DZKM  =  DZ(K-1)
                DZKMM =  DZ(K-2)
                DZKP  =  DZ(K+1)
               DDZK   =  DDZ(K)
              RDDZK   = RDDZ(K)

                AVX   =  DYJ*DDZK
                AVY   = DDXI*DDZK

C                                 KONVEKTIVE F-TERME
      FVE  =     AVX*(U(K,J,I)  + U(K,J+1,I)  )         *0.5
      FVW  =     AVX*(U(K,J,I-1)+ U(K,J+1,I-1))         *0.5
      FVN  =     AVY*(V(K,J,I)  +(V(K,J+1,I)-V(K,J,I) ) *0.5*DYJ /DDYJP)
      FVS  =     AVY*(V(K,J-1,I)+(V(K,J,I)  -V(K,J-1,I))*0.5*DYJM/DDYJ )
      FVT  =     AVZ*(W(K,J,I)  + W(K,J+1,I)  )         *0.5
      FVB  =     AVZ*(W(K-1,J,I)+ W(K-1,J+1,I))         *0.5





C




C                                 JETZT KOMPAKTE INTERPOLATION
C                                 V IN I,J UND K-RICHTUNG
      QKVE = 0.5*FVE*(V(K,J,I)   + V(K,J,I+1))
      QKVW = 0.5*FVW*(V(K,J,I-1) + V(K,J,I)  )

       QKVN =FVN*(V(K,J,I) +(V(K,J+1,I)-V(K,J,I) ) *0.5*DYJ /DDYJP)
       QKVS =FVS*(V(K,J-1,I)+(V(K,J,I) -V(K,J-1,I))*0.5*DYJM/DDYJ )

      QKVT = 0.5*FVT*(V(K,J,I)   + V(K+1,J,I))
      QKVB = 0.5*FVB*(V(K-1,J,I) + V(K,J,I)  )



C                                 ZEITSCHRITT


      VO(K,J,I) =  FKDTV * RDDZK * (QKVE-QKVW+QKVN-QKVS+QKVT-QKVB)



  110 CONTINUE

 2120 CONTINUE

C                                 WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
C                                 --------------------------------------

      IF ( I .EQ. ISTART ) GOTO 2130
      IF ( J .EQ. JSTART ) GOTO 2130

      KSTOPM = KSTOP - 1 + NTW
      KSTARTM= KSTART - NBW

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

                AWX   = DDYJ*DZK
                AWY   = DDXI*DZK

C                                 KONVEKTIVE F-TERME
      FWE  =     AWX*(U(K,J,I)  + U(K+1,J,I)  )         *0.5
      FWW  =     AWX*(U(K,J,I-1)+ U(K+1,J,I-1))         *0.5
      FWN  =     AWY*(V(K,J,I)  + V(K+1,J,I)  )         *0.5
      FWS  =     AWY*(V(K,J-1,I)+ V(K+1,J-1,I))         *0.5
      FWT  =     AWZ*(W(K,J,I)  +(W(K+1,J,I)-W(K,J,I)  )*0.5*DZK /DDZKP)
      FWB  =     AWZ*(W(K-1,J,I)+(W(K,J,I)  -W(K-1,J,I))*0.5*DZKM/DDZK )
C                                 JETZT KOMPAKTE INTERPOLATION
C                                 W IN I,J UND K-RICHTUNG
      QKWE = 0.5*FWE*(W(K,J,I)   + W(K,J,I+1))
      QKWW = 0.5*FWW*(W(K,J,I-1) + W(K,J,I)  )

      QKWN = 0.5*FWN*(W(K,J,I)   + W(K,J+1,I))
      QKWS = 0.5*FWS*(W(K,J-1,I) + W(K,J,I)  )

       QKWT =FWT*(W(K,J,I)+(W(K+1,J,I)-W(K,J,I)  )*0.5*DZK /DDZKP)
       QKWB =FWB*(W(K-1,J,I)+(W(K,J,I)-W(K-1,J,I))*0.5*DZKM/DDZK )


C                                 ZEITSCHRITT


      WO(K,J,I) =  FKDTW * RDZK  * (QKWE-QKWW+QKWN-QKWS+QKWT-QKWB)



  120 CONTINUE

 2130 CONTINUE

   20 CONTINUE
   10 CONTINUE
C                                 20 20 20 20 20 20 20 20 20 20 20 20 20
C                                 10 10 10 10 10 10 10 10 10 10 10 10 10

c 2200  IF(ABS(WDIF) .LE. SMALL) GOTO 2300
 2200  CONTINUE

C                                 **************************************
C                                 DIFFUSIVE TERME   DIFFUSIVE TERME   DI
C                                 **************************************


C                                  WALL KORRECTION FUER U,V,W
      IF(ABS(WDIF) .GT. SMALL) THEN

        CALL SWCLE3D  (KK,JJ,II,KMX,JMX,IMX,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,P,G,B,WCU,WCV,WCW,
     $                 GMOL,RHO,UGRID,MTURB,ZTOT,GRADPX,
     $                 IC1,IC2,JC1,JC2,KC1,KC2,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,
     $                 FUK,FVK,FUJ,FWJ,FVI,FWI,
     $                 RDDX,RDDY,RDDZ)
      ENDIF

      DO I = 1,II-1
         DO J = 1,JJ-1
            DO K = 1,KK-1
               GI(k,j,i) = G(k,j,i)*G(k,j,i+1)/(G(k,j,i)+G(k,j,i+1))
               GJ(k,j,i) = G(k,j,i)*G(k,j+1,i)/(G(k,j,i)+G(k,j+1,i))
               GK(k,j,i) = G(k,j,i)*G(k+1,j,i)/(G(k,j,i)+G(k+1,j,i))
            enddo
         enddo
      enddo

C***********************************************************************
C                                 11 11 11 11 11 11 11 11 11 11 11 11 11

      DO  11  I = ISTART, ISTOP

C                            BERECHNUNG DER DRUCKGRADIENTEN IN Z-RICHTUNG


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


C                                 22 22 22 22 22 22 22 22 22 22 22 22 22

      DO  22  J = JSTART, JSTOP
         
C                                  GEOMETRISCHE KONTANTEN (J)
                DYJ   =  DY(J)
                DYJM  =  DY(J-1)
                DYJP  =  DY(J+1)
               DDYJ   = DDY(J)
               DDYJM  = DDY(J-1)
               DDYJP  = DDY(J+1)

               RDYJ   =  RDY(J)
               RDYJM  =  RDY(J-1)
               RDYJP  =  RDY(J+1)
              RDDYJ   = RDDY(J)
              RDDYJM  = RDDY(J-1)
              RDDYJP  = RDDY(J+1)

                AUZ   =  DXI*DDYJ
                AVZ   = DDXI* DYJ
                AWZ   = DDXI*DDYJ

       IF(ABS(WDIF) .LE. SMALL) GOTO 2300

C                             KOEFFIZIENTEN BEIM DIFFUSIVEN ZEITSCHRITT


                FDDTU = -1.0/RHO*RDDY(J)* RDX(I)*WDIF
                FDDTV = -1.0/RHO*RDDX(I)* RDY(J)*WDIF
                FDDTW = -1.0/RHO*RDDX(I)*RDDY(J)*WDIF


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


      DO 200  K = KSTART, KSTOP
C                                  GEOMETRISCHE KONSTANTEN (K)
               RDZK   =  RDZ(K)
               RDZKM  =  RDZ(K-1)
               DDZK   =  DDZ(K)
              RDDZK   = RDDZ(K)

                AUX   = DDYJ*DDZK
                AUY   =  DXI*DDZK

C                                 DIFFUSIVE D-TERME MIT
C                                 TRANSPORTKOEFFIZIENTEN

      GUE  = G(K,J,I+1)
      GUW  = G(K,J,I)
      GUN  = GJ(k  ,j  ,i)+GJ(k  ,j  ,i+1)
      GUS  = GJ(k  ,j-1,i)+GJ(k  ,j-1,i+1)
      GUT  = GK(k  ,j  ,i)+GK(k  ,j  ,i+1)
      GUB  = GK(k-1,j  ,i)+GK(k-1,j  ,i+1)


C      GUE  = G(K,J,I+1)
C      GUW  = G(K,J,I)
C      GUN  = G(K,J,I)  *G(K,J+1,I)   / (G(K,J,I)  +G(K,J+1,I)  )
C     $     + G(K,J,I+1)*G(K,J+1,I+1) / (G(K,J,I+1)+G(K,J+1,I+1))
C      GUS  = G(K,J,I)  *G(K,J-1,I)   / (G(K,J,I)  +G(K,J-1,I)  )
C     $     + G(K,J,I+1)*G(K,J-1,I+1) / (G(K,J,I+1)+G(K,J-1,I+1))
C      GUT  = G(K,J,I)  *G(K+1,J,I)   / (G(K,J,I)  +G(K+1,J,I)  )
C     $     + G(K,J,I+1)*G(K+1,J,I+1) / (G(K,J,I+1)+G(K+1,J,I+1))
C      GUB  = G(K,J,I)  *G(K-1,J,I)   / (G(K,J,I)  +G(K-1,J,I)  )
C     $     + G(K,J,I+1)*G(K-1,J,I+1) / (G(K,J,I+1)+G(K-1,J,I+1))


C                                     HIER KOMPAKTE BERECHNUNG
C                                     DER ERSTEN ABLEITUNG U
      QDUE = -GUE * AUX*RDDXIP * (U(K,J,I+1) - U(K,J,I))
      QDUW = -GUW * AUX*RDDXI  * (U(K,J,I)   - U(K,J,I-1))

      QDUN = -GUN * AUY*RDYJ   * (U(K,J+1,I) - U(K,J,I))
      QDUS = -GUS * AUY*RDYJM  * (U(K,J,I)   - U(K,J-1,I))

      QDUT = -GUT * AUZ*RDZK   * (U(K+1,J,I) - U(K,J,I))
      QDUB = -GUB * AUZ*RDZKM  * (U(K,J,I)   - U(K-1,J,I))

      QDUC = WCU(K,J,I)



C                                 ZEITSCHRITT

      UO(K,J,I) = UO(K,J,I) 
     $            + FDDTU * RDDZK * (QDUE-QDUW+QDUN-QDUS+QDUT-QDUB-QDUC)

  200 CONTINUE

 2210 CONTINUE


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


      DO 210  K = KSTART, KSTOP
C                                  GEOMETRISCHE KONSTANTEN (K)
                RDZK   =  RDZ(K)
                RDZKM  =  RDZ(K-1)
               DDZK   = DDZ(K)
              RDDZK   =RDDZ(K)

                AVX   =  DYJ*DDZK
                AVY   = DDXI*DDZK

C                                 DIFFUSIVE D-TERME MIT
C                                 TRANSPORTKOEFFIZIENTEN

       GVE  = GI(k  ,j  ,i  )+GI(k  ,j+1,i  )
       GVW  = GI(k  ,j  ,i-1)+GI(k  ,j+1,i-1)
       GVN  = G (K  ,J+1,I  )
       GVS  = G (K  ,J  ,I  )
       GVT  = GK(k  ,j  ,i  )+GK(k  ,j+1,i  )
       GVB  = GK(k-1,j  ,i  )+GK(k-1,j+1,i  )

C      GVE  = G(K,J,I)  *G(K,J,I+1)   / (G(K,J,I)+  G(K,J,I+1)  )
C     $     + G(K,J+1,I)*G(K,J+1,I+1) / (G(K,J+1,I)+G(K,J+1,I+1))
C      GVW  = G(K,J,I)  *G(K,J,I-1)   / (G(K,J,I)  +G(K,J,I-1)  )
C     $     + G(K,J+1,I)*G(K,J+1,I-1) / (G(K,J+1,I)+G(K,J+1,I-1))
C      GVN  = G(K,J+1,I)
C      GVS  = G(K,J,I)
C      GVT  = G(K,J,I)  *G(K+1,J,I)   / (G(K,J,I)  +G(K+1,J,I)  )
C     $     + G(K,J+1,I)*G(K+1,J+1,I) / (G(K,J+1,I)+G(K+1,J+1,I))
C      GVB  = G(K,J,I)  *G(K-1,J,I)   / (G(K,J,I)  +G(K-1,J,I)  )
C     $     + G(K,J+1,I)*G(K-1,J+1,I) / (G(K,J+1,I)+G(K-1,J+1,I))


C                                     HIER KOMPAKTE BERECHNUNG
C                                     DER ERSTEN ABLEITUNG V
      QDVE = -GVE * AVX*RDXI   * (V(K,J,I+1) - V(K,J,I))
      QDVW = -GVW * AVX*RDXIM  * (V(K,J,I)   - V(K,J,I-1))

      QDVN = -GVN * AVY*RDDYJP * (V(K,J+1,I) - V(K,J,I))
      QDVS = -GVS * AVY*RDDYJ  * (V(K,J,I)   - V(K,J-1,I))

      QDVT = -GVT * AVZ*RDZK   * (V(K+1,J,I) - V(K,J,I))
      QDVB = -GVB * AVZ*RDZKM  * (V(K,J,I)   - V(K-1,J,I))

      QDVC =  WCV(K,J,I)


C                                  ZEITSCHRITT

      VO(K,J,I) = VO(K,J,I) 
     $            + FDDTV * RDDZK * (QDVE-QDVW+QDVN-QDVS+QDVT-QDVB-QDVC)

  210 CONTINUE

 2220 CONTINUE

C                                 WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
C                                 --------------------------------------

      IF ( I .EQ. ISTART ) GOTO 2230
      IF ( J .EQ. JSTART ) GOTO 2230

      KSTOPM = KSTOP - 1 + NTW
      KSTARTM= KSTART - NBW

      DO 220  K = KSTARTM, KSTOPM
C                                  GEOMETRISCHE KONSTANTEN (K)
                DZK   =   DZ(K)
               RDZK   =  RDZ(K)
              RDDZK   = RDDZ(K)
              RDDZKP  = RDDZ(K+1)

                AWX   = DDYJ*DZK
                AWY   = DDXI*DZK

C                                 DIFFUSIVE D-TERME MIT
C                                 TRANSPORTKOEFFIZIENTEN


       GWE  = GI(k,j  ,i  )+GI(k+1,j  ,i  )
       GWW  = GI(k,j  ,i-1)+GI(k+1,j  ,i-1)
       GWN  = GJ(k,j  ,i  )+GJ(k+1,j  ,i  )
       GWS  = GJ(k,j-1,i  )+GJ(k+1,j-1,i  )
       GWT  = G(K+1,J,I)
       GWB  = G(K,J,I)

C      GWE  = G(K,J,I)  *G(K,J,I+1)   / (G(K,J,I)  +G(K,J,I+1)  )
C     $     + G(K+1,J,I)*G(K+1,J,I+1) / (G(K+1,J,I)+G(K+1,J,I+1))
C      GWW  = G(K,J,I)  *G(K,J,I-1)   / (G(K,J,I)  +G(K,J,I-1)  )
C     $     + G(K+1,J,I)*G(K+1,J,I-1) / (G(K+1,J,I)+G(K+1,J,I-1))
C      GWN  = G(K,J,I)  *G(K,J+1,I)   / (G(K,J,I)  +G(K,J+1,I)  )
C     $     + G(K+1,J,I)*G(K+1,J+1,I) / (G(K+1,J,I)+G(K+1,J+1,I))
C      GWS  = G(K,J,I)  *G(K,J-1,I)   / (G(K,J,I)  +G(K,J-1,I)  )
C     $     + G(K+1,J,I)*G(K+1,J-1,I) / (G(K+1,J,I)+G(K+1,J-1,I))
C      GWT  = G(K+1,J,I)
C      GWB  = G(K,J,I)


C                                     HIER KOMPAKTE BERECHNUNG
C                                     DER ERSTEN ABLEITUNG W
      QDWE = -GWE * AWX*RDXI   * (W(K,J,I+1) - W(K,J,I))
      QDWW = -GWW * AWX*RDXIM  * (W(K,J,I)   - W(K,J,I-1))

      QDWN = -GWN * AWY*RDYJ   * (W(K,J+1,I) - W(K,J,I))
      QDWS = -GWS * AWY*RDYJM  * (W(K,J,I)   - W(K,J-1,I))

      QDWT = -GWT * AWZ*RDDZKP * (W(K+1,J,I) - W(K,J,I))
      QDWB = -GWB * AWZ*RDDZK  * (W(K,J,I)   - W(K-1,J,I))
      QDWC = WCW(K,J,I)



C                                 ZEITSCHRITT
      WO(K,J,I) = WO(K,J,I) 
     $            + FDDTW * RDZK  * (QDWE-QDWW+QDWN-QDWS+QDWT-QDWB-QDWC)

  220 CONTINUE

 2230 CONTINUE

 2300 IF(IDUZ .EQ. 1) GOTO 2400

C                                 **************************************
C                                 QUELL-TERME   QUELL-TERME   QUELL-TERM
C                                 **************************************
C
C                                 UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C     
C                                 NUR BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT

      IF ((NBU.EQ.0).AND.(I .EQ. ISTOP))  GOTO 2310
      IF ((NFU.EQ.0).AND.(I .EQ. ISTART)) GOTO 2310
      IF ( J .EQ. JSTART ) GOTO 2310


C                                 ACHTUNG: IN X-RICHTUNG KANN EIN
C                                 MITTLERER DRUCKGRADIENT VORGEGEBEN
C                                 WERDEN
C
C                                 ZEITSCHRITT
C


      DO K = KSTART, KSTOP
         UO(K,J,I) = UO(K,J,I) - 1.0/(RHO*DXI) * WSOR
     $        * (P(K,J,I+1)-P(K,J,I) + GRADPX*DXI)
      ENDDO






 2310 CONTINUE

C                                 VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
C                                 --------------------------------------
C                                 NUR BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT

      IF ((NLV.EQ.0).AND.(J .EQ. JSTOP))  GOTO 2320
      IF ((NRV.EQ.0).AND.(J .EQ. JSTART)) GOTO 2320
      IF (I .EQ. ISTART) GOTO 2320


C                                  ZEITSCHRITT


      DO K = KSTART, KSTOP
         VO(K,J,I) = VO(K,J,I) - 1.0/(RHO*DYJ) * WSOR
     $        * (P(K,J+1,I)-P(K,J,I))
      ENDDO



 2320 CONTINUE

C                                 WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
C                                 --------------------------------------

      IF (I .EQ. ISTART) GOTO 2400
      IF (J .EQ. JSTART) GOTO 2400

      KSTOPM = KSTOP - 1 + NTW
      KSTARTM= KSTART - NBW

C                                 ZEITSCHRITT

      DO K = KSTARTM,KSTOPM
         WO(K,J,I) = WO(K,J,I) - 1.0/(RHO*DZ(K)) * WSOR
     $             * (P(K+1,J,I)-P(K,J,I))
      ENDDO




 2400 CONTINUE
 


   22 CONTINUE
   11 CONTINUE

      RETURN
      END
