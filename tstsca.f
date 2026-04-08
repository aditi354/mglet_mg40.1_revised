










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
C#ifdef _TSCAL_
      SUBROUTINE TSTSCA  (KK,JJ,II,KMX,JMX,IMX,
     $                    X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                    RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                    U,V,W,T,UO,VO,WO,TO,TP,P,G,
     $                    B,DT,LTXSTART,TXSTART,TFR,
     $                    NBUF,ITSTEP,WPHISCA,WCONSCA,WDIFSCA,
     $                    WSORSCA,IDUZSCA,IC1,IC2,JC1,JC2,KC1,KC2,
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
C*STARLET***************************************************************
C        T S T S C A      TIMESTEP - SCALAR TRANSPORT - TURBULENT - LES
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAY DIMENSIONS
C        KMX, JMX, IMX  - GRID SIZE IN Z,Y,X DIRECTION
C        DX,DY,DZ       - GRID SPACING (CENTER) 
C        DDX,DDY,DDZ    - CONTROL VOLUME SIZE IN X,Y,Z DIRECTION
C        U(KK,JJ,II)    - NEW VELOCITY FIELD
C        V(KK,JJ,II)    - NEW VELOCITY FIELD
C        W(KK,JJ,II)    - NEW VELOCITY FIELD
C        T(KK,JJ,II)    - NEW TEMPERATURE FIELD
C        UO(KK,JJ,II)   + OLD VELOCITY FIELD (CHANGED FOR
C                         IDUZSCA = 1 AND IDUZSCA = 2 !!)
C        VO(KK,JJ,II)   + OLD VELOCITY FIELD (  "    "  "    "
C        WO(KK,JJ,II)   + OLD VELOCITY FIELD (  "    "  "    "
C        TO(KK,JJ,II)   + OLD TEMPERATURE FIELD (CHANGED FOR
C                         IDUZSCA = 1 AND IDUZSCA = 2 !!)
C        TP(KK,JJ,NBUF) + BUFFERARRAYS
C        G(KK,JJ,II)    - EFFECTIF DYNAMIC VISCOSITY (= MUE)
C        GMOL           - MOLEKULARE DYNAMISCHE VISKOSITAET
C        PRMOL          - MOLECULAR PRANDTL NUMBER
C        PRTURB         - TURBULENT PRANDTL NUMBER
C                         (IF 0 => KAYS FORMULATION!)
C        LAMDA          - SPECIFIC HEAT TRANSFER COEFFICIENT
C        DT             - TIME STEP SIZE
C        RHO            - DENSITY (= CONST)
C        WCT            - WALL CORRECTION FOR T (Diffusive transport by
C                                               wall heat flux)
C        MTURB          - 0 = LAMINAR,  1 = TURBULENT
C        NBUF           - NUMBER OF BUFFER LAYERS (F. QUICK: NBUF=3)
C        ITSTEP         - TIME STEP COUNTER
C        WCONSCA        - WEIGHT FACTOR FOR CONVECTIVE TRANSPORT
C        WDIFSCA        -                  DIFFUSIVE TRANSPORT
C        WSORSCA        -                  HEAT SOURCES
C        WPHISCA        -                  POINT VALUES
C        IDUZSCA        - FOR IDUZSCA= 1: OLD TIME LEVEL 
C                         FOR IDUZSCA= 2: NEW TIME LEVEL
C        IC1,JC1,KC1    - LEFT BORDER INDICES OF BOUNDING BOX 
C        IC2,JC2,KC2    - RIGHT BORDER INDICES OF BOUNDING BOX          
C                         (=1 LAYER AROUND BODY)
C                  X    - X-KOORDINATE DES ZELLMITTELPUNKTES            
C
C DEFINE DIREKTIVEN     : ZEN, QUD, UPW,EULERSCA, LEAPFSCA, ADBASCA, 
C                         BAPAR_FIXED, TOPAR_FIXED
C UPROG                 : SWCLESCA
C
C        26.09.02 (TB)  : TSTSCA : BASED ON TSTLE2
C                                  ZEITSCHRITT WAHLWEISE MITTELS EULER-,
C                                  LEAPFROG- ODER ADAMS/BASHFORTH-VER-
C                                  FAHREN
C*STARLET***************************************************************
C

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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
C     VERS:   28.09.95 (AO) JET VARIABLES INTRODUCED
C             12.12.02 (TB) SCALAR PARAMETERS (FIX VALUE, GRADIENT) ADDED
C
C     XMPOS1,XMPOS2,YMPOS1,YMPOS2,ZMPOS1,ZMPOS2: BEREICH DER EFFEKTMESSUNG 
C                                                AUFGRUND DER MANIPULATION

      COMMON /COBOUND/
     &                 NBOCD,
     &                NBOCONDS,     LARBOCONDS,     ITYPBOCONDS,
     &                LBOGRIDS,    LPOSBOGRIDS,
     &                 FRONT,    BACK,     RIGHT,     LEFT,
     &                 BOTTOM,   TOP,      CUBE,
     &                 IBPOS,   JBPOS,    KBPOS,
     &                 IBANF,   JBANF,    KBANF,
     &                 IBEND,   JBEND,    KBEND,
     &                 XBANF,   YBANF,    ZBANF,
     &                 XBEND,   YBEND,    ZBEND,
     &                 ANIVEAU,
     &                 AUB, AVB, AWB,
     &                 FREQB,  FLOWTYP, WAVENUMBER,
     &                 LOPT,NTOPT1,NTOPT2,
     &                 XMPOS1,YMPOS1,ZMPOS1,
     &                 XMPOS2,YMPOS2,ZMPOS2,
     &                 TIMEALT,XRTALT1,XRTALT2,FREQOPT,PERIODE,ITALT,
     &                 FREQALT1,FREQALT2,XRMIN,FXRMIN,NXRMIN,
     &                 RANNUM,PHASE


      INTEGER
     &       NBOCD   (                9, MAXGRIDS),
     &       NBOCONDS(                9, MAXGRIDS),
     &     LARBOCONDS( 6, MAXBOCONDS, 9, MAXGRIDS),
     &    ITYPBOCONDS(    MAXBOCONDS, 9, MAXGRIDS),
     &       LBOGRIDS(    MAXBOCONDS, 9, MAXGRIDS),
     &    LPOSBOGRIDS( 3, MAXBOCONDS, 9, MAXGRIDS),
     &   IBPOS(MAXBOCONDS,9,MAXGRIDS),  JBPOS(MAXBOCONDS,9,MAXGRIDS),
     &   KBPOS(MAXBOCONDS,9,MAXGRIDS),
     &   IBANF(MAXBOCONDS,9,MAXGRIDS),  IBEND(MAXBOCONDS,9,MAXGRIDS),
     &   JBANF(MAXBOCONDS,9,MAXGRIDS),  JBEND(MAXBOCONDS,9,MAXGRIDS),
     &   KBANF(MAXBOCONDS,9,MAXGRIDS),  KBEND(MAXBOCONDS,9,MAXGRIDS),
     &   LOPT(MAXBOCONDS,MAXGRIDS),
     &   NTOPT1(MAXBOCONDS,MAXGRIDS),NTOPT2(MAXBOCONDS,MAXGRIDS)



      CHARACTER (LEN=16)
     &      FRONT(MAXBOCONDS,MAXGRIDS),   BACK(MAXBOCONDS,MAXGRIDS),
     &      RIGHT(MAXBOCONDS,MAXGRIDS),   LEFT(MAXBOCONDS,MAXGRIDS),
     &     BOTTOM(MAXBOCONDS,MAXGRIDS),    TOP(MAXBOCONDS,MAXGRIDS),
     &       CUBE(MAXBOCONDS,MAXGRIDS),
     &      FLOWTYP(MAXBOCONDS,9,MAXGRIDS)

      REAL  ANIVEAU(MAXBOCONDS,9,MAXGRIDS), AUB(MAXBOCONDS,9,MAXGRIDS),
     &      AVB(MAXBOCONDS,9,MAXGRIDS), AWB(MAXBOCONDS,9,MAXGRIDS),
     &      FREQB(MAXBOCONDS,9,MAXGRIDS),
     &      XBANF(MAXBOCONDS,9,MAXGRIDS),XBEND(MAXBOCONDS,9,MAXGRIDS),
     &      YBANF(MAXBOCONDS,9,MAXGRIDS),YBEND(MAXBOCONDS,9,MAXGRIDS),
     &      ZBANF(MAXBOCONDS,9,MAXGRIDS),ZBEND(MAXBOCONDS,9,MAXGRIDS),
     &      XM1(MAXBOCONDS,9,MAXGRIDS),XM2(MAXBOCONDS,9,MAXGRIDS),
     &      XM3(MAXBOCONDS,9,MAXGRIDS),
     &      YM1(MAXBOCONDS,9,MAXGRIDS),YM2(MAXBOCONDS,9,MAXGRIDS),
     &      YM3(MAXBOCONDS,9,MAXGRIDS),
     &      ZM1(MAXBOCONDS,9,MAXGRIDS),ZM2(MAXBOCONDS,9,MAXGRIDS),
     &      ZM3(MAXBOCONDS,9,MAXGRIDS),
     &      WAVENUMBER(MAXBOCONDS,9,MAXGRIDS),
     &      XMPOS1(MAXBOCONDS,MAXGRIDS),
     &      XMPOS2(MAXBOCONDS,MAXGRIDS),
     &      YMPOS1(MAXBOCONDS,MAXGRIDS),
     &      YMPOS2(MAXBOCONDS,MAXGRIDS),
     &      ZMPOS1(MAXBOCONDS,MAXGRIDS),
     &      ZMPOS2(MAXBOCONDS,MAXGRIDS),
     &      RANNUM,
     &      PHASE(MAXBOCONDS,9,MAXGRIDS)

C
      REAL        X(II),         Y(JJ),         Z(KK),
     &           DX(II),        DY(JJ),        DZ(KK),
     &          DDX(II),       DDY(JJ),       DDZ(KK),
     $        RDX(II),RDY(JJ),RDZ(KK),RDDX(II),RDDY(JJ),RDDZ(KK)
C
      REAL   U(KK,JJ,II),   V(KK,JJ,II),   W(KK,JJ,II),T(KK,JJ,II),
     $       UO(KK,JJ,II),  VO(KK,JJ,II),  WO(KK,JJ,II),TO(KK,JJ,II),
     $       TP(KK,JJ,NBUF),TFR(KK,JJ,2)
C
      REAL        P(KK,JJ,II),   G(KK,JJ,II)
      REAL        B(KK,JJ,II)    

      LOGICAL LTXSTART
      REAL TXSTART  
         
C      REAL   QSUMD(KK,JJ,NBUF),QSUMC(KK,JJ,NBUF),TPTMP(KK,JJ,NBUF),
C     $       QCE(KK,JJ,NBUF),QCW(KK,JJ,NBUF),QCS(KK,JJ,NBUF),
C     $       QCN(KK,JJ,NBUF),QCB(KK,JJ,NBUF),QCT(KK,JJ,NBUF),
C     $       QDE(KK,JJ,NBUF),QDW(KK,JJ,NBUF),QDS(KK,JJ,NBUF),
C     $       QDN(KK,JJ,NBUF),QDB(KK,JJ,NBUF),QDT(KK,JJ,NBUF),
C     $       QDC(KK,JJ,NBUF)     
C
C-----------------------------------------------------------------------
C      
C-----------------------------------------------------------------------
      KSTART = 3
      JSTART = 2 
      ISTART = 2
      KSTOP = KMX - 2
      JSTOP = JMX - 2
      ISTOP = IMX - 2

C
C      HP    =  0.5
C      HN    = -0.5
C
      TMAX   = 0.0
      DTMAX  = 0.0
      WCTMAX = 0.0
      DTKMAX = 0
      DTJMAX = 0
      DTIMAX = 0 
C      WCTKMAX= 0
C      WCTJMAX= 0
C      WCTIMAX= 0 
C      TPMAX  = 0.0
C      TEMAX  = 0.0
C      TWMAX  = 0.0
C      TNMAX  = 0.0
C      TSMAX  = 0.0
C      TTMAX  = 0.0
C      TBMAX  = 0.0 
C      
C      QCMAX  = 0.0 
C      QDMAX  = 0.0 
C      QCEMAX = 0.0
C      QCWMAX = 0.0      
C      QCNMAX = 0.0
C      QCSMAX = 0.0   
C      QCTMAX = 0.0
C      QCBMAX = 0.0  
C      QDEMAX = 0.0
C      QDWMAX = 0.0      
C      QDNMAX = 0.0
C      QDSMAX = 0.0   
C      QDTMAX = 0.0
C      QDBMAX = 0.0 
C      QDCMAX = 0.0                                                  
C      

c      DO 1 I= ISTART,ISTOP
c       DO 2 J= JSTART,JSTOP
c        DO 3 K= KSTART,KSTOP
c          WRITE(6,*)'T(',K,',',J,',',I,')= ',T(K,J,I)  
c          WRITE(6,*)'TO(',K,',',J,',',I,')= ',TO(K,J,I)                
c    3   CONTINUE        
c    2  CONTINUE
c    1 CONTINUE 
CC-----------------------------------------------------------------------
CC                                 BEI PERIODISCHEN RANDBEDINGUNGEN
CC                                 UND BEI GITTERKOPPLUNG
CC                                 SOWIE BEI AUSFLUSSBEDINGUNG
CC                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
CC                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
CC                                 BOUNDMG GESETZT
C
C      NBU = 0
C      IF ((NBAC.EQ.1).OR.(NBAC.EQ.7))  NBU = 1
C      IF ((NBAC.EQ.3).OR.(NBAC.EQ.4).OR.(NBAC.EQ.13))  NBU = 1
C#ifndef _BAPAR_FIXED_
C      IF (NBAC.EQ.8)  NBU = 1
C#endif
C      NLV = 0
C      IF ((NLFT.EQ.1).OR.(NLFT.EQ.7).OR.(NLFT.EQ.8))  NLV = 1
C      IF ((NLFT.EQ.3).OR.(NLFT.EQ.4).OR.(NLFT.EQ.13)) NLV = 1
C      NTW = 0
C      IF ((NTOP.EQ.1).OR.(NTOP.EQ.7).OR.(NTOP.EQ.3).OR.(NTOP.EQ.4)
C     &               .OR.(NTOP.EQ.13)) NTW = 1
C#ifndef _TOPAR_FIXED_
C      IF (NTOP.EQ.8)  NTW = 1
C#endif
CC                                 BEI LOKAL VERFEINERTEM GITTER
CC                                 WIRD DIE NORMALKOMPONENTE DER ERSTEN
CC                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
CC                                 BOUNDMG GESETZT
CC
      NFT = 0
C      IF ((NFRO.EQ.8).OR.(NFRO.EQ.3).OR.(NFRO.EQ.4).OR.(NFRO.EQ.13)) 
C     &    NFU = 1
C      NRV = 0
C      IF ((NRGTO.EQ.8).OR.(NRGT.EQ.3).OR.(NRGT.EQ.4).OR.(NRGT.EQ.13))
C     &    NRV = 1
C      NBW = 0
C      IF ((NBOT.EQ.8).OR.(NBOT.EQ.3).OR.(NBOT.EQ.4).OR.(NBOT.EQ.13))
C     &    NBW = 1
CC
C=======================================================================
C                                 10 10 10 10 10 10 10 10 10 10 10 10 10
      DO  10  I = ISTART, ISTOP
C
      ICOM = 1 + MOD(I,NBUF)   
C      WRITE(6,*) 'ICOM= ',ICOM       
C
C                                  *************************************
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
   
C=======================================================================
C                                 20 20 20 20 20 20 20 20 20 20 20 20 20
      DO  20  J = JSTART, JSTOP
C                                  *************************************



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
                ATZ   = DDXI*DDYJ
C                                 **************************************
C                                 PUFFERFELDER WERDEN MIT NULL BELEGT
C
      DO  30 K = 1, KMX
          TP(K,J,ICOM) = 0.0 
         
C          QSUMC(K,J,ICOM) = 0.0 
C          QSUMD(K,J,ICOM) = 0.0 
C          TPTMP(K,J,ICOM) = 0.0
C          QCE(K,J,ICOM)   = 0.0
C          QCW(K,J,ICOM)   = 0.0   
C          QCS(K,J,ICOM)   = 0.0
C          QCN(K,J,ICOM)   = 0.0 
C          QCB(K,J,ICOM)   = 0.0
C          QCT(K,J,ICOM)   = 0.0  
C          QDE(K,J,ICOM)   = 0.0
C          QDW(K,J,ICOM)   = 0.0   
C          QDS(K,J,ICOM)   = 0.0
C          QDN(K,J,ICOM)   = 0.0 
C          QDB(K,J,ICOM)   = 0.0
C          QDT(K,J,ICOM)   = 0.0 
C          QDC(K,J,ICOM)   = 0.0                                                                    
   30 CONTINUE
    
      IF(ABS(WCONSCA) .LE. SMALL) GOTO 2200
C
C***********************************************************************
C                                 KONVEKTIVE TERME   KONVEKTIVE TERME
C***********************************************************************
C
C                             KOEFFIZIENTEN BEIM KONVEKTIVEN ZEITSCHRITT
C

		FKDTT = -1.0*DT*RDDX(I)*RDDY(J)*WCONSCA   
C      WRITE(6,*) 'FKDDT= ',FKDTT   
C
C                                 **************************************
C                                 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
C                                 NUR BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 UND BEI GITTERKOPPLUNG
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT
C
C      IF ((NBU.EQ.0).AND.(I .EQ. ISTOP))  GOTO 2110
C      IF ((NFU.EQ.0).AND.(I .EQ. ISTART)) GOTO 2110
C      IF ( J .EQ. JSTART ) GOTO 2110
C-----------------------------------------------------------------------
      DO 100  K = KSTART, KSTOP
C                                  GEOMETRISCHE KONSTANTEN (K)
                DZK   =  DZ(K)
                DZKM  =  DZ(K-1)
                DZKMM =  DZ(K-2)
                DZKP  =  DZ(K+1)
               DDZK   =  DDZ(K)
              RDDZK   = RDDZ(K)
C
                ATX   = DDYJ*DDZK
                ATY   = DDXI*DDZK
C                                ***************************************
C                                 KONVEKTIVE TERME
C
      FTE  =     ATX*U(K,J,I)
      FTW  =     ATX*U(K,J,I-1)
      FTN  =     ATY*V(K,J,I)  
      FTS  =     ATY*V(K,J-1,I)
      FTT  =     ATZ*W(K,J,I)  
      FTB  =     ATZ*W(K-1,J,I)
C
C                                 CENTRAL DIFFERENCING

      QCTE = 0.5*FTE*(T(K,J,I)   + T(K,J,I+1))
      QCTW = 0.5*FTW*(T(K,J,I-1) + T(K,J,I)  )
      QCTN = 0.5*FTN*(T(K,J,I)   + T(K,J+1,I))
      QCTS = 0.5*FTS*(T(K,J-1,I) + T(K,J,I)  )
      QCTT = 0.5*FTT*(T(K,J,I)   + T(K+1,J,I))
      QCTB = 0.5*FTB*(T(K-1,J,I) + T(K,J,I)  )

C
C                                 **************************************
C                                 TIME STEP
      QSUMC        = FKDTT * RDDZK * (QCTE-QCTW+QCTN-QCTS+QCTT-QCTB)
      TP(K,J,ICOM) = QSUMC
C      
C      QCE(K,J,ICOM) = QCTE*FKDTT*RDDZK
C      QCW(K,J,ICOM) = -QCTW*FKDTT*RDDZK
C      QCN(K,J,ICOM) = QCTN*FKDTT*RDDZK
C      QCS(K,J,ICOM) = -QCTS*FKDTT*RDDZK   
C      QCT(K,J,ICOM) = QCTT*FKDTT*RDDZK
C      QCB(K,J,ICOM) = -QCTB*FKDTT*RDDZK          

C      IF (ABS(QSUMC(J,I)) .GT. ABS(QCMAX)) THEN 
CC      	  WRITE(6,*) 'WCT= ',WCT                            
C          QCMAX  = QSUMC(K,J)
C          QCKMAX = K
C          QCJMAX = J 
C          QCIMAX = I
C      ENDIF      
C      WRITE(6,*)'TP(',K,',',J,',',ICOM,')= ',TP(K,J,ICOM)
C      WRITE(6,*)'QSUMC(',K,',',J,',',ICOM,')= ',QSUMC     
C
  100 CONTINUE
C-----------------------------------------------------------------------
 2110 CONTINUE
 
 2200 IF(ABS(WDIFSCA) .LE. SMALL) GOTO 2300
C
C***********************************************************************    
C                                 DIFFUSIVE TERME   DIFFUSIVE TERME   DI
C                                 **************************************
C
C                             KOEFFIZIENTEN BEIM DIFFUSIVEN ZEITSCHRITT
C
C      NEW: NO RHO IN FDDT       
       FDDTT = -1.0*DT*RDDX(I)*RDDY(J)*WDIFSCA  
C      WRITE(6,*) 'FDDDT= ',FDDTT        
C
C                                 **************************************
C                                 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
C                                 NUR BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT
C
C      IF ((NBU.EQ.0).AND.(I .EQ. ISTOP))  GOTO 2210
C      IF ((NFU.EQ.0).AND.(I .EQ. ISTART)) GOTO 2210
C      IF ( J .EQ. JSTART ) GOTO 2210
C-----------------------------------------------------------------------
      DO 200  K = KSTART, KSTOP
C                                  *************************************    
C                                  GEOMETRISCHE KONSTANTEN (K)
               RDZK   =  RDZ(K)
               RDZKM  =  RDZ(K-1)
               DDZK   =  DDZ(K)
              RDDZK   = RDDZ(K)
C
                ATX   = DDYJ*DDZK
                ATY   = DDXI*DDZK
                ATZ   = DDXI*DDYJ
C
C                                 **************************************
C                                 DIFFUSIVE D-TERME MIT
C                                 TRANSPORTKOEFFIZIENTEN
      GSCAMOL = GMOL / RHO / PRMOL  
      
      GSCAE  = GSCAMOL
      GSCAW  = GSCAMOL
      GSCAN  = GSCAMOL
      GSCAS  = GSCAMOL
      GSCAT  = GSCAMOL
      GSCAB  = GSCAMOL
   
      IF(MTURB .EQ. 1) THEN
      	 GTGMOLP = (G(K,J,I)-GMOL)   / GMOL  
      	 GTGMOLE = (G(K,J,I+1)-GMOL) / GMOL
      	 GTGMOLW = (G(K,J,I-1)-GMOL) / GMOL 
      	 GTGMOLN = (G(K,J+1,I)-GMOL) / GMOL
      	 GTGMOLS = (G(K,J-1,I)-GMOL) / GMOL      	 
      	 GTGMOLT = (G(K+1,J,I)-GMOL) / GMOL
      	 GTGMOLB = (G(K-1,J,I)-GMOL) / GMOL
C     1/RE * 1/PR + 1/REt * 1/PRt:                 	       	      
         GSCAE  = GSCAMOL + (G(K,J,I+1) + G(K,J,I) - 2.0*GMOL) / RHO
     $			  / (PRT(GTGMOLE) + PRT(GTGMOLP))
         GSCAW  = GSCAMOL + (G(K,J,I) + G(K,J,I-1) - 2.0*GMOL) / RHO
     $			  / (PRT(GTGMOLP) + PRT(GTGMOLW))   
         GSCAN  = GSCAMOL + (G(K,J+1,I) + G(K,J,I) - 2.0*GMOL) / RHO
     $			  / (PRT(GTGMOLN) + PRT(GTGMOLP))  
         GSCAS  = GSCAMOL + (G(K,J,I) + G(K,J-1,I) - 2.0*GMOL) / RHO
     $			  / (PRT(GTGMOLP) + PRT(GTGMOLS))  
         GSCAT  = GSCAMOL + (G(K+1,J,I) + G(K,J,I) - 2.0*GMOL) / RHO
     $			  / (PRT(GTGMOLT) + PRT(GTGMOLP))
         GSCAB  = GSCAMOL + (G(K,J,I) + G(K-1,J,I) - 2.0*GMOL) / RHO
     $			  / (PRT(GTGMOLP) + PRT(GTGMOLB))  
      ENDIF

C
C     ONLY "FREE" CELL FLUXES ARE CONSIDERED AS FOLLOWS. WALL CELL 
C     CONTRIBUTION IS TAKEN CARE OF IN THE WALL CORRECTION TERM 
C
C     BFIELD-FACS = 1 WHEN CELL FACE IS FREE (NO WALL), OTHERWISE = 0
      FACBE = MAX(0.0,SIGN(1.0,B(K,J,I)))* MAX(0.0,SIGN(1.0,B(K,J,I+1)))
      FACBW = MAX(0.0,SIGN(1.0,B(K,J,I)))* MAX(0.0,SIGN(1.0,B(K,J,I-1))) 
      FACBN = MAX(0.0,SIGN(1.0,B(K,J,I)))* MAX(0.0,SIGN(1.0,B(K,J+1,I)))           
      FACBS = MAX(0.0,SIGN(1.0,B(K,J,I)))* MAX(0.0,SIGN(1.0,B(K,J-1,I))) 
      FACBT = MAX(0.0,SIGN(1.0,B(K,J,I)))* MAX(0.0,SIGN(1.0,B(K+1,J,I)))    
      FACBB = MAX(0.0,SIGN(1.0,B(K,J,I)))* MAX(0.0,SIGN(1.0,B(K-1,J,I)))              
C
      QDTE = -GSCAE * ATX*RDXI   * (T(K,J,I+1) - T(K,J,I))   * FACBE
      QDTW = -GSCAW * ATX*RDXIM  * (T(K,J,I)   - T(K,J,I-1)) * FACBW    
      QDTN = -GSCAN * ATY*RDYJ   * (T(K,J+1,I) - T(K,J,I))   * FACBN 
      QDTS = -GSCAS * ATY*RDYJM  * (T(K,J,I)   - T(K,J-1,I)) * FACBS    
      QDTT = -GSCAT * ATZ*RDZK   * (T(K+1,J,I) - T(K,J,I))   * FACBT   
      QDTB = -GSCAB * ATZ*RDZKM  * (T(K,J,I)   - T(K-1,J,I)) * FACBB
C
C     CALCULATE WALL CORRECTION TERM 
C          
      CALL SWCLESCA  (KK,JJ,II,KMX,JMX,IMX,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,T,G,B,
     $                 WCT,I,J,K,
     $                 IC1,IC2,JC1,JC2,KC1,KC2,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,
     $                 RDDX,RDDY,RDDZ,IGRID,DDS,SM,SP)
      QDTC = WCT 
      
C      IF (ABS(WCT) .GT. ABS(WCTMAX)) THEN 
CC      	  WRITE(6,*) 'WCT= ',WCT                            
C          WCTMAX  = WCT
C          WCTKMAX = K
C          WCTJMAX = J 
C          WCTIMAX = I
C          TPMAX=T(K,J,I)
C          TEMAX=T(K,J,I+1)
C          TWMAX=T(K,J,I-1)
C          TNMAX=T(K,J+1,I)
C          TSMAX=T(K,J-1,I)
C          TTMAX=T(K+1,J,I)
C          TBMAX=T(K-1,J,I)
C      ENDIF      
                   
C
C                                 **************************************
C                                 TIME STEP
      QSUMD = FDDTT*RDDZK*(QDTE-QDTW+QDTN-QDTS+QDTT-QDTB-QDTC)
      TP(K,J,ICOM) = TP(K,J,ICOM) +  QSUMD
C
C      QDE(K,J,ICOM) = QDTE*FDDTT*RDDZK
C      QDW(K,J,ICOM) = -QDTW*FDDTT*RDDZK
C      QDN(K,J,ICOM) = QDTN*FDDTT*RDDZK
C      QDS(K,J,ICOM) = -QDTS*FDDTT*RDDZK   
C      QDT(K,J,ICOM) = QDTT*FDDTT*RDDZK
C      QDB(K,J,ICOM) = -QDTB*FDDTT*RDDZK 
C      QDC(K,J,ICOM) = -QDTC*FDDTT*RDDZK      
C      WRITE(6,*)'TP(',K,',',J,',',ICOM,')= ',TP(K,J,ICOM)
C      WRITE(6,*)'QSUMD(',K,',',J,',',ICOM,')= ',QSUMD 
C
C      IF (ABS(QSUMD(K,J)) .GT. ABS(QDMAX)) THEN 
CC      	  WRITE(6,*) 'WCT= ',WCT                            
C          QDMAX  = QSUMD(K,J)
C          QDKMAX = K
C          QDJMAX = J 
C          QDIMAX = I
C      ENDIF 
      
  200 CONTINUE
C-----------------------------------------------------------------------
 2210 CONTINUE
C
 2300 IF(ABS(WSORSCA) .LE. SMALL) GOTO 2400
C
C      CALL ERRR (502,' TSTSCA ')
C***********************************************************************
C                                 QUELL-TERME   QUELL-TERME   QUELL-TERM
C                                 **************************************
C
C                                 NUR BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT
C
C      IF ((NBU.EQ.0).AND.(I .EQ. ISTOP))  GOTO 2310
C      IF ((NFU.EQ.0).AND.(I .EQ. ISTART)) GOTO 2310
C      IF ( J .EQ. JSTART ) GOTO 2310
C-----------------------------------------------------------------------
C      DO 300  K = KSTART, KSTOP
C
C
C                                 **************************************
C                                 ZEITSCHRITT
C
C      TP(K,J,ICOM) = TP(K,J,ICOM) - DT/(RHO*LAMDA) * Q(K,J,I) * WSORSCA
C
C  300 CONTINUE
C-----------------------------------------------------------------------
C 2310 CONTINUE

 2400 IF(ABS(WPHISCA) .LE. SMALL) GOTO 2500
C
C***********************************************************************
C                                 ANTEIL DER PUNKTWERTE   ANTEIL DER PUN
C                                 **************************************
C
C                                 NUR BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT
C
C      IF ((NBU.EQ.0).AND.(I .EQ. ISTOP))  GOTO 2410
C      IF ((NFU.EQ.0).AND.(I .EQ. ISTART)) GOTO 2410
C      IF ( J .EQ. JSTART ) GOTO 2410
C-----------------------------------------------------------------------
      DO 400  K = KSTART, KSTOP
C                                 **************************************
C                                 ZEITSCHRITT
C
      TP(K,J,ICOM) = TP(K,J,ICOM) + WPHISCA * T(K,J,I)
  400 CONTINUE
C-----------------------------------------------------------------------
C 2410 CONTINUE
 2500 CONTINUE
C
C
   20 CONTINUE
C                                 20 20 20 20 20 20 20 20 20 20 20 20 20
C=======================================================================
C 
C                                 JETZT ZURUECKSPEICHERN DER TEMPORAEREN
C                                 PUFFER. NUR TATSAECHLICH BERECHNETE
C                                 PUFFERWERTE WERDEN ZURUECKGESPEICHERT
C      IF(I .LT. (ISTART+NBUF-1)) GOTO 3010
      IF(I .LT. (ISTART+NBUF)) GOTO 3010      
      INEW   = I
      IREPMX = 1
      IF(I .EQ. ISTOP) IREPMX = NBUF
C-----------------------------------------------------------------------
      DO 810 IREP=1,IREPMX
         IF(I .LT. ISTOP) GOTO 3020
            INEW = ISTOP + IREP - 1
 3020    IBUFF = 1 + MOD(INEW+1,NBUF)
         IBACK = INEW - NBUF + 1
C-----------------------------------------------------------------------
         DO 820 J=JSTART,JSTOP
C
C                                 NUR BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT
C
C      IF ((NBU.EQ.0).AND.(INEW .EQ. (ISTOP+NBUF-1)))  GOTO 3030
C      IF ((NFT.EQ.0).AND.(I .EQ. ISTART)) GOTO 3030
      IF (J .EQ. JSTART) GOTO 3030
C-----------------------------------------------------------------------
            DO 830 K=KSTART,KSTOP
C               IF (I .EQ. ISTART) GOTO 3030
C               IF (J .EQ. JSTART) GOTO 3030            
C               WRITE(6,*)'TP(',K,',',J,',',IBUFF,')= ',TP(K,J,IBUFF) 
               IF (LTXSTART) THEN
CTBC  ONLY WRITE T IN CASE X >= TXSTART (SPARE OUT AREA)
                 IF ( X(IBACK) .GE. TXSTART ) THEN                     
                    TO(K,J,IBACK) = TO(K,J,IBACK)*ZERONE + TP(K,J,IBUFF)             
                 ENDIF
               ELSE
                 TO(K,J,IBACK) = TO(K,J,IBACK)*ZERONE + TP(K,J,IBUFF)
               ENDIF                 
C 
C               IF (ABS(TO(K,J,IBACK)) .GE. TMAX) THEN                             
C                   TMAX  = ABS(TO(K,J,IBACK))
C                   TKMAX = K
C                   TJMAX = J 
C                   TIMAX = IBACK
C               ENDIF
              IF (IDUZSCA .EQ. 1) THEN
               IF (ABS(TP(K,J,IBUFF)) .GE. ABS(DTMAX)) THEN                             
              DTMAX   = ABS(TP(K,J,IBUFF))*SIGN(1.0,TP(K,J,IBUFF))
                   DTKMAX = K
                   DTJMAX = J 
                   DTIMAX = IBACK
                   TMAX  = ABS(TO(K,J,IBACK))
C                   QCMAX  = QSUMC(K,J,IBUFF)                   
C                   QDMAX  = QSUMD(K,J,IBUFF) 
C                   QCEMAX = QCE(K,J,IBUFF) 
C                   QCWMAX = QCW(K,J,IBUFF)  
C                   QCNMAX = QCN(K,J,IBUFF) 
C                   QCSMAX = QCS(K,J,IBUFF)  
C                   QCTMAX = QCT(K,J,IBUFF) 
C                   QCBMAX = QCB(K,J,IBUFF)   
C                   QDEMAX = QDE(K,J,IBUFF) 
C                   QDWMAX = QDW(K,J,IBUFF)  
C                   QDNMAX = QDN(K,J,IBUFF) 
C                   QDSMAX = QDS(K,J,IBUFF)  
C                   QDTMAX = QDT(K,J,IBUFF) 
C                   QDBMAX = QDB(K,J,IBUFF) 
C                   QDCMAX = QDC(K,J,IBUFF)                                                                                                                                            
               ENDIF            
              ENDIF                
C
  830       CONTINUE

C-----------------------------------------------------------------------
 3030       CONTINUE
C
  820    CONTINUE
C-----------------------------------------------------------------------
  810 CONTINUE
C-----------------------------------------------------------------------
 3010 CONTINUE
C
   10 CONTINUE  
C                                 10 10 10 10 10 10 10 10 10 10 10 10 10
C=======================================================================
C     

C      IF (IDUZSCA .EQ. 2) THEN 
C       WRITE(6,*) 
C       WRITE(6,6004 IGRID,INT(DTIMAX),INT(DTJMAX),INT(DTKMAX)        
C     $                    ,DTMAX,TMAX
C       WRITE(6,6005) QCMAX,QDMAX    
C       WRITE(6,6006) QCEMAX,QDEMAX  
C       WRITE(6,6007) QCWMAX,QDWMAX   
C       WRITE(6,6008) QCNMAX,QDNMAX  
C       WRITE(6,6009) QCSMAX,QDSMAX  
C       WRITE(6,6010) QCTMAX,QDTMAX   
C       WRITE(6,6011) QCBMAX,QDBMAX  
C       WRITE(6,6012) QDCMAX 

C      ENDIF                                                      
 6004 FORMAT(1X,'GRID!:',I3,' (I,J,K)= (',I4,',',I4,',',I4,'):'
     $                ,'     DT=',E14.6,'   T= ',F8.4)
 6005 FORMAT(41X,'QC=',E14.6,'  QD= ',E14.6)    
 6006 FORMAT(40X,'QCE=',E14.6,' QDE= ',E14.6)
 6007 FORMAT(40X,'QCW=',E14.6,' QDW= ',E14.6)
 6008 FORMAT(40X,'QCN=',E14.6,' QDN= ',E14.6)  
 6009 FORMAT(40X,'QCS=',E14.6,' QDS= ',E14.6)
 6010 FORMAT(40X,'QCT=',E14.6,' QDT= ',E14.6)
 6011 FORMAT(40X,'QCB=',E14.6,' QDB= ',E14.6) 
 6012 FORMAT(59X,'QDC= ',E14.6)   
       
 6013 FORMAT(1X,'  TP=',F8.4,'  TE=',F8.4,'  TW=',F8.4
     $         ,'  TN=',F8.4,'  TS=',F8.4,'  TT=',F8.4,'  TB=',F8.4)     
      RETURN
      END
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C#endif
