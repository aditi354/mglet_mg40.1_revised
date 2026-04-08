










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
      SUBROUTINE TSTSCA4(KK,JJ,II,KMX,JMX,IMX,
     $                  X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                  RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                  U,V,W,T,TO,P,G,B,DT,WCT,
     $                  BP,BU,BV,BW,
     $                  LTXSTART,TXSTART,NBUF,ITSTEP,
     $                  WKONSCA,WDIFSCA,WSORSCA,
     $                  IC1,IC2,JC1,JC2,KC1,KC2,
     $                  NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,
     $                  NFROS,NBACS,NRGTS,NLFTS,NBOTS,NTOPS,NCUBS,
     $                  RSGS3,FAKTOR,
     $                  FTI,FTJ,FTK,
     $                  COEFTX,COEFTY,COEFTZ,
     $                  LCOL,DIAG,RCOL,UZ,RSP,IGRID,
     $                  PRTMAX,XPRTMAX,YPRTMAX,ZPRTMAX,
     $                  PRTMIN,XPRTMIN,YPRTMIN,ZPRTMIN
     $                  ,GI,GJ,GK)
C*STARLET***************************************************************
C        T S T S C A 4  RK TIMESTEP - SCALAR TRANSPORT - TURBULENT - LES
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAY DIMENSIONS
C        KMX, JMX, IMX  - TOTAL GRID SIZE IN Z,Y,X DIRECTION
C        DX,DY,DZ       - GRID SPACING (CENTER) 
C        DDX,DDY,DDZ    - CONTROL VOLUME SIZE IN X,Y,Z DIRECTION
C        U(KK,JJ,II)    - NEW VELOCITY FIELD
C        V(KK,JJ,II)    - NEW VELOCITY FIELD
C        W(KK,JJ,II)    - NEW VELOCITY FIELD
C        T(KK,JJ,II)    - NEW TEMPERATURE FIELD
C        TO(KK,JJ,II)   + OLD TEMPERATURE FIELD
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
C
C        IC1,JC1,KC1    - LEFT BORDER INDICES OF BOUNDING BOX 
C        IC2,JC2,KC2    - RIGHT BORDER INDICES OF BOUNDING BOX          
C                         (=1 LAYER AROUND BODY)
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
C DEFINE DIREKTIVEN     : ZEN, QUD, UPW,EULERSCA, LEAPFSCA, ADBASCA, 
C                         BAPAR_FIXED, TOPAR_FIXED
C UPROG                 : SWCLESCA
C
C        26.05.03 (TB)  : TSTSCA : BASED ON TSTLE4
C                         RUNGE KUTTA TIMESTEP FOR SCALAR TRANSPORT
C                         
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
     $       TO(KK,JJ,II), WCT(KK,JJ,II)
C
      REAL        P(KK,JJ,II),   G(KK,JJ,II)
      REAL        B(KK,JJ,II)
      REAL        GI(KK,JJ,II),GJ(KK,JJ,II),GK(KK,JJ,II)
      REAL  BU(KK,JJ,II),BV(KK,JJ,II),BW(KK,JJ,II),BP(KK,JJ,II)

      LOGICAL LTXSTART
      REAL TXSTART 
 
      REAL      RSGS3(KK,JJ,II)

C----------------DIMENSIONIERUNG NICHT ALLGEMEIN GUELTIG---------
C      REAL      COEFFX(II,12),COEFFY(JJ,12),COEFFZ(KK,12),
C     $          COEFDX(II,12),COEFDY(JJ,12),COEFDZ(KK,12)
C
C      REAL      LCOL(II),DIAG(II),RCOL(II),UZ(II),RSP(II),
C     $          FAKTOR(II)
C-----------------------------------------------------------------

      REAL        FTI(1)      
      REAL        FTJ(1)
      REAL        FTK(1)
C
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

      HP    =  0.5
      HN    = -0.5


C
C                                 BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 UND BEI GITTERKOPPLUNG
C                                 SOWIE BEI AUSFLUSSBEDINGUNG
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT

C      NBU = 0
C      IF ((NBAC.EQ.1).OR.(NBAC.EQ.7))  NBU = 1
C      IF ((NBAC.EQ.3).OR.(NBAC.EQ.4))  NBU = 1
C#ifndef _BAPAR_FIXED_
C      IF (NBAC.EQ.8)  NBU = 1
C#endif
C      NLV = 0
C      IF ((NLFT.EQ.1).OR.(NLFT.EQ.7).OR.(NLFT.EQ.8))  NLV = 1
C      IF ((NLFT.EQ.3).OR.(NLFT.EQ.4))  NLV = 1
C      NTW = 0
C      IF ((NTOP.EQ.1).OR.(NTOP.EQ.7).OR.(NTOP.EQ.3))  NTW = 1
C#ifndef _TOPAR_FIXED_
C      IF (NTOP.EQ.8)  NTW = 1
C#endif
C
C                                 BEI LOKAL VERFEINERTEM GITTER
C                                 WIRD DIE NORMALKOMPONENTE DER ERSTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT

C      NFU = 0
C      IF ( NFRO .EQ. 8 .OR. NFRO .EQ. 3 ) NFU = 1
C      NRV = 0
C      IF ( NRGT .EQ. 8 .OR. NRGT .EQ. 3 ) NRV = 1
C      NBW = 0
C      IF ( NBOT .EQ. 8 .OR. NBOT .EQ. 3 ) NBW = 1
C***********************************************************************
C*     DECONVOLUTION OF FIELD FOR IMPROVED RESOLUTION KONVEKTIVE TERM  *
C***********************************************************************
C
C***********************************************************************
C           HIER KOMPAKTVERFAHREN VIERTER ORDNUNG IN X-RICHTUNG        *
C***********************************************************************
CC                                PERIODISCHE RANDBEDINGUNGEN 
CC                                NUR AEQUIDISTANTE GITTER !! 
      IF (NFRO .EQ.1) THEN
C

CC                                NICHT-PERIODISCHE RANDBED.
      ELSE


      ENDIF


C     INITIALIZE PRTMAX,PRTMIN AND XYZPRTMINMAX
      PRTMAX = 0.0
      PRTMIN = GREAT
      XPRTMAX= 0.0
      YPRTMAX= 0.0
      ZPRTMAX= 0.0  
      XPRTMIN= 0.0
      YPRTMIN= 0.0
      ZPRTMIN= 0.0
      IF(ABS(WKONSCA) .LE. SMALL) GOTO 2200
C
C***********************************************************************
C                                 10 10 10 10 10 10 10 10 10 10 10 10 10
      DO  10  I = ISTART, ISTOP
C
C***********************************************************************
C           HIER KOMPAKTVERFAHREN VIERTER ORDNUNG IN Y-RICHTUNG        *
C***********************************************************************
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

              ATZ   =  DDXI*DDYJ
                

C                                 **************************************
C                                 KONVEKTIVE TERME   KONVEKTIVE TERME
C                                 **************************************
C
C                             KOEFFIZIENTEN BEIM KONVEKTIVEN ZEITSCHRITT
C


		FKDTT = -1.0*RDDY(J)* RDDX(I)* WKONSCA

C
C                                 **************************************
C                                 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
C                                 NUR BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 UND BEI GITTERKOPPLUNG
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT

C      IF (I .EQ. ISTOP)  GOTO 2110
      IF (I .EQ. ISTART) GOTO 2110
      IF (J .EQ. JSTART) GOTO 2110

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
                ATX   = DDYJ*DDZK
                ATY   = DDXI*DDZK
C                                 KONVEKTIVE TERME
C---------------------- Fred Body
C#ifdef 
C      FACBE = MAX(BP(K,J,I),BP(K,J,I+1))
C      FACBW = MAX(BP(K,J,I),BP(K,J,I-1))
C      FACBN = MAX(BP(K,J,I),BP(K,J+1,I))
C      FACBS = MAX(BP(K,J,I),BP(K,J-1,I))
C      FACBT = MAX(BP(K,J,I),BP(K+1,J,I))
C      FACBB = MAX(BP(K,J,I),BP(K-1,J,I))
C#else
      FACBE = 1.0 
      FACBW = 1.0
      FACBN = 1.0
      FACBS = 1.0
      FACBT = 1.0
      FACBB = 1.0
C#endif

C                         CENTRAL DIFFERENCING / KOMPAKTE INTERPOLATION 
C                                 T IN I,J UND K DIRECTION
      FTE  =     ATX*U(K,J,I)
      FTW  =     ATX*U(K,J,I-1)
      FTN  =     ATY*V(K,J,I)  
      FTS  =     ATY*V(K,J-1,I)
      FTT  =     ATZ*W(K,J,I)  
      FTB  =     ATZ*W(K-1,J,I)



      QCTE = 0.5*FTE*(T(K,J,I)   + T(K,J,I+1))*FACBE
      QCTW = 0.5*FTW*(T(K,J,I-1) + T(K,J,I)  )*FACBW


      QCTN = 0.5*FTN*(T(K,J,I)   + T(K,J+1,I))*FACBN
      QCTS = 0.5*FTS*(T(K,J-1,I) + T(K,J,I)  )*FACBS

      QCTT = 0.5*FTT*(T(K,J,I)   + T(K+1,J,I)) *FACBT
      QCTB = 0.5*FTB*(T(K-1,J,I) + T(K,J,I)  )*FACBB
C
C
C                                 ZEITSCHRITT
C
      TO(K,J,I)    = FKDTT * RDDZK * (QCTE-QCTW+QCTN-QCTS+QCTT-QCTB)

C
  100 CONTINUE
C-----------------------------------------------------------------------
 2110 CONTINUE
   20 CONTINUE
   10 CONTINUE
 2200 CONTINUE



      CALL SWCLESCA  (KK,JJ,II,KMX,JMX,IMX,DX,DY,DZ,DDX,DDY,DDZ,
     $                 U,V,W,T,G,B,
     $                 WCT,
     $                 IC1,IC2,JC1,JC2,KC1,KC2,
     $                 NFROS,NBACS,NRGTS,NLFTS,NBOTS,NTOPS,NCUBS,
     $                 RDDX,RDDY,RDDZ,IGRID,DDS,SM,SP,FTI,FTJ,FTK)

      DO  11  I = ISTART, ISTOP

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
      DO  22  J = JSTART, JSTOP
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

              ATZ   =  DDXI*DDYJ
C                                 **************************************
C                                 DIFFUSIVE TERME   DIFFUSIVE TERME   DI
C                                 **************************************
C
C                             KOEFFIZIENTEN BEIM DIFFUSIVEN ZEITSCHRITT
C
		FDDTT = -1.0 *RDDY(J)* RDDX(I)*WDIFSCA
C
C                                 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C                                 --------------------------------------
C
C
C                                 NUR BEI PERIODISCHEN RANDBEDINGUNGEN
C                                 WIRD DIE NORMALKOMPONENTE DER LETZTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT
C
C      IF (I .EQ. ISTOP)  GOTO 2210
      IF (I .EQ. ISTART) GOTO 2210
      IF (J .EQ. JSTART) GOTO 2210


C
      DO 200  K = KSTART, KSTOP
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

C
C     ONLY "FREE" CELL FLUXES ARE CONSIDERED AS FOLLOWS. WALL CELL 
C     CONTRIBUTION IS TAKEN CARE OF IN THE WALL CORRECTION TERM 
C
C     FACS = 1 WHEN CELL FACE IS FREE (NO BOUNDING WALL), OTHERWISE = 0
c      FACBE = MAX(BP(K,J,I),BP(K,J,I+1))
c     FACBW = MAX(BP(K,J,I),BP(K,J,I-1))
c     FACBN = MAX(BP(K,J,I),BP(K,J+1,I))
c     FACBS = MAX(BP(K,J,I),BP(K,J-1,I))
c     FACBT = MAX(BP(K,J,I),BP(K+1,J,I))
c     FACBB = MAX(BP(K,J,I),BP(K-1,J,I))

C      FACBE = BP(K,J,I)*BP(K,J,I+1)
C      FACBW = BP(K,J,I)*BP(K,J,I-1)
CC      FACBN = BP(K,J,I)*BP(K,J+1,I)
C      FACBS = BP(K,J,I)*BP(K,J-1,I)
C      FACBT = BP(K,J,I)*BP(K+1,J,I)
C      FACBB = BP(K,J,I)*BP(K-1,J,I)
C
C                                     HIER KOMPAKTE BERECHNUNG
C                                     DER ERSTEN ABLEITUNG U
      QDTE = -GSCAE * ATX*RDXI   * (T(K,J,I+1) - T(K,J,I))   * FACBE
      QDTW = -GSCAW * ATX*RDXIM  * (T(K,J,I)   - T(K,J,I-1)) * FACBW 
      QDTN = -GSCAN * ATY*RDYJ   * (T(K,J+1,I) - T(K,J,I))   * FACBN
      QDTS = -GSCAS * ATY*RDYJM  * (T(K,J,I)   - T(K,J-1,I)) * FACBS

      QDTT = -GSCAT * ATZ*RDZK   * (T(K+1,J,I) - T(K,J,I))   * FACBT   
      QDTB = -GSCAB * ATZ*RDZKM  * (T(K,J,I)   - T(K-1,J,I)) * FACBB
C
C                                 ZEITSCHRITT
C
       QSUMD = FDDTT*RDDZK*(QDTE-QDTW+QDTN-QDTS+QDTT-QDTB-WCT(K,J,I))      
      TO(K,J,I) = TO(K,J,I) + QSUMD
C
  200 CONTINUE
C
 2210 CONTINUE
C
C 2300 IF(ABS(WSORSCA) .LE. SMALL) GOTO 2400
 2300 CONTINUE
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
CC      IF (I .EQ. ISTOP)  GOTO 2310
C      IF (I .EQ. ISTART) GOTO 2310
C      IF (J .EQ. JSTART) GOTO 2310
C-----------------------------------------------------------------------
C      DO 300  K = KSTART, KSTOP
C
C
C                                 **************************************
C                                 ZEITSCHRITT
C
C      TO(K,J,I) = TO(K,J,I) - 1.0/(RHO*LAMDA) * Q(K,J,I) * WSORSCA
C
C  300 CONTINUE
C-----------------------------------------------------------------------
C 2310 CONTINUE
C
 2400 CONTINUE
 

   22 CONTINUE   
C              
C=======================================================================

   11 CONTINUE   
C              
C=======================================================================
C
      RETURN
      END
C#endif
