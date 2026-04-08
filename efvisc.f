










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

      SUBROUTINE EFVISC  (KK,JJ,II,KMX,JMX,IMX,X,Y,Z,
     $                    DX,DY,DZ,DDX,DDY,DDZ,
     $                    RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                    U,V,W,G,B,
     $                    DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY,
     $                    GMOL,RHO,UGRID,
     $                    IC1,IC2,JC1,JC2,KC1,KC2,
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,CONV1S, 
     $   HILF0,HILF1,HILF2,HILF3,HILF4,HILF5,HILF6,
     $   LUX,LUY,LUZ,LVY,LVZ,LWZ,MUX,MUY,MUZ,MVY,MVZ,MWZ,
     $   FA11,FA12,FA13,FA21,FA22,FA23,FA31,FA32,FA33,
     $   DUDX3,DUDY3,DUDZ3,DVDX3,DVDY3,DVDZ3,
     $   DWDX3,DWDY3,DWDZ3,NORMS,SDUDX3,SDUDY3,SDUDZ3,
     $   SDVDY3,SDVDZ3,SDWDZ3,
     $   UC,VC,WC,UUC,UVC,UWC,VVC,VWC,WWC,CDELTA)
C*STARLET***************************************************************
C        E F V I S C      BERECHNUNG DER EFFEKTIVEN DYNAMISCHEN
C                         VISKOSITAET
C                         EINBAUTEN, DIE AUS SAEULEN AUFGEBAUT WERDEN
C                         KOENNEN, WERDEN BERUECKSICHTIGT. DIE GRADIEN-
C                         TEN AN FESTEN (NOSLIP-) WAENDEN WERDEN MITTELS
C                         DES 1/7 POTENZGESETZES BESTIMMT.
C                         ZWISCHEN ZWEI NOSLIP WAENDEN MUESSEN MINDES-
C                         TENS  Z W E I  ZELLEN LIEGEN!
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        X, Y, Z        - KOORDINATEN DER ZELLDEFINITIONSPUNKTE
C        DX,DY,DZ       - ABSTAND DER GITTERPUNKTE
C        DDX,DDY,DDZ    - KANTENLAENGEN DER KONTROLLVOLUMINA
C        U(KK,JJ,II)    - GESCHWINDIGKEITSFELD
C        V(KK,JJ,II)    - GESCHWINDIGKEITSFELD
C        W(KK,JJ,II)    - GESCHWINDIGKEITSFELD
C        G(KK,JJ,II)    + EFFEKTIVE DYNAMISCHE VISKOSITAET (= MUE)
C        GMOL           - MOLEKULARE DYNAMISCHE VISKOSITAET
C        RHO            - DICHTE (= CONST)
C        UGRID          - GESCHWINDIGKEIT DES KOORDINATENSYSTEMS IN
C                         XeRICHTUNG (GALILEI-TRANSFORMATION)
C        IC1,JC1,KC1    - LINKE  RAENDER DER BOUNDING BOX
C        IC2,JC2,KC2    - RECHTE RAENDER DER BOUNDING BOX
C
C UPROG                 : ERRR, BVISK
C
C VERS:  08.01.86 (HW)  : ORIGINAL
C        26.03.86 (HW)  : BVISK EINGEFUEHRT
C        20.08.86 (HW)  : UEBERGABE WEGEN EINFUEHRUNG VON BVISKK UM
C                         IB1,IB2,... ERWEITERT.
C        17.09.86 (HW)  : DDMIX4 KANN DEFINIERT WERDEN. DER MISCHUNGS-
C                         WEG WIRD AUS DEM VERHAELTNIS VON ZELLVOLUMEN
C                         ZU ZELLOBERFLAECHE GEBILDET.
C        12.07.88 (HW)  : VORZEICHENKORREKTUR DER GRADIENTEN IN DER
C                         NAEHE FESTER WAENDE (WEGEN AUSWERTUNG VON
C                         OMEGA-X, -Y UND -Z NOETIG GEWORDEN)
C        17.03.92 (MM)  : LAENGENMASS DMIX UND DMIXG KOMMT UEBER
C                         B-FELD HEREIN
C        18.03.92 (MM)  : FELDER DUDY.....DWDY EINGEFUEHRT
C         4. 4.92 (MM)  : SCHLEIFEN ZUR WANDKORREKTUR AM KOERPER
C                         LAUFEN NUR INNERHALB DER
C                         BOUNDING BOX (IC1,IC2,JC1,JC2,KC1,KC2)
C        17. 4.93 (MM)  : FEHLER IN _LENOS_ BESEITIGT
C                 
C        08.02.97 (AO)  : DYNAMISCHES SUBGRIDSCALE-MODELL EINGEBAUT
C                         NACH GERMANO. ROUTINEN UEBERNOMMEN VON
C                         G: BAERWOLFF (HFI).(DYNSGSM)
C                         ORIGINALROUTINEN AUS STANFORD
C
C DEFINE-DIREKTIVEN     : DDMIX2, DDMIX3, DDMIX4, DYNSGSM
C
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

      COMMON /CONLES/  CONV2S,CMUE,CAPPA,ECONST
      SAVE   /CONLES/
C
      CHARACTER (LEN=1)  ITYP
C
      REAL        X(II),         Y(JJ),         Z(KK),
     $           DX(II),        DY(JJ),        DZ(KK),
     $          DDX(II),       DDY(JJ),       DDZ(KK),
     $          RDX(II),       RDY(JJ),       RDZ(KK),
     $         RDDX(II),      RDDY(JJ),      RDDZ(KK),
     $            U(KK,JJ,II),   V(KK,JJ,II),   W(KK,JJ,II),
     $            G(KK,JJ,II),   B(KK,JJ,II),rho,gmol
C
      REAL        DUDY(KK,JJ),  DUDZ(KK,JJ),
     $            DVDX(KK,JJ),  DVDZ(KK,JJ),
     $            DWDX(KK,JJ),  DWDY(KK,JJ),ugrid,
     $            mineddyv,maxeddyv,maxs,maxcon,
     $  third,smal2,econst,cap05,rindef,smaone,small,great,cmue,
     $  cappa,conv2s,preset

      REAL CONV1S(II)
C
      integer KK,JJ,II,KMX,JMX,IMX,Ic1,ic2,jc1,jc2,kc1,kc2,
     $        km2,jm2,im2,i,j,k,i0,j0,k0

      integer  istart,jstart,kstart,iend,jend,kend

      real     gradpp,ddzk,rddzpl,dzf,dyf2,ddyj,
     $         rddypl,dyf,dwdz,sabs,dmixg,dvdy,dzf2,vorz,dudx,
     $         dxf2,cmaxr,ddxi,rddxpl,dxf,cmaxs      
      external gradpp 


      real 
     $  cdelta(1)
                           
      real
     $  hilf0(1),hilf1(1),hilf2(1),
     $  hilf3(1),hilf4(1),hilf5(1),
     $  hilf6(1),
     $  lux(1),luy(1),luz(1),
     $  lvy(1),lvz(1),lwz(1),
     $  mux(1),muy(1),muz(1),
     $  mvy(1),mvz(1),mwz(1),
     $  fa11(1),fa12(1),fa13(1),
     $  fa21(1),fa22(1),fa23(1),
     $  fa31(1),fa32(1),fa33(1)

      real 
     $ dudx3(1),dudy3(1),dudz3(1),
     $ dvdx3(1),dvdy3(1),dvdz3(1),
     $ dwdx3(1),dwdy3(1),dwdz3(1),
     $ norms(1),
     $ sdudx3(1),sdudy3(1),sdudz3(1),
     $ sdvdy3(1),sdvdz3(1),sdwdz3(1),
     $ uc(1),vc(1),wc(1),
     $ uuc(1),uvc(1),uwc(1),
     $ vvc(1),vwc(1),wwc(1)

      real
     $        fac1,fac2,fac3,asq,  ! related to filtering
     $        sux,suy,suz,svy,svz,swz,r,s 
      integer nf_ratio,nf_type,i1,j1,k1,i2,j2,k2,i3,j3,k3,
     $        i4,j4,k4



      INTEGER  ie0,je0,ke0
C
      NBND = 2

      KM2    = KMX-NBND
      JM2    = JMX-NBND
      IM2    = IMX-NBND

      kstart = 1+NBND
      jstart = 1+NBND
      istart = 1+NBND       

       iend  = im2
       jend  = jm2
       kend  = km2  
C
      KM2    = KMX-2
      JM2    = JMX-2
      IM2    = IMX-2
C
      CAP05  = 0.5 * CAPPA
      SMAL2  = 2.0 * SMALL
      THIRD  = 1.0 / 3.0

      maxeddyv=-10000000000.0
      mineddyv= 10000000000.0
      maxs    =0.0
      maxcon  =0.0 
      cmaxr   =0.0
      cmaxs   =0.0 

      i0=0
      j0=0
      k0=0  

C
C                                   AUF G-FELD LIEGT NUR MOLEKULARE
C                                   VISKOSITAET, FUER DIREKTE SIMULATION
C
      DO I=3,IM2
      DO J=3,JM2
      DO K=3,KM2

             G(K,J,I) = GMOL

      ENDDO
      ENDDO
      ENDDO
C

C     
C
      RETURN
      END

      FUNCTION GRADP2  (UQUER,DDS)
C*STARLET***************************************************************
C        G R A D P P    BERECHNUNG DES GRADIENTEN D(UQUER)/D(DDS)
C                       DER GRADIENT WIRD ALS PUNKTWERT AM PUNKT Z=DDS/2
C                       BESTIMMT. DAS VORZEICHEN ENTSPRICHT DEM VON
C                       UQUER !
C*STARLET***************************************************************
C
C PARAM: UQUER          - GESCHWINDIGKEITSKOMPONENTE AM WANDNAECHSTEN
C                         GITTERPUNKT (UQUER IST EIN MITTELWERT UEBER
C                         DIE MASCHENFLAECHE SENKRECHT ZU UQUER)
C        DDS            - KANTENLAENGE DER WANDNAECHSTEN ZELLE
C                         (SENKRECHT ZU UQUER UND ZUR WAND )
C
C VERS:  10.09.85 (HW)  : ORIGINAL
C        14.01.86 (HW)  : VEREINFACHUNG DES AUSDRUCKES
C        11.07.88 (HW)  : DAS VORZEICHEN WIRD VON UQUER UEBERNOMMEN
C                         (NOTWENDIG WEGEN DER AUSWERTEROUTINEN FUER
C                          DIE KOMPONENTEN DER VORTICITY)
C
C DEFINE-DIREKTIVEN     : NATWBC
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


      COMMON /CONPOT/  CWA,CWB,CPO1,CPO2,CPO3,CPO4,CPO5,CPO6,CPO7,CPO8,
     $                 CPO9,CPO10,CPO11,CPO12
      SAVE   /CONPOT/
C
C
      UQUERN = ABS(UQUER)
      VZ     = SIGN(1.0,UQUER)
C
C
C                                 HIER: UQUER LIEGT IN DER VISKOSEN
C                                 UNTERSCHICHT
 2010 GRADP2 = 2.0*UQUERN/DDS * VZ
      RETURN
C
      END

















