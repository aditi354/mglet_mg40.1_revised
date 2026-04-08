










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
      SUBROUTINE BFROMG
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,    
     $              U,V,W,P,G,B,UFR,VFR,WFR,HILF,
     $              ITYP,IGRID,IRB,JRB,KRB,TIMEPH,PFR,GFR
     $             )     
C*STARLET***************************************************************
C        B F R O M G      SETZEN DER RANDBEDINGUNGEN FUER DIE
C                         GESCHWINDIGKEITSFELDER.
C                         (LARGE-EDDY-SIMULATION)
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        X, Y, Z        - KOORDINATEN DER ZELLDEFINITIONSPUNKTE
C        DX,DY,DZ       - ABSTAND DER GITTERPUNKTE
C        DDX,DDY,DDZ    - KANTENLAENGEN DER KONTROLLVOLUMINA
C        U(KK,JJ,II)    + GESCHWINDIGKEITSFELD
C        V(KK,JJ,II)    + GESCHWINDIGKEITSFELD
C        W(KK,JJ,II)    + GESCHWINDIGKEITSFELD
C        T(KK,JJ,II)    + SCALAR T FIELD
C        P(KK,JJ,II)    - DRUCKFELD
C        G(KK,JJ,II)    - EFFEKTIVE DYNAMISCHE VISKOSITAET (= MUE)
C        TK(KK,JJ,II)   - TURBULENZENERGIE
C        UFR(KK,JJ, 2)  - ENTHAELT DIE I=1 U. I=2 EBENEN DES U-PROFILS
C        VFR(KK, 2,II)  -    "      "  J=1 U. J=2    "    "  V-PROFILS
C        WFR( 2,JJ,II)  -    "      "  K=1 U. K=2    "    "  W-PROFILS
C        TFR(KK,JJ, 2)  - ENTHAELT DIE I=1 U. I=2 EBENEN DES T-PROFILS
C        GMOL           - MOLEKULARE DYNAMISCHE VISKOSITAET
C        RHO            - DICHTE (= CONST)
C        UGRID          - GESCHWINDIGKEIT DES KOORDINATENSYSTEMS IN
C                         X-RICHTUNG (GALILEI-TRANSFORMATION)
C        XTOT           - LAENGE DES BERECHNUNGSGEBIETES IN X-RICHTUNG
C        YTOT           -    "                "          IN Y-RICHTUNG
C        ZTOT           -    "                "          IN Z-RICHTUNG
C        MTURB          - SCHALTER: 0 = LAMINAR; 1 = TURBULENT
C        ITYP           - HOLLERITH-KONSTANTE:
C                         'P'  VOR DER DRUCKKORREKTUR
C                         'T'  VOR DEM ZEITSCHRITT
C
C VERS:  18.03.86 (HW)  : ORIGINAL
C VERS:   5. 3.92 (MM)  : BEINHALTET BFROL UND BFOLI, 
C                         VOELLIG NEUE VERSION
C                         SCHLEIFEN LAUFEN GRUNDSAETZLICH UEBER GANZES
C                         BERECHNUNGSGEBIET
C        19.12.93 (MM)  : NACHBAR EINGEFUEHRT
C        09.12.98 (AM)  : FRONT-FLAECHE EINGEFUEHRT (BFR,GFR)
C        28.02.03 (TB)  : SCALAR BOUNDARY TREATMENT IMPLEMENTED
C
C DEFINE-DIREKTIVEN     : QUD, 
C
C*STARLET***************************************************************
C

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


      COMMON /COMGRID/
     &               NGRID,    NGRDOLD,   NGRDSET, NGRDDFD,
     &               NAUFP,    NAUFOLD,
     &                 KMX,  JMX,  IMX,
     &                 KMXA, JMXA, IMXA,
     &                IP3D, IP2D, IP1D, IPBB,IPB3, IPBU,
     &                NOF3D,NOF2D,NOF1D,NOFBB,NOFB3,NOFBU,
     &                IPA,   IP1L,  IP2L,
     &                NOFA,  NOF1L, NOF2L,  
     &                 IC1,  IC2,  JC1,  JC2,  KC1,  KC2

      INTEGER
     &         KMX(MAXGRIDS),     JMX(MAXGRIDS),       IMX(MAXGRIDS),
     &        KMXA(MAXGRIDS),    JMXA(MAXGRIDS),      IMXA(MAXGRIDS),
     &        IP3D(MAXGRIDS),    IP2D(MAXGRIDS),      IP1D(MAXGRIDS),
     &        IPBB(MAXGRIDS),    IPB3(MAXGRIDS),      IPBU(MAXGRIDS),
     &        NOF3D,   NOF2D,   NOF1D,    NOFBB,   NOFBU,
     &         IPA(MAXGRIDS),    IP1L(MAXGRIDS),      IP2L(MAXGRIDS),
     &         NOFA,   NOF1L,   NOF2L,    NAUFP,
     &         IC1(MAXGRIDS),     IC2(MAXGRIDS),
     &         JC1(MAXGRIDS),     JC2(MAXGRIDS),
     &         KC1(MAXGRIDS),     KC2(MAXGRIDS)
 

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
      CHARACTER (LEN=1)  ITYP
      CHARACTER (LEN=16) FLOWT
C
      REAL        X(IDIM1D),         Y(IDIM1D),         Z(IDIM1D),
     $           DX(IDIM1D),        DY(IDIM1D),        DZ(IDIM1D),
     $          DDX(IDIM1D),       DDY(IDIM1D),       DDZ(IDIM1D),
     $          U( IDIM3D ),       V( IDIM3D ),       W( IDIM3D ),
     $          P( IDIM3D ),       G( IDIM3D ),       B( IDIM3D ),
     $        UFR(IDIM2D*2),     VFR(IDIM2D*2),     WFR(IDIM2D*2),
     $        PFR(IDIM2D*2),     GFR(IDIM2D*2),
     $       HILF( IDIM3D ),     TIMEPH
C
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

       CALL MGDIMS  (KK,JJ,II,IGRID)
       CALL MGPOINT (IP3,IP2,IP1,IBB,IB3,IBU,IGRID)


CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      DO IBOCOND = 1,NBOCONDS( 1 , IGRID )
CTBC     SCALAR START/STOP INDICES SET IN MGBLOINDSET
         JSTART = LARBOCONDS( 3 , IBOCOND , 1 , IGRID )
         JSTOP  = LARBOCONDS( 4 , IBOCOND , 1 , IGRID )
         KSTART = LARBOCONDS( 5 , IBOCOND , 1 , IGRID )
         KSTOP  = LARBOCONDS( 6 , IBOCOND , 1 , IGRID )


         ITYPBC = ITYPBOCONDS( IBOCOND , 1 , IGRID )

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72 

         IF ( ITYPBC .EQ. 1 ) THEN

              CALL BFRPER (KK,JJ,II,             
     $                      U(IP3),V(IP3),W(IP3),P(IP3),G(IP3),
     $                      JSTART,JSTOP,KSTART,KSTOP,ITYP
     $                    )
     
         ELSEIF ( ITYPBC .EQ.  2 .OR. ITYPBC .EQ. 11  
     $       .OR. ITYPBC .EQ. 12                     ) THEN

              CALL BFRFIX (KK,JJ,II,              
     $                     U(IP3),V(IP3),W(IP3),P(IP3),G(IP3),UFR(IBB),
     $                VFR(IBB),WFR(IBB),JSTART,JSTOP,KSTART,KSTOP,ITYP
     $                    )
     
C        ELSEIF ( ITYPBC .EQ. 3 ) THEN
C
C             CALL BFROP1 (KK,JJ,II,U(IP3),V(IP3),
C    $                     W(IP3),P(IP3),G(IP3),
C    $                      JSTART,JSTOP,KSTART,KSTOP,ITYP)
         ELSEIF ( ITYPBC .EQ. 3 ) THEN

              CALL BFROP1 (KK,JJ,II,
     $                      U(IP3),V(IP3),W(IP3),P(IP3),G(IP3),
     $                      JSTART,JSTOP,KSTART,KSTOP,ITYP
     $                    )       




         ELSEIF ( ITYPBC .EQ. 9 ) THEN
                  FREQ =   FREQB(IBOCOND,1,IGRID)
                    AU =   AUB(IBOCOND,1,IGRID)
                    AV =   AVB(IBOCOND,1,IGRID)
                    AW =   AWB(IBOCOND,1,IGRID)
                  ANIV =   ANIVEAU(IBOCOND,1,IGRID)
		  FLOWT =   FLOWTYP(IBOCOND,1,IGRID)
C
C                   FREQ =   0.0
C                     AU =   0.0                       
C                     AV =   0.0                 
C                     AW =   0.0                   
C                   ANIV =   ANIVEAU(IBOCOND,1,IGRID)
C 		  FLOWT =   'UNIFORM'                  
              CALL BFRBLO  (KK,JJ,II,              
     $                      U(IP3),V(IP3),W(IP3),P(IP3),G(IP3),
     $                      JSTART,JSTOP,KSTART,KSTOP,
     $                      FREQ,AU,AV,AW,ANIV,TIMEPH,
     $                      X(IP1),Y(IP1),Z(IP1),FLOWT,     
     $                      UFR(IBB),VFR(IBB),WFR(IBB)
     $                     )     
C 
         ELSEIF ( ITYPBC .EQ. 5 ) THEN
                                    
              CALL BFRNOS (KK,JJ,II,
     $                      U(IP3),V(IP3),W(IP3),P(IP3),G(IP3),
     $                      JSTART,JSTOP,KSTART,KSTOP,ITYP)


         ELSEIF ( ITYPBC .EQ. 7 ) THEN
C
C                                     BFRCON IS CALLED BY CONFROMG
C
         ELSEIF ( ITYPBC .EQ. 8 ) THEN
C
              CALL BFRPAR (KK,JJ,II,              
     $                      U(IP3),V(IP3),W(IP3),P(IP3),G(IP3),
     $                 UFR(IBB),VFR(IBB),WFR(IBB),PFR(IBB),GFR(IBB),
     $                      JSTART,JSTOP,KSTART,KSTOP,ITYP,IRB,JRB,KRB
     $                     )     
C
C         TURBULENT INFLOW CONDITION (NOT YET IMPLEMENTED A.O)
C        ELSEIF ( ITYPBC .EQ. 10 ) THEN

C             CALL BFRTUR (KK,JJ,II,
C    $                     U(IP3),V(IP3),W(IP3),P(IP3),G(IP3),UFR(IBB),
C    $                     JSTART,JSTOP,KSTART,KSTOP,ITYP)

         ELSE

            CALL ERRR (711,' BFROMG ')

         ENDIF


      ENDDO

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72


      RETURN
      END

