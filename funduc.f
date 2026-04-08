










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
      FUNCTION CH1680  (TEXTIN)
C*STARLET***************************************************************
C        C H 1 6 8 0    DIE ERSTEN 16 CHARACTERS DER CHARACTER (LEN=80)
C                       VARIABLEN CH1680 WERDEN MIT 'TEXTIN' BELEGT;
C                       CH1680(17:80) WERDEN MIT BLANKS GEFUELLT.
C*STARLET***************************************************************
C
C PARAM: TEXTIN         - CHARACTER (LEN=16)  VARIABLE
C
C VERS:  17.12.86 (HW)  : ORIGINAL
C
C DEFINE-DIREKTIVEN     : KEINE
C
C*STARLET***************************************************************
C
      CHARACTER (LEN=1)   BLANK
      CHARACTER (LEN=16)  TEXTIN
      CHARACTER (LEN=80)  CH1680
C
      DATA           BLANK  /' '/
C
      DO 100 N = 1,80
  100    CH1680(N:N) = BLANK
C
      CH1680(1:16) = TEXTIN(1:16)
C
      RETURN
      END
      FUNCTION GRADPP  (UQUER,DDS)
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
 2010 GRADPP = 2.0*UQUERN/DDS * VZ
      RETURN
C
      END
      FUNCTION GRAPIO  (UQUER,DDS)
C*STARLET***************************************************************
C        G R A P I O    BERECHNUNG DES GRADIENTEN D(UQUER)/D(DDS)
C                       DAS VORZEICHEN ENTSPRICHT DEM VON UQUER !
C                       DER GRADIENT WIRD ALS INTEGRALER MITTELWERT VON
C                       Z=ZMATCH (GRENZE DER VISK. UNTERSCHICHT) BIS
C                       Z=DDS ERMITTELT. DER STEILE GRADIENT IN DER VIS-
C                       KOSEN UNTERSCHICHT BLEIBT DAHER UNBERUECKSICH-
C                       TIGT.
C*STARLET***************************************************************
C
C PARAM: UQUER          - GESCHWINDIGKEITSKOMPONENTE AM WANDNAECHSTEN
C                         GITTERPUNKT (UQUER IST EIN MITTELWERT UEBER
C                         DIE MASCHENFLAECHE SENKRECHT ZU UQUER)
C        DDS            - KANTENLAENGE DER WANDNAECHSTEN ZELLE
C                         (SENKRECHT ZU UQUER UND ZUR WAND )
C
C VERS:  13.09.85 (HW)  : ORIGINAL
C        12.07.88 (HW)  : DAS VORZEICHEN WIRD VON UQUER UEBERNOMMEN
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
 2010 GRAPIO = 2.0*UQUERN/DDS * VZ
      RETURN
C
      END
      FUNCTION UTAUP   (UQUER,DDS)
C*STARLET***************************************************************
C        U T A U P      BERECHNUNG DES BETRAGES DER SCHUBSPANNUNGS-
C                       GESCHWINDIGKEIT MITTELS DES 1/7 POTENZGESETZES
C*STARLET***************************************************************
C
C PARAM: UQUER          - GESCHWINDIGKEITSKOMPONENTE AM WANDNAECHSTEN
C                         GITTERPUNKT (UQUER IST EIN MITTELWERT UEBER
C                         DIE MASCHENFLAECHE SENKRECHT ZU UQUER)
C        DDS            - KANTENLAENGE DER WANDNAECHSTEN ZELLE
C                         (SENKRECHT ZU UQUER UND ZUR WAND )
C
C VERS:  10.09.85 (HW)  : ORIGINAL
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
C
C
C                                 HIER: UQUER LIEGT IN DER VISKOSEN
C                                 UNTERSCHICHT
C
 2010 UTAUP  = SQRT(2.0*GMOL*UQUERN/(RHO*DDS))
      RETURN
C
      END
      FUNCTION TAUWP   (UQUER,DDS)
C*STARLET***************************************************************
C        T A U W P      BERECHNUNG DER WANDSCHUBSPANNUNG MITTELS DES
C                       1/7 POTENZGESETZES. DAS VORZEICHEN DER WAND-
C                       SCHUBSPG. ENTSPRICHT DEM VON UQUER !
C*STARLET***************************************************************
C
C PARAM: UQUER          - GESCHWINDIGKEITSKOMPONENTE AM WANDNAECHSTEN
C                         GITTERPUNKT (UQUER IST EIN MITTELWERT UEBER
C                         DIE MASCHENFLAECHE SENKRECHT ZU UQUER)
C        DDS            - KANTENLAENGE DER WANDNAECHSTEN ZELLE
C                         (SENKRECHT ZU UQUER UND ZUR WAND )
C
C VERS:  10.09.85 (HW)  : ORIGINAL
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
      VZ     = SIGN(1.0,UQUER)
      UQUERN = ABS(UQUER)
C
C
C                                 HIER: UQUER LIEGT IN DER VISKOSEN
C                                 UNTERSCHICHT
C
 2010 TAUWP  = 2.0*GMOL*UQUERN/DDS * VZ
      RETURN
C
      END
      FUNCTION UQUEP   (UTAU,DDS)
C*STARLET***************************************************************
C        U Q U E P      BERECHNUNG DES BETRAGES DER UEBER DDS INTEGRAL
C                       GEMITTELTEN GESCHWINDIGKEIT 'UQUEP' BEI
C                       BEKANNTEM UTAU (1/7 POTENZGESETZ)
C*STARLET***************************************************************
C
C PARAM: UTAU           - WANDSCHUBSPANNUNGSGESCHWINDIGKEIT
C        DDS            - KANTENLAENGE DER WANDNAECHSTEN ZELLE
C                         (SENKRECHT ZU UQUEP UND ZUR WAND )
C
C VERS:  11.09.85 (HW)  : ORIGINAL
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
      UTAUN  = ABS(UTAU)
C
C
C                                 HIER: UQUEP LIEGT IN DER VISKOSEN
C                                 UNTERSCHICHT
 2010 UQUEP  = UTAUN**2 * 0.5*DDS * RHO/GMOL
      RETURN
C
      END
      FUNCTION DPRVD   (DWALL,UTAU)
C*STARLET***************************************************************
C        D P R V D      BERECHNUNG DER MISCHUNGSWEGLAENGE NACH
C                       PRANDTL (L = CAPPA*Y) MIT BERUECKSICHTIGUNG
C                       DER VAN DRIEST''SCHEN DAEMPFUNGSFUNKTION
C*STARLET***************************************************************
C
C PARAM: DWALL          - ABSTAND DES BETRACHTETEN PUNKTES NORMAL
C                         ZUR WAND
C        UTAU           - WANDSCHUBSPANNUNGSGESCHWINDIGKEIT
C                         (GEBILDET AUS DER VEKTORIELLEN SUMME BEIDER
C                         KOMPONENTEN PARALLEL ZUR WAND)
C
C VERS:  17.10.85 (HW)  : ORIGINAL
C        22.10.85 (HW)  : AENDERUNGEN
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


      COMMON /CONLES/  CONV2S,CMUE,CAPPA,ECONST
      SAVE   /CONLES/
C
C
      UTAUN = ABS(UTAU)
C
      DPRVD = CAPPA*DWALL
     $      * (1.0 - EXP(-DWALL*UTAUN*RHO/(GMOL*26.0)))
C
      RETURN
      END
      FUNCTION I2DO2 (X,Y,X1,X2,Y1,Y2,U1,U2,U3,U4)
C*MGLET***************************************************************
C      I 2 D O 2    INTERPOLATION 2D, ORDNUNG 2
C
C                X,Y:    POSITION IN DIE HINEININTERPOLIERT WIRD
C
C            U1,X1,Y1....:   FUNKTIONSWERT UND POSITIONEN:
C                                                         
C                                                         
C              U3(X1,Y2)                    U4(X2,Y2)                                           
C                               INT(X,Y)                          
C                                                         
C                                                         
C                                                         
C              U1(X1,Y1)                    U2(X2,Y1)                                           
C                                                         
C
C     VERS:  02.02.96 (MM)  :ORIGINAL
C
C*MGLET***************************************************************
C
          REAL I2DO2,X,Y,X1,X2,Y1,Y2,U1,U2,U3,U4
C

        DXE = ABS(X2 - X )
        DXW = ABS(X  - X1)
        DX  = DXE + DXW

        DYN = ABS(Y2 - Y )
        DYS = ABS(Y  - Y1)
        DY  = DYN + DYS

          I2DO2 =  1./(DX*DY) *
     +                              (  U1*DXE*DYN
     +                               + U2*DXW*DYN
     +                               + U3*DXE*DYS
     +                               + U4*DXW*DYS )
C
C*MGLET***************************************************************
C
      RETURN
      END

      FUNCTION I3DO2 (X,Y,Z,X1,X2,Y1,Y2,Z1,Z2,
     $                     U1,U2,U3,U4,U5,U6,U7,U8)
C*MGLET***************************************************************
C      I 3 D O 2    INTERPOLATION 3D  ORDNUNG 2
C
C                X,Y,Z:    POSITION IN DIE HINEININTERPOLIERT WIRD
C
C            U1,X1,Y1,Z1....:   FUNKTIONSWERTE UND POSITIONEN:
C                                                         
C                 U7(X1,Y2,Z2)--------------------U8(X2,Y2,Z2)                    
C                /|                              /|          
C              U5(X1,Y1,Z2)--------------------U6(X2,Y1,Z2)                                           
C              |  |              INT(X,Y,Z)    |  |                    
C              |  |                            |  |          
C              |  |                            |  |          
C              |  |                            |  |          
C              |  U3(X1,Y2,Z1)-----------------|--U4(X2,Y2,Z1)                    
C              | /                             | /          
C              U1(X1,Y1,Z1)--------------------U2(X2,Y1,Z1)                                           
C                                                         
C
C     VERS:  02.02.96 (MM)  :ORIGINAL
C
C*MGLET***************************************************************
C
          REAL I3DO2,X,Y,Z,X1,X2,Y1,Y2,Z1,Z2,
     $                     U1,U2,U3,U4,U5,U6,U7,U8
C

        DXE = ABS(X2 - X )
        DXW = ABS(X  - X1)
        DX  = DXE + DXW

        DYN = ABS(Y2 - Y )
        DYS = ABS(Y  - Y1)
        DY  = DYN + DYS

        DZT = ABS(Z2 - Z )
        DZB = ABS(Z  - Z1)
        DZ  = DZT + DZB

          I3DO2 =  1./(DX*DY*DZ) *
     +                      ( U1*DXE*DYN*DZT
     +                       +U2*DXW*DYN*DZT
     +                       +U3*DXE*DYS*DZT
     +                       +U4*DXW*DYS*DZT
     +                       +U5*DXE*DYN*DZB
     +                       +U6*DXW*DYN*DZB
     +                       +U7*DXE*DYS*DZB
     +                       +U8*DXW*DYS*DZB )

C
C*MGLET***************************************************************
C
      RETURN
      END

      FUNCTION I1DO3 (C,A,B,D,E,Q,R,S,T)
C**MGLET********************************************************
C         I 1 D O 3   FUNCTION FOR THIRD-ORDER-INTERPOLATION
C                     DERIVED BY TAYLOR-SERIES EXPANSION IN
C                     PROGRAMM "DERIVE" ON A PC
C
C        C:           INTERPOLATED VALUE F AT (X)
C        A:           F AT X - (R+Q)
C        B:           F AT X - R
C        D:           F AT X + S
C        E:           F AT X + (S+T)
C
C   VERSION: 08.02.96 (MM)    ORIGINAL
C
C**MGLET********************************************************


      C = 
     $ -R*S*A*(S+T)/((Q+R+S+T)*(Q*(S+T)+R*T))
     $ +  S*B*(S+T)*(Q+R)/((Q*(S+T)+R*T)*(R+S))
     $ +  R*D*(S+T)*(Q+R)/((Q*(S+T)+R*T)*(R+S))
     $ -R*S*E*(Q+R)/((Q+R+S+T)*(Q*(S+T)+R*T))
      I1DO3 = C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION QWL   (DT,UQ,DDS,DS)
C*STARLET***************************************************************
C        Q W           CALCULATION OF SCALAR WALL FLUX 
C*STARLET***************************************************************
C
C PARAM: 
C        DT          - TEMPERATURE DIFFERENCE T-TWALL
C        UQ          - TANGENTIAL VELOCITY AT WALL NEAREST CELL CENTER
C        DDS         - DISTANCE Of WALL NEARES CELL CENTER FROM THE WALL
C
C VERS:  04.11.02 (TB)  : ORIGINAL
C
C DEFINE-DIREKTIVEN     : DNS
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
      UQN = ABS(UQ)
      RPRMOL = 1.0/PRMOL
      RHOTAU= (RHO * TAUWIN(UQN,DDS))**0.5
      SM = GMOL / RHOTAU**0.5 * CWA**(1.0/(1.0-CWB))
      FSM = 1/GMOL *( 0.5*SM**2.0*PRMOL + 
     $                  (RHOTAU/GMOL)**(-CPO1) * CWA/CPO2 *
     $                  (DS**CPO2 - SM**CPO2) )      
C
C
C                                 UQ IN VISCOUS SUBLAYER
C                                 TW AT WALL BOUNDARY FACE => DS/2
C
 3011 QWL     = GMOL / PRMOL * DT / (DS*0.5)
      RETURN
C
 3030 CALL ERRR (502,'FUNDUC QWL')
C
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    
      FUNCTION PRT(GTGMOL)
C*STARLET***************************************************************
C        P R T          CALCULATION OF LOCAL TURBULENT PRANDTL NUMBER
C                       IF PRTURB > 0 THEN PRTURB=CONST IS USED
C                       IF PRTURB <== THEN KAYS/CRAWFORD FORMULA IS USED
C
C                       LIMITATIONS OF KAYS/CRAWFORD: 0.5<PRMOL<7, 
C                                                     ANY RE
C                                                     ANY DPDX
C
C                      KAYS/CRAWFORD PROVIDES A RELATIVELY HIGH PRTURB 
C                      NEAR THE WALL (IN THE SUBLAYER) BUT APPROACHES 
C                      0.85 AS Y+ INCREASES INTO THE LOG LAYER 
C 
C*STARLET***************************************************************
C
C PARAM: GTGMOL         = GTURB/GMOL
C
C VERS:  24.10.2002 (TB) : ORIGINAL
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

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
      IF(PRTURB .LE. SMALL) CALL ERRR (503,'FUNDUC PRT')
C
C                              CONSTANT PRTURB
C
      PRT = PRTURB
C
      RETURN
C
      END
