










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
      SUBROUTINE TAUNN  (KK,JJ,II,X,Y,Z,
     $                    DX,DY,DZ,DDX,DDY,DDZ,
     $                    RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                    U,V,W,G,B,
     $                    DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY,
     $                    TAU11,TAU12,TAU13,
     $                       UP,TAU22,TAU23,
     $                       VP,   WP,TAU33)
C***MGLET***************************************************************
C        T A U N N        BERECHNUNG DES NICHT-NEWTONSCHEN
C                         SPANNUNGSTENSORS
C***MGLET***************************************************************
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
C        TAU11..TAU33   - ELEMENTE DES SPANNUNGSTENSORS
C        RHO            - DICHTE (= CONST)
C        UGRID          - GESCHWINDIGKEIT DES KOORDINATENSYSTEMS IN
C                         XeRICHTUNG (GALILEI-TRANSFORMATION)
C
C UPROG                 : ERRR,
C
C VERS:  13.04.99 (MM)  : ORIGINAL
C
C DEFINE-DIREKTIVEN     : _VA_TOONDER_
C
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
C
      REAL        X(II),         Y(JJ),         Z(KK),
     $           DX(II),        DY(JJ),        DZ(KK),
     $          DDX(II),       DDY(JJ),       DDZ(KK),
     $          RDX(II),       RDY(JJ),       RDZ(KK),
     $         RDDX(II),      RDDY(JJ),      RDDZ(KK),
     $            U(KK,JJ,II),   V(KK,JJ,II),   W(KK,JJ,II),
     $            G(KK,JJ,II),   B(KK,JJ,II)
C
      REAL 	TAU11(KK,JJ,II), TAU12(KK,JJ,II), TAU13(KK,JJ,II),
     $                           TAU22(KK,JJ,II), TAU23(KK,JJ,II),
     $                                            TAU33(KK,JJ,II)
C
      REAL UP(KK,JJ,II), VP(KK,JJ,II), WP(KK,JJ,II)
C
      REAL        DUDY(KK,JJ),  DUDZ(KK,JJ),
     $            DVDX(KK,JJ),  DVDZ(KK,JJ),
     $            DWDX(KK,JJ),  DWDY(KK,JJ)

C
      integer KK,JJ,II
C
      NBND = 2

      KM2    = KMX-NBND
      JM2    = JMX-NBND
      IM2    = IMX-NBND

C     
C
C
      RETURN
      END

