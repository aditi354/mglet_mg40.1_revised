










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
       SUBROUTINE ICOPHYSPAR ()


C*STARLET***************************************************************
C        I C O P H Y S P A R     INITIALISIEREN DES COMMON-BLOCKES
C                          COPHYSPAR
C                          ZUM STEUERN DER PHYSIKALISCHEN PARAMETER
C*STARLET***************************************************************
C
C PARAM: MTURB
C        TU_LEVEL
C        RHO
C        GMOL
C        UGRID
C        VREF
C        EXPON
C        CIDUFR
C        UFRCON
C        UFRFREQ
C        DELTA
C        XREF
C        IPP, JPP, KPP
C        XPER, YPER, ZPER
C        NXPER, NYPER, NZPER
C        NPPHYS
C        PPHYS
C        XPPHYS
C        
C
C VERS:   3.7.97 (MM)  : ORIGINAL
C         12.12.02 (TB)  : PARAMETERS FOR SCALAR(HEAT) TRANSPORT ADDED
C
C
C*STARLET***************************************************************
C

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
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  VORBELEGUNG 
C
 	 MTURB         = 0
         TU_LEVEL      = 0.1
         RHO           = 1.0
         GMOL          = 0.1 
         UGRID         = 0.0 
         VREF          = 1.0
         EXPON         = 0.14286
         CIDUFR        = 'UNIFORM'
         UFRCON        = 1.0
         UFRFREQ       = 0.0
         DELTA         = 1.0
         XREF          = 1.0
         IPP           = 3
         JPP           = 3
         KPP           = 3 
         XPER          = 1.0
         YPER          = 1.0
         ZPER          = 1.0
         NXPER         = 0 
         NYPER         = 0
         NZPER         = 0
         NPPHYS        = 0
         DO I=1,NPPHYS_MAX

            PPHYS (I)   = 0.0
            XPPHYS(I)   = FLOAT(I)

        ENDDO
         LAMDA	       = 1.0
         PRMOL         = 0.7
         PRTURB        = 0.7 
         TREF          = 1.0
         EXPONT        = 0.14286
         CIDUFR        = 'UNIFORM'
         TFRCON        = 1.0
         DELTAT        = 1.0        
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
C
      RETURN
      END


