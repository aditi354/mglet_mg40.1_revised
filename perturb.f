










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
      SUBROUTINE PERTURB (KK,JJ,II,KMX,JMX,IMX,
     $                    DX,DY,DZ,
     $                    U,V,W,UOLD,VOLD,WOLD,X,Y,Z,VCON,WCON,
     $                    YBANF,YBEND,ZBANF,ZBEND,HI,GRADPX)
C
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


      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

      COMMON /CONLES/  CONV2S,CMUE,CAPPA,ECONST
      SAVE   /CONLES/
C
      PARAMETER    (m=1024)
      REAL         U(KK,JJ,II),   V(KK,JJ,II),   W(KK,JJ,II),
     $             UOLD(KK,JJ,II),VOLD(KK,JJ,II),WOLD(KK,JJ,II),
     $             Z(KK),
     $             DX(II),DY(JJ),DZ(KK),
     $             X(II),   Y(JJ)
      REAL         HI(KK*JJ*II)
C
C
      IM2 = IMX-2
      IM1 = IMX-1
      JM1 = JMX-1
      JM2 = JMX-2
      JM3 = JMX-3
      KM1 = KMX-1
      KM2 = KMX-2
      KM3 = KMX-3

      DO I=1,II
       DO J=1,JJ
        DO K=1,KK
         UOLD(K,J,I) = U(K,J,I)
         VOLD(K,J,I) = 0.0
         WOLD(K,J,I) = 0.0
        ENDDO
       ENDDO
      ENDDO
C.... Hinzugeben der Störung
      IF(CIDUFR(1:9).EQ.'ORRSOMMER') THEN
      CALL ORRSOMMER(KK,JJ,II,U,W,Z,X,DZ,HI(1),HI(1*m+1),HI(2*m+1),
     $               HI(3*m+1),HI(4*m+1),HI(5*m+1),HI(6*m+1),
     $               HI(7*m+1),HI(8*m+1),HI(9*m+1),HI(10*m+1),
     $               HI(11*m+1),HI(12*m+1),HI(13*m+1),HI(14*m+1),
     $               HI(15*m+1),1.0)

      ENDIF

      RETURN
      END
