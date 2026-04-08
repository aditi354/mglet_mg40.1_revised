










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
      SUBROUTINE CONNECTMG 
     $                   (IDIM3D,IDIM2D,IDIM1D,
     $                    X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                    U,V,W,P,G,B,UFR,VRI,WBO,HILF,
     $                    UTO,VTO,WTO,PTO,GTO,
     $                    VFR,WFR,PFR,GFR,
     $                    ITYP,ILEVEL,IRB,JRB,KRB,NGRIDS,LOFGRIDS)
C*MGLET***************************************************************
C    C O N N E C T M G    KOMMUNIKATION ZWISCHEN VERSCHIEDENEN GITTERN
C                         FUER ALLE GITTER IN MULTIGRID-UMGEBUNG
C*MGLET***************************************************************
C
C VERS:  15.05.95 (MM)  : ORIGINAL AUS BOUNDMG ABGELEITET
C VERS:  09.12.98 (AM)  : aenderung wegen der front flaeche
C
C*MGLET***************************************************************
C
C
      CHARACTER (LEN=1)  ITYP
C
      REAL        X(IDIM1D),         Y(IDIM1D),         Z(IDIM1D),
     $           DX(IDIM1D),        DY(IDIM1D),        DZ(IDIM1D),
     $          DDX(IDIM1D),       DDY(IDIM1D),       DDZ(IDIM1D),
     $          U( IDIM3D ),       V( IDIM3D ),       W( IDIM3D ),
     $          P( IDIM3D ),       G( IDIM3D ),       B( IDIM3D ),
     $        UFR(IDIM2D,2),     VRI(IDIM2D,2),     WBO(IDIM2D,2),
     $        VFR(IDIM2D,2),     WFR(IDIM2D,2),
     $        PFR(IDIM2D,2),     GFR(IDIM2D,2),
     $        UTO(IDIM2D,2),     VTO(IDIM2D,2),     WTO(IDIM2D,2),
     $        PTO(IDIM2D,2),     GTO(IDIM2D,2),
     $       HILF( IDIM3D )
C
C
      INTEGER NGRIDS,LOFGRIDS(NGRIDS)
C
      IF (IRB .NE. 2 .AND .KRB .NE. 2 .AND. JRB .NE. 2) THEN
         JJRB = 0
         KKRB = 0
      CALL CONFROMG  
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,PFR,GFR,VRI,WBO,HILF,
     $              ITYP,ILEVEL,IRB,JJRB,KKRB,NGRIDS,LOFGRIDS)
      ENDIF


      IF (IRB .NE. 1 .AND .KRB .NE. 2 .AND. JRB .NE. 2) THEN
         JJRB = 0
         KKRB = 0
      CALL CONBACMG  
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,
     $              ITYP,ILEVEL,IRB,JJRB,KKRB,NGRIDS,LOFGRIDS)
      ENDIF


      IF (JRB .NE. 2) THEN
      CALL CONRGTMG  
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,
     $              ITYP,ILEVEL,IRB,JRB,KRB,NGRIDS,LOFGRIDS)
      ENDIF

      IF (JRB .NE. 1) THEN
      CALL CONLFTMG  
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,
     $              ITYP,ILEVEL,IRB,JRB,KRB,NGRIDS,LOFGRIDS)
      ENDIF


      IF (KRB .NE. 2) THEN
      CALL CONBOTMG  
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,
     $              ITYP,ILEVEL,IRB,JRB,KRB,NGRIDS,LOFGRIDS)
      ENDIF

      IF (KRB .NE. 1) THEN
      CALL CONTOPMG  
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,UTO,VTO,WTO,PTO,GTO,VRI,WBO,HILF,
     $              ITYP,ILEVEL,IRB,JRB,KRB,NGRIDS,LOFGRIDS)
      ENDIF
CTEST
CTEST      CALL BCUBMG  
CTEST     $             (IDIM3D,IDIM2D,IDIM1D,
CTEST     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
CTEST     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,
CTEST     $              ITYP,IGRID,IRB,JRB,KRB)
CTEST
CTESTC
CTEST      RETURN
      END
