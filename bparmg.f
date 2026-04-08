










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
      SUBROUTINE BPARMG 
     $                   (IDIM3D,IDIM2D,IDIM1D,
     $                    X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                    U,V,W,P,G,B,UFR,VRI,WBO,H2D1,H2D2,H2D3,
     $                    UTO,VTO,WTO,PTO,GTO,UBA,VBA,WBA,PBA,GBA,
     $                    VFR,WFR,PFR,GFR,
     $                    ITYP,ILEVEL,IRB,JRB,KRB,NGRIDS,LOFGRIDS
     $                   )     
C*MGLET***************************************************************
C    B P A R M G       KOMMUNIKATION ZWISCHEN Parent-GITTER
C                      UND CHILD-GITTER, BELEGT DIE RANDBED.-BUFFER
C                      DER LOKALEN GITTER (CHILDS)
C*MGLET***************************************************************
C
C VERS:  15.09.96 (MM)  : ORIGINAL AUS CONNECTMG ABGELEITET
C VERS:  09.12.98 (AM)  : PARALL. DER FRONT-FLAECHE
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
     $        UTO(IDIM2D,2),     VTO(IDIM2D,2),     WTO(IDIM2D,2),
     $        PTO(IDIM2D,2),     GTO(IDIM2D,2),
     $        VFR(IDIM2D,2),     WFR(IDIM2D,2),
     $        PFR(IDIM2D,2),     GFR(IDIM2D,2),
     $       H2D1( IDIM2D ),    H2D2( IDIM2D )
C
C
      INTEGER NGRIDS,LOFGRIDS(NGRIDS)
C
C      IF (IRB .NE. 2) THEN
      CALL PARFROMG  
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,U,V,W,P,G,B,
     $              UFR,VFR,WFR,PFR,GFR,H2D1,H2D2,H2D3,
     $              ITYP,ILEVEL,IRB,JRB,KRB,NGRIDS,LOFGRIDS
     $             )      
C      ENDIF


C      IF (IRB .NE. 1) THEN
      CALL PARBACMG  
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,U,V,W,P,G,B,
     $              UBA,VBA,WBA,PBA,GBA,H2D1,H2D2,H2D3,
     $              ITYP,ILEVEL,IRB,JRB,KRB,NGRIDS,LOFGRIDS
     $             )      
C      ENDIF


C      IF (JRB .NE. 2) THEN
      CALL PARRGTMG  
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,H2D1,
     $              ITYP,ILEVEL,IRB,JRB,KRB,NGRIDS,LOFGRIDS
     $             )      
C      ENDIF

C      IF (JRB .NE. 1) THEN
      CALL PARLFTMG  
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,H2D1,
     $              ITYP,ILEVEL,IRB,JRB,KRB,NGRIDS,LOFGRIDS
     $             )      
C      ENDIF


C      IF (KRB .NE. 2) THEN
      CALL PARBOTMG  
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,H2D1,
     $              ITYP,ILEVEL,IRB,JRB,KRB,NGRIDS,LOFGRIDS
     $             )      
C      ENDIF

C      IF (KRB .NE. 1) THEN
      CALL PARTOPMG  
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,U,V,W,P,G,B,
     $              UTO,VTO,WTO,PTO,GTO,H2D1,H2D2,H2D3,
     $              ITYP,ILEVEL,IRB,JRB,KRB,NGRIDS,LOFGRIDS
     $             )      
C      ENDIF


      END
