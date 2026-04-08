










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
      SUBROUTINE BOUNDMG 
     $                   (IDIM3D,IDIM2D,IDIM1D,
     $                    X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,    
     $                    U,V,W,P,G,B,UFR,VFR,WFR,VRI,HILF,
     $                    UTO,VTO,WTO,PTO,GTO,
     $                    ITYP,IGRID,IRB,JRB,KRB,
     $                    TIMEPHYS,UBA,VBA,WBA,PBA,GBA,UBO,VBO,WBO,
     $                    PFR,GFR
     $                   )     
C*STARLET***************************************************************
C        B O U N D M G    SETZEN DER RANDBEDINGUNGEN VOR DER
C                         DRUCKKORREKTUR UND VOR DEM ZEITSCHRITT.
C                         (LARGE-EDDY-SIMULATION)
C                         FUER EIN GITTER IN MULTIGRID-UMGEBUNG
C*STARLET***************************************************************
C
C VERS:  10. 3.93 (MM)  : ORIGINAL AUS BOUND ABGELEITET
C                         AENDERUNGEN: KOPF, DIMENSIONIERUNGEN
C VERS:  10.09.95 (AO)  : BLOWING/SUCTION AND CONVECTIVE BOUNDARY CONDITION
C                         INTRODUCED 
C VERS:  21.04.98 (AO)  : BBOTMG WIRD NUN NACH BCUB AUFGERUFEN, DA
C                         MANIPULATION AN KOERPERWAND MOEGLICH SEIN SOLL.
C                         ZANF UND ZEND WERDEN IN MGBLOINDEX BELEGT AUCH
C                         BEI BOTTOM-RB.
C 
C VERS:  12.12.98 (AM)  : AENDERUNG WEGEN DER FRONT-FLAECHE
C        28.02.03 (TB)  : SCALAR BOUNDARY TREATMENT IMPLEMENTED
C
C*STARLET***************************************************************
C
C
      CHARACTER (LEN=1)  ITYP
C
      REAL        X(IDIM1D),         Y(IDIM1D),         Z(IDIM1D),
     $           DX(IDIM1D),        DY(IDIM1D),        DZ(IDIM1D),
     $          DDX(IDIM1D),       DDY(IDIM1D),       DDZ(IDIM1D),
     $          U( IDIM3D ),       V( IDIM3D ),       W( IDIM3D ),
     $          P( IDIM3D ),       G( IDIM3D ),       B( IDIM3D ),
     $        UFR(IDIM2D,2),     VFR(IDIM2D,2),     WFR(IDIM2D,2),
     $        PFR(IDIM2D*2),     GFR(IDIM2D*2),
     $        VRI(IDIM2D,2),    
     $        UBO(IDIM2D,2),     VBO(IDIM2D,2),     WBO(IDIM2D*2),
     $        UBA(IDIM2D*2),     VBA(IDIM2D*2),     WBA(IDIM2D*2),
     $        PBA(IDIM2D*2),     GBA(IDIM2D*2),
     $        UTO(IDIM2D,2),     VTO(IDIM2D,2),     WTO(IDIM2D,2),
     $        PTO(IDIM2D,2),     GTO(IDIM2D,2),
     $       HILF( IDIM3D ), TIMEPHYS
C
C
      IF((ITYP.NE.'X').AND.(ITYP.NE.'Y').AND.(ITYP.NE.'Z')
     $    .AND.(ITYP.NE.'A')) THEN
      IF (IRB .NE. 2) THEN
       CALL BFROMG  
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,     
     $              U,V,W,P,G,B,UFR,VFR,WFR,HILF,
     $              ITYP,IGRID,IRB,JRB,KRB,TIMEPHYS,
     $              PFR,GFR
     $             )
      ENDIF

      IF (IRB .NE. 1) THEN
      CALL BBACMG  
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VFR,WFR,HILF,
     $              ITYP,IGRID,IRB,JRB,KRB,
     $              UBA,VBA,WBA,PBA,GBA
     $             )
      ENDIF

      IF (JRB .NE. 2) THEN
      CALL BRGTMG  
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WFR,HILF,
     $              ITYP,IGRID,IRB,JRB,KRB
     $             )     
      ENDIF

      IF (JRB .NE. 1) THEN
      CALL BLFTMG  
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WFR,HILF,
     $              ITYP,IGRID,IRB,JRB,KRB
     $             )     
      ENDIF

      IF (KRB .NE. 1) THEN
      CALL BTOPMG  
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UTO,VTO,WTO,PTO,GTO,HILF,
     $              ITYP,IGRID,IRB,JRB,KRB
     $             ) 
      ENDIF
      ENDIF

      IF ((ITYP.EQ.'X').OR.(ITYP.EQ.'Y').OR.(ITYP.EQ.'Z')
     $    .OR.(ITYP.EQ.'A')) THEN
      CALL BCUBMG  
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,WBO,HILF,
     $              ITYP,IGRID,IRB,JRB,KRB
     $             )     
      END IF
      IF ((ITYP.NE.'X').AND.(ITYP.NE.'Y').AND.(ITYP.NE.'Z')
     $    .AND.(ITYP.NE.'A')) THEN

      IF (KRB .NE. 2) THEN
      CALL BBOTMG  
     $             (IDIM3D,IDIM2D,IDIM1D,
     $              X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $              U,V,W,P,G,B,UFR,VRI,HILF,
     $              ITYP,IGRID,IRB,JRB,KRB,
     $              TIMEPHYS,UBO,VBO,WBO
     $             )      
      ENDIF
      ENDIF

C
      RETURN
      END
