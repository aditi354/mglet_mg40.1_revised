










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
      SUBROUTINE DOGRID (CIDENT,IIDENT,RIDENT,
     $                   DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,
     $              U,V,W,P,G,B,HILF,	
     $              IDIM3D,IDIM2D,IDIM1D,NBUF,NBND,IDIMA,IDIM1L,IDIM2L,
     $              IGRID)
C**** MGLET ************************************************************
C           D O G R I D
C                           SCHREIBT EIN GITTER RAUS
C
C         7. 6.93 (MM) :   ORIGINAL  (AUS DIGRID ABGELEITET)
C        26. 5.95 (MM) :   KOPF UM HILFSFELD ERWEITERT, FUER MPI
C        04.02.03 (TB) :   SCALAR TRANSPORT VARIABLE T ADDED
C
C**** MGLET ************************************************************


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
 


      COMMON /COSTRLES/
     $                  VERS,   NRRUN,
     $                  DREAD,  DWRITE, DCONT,  DCUB,   
     $                  LBODYINSERT,
     $                  LGRIDNEW,    LGRIDREMOVE,
     $                  NPRNEU, FPRNEU, MTSTEP,
     $                  DT,     ITPRIN, IPINF,  ITINT,
     $                  IPRINT_WSS, IPRINT_WNS, ITFLUC, ITMIT,  
     $                  MPCORR, EPCORR, EPFAK,  MPCVOR, MPCNACH,
     $                  IVPINF, 
     $                  OMG,    LDIMLO, ISETRE,
     $                  LINPRN, IWRB,   LREC

      INTEGER 
     $                  NRRUN,
     $                  MTURB,  NPRNEU, MTSTEP,
     $                          ITPRIN, IPINF,  ITINT,
     $                  ITFLUC, ITMIT,  MPCORR,        
     $                                          ISETRE,
     $                          MSLIN,  IWRB          

      LOGICAL
     $                  DREAD,  DWRITE, DCONT,  DCUB,   
     $                  LBODYINSERT,
     $                  LGRIDNEW,       LGRIDREMOVE,
     $                  LDIMLO, LINPRN, LREC

      REAL
     $                  DT,     EPCORR,  OMG, FPRNEU

      CHARACTER (LEN=8)       VERS


      COMMON /CLINOU/

     $                 LDIB,   LDIC,   LDIEIB, LDIEIC,
     $                 LDOB,   LDOC,   LDOEIB, LDOEIC

      SAVE   /CLINOU/
      LOGICAL

     $                 LDIB,   LDIC,   LDIEIB, LDIEIC,
     $                 LDOB,   LDOC,   LDOEIB, LDOEIC

      CHARACTER (LEN=8)    CIDENT(10)
C
      INTEGER IIDENT(100)

      REAL    RIDENT(100),
     $            X(IDIM1D),     Y(IDIM1D),     Z(IDIM1D),
     $           DX(IDIM1D),    DY(IDIM1D),    DZ(IDIM1D),
     $          DDX(IDIM1D),   DDY(IDIM1D),   DDZ(IDIM1D),
     $        U( IDIM3D ), V( IDIM3D ), W( IDIM3D ),
     $        P( IDIM3D ), G( IDIM3D ), B( IDIM3D ), HILF( IDIM3D )	


C
C                                 BEREITSTELLUNG DER POINTER
C                                 BASIS-RANDBEDINGUNGEN
C                                 DIMENSIONIERUNGEN
C
       CALL MGDIMS  (KK,JJ,II,IGRID)
       CALL MGPOINT (IP3,IP2,IP1,IBB,IB3,IBU,IGRID)
       CALL MGBASB  (NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)




         IF(LDOB  ) THEN
         CALL DOB    (CIDENT,IIDENT,RIDENT,KK,JJ,II,KK,JJ,II,NBND,
     $              DDX(IP1),DDY(IP1),DDZ(IP1),
     $              DX(IP1),DY(IP1),DZ(IP1),X(IP1),Y(IP1),Z(IP1),
     $              5,U,V,W,
     $              P,G,B,B,B,HILF,IGRID,IDIM3D)
         END IF
         IF(LDOC  ) THEN
         CALL DOC    (CIDENT,IIDENT,RIDENT,KK,JJ,II,KK,JJ,II,NBND,
     $              DDX(IP1),DDY(IP1),DDZ(IP1),
     $              DX(IP1),DY(IP1),DZ(IP1),X(IP1),Y(IP1),Z(IP1),
     $              5,U,V,W,
     $              P,G,B,TK,TE,HILF,IGRID,IDIM3D)
         END IF

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C

      RETURN
      END


