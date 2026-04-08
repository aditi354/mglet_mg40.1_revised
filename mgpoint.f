










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
       SUBROUTINE MGPOINT (IP3,IP2,IP1,IBB,IB3,IBU,IGRID)

C*MGLET*****************************************************************
C      M G P O I N T         MULTIGRID-VERWALTUNG
C                            LIEFERT ANFANGSPOSITIONEN DES 
C                            GITTERS IGRID IN EINDIMENSIONALEM FELD
C*MGLET*****************************************************************
C
C PARAM: IP3      : + POINTER FUER 3D-FELDER
C        IP2      : + POINTER FUER 2D-FELDER
C        IP1      : + POINTER FUER 1D-FELDER
C        IPB      : + POINTER FUER BUFFER-FELDER
C        IB3      : + POINTER FUER BUFFER-FELDER 3*2D
C        IGRID    : - GITTER, FUER DAS DIE GROESSEN GESETZT WERDEN
C
C
C VERS:   8. 4.93 (MM)  : ORIGINAL
C         22.04.98 (AO) : IB3 FUER DYNAMISCHES MODELL BUFFER-FELD
C
C
C*MGLET*****************************************************************
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
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                           UEBERPRUEFUNGEN

      IF (IGRID .GT. MAXGRIDS) CALL ERRR (  555 ,'MGPOINT')

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
        IP3   =    IP3D(IGRID)
        IP2   =    IP2D(IGRID)
        IP1   =    IP1D(IGRID)
        IBB   =    IPBB(IGRID)
        IB3   =    IPB3(IGRID)
        IBU   =    IPBU(IGRID)
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      IF ( IP3 .EQ. 0 ) then
         write(6,*)'mgpoint, igrid, ip3, nof3d:',igrid, ip3, nof3d
         CALL ERRR (  556 ,'MGPOINT')
      ENDIF
      IF ( IP2 .EQ. 0 ) CALL ERRR (  557 ,'MGPOINT')
      IF ( IP1 .EQ. 0 ) CALL ERRR (  558 ,'MGPOINT')
      IF ( IBB .EQ. 0 ) CALL ERRR (  559 ,'MGPOINT')
      IF ( IB3 .EQ. 0 ) CALL ERRR (  571 ,'MGPOINT')
      IF ( IBU .EQ. 0 ) CALL ERRR (  560 ,'MGPOINT')
C
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      IF ( IP3 .GT. NOF3D ) then
          write(6,*)'mgpoint, igrid, ip3, nof3d:',igrid, ip3, nof3d
         CALL ERRR (  566 ,'MGPOINT')
      endif
      IF ( IP2 .GT. NOF2D ) CALL ERRR (  567 ,'MGPOINT')
      IF ( IP1 .GT. NOF1D ) CALL ERRR (  568 ,'MGPOINT')
      IF ( IBB .GT. NOFBB ) CALL ERRR (  569 ,'MGPOINT')
      IF ( IBU .GT. NOFBU ) CALL ERRR (  570 ,'MGPOINT')
      IF ( IB3 .GT. NOFB3 ) CALL ERRR (  572 ,'MGPOINT')
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C

      RETURN
      END
