










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
      SUBROUTINE SETCOBONE (IGRID,IDIR,IBOCOND,CTYP,
     &                      IA,IE,JA,JE,KA,KE,INBR)
C*MGLET*****************************************************************
C    S E T C O B O N E     INITIALISIEREN DER COMMON-BLOCKES
C                          COBOUND
C                          HIER: INITIALISIEREN EINER RANDBEDINGUNG
C*MGLET****************************************************************
C
C PARAM:  
C       IGRID     :   GITTER, FUER DAS DIE RB. GESETZT WERDEN
C     CTYP(1:2),IDIR     : ZEIGER AUF DIE RICHTUNG:
C
C                         = "FR",1 : FRONT
C                         = "BA",2 : BACK 
C                         = "RI",3 : RIGHT
C                         = "LE",4 : LEFT
C                         = "BO",5 : BOTTOM
C                         = "TO",6 : TOP
C                         = "CU",7 : CUBUS
C
C         IBOCOND : + NUMMER DER SPEZIFIZIERTEN RANDBEDINGUNGEN
C                     PRO SEITENFLAECHE JEDES GITTERS
C
C   CTYP(3:5),ITYP: TYP DER RANDBEDINGUNG:
C
C                  = "PER",1 : PERIODISCH
C                  = "FIX",2 : FIXED-CONDITION
C                  = "OP1",3 : OPEN-WALL GRADIENT=0
C                  = "OP2",4 : OPEN-WALL ZWEITE ABLEITUNG=0
C                  = "NOS",5 : NOSLIP
C                  = "SLI",6 : SLIP
C                  = "CON",7 : CONNECT (GITTERKOPPLUNG)
C                  = "PAR",8 : PARENT  (FOR LOCALLY REFINED GRIDS)
C                  = "BLO",9 : BLOWING/SUCTION
C                  = "KON",10 : KONVECTIVE OUTLET
C                  = "FLU",11 : FLUCTUATIONS ARE TAKEN FROM DOWNSTREAM
C                  = "REC",12 : "Recycling", Velocity profile from	
C                               downstream is rescaled according to 
C                               Spalart & Lund
C                  = "OP3",13: Extrapolation dritter Ordnung
C                  = "GRA",15: FXGRADSCA, FIXED-GRADIENT-CONDITION
C                        (-> e.g. const. heat flux)
C                  = "SCA",16:: FXSCA, FIXED-SCALAR-CONDITION
C                        (-> e.g. const. wall temperature flux)
C
C   IA,IE,JA,JE:      SCHLEIFENGRENZEN IN RANDBEDINGUNGSROUTINEN:
C                     "FR","BA":JMX,KMX
C                     "RI","LE":IMX,KMX
C                     "BO","TO":JMX,IMX
C
C      LARBOCONDS : + BEGRENZUNGSINDIZES DER FLAECHE, AUF DER
C                     JEWEILIGE RANDBEDINGUNG GILT
C             INBR:   NACHBARGITTER
C        LBOGRIDS : + GITTER, AUS DEM RANDBEDINGUNG KOMMT
C     IPOS,JPOS,KPOS: POSITIONEN IM GITTER, ASU DEM RANDB. KOMMT
C     LPOSBOGRIDS : + POSITION IM GITTER, AUS DEM RANDBEDINGUNG KOMMT
C
C       IBOUNDGRID:   GITTER, VON DEM DER TYP DER RANDBEDINGUNGEN
C                     UEBERNOMMEN WIRD. MUSS BEI REGULAEREM GITTER
C                     GLEICH IGRID SEIN
C
C
C VERS:  19.10.93 (MM)  : ORIGINAL (AUS SETCOBOUND ABGELEITET
C VERS:  24.09.97 (AM)  : NEUE RANDBEDINGUNG (EXTRAPOLATION HOEHERER ORDNUNG)
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
 

      
      COMMON /COGRDCON/
     &                 IVPCHILD,
     &                 IPARENT,ISLPAR,
     &                 IPOSITION, JPOSITION, KPOSITION,
     &                 NOFSLCHILDS,IGRDOFSLCHILD,
     &                 ISLPOS, JSLPOS, KSLPOS,
     &                 IFRNBR, IBANBR, IRINBR, ILENBR,
     &                 IBONBR, ITONBR


      INTEGER
     &       IVPCHILD(MAXGRIDS),
     &       IPARENT(MAXGRIDS),ISLPAR(MAXGRIDS),
     & IPOSITION(MAXGRIDS), JPOSITION(MAXGRIDS), KPOSITION(MAXGRIDS),
     & NOFSLCHILDS(MAXGRIDS),IGRDOFSLCHILD(MAXGRIDS,MAXGRIDS),
     &    ISLPOS(MAXGRIDS), JSLPOS(MAXGRIDS), KSLPOS(MAXGRIDS),
     &       IFRNBR(MAXBOCONDS,MAXGRIDS), IBANBR(MAXBOCONDS,MAXGRIDS),
     &       IRINBR(MAXBOCONDS,MAXGRIDS), ILENBR(MAXBOCONDS,MAXGRIDS),
     &       IBONBR(MAXBOCONDS,MAXGRIDS), ITONBR(MAXBOCONDS,MAXGRIDS)
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
C
CCCCC      REAL X(IDIM1D), Y(IDIM1D), Z(IDIM1D)
C
      CHARACTER (LEN=8) CTYP
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  FUER ALLE RANDB. GLEICH
C

      NBOCONDS(            IDIR,IGRID) = IBOCOND


      LARBOCONDS(  1,IBOCOND,IDIR,IGRID) = IA
      LARBOCONDS(  2,IBOCOND,IDIR,IGRID) = IE
      LARBOCONDS(  3,IBOCOND,IDIR,IGRID) = JA
      LARBOCONDS(  4,IBOCOND,IDIR,IGRID) = JE
      LARBOCONDS(  5,IBOCOND,IDIR,IGRID) = KA
      LARBOCONDS(  6,IBOCOND,IDIR,IGRID) = KE
      
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  PERIODISCH
C              
      IF (CTYP(3:5) .EQ. 'PER') THEN

         ITYPBOCONDS(IBOCOND,IDIR,IGRID) = 1

      ELSEIF (CTYP(3:5) .EQ. 'FIX') THEN

         ITYPBOCONDS(IBOCOND,IDIR,IGRID) = 2

      ELSEIF (CTYP(3:5) .EQ. 'OP1') THEN

         ITYPBOCONDS(IBOCOND,IDIR,IGRID) = 3

      ELSEIF (CTYP(3:5) .EQ. 'OP2') THEN

         ITYPBOCONDS(IBOCOND,IDIR,IGRID) = 4

      ELSEIF (CTYP(3:5) .EQ. 'NOS') THEN

         ITYPBOCONDS(IBOCOND,IDIR,IGRID) = 5

      ELSEIF (CTYP(3:5) .EQ. 'SLI') THEN

         ITYPBOCONDS(IBOCOND,IDIR,IGRID) = 6

      ELSEIF (CTYP(3:5) .EQ. 'CON') THEN

         ITYPBOCONDS   (IBOCOND,IDIR,IGRID) = 7
         LBOGRIDS      (IBOCOND,IDIR,IGRID) = INBR

         CALL SETPOSBOGRIDS (LPOSBOGRIDS(1,IBOCOND,IDIR,IGRID),
     $                       IBOCOND,IDIR,IGRID,INBR,'CON')

      ELSEIF (CTYP(3:5) .EQ. 'PAR') THEN

         IF (IPARENT(IGRID) .GE. 1) THEN
         ITYPBOCONDS   (IBOCOND,IDIR,IGRID) = 8
         LBOGRIDS      (IBOCOND,IDIR,IGRID) = IPARENT(IGRID)

C         write (6,*)'calling setposbogrids', igrid,iparent(igrid)
         CALL SETPOSBOGRIDS (LPOSBOGRIDS(1,IBOCOND,IDIR,IGRID),
     $                       IBOCOND,IDIR,IGRID,IPARENT(IGRID),'PAR')
         ENDIF

        ELSEIF (CTYP(3:5) .EQ. 'BLO') THEN
C         BLOWING/SUCTION INDICES ARE DEFINED LATER IN MGBLOINDSET (AO)
C         GRID FIELD HAS TO BE BUILD FIRST !
C         DEFAULT FOR BLOWING/SUCTION = 0
        ISTART  = 0
        ISTOP   = 0
        JSTART  = 0
        JSTOP   = 0
        KSTART  = 0
        KSTOP   = 0
         ITYPBOCONDS   (IBOCOND,IDIR,IGRID) = 9
      ELSEIF (CTYP(3:5) .EQ. 'KON') THEN

         ITYPBOCONDS(IBOCOND,IDIR,IGRID) = 10

      ELSEIF (CTYP(3:5) .EQ. 'FLU') THEN

         ITYPBOCONDS(IBOCOND,IDIR,IGRID) = 11

      ELSEIF (CTYP(3:5) .EQ. 'REC') THEN

         ITYPBOCONDS(IBOCOND,IDIR,IGRID) = 12

      ELSEIF (CTYP(3:5) .EQ. 'OP3') THEN

         ITYPBOCONDS(IBOCOND,IDIR,IGRID) = 13
      ELSE
          WRITE(6,*)'SETCOBONE: ',CTYP(3:5) ,IBOCOND
          CALL ERRR (801,' SETCOBONE   ')

      ENDIF

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
C
      RETURN
      END

      SUBROUTINE SETPOSBOGRIDS(LPOSBOGRIDS,IBOCOND,IDIR,IGRID,INBR,CID)

C*MGLET*****************************************************
C
C                   SETZT DIE POSITIONEN, IN DEN NACHBARGITTER
C                   FUER DIR RANDBEDINGUNG 'CONNECT'
C
C
C     20.02.94 (MM):   ORIGINAL
C                      ACHTUNG!!!   FALLS TEILUEBERDECKTE GITTER HINZU
C*MGLET*****************************************************


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
 

      
      COMMON /COGRDCON/
     &                 IVPCHILD,
     &                 IPARENT,ISLPAR,
     &                 IPOSITION, JPOSITION, KPOSITION,
     &                 NOFSLCHILDS,IGRDOFSLCHILD,
     &                 ISLPOS, JSLPOS, KSLPOS,
     &                 IFRNBR, IBANBR, IRINBR, ILENBR,
     &                 IBONBR, ITONBR


      INTEGER
     &       IVPCHILD(MAXGRIDS),
     &       IPARENT(MAXGRIDS),ISLPAR(MAXGRIDS),
     & IPOSITION(MAXGRIDS), JPOSITION(MAXGRIDS), KPOSITION(MAXGRIDS),
     & NOFSLCHILDS(MAXGRIDS),IGRDOFSLCHILD(MAXGRIDS,MAXGRIDS),
     &    ISLPOS(MAXGRIDS), JSLPOS(MAXGRIDS), KSLPOS(MAXGRIDS),
     &       IFRNBR(MAXBOCONDS,MAXGRIDS), IBANBR(MAXBOCONDS,MAXGRIDS),
     &       IRINBR(MAXBOCONDS,MAXGRIDS), ILENBR(MAXBOCONDS,MAXGRIDS),
     &       IBONBR(MAXBOCONDS,MAXGRIDS), ITONBR(MAXBOCONDS,MAXGRIDS)

      INTEGER LPOSBOGRIDS(3)
      CHARACTER (LEN=3) CID

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

      IF (INBR.EQ.-1) RETURN
      IF (INBR.LE.0) CALL ERRR(502,'SETPOSBOGRIDS')

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                             POSITIONS FOR THE CONNECTED GRIDS ON
C                             THE SAME LEVEL
      IF (CID(1:3) .EQ. 'CON') THEN

      IF (IDIR.EQ.1) THEN

         LPOSBOGRIDS(1) = IMX(INBR)-1
         LPOSBOGRIDS(2) = 3
         LPOSBOGRIDS(3) = 3

      ELSEIF (IDIR.EQ.2) THEN

         LPOSBOGRIDS(1) = 2
         LPOSBOGRIDS(2) = 3
         LPOSBOGRIDS(3) = 3

      ELSEIF (IDIR.EQ.3) THEN

         LPOSBOGRIDS(1) = 3
         LPOSBOGRIDS(2) = JMX(INBR)-1
         LPOSBOGRIDS(3) = 3

      ELSEIF (IDIR.EQ.4) THEN

         LPOSBOGRIDS(1) = 3
         LPOSBOGRIDS(2) = 2
         LPOSBOGRIDS(3) = 3

      ELSEIF (IDIR.EQ.5) THEN

         LPOSBOGRIDS(1) = 3
         LPOSBOGRIDS(2) = 3
         LPOSBOGRIDS(3) = KMX(INBR)-1

      ELSEIF (IDIR.EQ.6) THEN

         LPOSBOGRIDS(1) = 3
         LPOSBOGRIDS(2) = 3
         LPOSBOGRIDS(3) = 2

      ELSE

         CALL ERRR (501,'SETPOSBOGRIDS')

         ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                             POSITIONS FOR THE LOCALLY REFINED GRIDS
C                             PARENT (=INBR) IS ON DIFFERENT LEVEL
      ELSE
C
C                             CHILD-GRIDS HAVE TO HAVE AN EVEN NUMBER
C                             OF GRID-POINTS
C
      IF ( IMX(IGRID)/2 .NE. (IMX(IGRID)+1)/2 ) 
     $                                 CALL ERRR (510,'SETPOSBOGRIDS')
      IF ( JMX(IGRID)/2 .NE. (JMX(IGRID)+1)/2 ) 
     $                                 CALL ERRR (511,'SETPOSBOGRIDS')
      IF ( KMX(IGRID)/2 .NE. (KMX(IGRID)+1)/2 ) 
     $                                 CALL ERRR (512,'SETPOSBOGRIDS')

      IF (IDIR.EQ.1) THEN

         LPOSBOGRIDS(1) = IPOSITION(IGRID)
         LPOSBOGRIDS(2) = JPOSITION(IGRID)
         LPOSBOGRIDS(3) = KPOSITION(IGRID)

      ELSEIF (IDIR.EQ.2) THEN

CTEST         LPOSBOGRIDS(1) = IPOSITION(IGRID) + (IMX(IGRID)-4)/2
         LPOSBOGRIDS(1) = IPOSITION(IGRID)
         LPOSBOGRIDS(2) = JPOSITION(IGRID)
         LPOSBOGRIDS(3) = KPOSITION(IGRID)

      ELSEIF (IDIR.EQ.3) THEN

         LPOSBOGRIDS(1) = IPOSITION(IGRID)
         LPOSBOGRIDS(2) = JPOSITION(IGRID)
         LPOSBOGRIDS(3) = KPOSITION(IGRID)

      ELSEIF (IDIR.EQ.4) THEN

         LPOSBOGRIDS(1) = IPOSITION(IGRID)
CTEST         LPOSBOGRIDS(2) = JPOSITION(IGRID) + (JMX(IGRID)-4)/2
         LPOSBOGRIDS(2) = JPOSITION(IGRID)
         LPOSBOGRIDS(3) = KPOSITION(IGRID)

      ELSEIF (IDIR.EQ.5) THEN

         LPOSBOGRIDS(1) = IPOSITION(IGRID)
         LPOSBOGRIDS(2) = JPOSITION(IGRID)
         LPOSBOGRIDS(3) = KPOSITION(IGRID)

      ELSEIF (IDIR.EQ.6) THEN

         LPOSBOGRIDS(1) = IPOSITION(IGRID)
         LPOSBOGRIDS(2) = JPOSITION(IGRID)
CTEST         LPOSBOGRIDS(3) = KPOSITION(IGRID) + (KMX(IGRID)-4)/2
         LPOSBOGRIDS(3) = KPOSITION(IGRID)

      ELSE

         CALL ERRR (501,'SETPOSBOGRIDS')

      ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      ENDIF

      RETURN
      END
