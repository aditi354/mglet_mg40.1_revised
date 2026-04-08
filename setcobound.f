










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
       SUBROUTINE SETCOBOUND (IGRID,IBOUNDGRID)


C*MGLET*****************************************************************
C    S E T C O B O U N D   INITIALISIEREN DER COMMON-BLOCKES
C                          COBOUND
C                          ZUM STEUERN DER RANDBEDINGUNGSROUTINEN
C*MGLET*****************************************************************
C
C PARAM: NBOCONDS : + ANZAHL DER SPEZIFIZIERTEN RANDBEDINGUNGEN
C                     PRO SEITENFLAECHE JEDES GITTERS
C      LARBOCONDS : + BEGRENZUNGSINDIZES DER FLAECHE, AUF DER
C                     JEWEILIGE RANDBEDINGUNG GILT
C     ITYPBOCONDS : + TYP DER JEWEILIGEN RANDBEDINGUNG
C        LBOGRIDS : + GITTER, AUS DEM RANDBEDINGUNG KOMMT
C     LPOSBOGRIDS : + POSITION IM GITTER, AUS DEM RANDBEDINGUNG KOMMT
C
C       IGRID     :   GITTER, FUER DAS DIE RB. GESETZT WERDEN
C       IBOUNDGRID:   GITTER, VON DEM DER TYP DER RANDBEDINGUNGEN
C                     UEBERNOMMEN WIRD. MUSS BEI REGULAEREM GITTER
C                     GLEICH IGRID SEIN
C
C
C        IDIR     : ZEIGER AUF DIE RICHTUNG:
C
C                         = 1 : FRONT
C                         = 2 : BACK 
C                         = 3 : RIGHT
C                         = 4 : LEFT
C                         = 5 : BOTTOM
C                         = 6 : TOP
C                         = 7 : CUBUS
C
C        ITYP     : TYP DER RANDBEDINGUNG:
C
C                       = 1 : PERIODIC
C                       = 2 : FIXED-CONDITION
C                       = 3 : OPEN-WALL GRADIENT=0
C                       = 4 : OPEN-WALL ZWEITE ABLEITUNG=0
C                       = 5 : NOSLIP
C                       = 6 : SLIP
C                       = 7 : CONNECT (GITTERKOPPLUNG)
C                       = 8 : PARENT  (FOR LOCALLY REFINED GRIDS)
C                       = 9 : BLOWING/SUCTION
C                       = 10 : KONVECTIVE OUTLET
C                       = 11 : FLUCTUATIONS ARE TAKEN FROM DOWNSTREAM
C                       = 12 : "Recycling", Velocity profile from	
C                              downstream is rescaled according to 
C                              Spalart & Lund
C                       = 13: OP3: OPEN-WALL 3.DERIVATIVE=0 
C                       = 15: FXGRADSCA: FIXED-GRADIENT-CONDITION
C                             (-> e.g. const. heat flux)
C                       = 16: FXSCA: FIXED-SCALAR-CONDITION
C                             (-> e.g. const. wall temperature flux)
C
C
C VERS:  10. 3.93 (MM)  : ORIGINAL
C        17. 6.93 (MM)  : IBOUNDGRID EINGEFUEHRT
C        19.10.93 (MM)  : SETCOBONE WIRD AUFGERUFEN
C        20. 2.94 (MM)  : UEBERGABE AN SETCOBONE GEAENDERT
C VERS:  05.12.01 (TB)  : NEW FIXED-GRADIENT-CONDITION (Heat Flux BC)
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
C
C                                VORBELEGUNG  GESCHIEHT IN
C                                ICOBOUND
C
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  GITTER

      IG = IGRID
      IB = IBOUNDGRID

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  FRONT 


      DO I = 1,NBOCD( 1 , IG )
      FRONT(I,IG) = FRONT(I,IB)
      CALL SETCOBONE (IG,  1 , I ,FRONT(I,IB),0,0,1,JMX(IG),1,KMX(IG),
     $                IFRNBR(I,IG))
      ENDDO

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  BACK

      DO I = 1,NBOCD( 2 , IG )
      BACK(I,IG) = BACK(I,IB)
      CALL SETCOBONE (IG,  2 , I ,BACK(I,IB),0,0,1,JMX(IG),1,KMX(IG),
     $                IBANBR(I,IG))
      ENDDO

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  RIGHT

      DO I = 1,NBOCD( 3 , IG )
      RIGHT(I,IG) = RIGHT(I,IB)
      CALL SETCOBONE (IG,  3 , I ,RIGHT(I,IB),1,IMX(IG),0,0,1,KMX(IG),
     $                IRINBR(I,IG))
      ENDDO

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  LEFT

      DO I = 1,NBOCD( 4 , IG )
      LEFT(I,IG) = LEFT(I,IB)
      CALL SETCOBONE (IG,  4 , I ,LEFT(I,IB),1,IMX(IG),0,0,1,KMX(IG),
     $                ILENBR(I,IG))
      ENDDO

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  BOTTOM

      DO I = 1,NBOCD( 5 , IG )
      BOTTOM(I,IG) = BOTTOM(I,IB)
      CALL SETCOBONE (IG, 5 ,I,BOTTOM(I,IB),1,IMX(IG),1,JMX(IG),0,0,
     $                IBONBR(I,IG))
      ENDDO

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  TOP

      DO I = 1,NBOCD( 6 , IG )
      TOP(I,IG) = TOP(I,IB)
      CALL SETCOBONE (IG,  6 , I ,TOP(I,IB),1,IMX(IG),1,JMX(IG),0,0,
     $                ITONBR(I,IG))
      ENDDO

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  KUBUS

      DO I = 1,NBOCD( 7 , IG )
      CUBE(I,IG) = CUBE(I,IB)
      CALL SETCOBONE (IG,  7 , I ,CUBE(I,IB),1,0 ,1,0 ,1,0,
     $                0)
      ENDDO

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
C
      RETURN
      END


