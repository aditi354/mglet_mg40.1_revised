










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
         SUBROUTINE MGCONINF (DX,DY,DZ,X,Y,Z,IDIM1D,NBND)
C
C--MGLET----------------------------------------------------------------
C
C                  DRUCKT INFORMATION UEBER DIE GITTERVERNETZUNG
C                  
C
C
C        26. 2.94 (MM)  : ORIGINAL
C        21.03.03 (TB)  : SCALAR BOUNDARY TREATMENT IMPLEMENTED
C
C--MGLET----------------------------------------------------------------



      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


      INTEGER MAXLVLMAX,MINLVLMIN,MAXGRDSOFLVL

      PARAMETER ( MAXLVLMAX = 4,  MINLVLMIN = -4 )
      PARAMETER ( MAXGRDSOFLVL =260)

      COMMON /COLEVEL/ MAXLEVEL,MINLEVEL
      COMMON /COLEVEL/ NOFLEVEL,IGRDOFLEVEL,NOFPSDIR,IGRDOFPSDIR
      COMMON /COLEVEL/ NOFVPIT ,IGRDOFVPIT ,NOFTST,  IGRDOFTST
      COMMON /COLEVEL/ NOFSLCED,IGRDOFSLCED,NOFSLICE,IGRDOFSLICE
      COMMON /COLEVEL/ NOFSCAI, IGRDOFSCAI,NGRIDPERLEVEL,NOFTHISGRID

      INTEGER MAXLEVEL,MINLEVEL

      INTEGER NOFLEVEL   (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFLEVEL(MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFVPIT    (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFVPIT (MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFTST     (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFTST  (MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFSCAI    (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFSCAI (MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFSLCED   (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFSLCED(MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFSLICE   (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFSLICE(MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFPSDIR   (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFPSDIR(MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NGRIDPERLEVEL    
      INTEGER NOFTHISGRID(MAXGRDSOFLVL)
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
      REAL    
     $            X(IDIM1D),     Y(IDIM1D),     Z(IDIM1D),
     $           DX(IDIM1D),    DY(IDIM1D),    DZ(IDIM1D)

      CHARACTER (LEN=8) CGRID,CDIR,CBOC
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                                    ALLE GITTER JEDES LEVELS WERDEN
C                                    AUSGEDRUCKT
          WRITE (6,*) '************************************************'
     $,'*****************'
          WRITE (6,*) '***'
          WRITE (6,*) '***'
          WRITE (6,*) '***     INFORMATION UEBER DIE RANDBEDINGUNGEN'
          WRITE (6,*) '***     UND DIE GITTERVERNETZUNG'
          WRITE (6,*) '***'
          WRITE (6,*) '***'
          WRITE (6,*) '************************************************'
     $,'*****************'
C
      DO ILEVEL = MINLEVEL,MAXLEVEL

          WRITE (6,*) '************************************************'
     $,'*****************'
          WRITE (6,*) '***'
          WRITE (6,*) '***      LEVEL:',ILEVEL
          WRITE (6,*) '***'
          WRITE (6,*) '***'

        DO I = 1, NOFLEVEL(ILEVEL)
          IGRID = IGRDOFLEVEL(I,ILEVEL)

          WRITE (CGRID,'("GRD ",I3,":")') IGRID

          WRITE (6,*) '************************************************'
     $,'*****************'
          WRITE (6,*) CGRID
C
C                                    ALLE RICHTUNGEN WERDEN UNTERSUCHT
C

          DO IDIR = 1,6

            IF (IDIR .EQ. 1) THEN
               WRITE (CDIR,'(" FRONT  ")')
            ELSEIF (IDIR .EQ. 2) THEN
               WRITE (CDIR,'(" BACK   ")')
            ELSEIF (IDIR .EQ. 3) THEN
               WRITE (CDIR,'(" RIGHT  ")')
            ELSEIF (IDIR .EQ. 4) THEN
               WRITE (CDIR,'(" LEFT   ")')
            ELSEIF (IDIR .EQ. 5) THEN
               WRITE (CDIR,'(" BOTTOM ")')
            ELSEIF (IDIR .EQ. 6) THEN
               WRITE (CDIR,'(" TOP    ")')
            ENDIF



C                       = "PER",1 : PERIODISCH
C                       = "FIX",2 : FIXED-CONDITION
C                       = "OP1",3 : OPEN-WALL GRADIENT=0
C                       = "OP2",4 : OPEN-WALL ZWEITE ABLEITUNG=0
C                       = "NOS",5 : NOSLIP
C                       = "SLI",6 : SLIP
C                       = "CON",7 : CONNECT (GITTERKOPPLUNG)
C                       = "PAR",8 : PARENT  (FOR LOCALLY REFINED GRIDS)
C                       = "BLO",9 : BLOWING/SUCTION
C                       = "KON",10 : KONVECTIVE OUTLET
C                       = "FLU",11 : FLUCTUATIONS ARE TAKEN FROM DOWNSTREAM
C                       = "GRA",15 : FIXED SCALAR GRADIENT
C                       = "SCA",16 : FIXED SCALAR CONDITION

            DO IBOC = 1,NBOCONDS(IDIR,IGRID)
              IF (ITYPBOCONDS(IBOC,IDIR,IGRID) .EQ. 1) THEN
               WRITE (6,*) CGRID,CDIR,'PERIODIC'

              ELSEIF (ITYPBOCONDS(IBOC,IDIR,IGRID) .EQ. 2) THEN
               WRITE (6,*) CGRID,CDIR,'FIXED'

              ELSEIF (ITYPBOCONDS(IBOC,IDIR,IGRID) .EQ. 3) THEN
               WRITE (6,*) CGRID,CDIR,'OPEN-WALL 1'

              ELSEIF (ITYPBOCONDS(IBOC,IDIR,IGRID) .EQ. 4) THEN
               WRITE (6,*) CGRID,CDIR,'OPEN-WALL 2'

              ELSEIF (ITYPBOCONDS(IBOC,IDIR,IGRID) .EQ. 5) THEN
               WRITE (6,*) CGRID,CDIR,'NO-SLIP'

              ELSEIF (ITYPBOCONDS(IBOC,IDIR,IGRID) .EQ. 6) THEN
               WRITE (6,*) CGRID,CDIR,'SLIP'

              ELSEIF (ITYPBOCONDS(IBOC,IDIR,IGRID) .EQ. 7) THEN
               WRITE (6,*) CGRID,CDIR,'CONNECT WITH: ',
     $                     LBOGRIDS(IBOC,IDIR,IGRID)

              ELSEIF (ITYPBOCONDS(IBOC,IDIR,IGRID) .EQ. 8) THEN
               WRITE (6,*) CGRID,CDIR,'PARENT IS: ',
     $                     LBOGRIDS(IBOC,IDIR,IGRID)

              ELSEIF (ITYPBOCONDS(IBOC,IDIR,IGRID) .EQ. 9) THEN
               WRITE (6,*) CGRID,CDIR,'BLOWING/SUCTION'

              ELSEIF (ITYPBOCONDS(IBOC,IDIR,IGRID) .EQ. 10) THEN
               WRITE (6,*) CGRID,CDIR,'KONVECTIVE OUTLET'

              ELSEIF (ITYPBOCONDS(IBOC,IDIR,IGRID) .EQ. 11) THEN
               WRITE (6,*) CGRID,CDIR,'FLUCTUATIONS FROM DOWNSTREAM'
              ENDIF
C
C                                    
            ENDDO

          ENDDO

          WRITE (6,*) CGRID
          WRITE (6,*) '************************************************'
     $,'*****************'

        ENDDO
      ENDDO
C
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C

        DO I = 1, NOFLEVEL(ILEVEL)
          IGRID = IGRDOFLEVEL(I,ILEVEL)
          DO IDIR = 1,6
            DO IBOC = 1,NBOCONDS(IDIR,IGRID)
              IF (ITYPBOCONDS(IBOC,IDIR,IGRID) .EQ. 7) THEN
               IF(LBOGRIDS(IBOC,IDIR,IGRID) .LE. 0) THEN
                 CALL ERRR (700,'MGCONINF')
               ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C

      RETURN
      END


