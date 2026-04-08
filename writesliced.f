










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
      SUBROUTINE WRITESLICED (KKA,JJA,IIA,PHI,HILF,IGRID,
     $     KANAL,MODUS,IDIM3D)
C-MGLET---------------------------------------------------------------
C
C-MGLET---------------------------------------------------------------
      CHARACTER (LEN=8)   MODUS
      REAL PHI(IDIM3D),HILF(IDIM3D)

C

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
      
      COMMON /COGRDPRO/ LEVEL,LCHILD,XMIN,YMIN,ZMIN,XTOT,YTOT,ZTOT
      COMMON /COGRDPRO/ XHOMOG,YHOMOG,ZHOMOG,NXGRAE,NYGRAE,NZGRAE
      COMMON /COGRDPRO/ GRADPX,UBULKX,LTST,LVP,LSCAI,LPLEVEL,LPOISSONDIR
      COMMON /COGRDPRO/ LSLICE,NXSLICE,NYSLICE,NZSLICE,NVPGRIDS
      COMMON /COGRDPRO/ CONV1SANF,CONV1SEND,TRANSLES1,TRANSLES2
      COMMON /COGRDPRO/ GRADPXOLD

      INTEGER LEVEL(MAXGRIDS),NVPGRIDS(MAXGRIDS)
      INTEGER NXGRAE(MAXGRIDS),NYGRAE(MAXGRIDS),NZGRAE(MAXGRIDS)
      INTEGER NXSLICE(MAXGRIDS),NYSLICE(MAXGRIDS),NZSLICE(MAXGRIDS)

      REAL XTOT(MAXGRIDS),YTOT(MAXGRIDS),ZTOT(MAXGRIDS)
      REAL XMIN(MAXGRIDS),YMIN(MAXGRIDS),ZMIN(MAXGRIDS)
      REAL GRADPX(MAXGRIDS),UBULKX(MAXGRIDS)
      REAL CONV1SANF(MAXGRIDS),CONV1SEND(MAXGRIDS)
      REAL TRANSLES1(MAXGRIDS),TRANSLES2(MAXGRIDS),GRADPXOLD(MAXGRIDS)

      LOGICAL LCHILD(MAXGRIDS),LSLICE(MAXGRIDS),LPOISSONDIR(MAXGRIDS)
      LOGICAL XHOMOG(MAXGRIDS),YHOMOG(MAXGRIDS),ZHOMOG(MAXGRIDS)
      LOGICAL LTST(MAXGRIDS),LVP(MAXGRIDS),LSCAI(MAXGRIDS)
      LOGICAL LPLEVEL(MAXGRIDS)
      
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

      INTEGER MAXPROCS,NUMPROCS,MYID,IDIM_MG_BUFS,IDIM_MG_BIG
      PARAMETER (MAXPROCS = 64)
      PARAMETER ( IDIM_MG_BUFS=  9000000 , IDIM_MG_BIG= 19500000 )
      COMMON /COIDTHR/ NUMPROCS,MYID,IDPROCOFGRD
      COMMON /COIDTHR/ NVPOFPROC,IVPOFPROC,NTSOFPROC,ITSOFPROC
      COMMON /COIDTHR/ NSIOFPROC,ISIOFPROC
      INTEGER IDPROCOFGRD(MAXGRIDS)
      INTEGER IVPOFPROC  (MAXGRDSOFLVL,0:MAXPROCS)
      INTEGER NVPOFPROC  (             0:MAXPROCS,MINLVLMIN:MAXLVLMAX)
      INTEGER ITSOFPROC  (MAXGRDSOFLVL,0:MAXPROCS)
      INTEGER NTSOFPROC  (             0:MAXPROCS,MINLVLMIN:MAXLVLMAX)
      INTEGER ISIOFPROC  (MAXGRDSOFLVL,0:MAXPROCS)
      INTEGER NSIOFPROC  (             0:MAXPROCS,MINLVLMIN:MAXLVLMAX)
      COMMON /MG_BUFS/ SENDBUF,RECVBUF,BIGBUF
      REAL SENDBUF(IDIM_MG_BUFS),RECVBUF(IDIM_MG_BUFS)
      REAL BIGBUF(IDIM_MG_BIG)
C------------------------------------------------------------ POINTER-  
         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $               NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)

C---------------------------------------------------------- CONSISTENCY?

               IF (ZHOMOG(IGRID)) THEN
                  CALL ERRR (501,' WRITESLICED ')
               ELSE
                  IF (KKA .NE. KK)  CALL ERRR (601,' WRITESLICED ') 
               ENDIF
               
               IF (YHOMOG(IGRID)) THEN
                  IF (JJA .NE.  1 .AND. JJA .NE. JJ)  
     $                 CALL ERRR (502,' WRITESLICED ')
               ELSE
                  IF (JJA .NE. JJ)  CALL ERRR (602,' WRITESLICED ') 
               ENDIF
               
               IF (XHOMOG(IGRID)) THEN
                  IF (IIA .NE. 1 .AND. IIA .NE. II) THEN
                  CALL ERRR (503,' WRITESLICED ')
                  ENDIF
               ELSE
                  IF (IIA .NE. II)  CALL ERRR (603,' WRITESLICED ') 
               ENDIF

C---------------------------------------------------------------------
C                           BLOCK-WISE READING OF FIELD
C                           BLOCK MUST BE SMALLER THAN ONE SUBGRID

      IF (LSLICE(IGRID)) THEN
CC           IMBLOCK = IIA / (NXSLICE(IGRID)*NYSLICE(IGRID)) + 1
         IMBLOCK = IDIM_MG_BIG / (KKA*JJA)
      ELSE
         IMBLOCK = IIA
      ENDIF
      NBLOCK  = IIA/IMBLOCK + 1
C      NBLOCK = 1
      ISTOP   = 0
      ILENGTH = 0
      IF(IIA .EQ. 1) NBLOCK = 1 

C-MGLET---------------------------------------------------------------

CCC         DO 700 I = 0,IIA-1

      DO 700 IBLOCK = 1,NBLOCK

         ISTART  = ISTOP+1
         ILENGTH = MIN ( IMBLOCK , (IIA-ISTART+1))
         ISTOP   = ISTART + ILENGTH - 1
C-MGLET---------------------------- INITIALIZING BIGBUF --------------
      CALL SETS   (IDIM_MG_BIG,1,1,IDIM_MG_BIG,1,1,BIGBUF,0.0)                  

C-MGLET-------------------------- SAMPLING DATA FROM PROCESSORS ------
C


         DO ICHILD=1,MAX(1,NOFSLCHILDS(IGRID))
            ISUB = ICHILD
            IF (LSLICE(IGRID)) THEN
               ISUB = IGRDOFSLCHILD(ICHILD,IGRID)
            ELSE
               ISUB = IGRID
            ENDIF
               
C------------------------------------------------------------ POINTER  
         CALL MGDPB (KKS,JJS,IIS,IP3S,IP2S,IP1S,IBBS,IBUS,
     $               NFROS,NBACS,NRGTS,NLFTS,NBOTS,NTOPS,NCUBS,ISUB)


C--------------------------------------------- J-DIRECTION HOMOGENEOUS?
                     IF (JJA .EQ. 1 .OR. IIA .EQ. 1) THEN
                        FAVERAGE = 1.0/(FLOAT(NYSLICE(IGRID))
     $                                 *FLOAT(NXSLICE(IGRID)))
                        IF(JJA .EQ. 1) JJS = JJA
                        IF(IIA .EQ. 1) IIS = IIA
                        CALL MGPOINA (IP3S,I1L,I2L,ISUB)
                     ELSE
                        FAVERAGE = 1.0
                     ENDIF

C------------------------------------------------ POSITIONS      -----

               IF (LSLICE(IGRID)) THEN
                  IPOS = ISLPOS(ISUB)
                  JPOS = JSLPOS(ISUB)
                  KPOS = KSLPOS(ISUB)
               ELSE
                  IPOS = 3
                  JPOS = 3
                  KPOS = 3
               ENDIF

                     IF (JJA .EQ. 1) JPOS=3
                     IF (IIA .EQ. 1) IPOS=3

C----------------------------------------------------- OVERLAPPING ?
               IF ( ((ISTART .LE. IPOS-3+IIS) .AND.  
     $              (ISTOP  .GT. IPOS-1)) .OR.
     $               (IIA .EQ. 1)) THEN

                     IPROC = 0
                     MYID  = 0

C------------------------------------------ OVERLAPPING INDICES ------

                     KA=1
                     KE=KKA
                     JA=3
                     JE=JJS-2

C------------------------------------------ CORRECTION, IF AT END-----


                     IF (JPOS       .EQ.  3) JA=2
                     IF (JPOS-3+JJS .EQ. JJ) JE=JJS-1

C------------------------------------- J-DIRECTION HOMOGENEOUS?

                     IF (JJA .EQ. 1) THEN
                        JA=1
                        JE=1
                     ENDIF

C-------------------------------------------- I-DIRECTION ------------
                     
                     IGLOBAL_START = MAX(IPOS,ISTART)
                     IF (IPOS .EQ. 3) IGLOBAL_START = MAX(2,ISTART)
C         write (6,*) 'writesliced',istart,ipos,iglobal_start

                     IGLOBAL_STOP = MIN(IPOS-3+IIS-2,ISTOP)
                     IF (IPOS-3+IIS .EQ. II) 
     $                          IGLOBAL_STOP = MIN(ISTOP,II-1)

                     IA = IGLOBAL_START - IPOS + 3
                     IE = IGLOBAL_STOP  - IPOS + 3
                     
                     IF (IIA .EQ. 1) THEN
                         IA = 1
                         IE = 1
                     ENDIF
 
                     ISHIFT = IPOS - ISTART - 2
                     JSHIFT=JPOS-3
                     KSHIFT=KPOS-3

C---------------------------------------------------------------------

                     IF (IE .GE. IA ) THEN
C           write (6,*) myid,'calling addregion',iproc,ip3s

               IP3S = MAX (IP3S,1)

                     CALL ADDREGION(KKS,JJS,IIS,PHI(IP3S),
     $                    IPROC,KA,KE,JA,JE,IA,IE,
     $                    KKA,JJA,ILENGTH,BIGBUF,0,
     $                    KSHIFT,JSHIFT,ISHIFT,
     $                    HILF,MYID,FAVERAGE)

C           write (6,*) myid,'called addregion'


                     ENDIF
                  

               ENDIF
            ENDDO

C-MGLET------------------------- DATA SAMPLED ------------------------
C-MGLET------------------------- WRITING DATA BLOCK ------------------
               DO IGLOBAL=ISTART,ISTOP
               
                  I = IGLOBAL - ISTART 
C                  I = IGLOBAL-1

C-MGLET---------------------------------------------------------------
                  IF(MODUS .EQ. 'BINAER  ') THEN

            WRITE (KANAL) 
     $           ((BIGBUF( K + J*KKA + I*JJA*KKA),K=1,KKA),
     $           J=0,JJA-1)
                  ELSEIF(MODUS .EQ. 'CODIERT ') THEN

                     WRITE (KANAL,6040) 
     $                    ((BIGBUF( K + J*KKA + I*JJA*KKA),
     $                    K=1,KKA),J=0,JJA-1)
                  
                  ELSE
                     IFAIL = 501
                  ENDIF

               ENDDO
            


  700 CONTINUE

      RETURN
 6040 FORMAT (6(E12.5E3,1X))
      END
