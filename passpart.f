










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
         SUBROUTINE PASSPART (IDIM3D,IDIM2D,IDIM1D,IDIMA,NBND,
     $                      X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                      RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                      U,V,W,P,G,B,HILF1,HILF2,
     $                      H2D1,H2D2,H2D3,H2D4,H2D5,H2D6,
     $     AU,AV,AW,
     $     FELD1,FELD2,FELD3,
     $      BP,
     $                      TIMEPH,ITSTEP,ITTOT)


C**MGLET****************************************************************
C
C       P A S S P A R T      
C                     TRANSPORTIERT POSITIONEN UND ORIENTIERUNGSVEKTOREN
C                     DER PARTIKEL
C                     VERTEILT PARTIKEL AUF PROZESSOREN UND GITTER
C
C     1.9.99 (MM):  ORIGINAL
C   10.06.03 (FS):  PASSIVER PARTIKELTRANSPORT FUER T-MISCHER
C
C**MGLET****************************************************************
C

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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

C
C                            DARF NUR NACH COMGRID STEHEN !!!!
C
      COMMON /COMGVP/
     &               IPCORR, IPCGES,   NDIVLEPS,  LDIVLEPS, DIVGMX

      INTEGER
     &         IPCORR(MAXGRIDS),   IPCGES(MAXGRIDS),   
     &       NDIVLEPS(MAXGRIDS), LDIVLEPS(MAXGRIDS)

      REAL   DIVGMX(MAXGRIDS)
 

      PARAMETER ( NPART_MAX = 20000, NFAM_MAX = 20000 )

      COMMON /COPART2/  NPART, NFAM, RUN,WLOP,
     $                 IPART_OUT,NXNPART,PARREAD,
     $                 IPART, UPART, DUPART,
     $                 XPART,FXPART,FYPART,FZPART,
     $                 START_VEC,DS,PFLOWCASE
      LOGICAL PARREAD

      INTEGER IPART_OUT,WLOP,PFLOWCASE
      INTEGER NPART, NPART_MAX, NFAM, NFAM_MAX,
     $     NXNPART,
     $        IPART (    NPART_MAX ),
     $     IINDEX( NFAM_MAX ),IINDEX_N(NFAM_MAX),
     $     JINDEX( NFAM_MAX ),JINDEX_N(NFAM_MAX),
     $     KINDEX( NFAM_MAX ),KINDEX_N(NFAM_MAX),
     $     RUN   ( NFAM_MAX )

      REAL   SHAPE_PART
      REAL   DISSI(NFAM_MAX),
     $     DISTURB(NFAM_MAX),
     $     PARTSCA(NFAM_MAX),
     $     START_VEC ( NFAM_MAX*3 ),
     $     XPART (NFAM_MAX*3),
     $     UPART (NFAM_MAX*3),
     $     DUPART (NFAM_MAX*9),
     $     DS(NFAM_MAX)
      REAL    FXPART(NFAM_MAX,4),FYPART(NFAM_MAX,4),FZPART(NFAM_MAX,4)

C--------------------------------------------------------
C       X,Y,Z DER FAMILIES BEFINDEN SICH DAMIT AUF:
C     XPART(       1), XPART(NFAM + 1) UND XPART(2*NFAM + 1)
C
C
      REAL        RDDX(IDIM1D),     RDDY(IDIM1D),     RDDZ(IDIM1D),
     $             RDX(IDIM1D),      RDY(IDIM1D),      RDZ(IDIM1D),
     $             DDX(IDIM1D),      DDY(IDIM1D),      DDZ(IDIM1D),
     $              DX(IDIM1D),       DY(IDIM1D),       DZ(IDIM1D),
     $               X(IDIM1D),        Y(IDIM1D),        Z(IDIM1D)


      REAL        U (   IDIM3D  ), V (   IDIM3D  ), W (   IDIM3D  ),
     $            P (   IDIM3D  ), G (   IDIM3D  ), B (   IDIM3D  )

      INTEGER KKA,JJA,IIA
      REAL AU( IDIMA ),AV ( IDIMA ),AW ( IDIMA ),
     $     FELD1(IDIM3D),FELD2(IDIM3D),FELD3(IDIM3D)

      REAL BP( IDIM3D )

      REAL     H2D1( IDIM2D ), H2D2( IDIM2D ), H2D3( IDIM2D ),
     $         H2D4( IDIM2D ), H2D5( IDIM2D ), H2D6( IDIM2D )
C
C
      REAL     HILF1( IDIM3D ),HILF2( IDIM3D )

            tsearch = 0.0
            tgetvel = 0.0
            tgetder = 0.0
            ttstorient = 0.0
            ttstpos = 0.0
            tbpart = 0.0
            tges   = 0.0
C
C------------------------------------------  GERECHNET WIRD ERST NUR
C                                            AUF DEM MAXIMALEN LEVEL
C                                            (FEINSTES GITTER)

      ILEVEL = MAXLEVEL

      DO I = 1,NOFTST(ILEVEL)
         IGRID = IGRDOFTST(I,ILEVEL)

         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $        NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
         CALL MGDIMA  (KKA,JJA,IIA,IGRID)
         CALL MGPOINA (IAV,I2L,I1L,IGRID)       
C------------------------------------------- SEARCHING INDEX

            t1 = GETSEC(0)
            CALL SEARCHINDEX(KK,Z(IP1),DZ(IP1),DDZ(IP1),
     $           NZGRAE(IGRID),NFAM,XPART(2*NFAM_MAX + 1),KINDEX,
     $           HILF1(IP3))

            CALL SEARCHINDEX(JJ,Y(IP1),DY(IP1),DDY(IP1),
     $           NYGRAE(IGRID),NFAM,XPART(  NFAM_MAX + 1),JINDEX,
     $           HILF1(IP3))

            CALL SEARCHINDEX(II,X(IP1),DX(IP1),DDX(IP1),
     $           NXGRAE(IGRID),NFAM,XPART(         1),IINDEX,
     $           HILF1(IP3))

            t2 = GETSEC(0)

C------------------------------------------- GETTING VELOCITIES

            CALL GETVELOCITIES (KK,JJ,II,  X(IP1),   Y(IP1),   Z(IP1),
     $                                    DX(IP1),  DY(IP1),  DZ(IP1),
     $                                  RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                          U(IP3),V(IP3),W(IP3),
     $                          HILF1(IP3),HILF2(IP3),
     $                          H2D1(IP1),H2D2(IP1),H2D3(IP1),
     $                          NFAM,XPART(1),
     $                          XPART(NFAM_MAX+1),XPART(2*NFAM_MAX+1),
     $                          IINDEX,JINDEX,KINDEX,UPART,
     $                          FXPART,FYPART,FZPART)

            t3 = GETSEC(0)
C------------------------------------------ GETTING VELOCITY DERIVATIVES

            CALL GETDERIVATIVES(KK,JJ,II,  X(IP1),   Y(IP1),   Z(IP1),
     $                                    DX(IP1),  DY(IP1),  DZ(IP1),
     $                                   RDX(IP1), RDY(IP1), RDZ(IP1),
     $                                  RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                          U(IP3),V(IP3),W(IP3),
     $                          HILF1(IP3),HILF2(IP3),
     $                          H2D1(IP1),H2D2(IP1),H2D3(IP1),
     $                          NFAM,XPART(1),
     $                          XPART(NFAM_MAX+1),XPART(2*NFAM_MAX+1),
     $                          IINDEX,JINDEX,KINDEX,DUPART,
     $                          FXPART,FYPART,FZPART)

            t4 = GETSEC(0)

C------------------------------ COMPUTING RANDOM NUMBERS

C------------------------------------------- NEW ORIENTATION
C
CC------------------------------------------- PARTICLE STRESSES

CC------------------------------------------- SETTING STRESSES FOR GRID
CC------------------------------------------- OLD STRESSES ARE WEIGHTED BY
CC------------------------------------------- FOLD = 0.1
C

CC------------------------------------------ HILF1 CONTAINS NUMBER OF
CC------------------------------------------ PARTICLES IN CELL
C

C            t5c = GETSEC(0)
C 1000  continue
CC------------------------------------------- NEW POSITION


            CALL TSTPOS (NFAM,UPART,
     $                   XPART(1),XPART(NFAM_MAX+1),XPART(2*NFAM_MAX+1),
     $                   DT,DS)

            t6 = GETSEC(0)
C------------------------------------------- Boundary Conditions
            CALL SEARCHINDEX(KK,Z(IP1),DZ(IP1),DDZ(IP1),
     $           NZGRAE(IGRID),NFAM,XPART(2*NFAM_MAX + 1),KINDEX_N,
     $           HILF1(IP3))

            CALL SEARCHINDEX(JJ,Y(IP1),DY(IP1),DDY(IP1),
     $           NYGRAE(IGRID),NFAM,XPART(  NFAM_MAX + 1),JINDEX_N,
     $           HILF1(IP3))

            CALL SEARCHINDEX(II,X(IP1),DX(IP1),DDX(IP1),
     $           NXGRAE(IGRID),NFAM,XPART(         1),IINDEX_N,
     $           HILF1(IP3))

            CALL BPART_BODY (KK,JJ,II,  X(IP1),   Y(IP1),   Z(IP1),
     $                            DX(IP1),  DY(IP1),  DZ(IP1),
     $                   NFAM,
     $                   XPART(1),XPART(NFAM_MAX+1),XPART(2*NFAM_MAX+1),
     $                   NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,RUN(1),
     $       START_VEC(1),START_VEC(NFAM_MAX+1),START_VEC(2*NFAM_MAX+1)
     $       ,IINDEX,JINDEX,KINDEX,IINDEX_N,JINDEX_N,KINDEX_N,BP(IP3))
            CALL BPART (KK,JJ,II,  X(IP1),   Y(IP1),   Z(IP1),
     $                            DX(IP1),  DY(IP1),  DZ(IP1),
     $                   NFAM,
     $                   XPART(1),XPART(NFAM_MAX+1),XPART(2*NFAM_MAX+1),
     $                   NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,RUN(1),
     $       START_VEC(1),START_VEC(NFAM_MAX+1),START_VEC(2*NFAM_MAX+1))


C------------------------------------------- CALCULATION OF DISSIPATION

            CALL PARTDISS (KK,JJ,II,X,Y,Z,DX(IP1),DY(IP1),DZ(IP1),
     $           RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $           DDX(IP1),DDY(IP1),DDZ(IP1),
     $           U(IP3),V(IP3),W(IP3),B(IP3),HILF1(IP3),HILF2(IP3),
     $           H2D1(IP2),H2D2(IP2),H2D3(IP2),
     $           NFAM,XPART(1),XPART(NFAM_MAX+1),XPART(2*NFAM_MAX+1),
     $           IINDEX,JINDEX,KINDEX,DISSI,DISTURB,
     $           FXPART,FYPART,FZPART,
     $           AU(IAV),AV(IAV),AW(IAV),
     $           KKA,JJA,IIA,
     $           FELD1(IP3),FELD2(IP3),FELD3(IP3),
     $           BP(IP3),
     $           NFRO,NBAC,NRGT,NLFT,NBOT,NTOP)

C------------------------------------------- LOOP OVER CLUSTER FINISHED
            tsearch    = tsearch    + t2-t1
            tgetvel    = tgetvel    + t3-t2
            tgetder    = tgetder    + t4-t3
            ttstorient = ttstorient + t5-t4
            ttstpos    = ttstpos    + t6-t5
            tbpart     = tbpart     + t7-t6
            tges       = tges       + t7-t1
C------------------------------------------- CALCULATION OF CONCENTRATION AT PARTICLE POSITION
C------------------------------------------- OUTPUT?


         IF (MOD (ITTOT,IPART_OUT) .EQ. 0) THEN
            CALL WRITE_PASSPART (KK,JJ,II,23,RUN,NFAM,
     $           XPART,DISSI,DISTURB,KINDEX,JINDEX,IINDEX,NFAM_MAX,WLOP,
     $           DS,TIMEPH,PARTSCA)
            
         ENDIF
      ENDDO

C      IF (ITSTEP .EQ. 1) THEN
C
C            tsearch = t2-t1
C            tgetvel = t3-t2
C            tgetder = t4-t3
C            ttstorient = t5-t4
C            trr = t5a1-t5
C            trrrr = t5a2-t5a1
C            tstress = t5a-t5a2
C            tparticle_stress = t5b-t5a
C            tsmooth_stress = t5c-t5b
C            ttstpos = t6-t5c
C            tbpart = t7-t6
C
C        if (myid .eq. 0) then
C            write (6,*)'   search index:',tsearch
C            write (6,*)' get velocities:',tgetvel
C            write (6,*)'get derivatives:',tgetder
C            write (6,*)'   orientations:',ttstorient
C            write (6,*)'             rr:',trr
C            write (6,*)'           rrrr:',trrrr
C            write (6,*)'         stress:',tstress
C            write (6,*)'particle_stress:',tparticle_stress
C            write (6,*)'  smooth_stress:',tsmooth_stress
C            write (6,*)'  new positions:',ttstpos
C            write (6,*)' boundary cond.:',tbpart
C            write (6,*)'            all:',tges
C        endif
C
C        ENDIF

         RETURN
         END
