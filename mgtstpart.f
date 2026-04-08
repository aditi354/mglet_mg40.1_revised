










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
         SUBROUTINE MGTSTPART (IDIM3D,IDIM2D,IDIM1D,NBND,
     $                      X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                      RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                      U,V,W,P,G,B,HILF1,HILF2,
     $                      H2D1,H2D2,H2D3,H2D4,H2D5,H2D6,
     $                      TIMEPH,ITSTEP,ITTOT,
     $                      TAU11,TAU12,TAU13,
     $                      TAU21,TAU22,TAU23,
     $                      TAU31,TAU32,TAU33)


C**MGLET****************************************************************
C
C   M G T S T P A R T      
C                     TRANSPORTIERT POSITIONEN UND ORIENTIERUNGSVEKTOREN
C                     DER PARTIKEL
C                     VERTEILT PARTIKEL AUF PROZESSOREN UND GITTER
C
C     1.9.99 (MM):  ORIGINAL
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
 

	PARAMETER ( NPART_MAX =  1, NFAM_MAX = 1 )
C      PARAMETER ( NPART_MAX =   1, NFAM_MAX = 1 )

      COMMON /COPART/  NPART, NFAM, SHAPE_PART,
     $                 ASPECT, VOL_FRAC, IRANDOM, FRANDOM, 
     $                 ITSKIP, 
     $                 IPART_OUT,
     $                 CONF, UPART, DUPART, XPART, RRANDOM,
     $                 RR, RRRREE, TAU,
     $                 FXPART,FYPART,FZPART

      INTEGER IPART_OUT, IRANDOM, ITSKIP
      INTEGER NPART, NPART_MAX, NFAM, NFAM_MAX,
     $     IINDEX( NFAM_MAX ),
     $     JINDEX( NFAM_MAX ),
     $     KINDEX( NFAM_MAX )

      REAL   SHAPE_PART, ASPECT, VOL_FRAC, FRANDOM
      REAL    CONF  ( NFAM_MAX * NPART_MAX * 3 )
      REAL    UPART ( NFAM_MAX * 3 ),
     $       DUPART ( NFAM_MAX * 9 ),
     $        XPART ( NFAM_MAX * 3 ),
     $      RRANDOM ( NPART_MAX* 3 ),
     $           RR ( NFAM_MAX*9 ),
     $       RRRREE ( NFAM_MAX*9 ),
     $          TAU ( NFAM_MAX*9 )

      REAL FXPART(NFAM_MAX,4),FYPART(NFAM_MAX,4),FZPART(NFAM_MAX,4)

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

      REAL     H2D1( IDIM2D ), H2D2( IDIM2D ), H2D3( IDIM2D ),
     $         H2D4( IDIM2D ), H2D5( IDIM2D ), H2D6( IDIM2D )
C
C
      REAL     HILF1( IDIM3D ),HILF2( IDIM3D )

      LOGICAL LOUT

c        write (6,*) 'MGTSTPART',nfam,npart,ittot,itstep
C        write (6,*) 'MGTSTPART',conf(1)
C        write (6,*) 'MGTSTPART',xpart(1)
C          write (6,*) "DUPART..",dupart(1),upart(1),xpart(1)
C
C       write (6,*) 'MGTSTPART',nfam,npart
C        write (6,*) 'MGTSTPART',conf(1),conf(nfam*npart + 1),
C     $                           conf(nfam*npart*2 + 1)
C        write (6,*) 'MGTSTPART',xpart(1),xpart(nfam + 1),
C     $                           xpart(nfam*2 + 1)

         CALL ERRR (501,' MGTSTPART! ')
            tsearch = 0.0
            tgetvel = 0.0
            tgetder = 0.0
            ttstorient = 0.0
            ttstpos = 0.0
            tbpart = 0.0
            tges   = 0.0

C
C
C
C------------------------------------------ SETUP
C
      SHAPE_PART = (ASPECT**2 - 1.0)/(ASPECT**2 + 1.0)
      DR = FRANDOM**2 / 2.0

      IF (ITSTEP .EQ. 1) THEN

      WRITE (6,*) 'BROWNIAN DIFFUSIVITY IS: ',DR
      WRITE (6,*) 'THE SHAPE FACTOR IS: ',SHAPE_PART
      WRITE (6,*) 'PECLET-NUMBER IS: ',SHEAR / DR

      ENDIF

      call brenner(rmu0,rmu1,rmu2,rmu3,rmu4,DR,VOL_FRAC,ASPECT)
C      write (6,*)' rmu0,... ',rmu0,rmu1,rmu2,rmu3,rmu4



C
C------------------------------------------  GERECHNET WIRD ERST NUR
C                                            AUF DEM MAXIMALEN LEVEL
C                                            (FEINSTES GITTER)
      ILEVEL = MAXLEVEL

      DO I = 1,NOFTST(ILEVEL)
         IGRID = IGRDOFTST(I,ILEVEL)

         CALL MGDPB (KK,JJ,II,IP3,IP2,IP1,IBB,IBU,
     $        NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
         
C------------------------------------------- SEARCHING INDEX

            t1 = GETSEC(0)

            CALL SEARCHINDEX(KK,Z(IP1),DZ(IP1),DDZ(IP1),
     $           NZGRAE(IGRID),NFAM,XPART(2*NFAM + 1),KINDEX,
     $           HILF1(IP3))

            CALL SEARCHINDEX(JJ,Y(IP1),DY(IP1),DDY(IP1),
     $           NYGRAE(IGRID),NFAM,XPART(  NFAM + 1),JINDEX,
     $           HILF1(IP3))

            CALL SEARCHINDEX(II,X(IP1),DX(IP1),DDX(IP1),
     $           NXGRAE(IGRID),NFAM,XPART(         1),IINDEX,
     $           HILF1(IP3))

            t2 = GETSEC(0)

C       write (6,*)'mgtstpart kindex',myid,kindex(100)
C       write (6,*)'mgtstpart jindex',myid,jindex(100)
C       write (6,*)'mgtstpart iindex',myid,iindex(100)

C------------------------------------------- GETTING VELOCITIES

            CALL GETVELOCITIES (KK,JJ,II,  X(IP1),   Y(IP1),   Z(IP1),
     $                                    DX(IP1),  DY(IP1),  DZ(IP1),
     $                                  RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                          U(IP3),V(IP3),W(IP3),
     $                          HILF1(IP3),HILF2(IP3),
     $                          H2D1(IP1),H2D2(IP1),H2D3(IP1),
     $                          NFAM,XPART(1),
     $                          XPART(  NFAM + 1),XPART(2*NFAM + 1),
     $                          IINDEX,JINDEX,KINDEX,UPART,
     $                          FXPART,FYPART,FZPART)

C       write (6,*)'mgtstpart upart',myid,upart(100)
            t3 = GETSEC(0)

C         goto 1000
C------------------------------------------ GETTING VELOCITY DERIVATIVES

            CALL GETDERIVATIVES(KK,JJ,II,  X(IP1),   Y(IP1),   Z(IP1),
     $                                    DX(IP1),  DY(IP1),  DZ(IP1),
     $                                   RDX(IP1), RDY(IP1), RDZ(IP1),
     $                                  RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                          U(IP3),V(IP3),W(IP3),
     $                          HILF1(IP3),HILF2(IP3),
     $                          H2D1(IP1),H2D2(IP1),H2D3(IP1),
     $                          NFAM,XPART(1),
     $                          XPART(  NFAM + 1),XPART(2*NFAM + 1),
     $                          IINDEX,JINDEX,KINDEX,DUPART,
     $                          FXPART,FYPART,FZPART)

C       write (6,*)'mgtstpart dupart',myid,dupart(100)
            t4 = GETSEC(0)



C------------------------------ COMPUTING RANDOM NUMBERS

            CALL RANDOM_NUMBERS (NPART*3,RRANDOM,IRANDOM, DT, FRANDOM)

C------------------------------------------- NEW ORIENTATION

               CALL TSTORIENT (NPART,NFAM,DUPART,CONF,
     $                      DT,SHAPE_PART,IRANDOM,FRANDOM,RRANDOM)

C       write (6,*)'mgtstpart conf',myid,conf(100)
            t5 = GETSEC(0)


            CALL CAL_RR (NPART, NFAM, CONF, LOUT,
     $                    IT, DT, RR)
            t5a1 = GETSEC(0)
            CALL CAL_RRRR (NPART, NFAM, CONF, LOUT,
     $                     IT, DT, DUPART,RRRREE)

            t5a2 = GETSEC(0)
C       write (6,*)'mgtstpart rrrr',myid,rr(100),rrrree(100)
C------------------------------------------- PARTICLE STRESSES

            CALL STRESS (NFAM, RR, RRRREE, DUPART,
     $                   TAU, ITTOT, DT,
     $                   RMU0,RMU1,RMU2,RMU3,RMU4,GMOL,VOL_FRAC,LOUT)

            t5a = GETSEC(0)

C       write (6,*)'mgtstpart tau',myid,tau(100)

C------------------------------------------- SETTING STRESSES FOR GRID
C------------------------------------------- OLD STRESSES ARE WEIGHTED BY
C------------------------------------------- FOLD = 0.1

            FOLD = 0.0001

C         IF (ITSTEP .EQ. 0) THEN
            CALL SETS (KK,JJ,II,KK,JJ,II,TAU11(IP3),0.0)
            CALL SETS (KK,JJ,II,KK,JJ,II,TAU12(IP3),0.0)
            CALL SETS (KK,JJ,II,KK,JJ,II,TAU13(IP3),0.0)
            CALL SETS (KK,JJ,II,KK,JJ,II,TAU21(IP3),0.0)
            CALL SETS (KK,JJ,II,KK,JJ,II,TAU22(IP3),0.0)
            CALL SETS (KK,JJ,II,KK,JJ,II,TAU23(IP3),0.0)
            CALL SETS (KK,JJ,II,KK,JJ,II,TAU31(IP3),0.0)
            CALL SETS (KK,JJ,II,KK,JJ,II,TAU32(IP3),0.0)
            CALL SETS (KK,JJ,II,KK,JJ,II,TAU33(IP3),0.0)
C         ENDIF

            t5b2 = GETSEC(0)
C------------------------------------------ HILF1 CONTAINS NUMBER OF
C------------------------------------------ PARTICLES IN CELL

            CALL SETS (KK,JJ,II,KK,JJ,II,HILF1(IP3),FOLD)

            CALL PARTICLE_STRESS (KK,JJ,II,  X(IP1),   Y(IP1),   Z(IP1),
     $                                    DX(IP1),  DY(IP1),  DZ(IP1),
     $                                  RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                          TAU11(IP3),TAU12(IP3),TAU13(IP3),
     $                          TAU21(IP3),TAU22(IP3),TAU23(IP3),
     $                          TAU31(IP3),TAU32(IP3),TAU33(IP3),
     $                          HILF1(IP3),
     $                          H2D1(IP1),H2D2(IP1),H2D3(IP1),
     $                          NFAM,XPART(1),
     $                          XPART(  NFAM + 1),XPART(2*NFAM + 1),
     $                          IINDEX,JINDEX,KINDEX,TAU)

            t5b = GETSEC(0)


C         goto 1000
          CALL SMOOTH_STRESS (KK,JJ,II,  X(IP1),   Y(IP1),   Z(IP1),
     $                                    DX(IP1),  DY(IP1),  DZ(IP1),
     $                                   DDX(IP1), DDY(IP1), DDZ(IP1),
     $                                   RDX(IP1), RDY(IP1), RDZ(IP1),
     $                                  RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                          TAU11(IP3),
     $                          HILF2(IP3),HILF1(IP3))

          CALL SMOOTH_STRESS (KK,JJ,II,  X(IP1),   Y(IP1),   Z(IP1),
     $                                    DX(IP1),  DY(IP1),  DZ(IP1),
     $                                   DDX(IP1), DDY(IP1), DDZ(IP1),
     $                                   RDX(IP1), RDY(IP1), RDZ(IP1),
     $                                  RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                          TAU12(IP3),
     $                          HILF2(IP3),HILF1(IP3))

          CALL SMOOTH_STRESS (KK,JJ,II,  X(IP1),   Y(IP1),   Z(IP1),
     $                                    DX(IP1),  DY(IP1),  DZ(IP1),
     $                                   DDX(IP1), DDY(IP1), DDZ(IP1),
     $                                   RDX(IP1), RDY(IP1), RDZ(IP1),
     $                                  RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                          TAU13(IP3),
     $                          HILF2(IP3),HILF1(IP3))

C          CALL SMOOTH_STRESS (KK,JJ,II,  X(IP1),   Y(IP1),   Z(IP1),
C     $                                    DX(IP1),  DY(IP1),  DZ(IP1),
C     $                                   DDX(IP1), DDY(IP1), DDZ(IP1),
C     $                                   RDX(IP1), RDY(IP1), RDZ(IP1),
C     $                                  RDDX(IP1),RDDY(IP1),RDDZ(IP1),
C     $                          TAU21(IP3),
C     $                          HILF2(IP3),HILF1(IP3))

          CALL SMOOTH_STRESS (KK,JJ,II,  X(IP1),   Y(IP1),   Z(IP1),
     $                                    DX(IP1),  DY(IP1),  DZ(IP1),
     $                                   DDX(IP1), DDY(IP1), DDZ(IP1),
     $                                   RDX(IP1), RDY(IP1), RDZ(IP1),
     $                                  RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                          TAU22(IP3),
     $                          HILF2(IP3),HILF1(IP3))

          CALL SMOOTH_STRESS (KK,JJ,II,  X(IP1),   Y(IP1),   Z(IP1),
     $                                    DX(IP1),  DY(IP1),  DZ(IP1),
     $                                   DDX(IP1), DDY(IP1), DDZ(IP1),
     $                                   RDX(IP1), RDY(IP1), RDZ(IP1),
     $                                  RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                          TAU23(IP3),
     $                          HILF2(IP3),HILF1(IP3))

C          CALL SMOOTH_STRESS (KK,JJ,II,  X(IP1),   Y(IP1),   Z(IP1),
C     $                                    DX(IP1),  DY(IP1),  DZ(IP1),
C     $                                   DDX(IP1), DDY(IP1), DDZ(IP1),
C     $                                   RDX(IP1), RDY(IP1), RDZ(IP1),
C     $                                  RDDX(IP1),RDDY(IP1),RDDZ(IP1),
C     $                          TAU31(IP3),
C     $                          HILF2(IP3),HILF1(IP3))

C          CALL SMOOTH_STRESS (KK,JJ,II,  X(IP1),   Y(IP1),   Z(IP1),
C     $                                    DX(IP1),  DY(IP1),  DZ(IP1),
C     $                                   DDX(IP1), DDY(IP1), DDZ(IP1),
C     $                                   RDX(IP1), RDY(IP1), RDZ(IP1),
C     $                                  RDDX(IP1),RDDY(IP1),RDDZ(IP1),
C     $                          TAU32(IP3),
C     $                          HILF2(IP3),HILF1(IP3))

          CALL SMOOTH_STRESS (KK,JJ,II,  X(IP1),   Y(IP1),   Z(IP1),
     $                                    DX(IP1),  DY(IP1),  DZ(IP1),
     $                                   DDX(IP1), DDY(IP1), DDZ(IP1),
     $                                   RDX(IP1), RDY(IP1), RDZ(IP1),
     $                                  RDDX(IP1),RDDY(IP1),RDDZ(IP1),
     $                          TAU33(IP3),
     $                          HILF2(IP3),HILF1(IP3))

C            write (6,*) 'smooth_stress finished'
C       write (6,*)'mgtstpart tau11',myid,tau11(395367)
C       write (6,*)'mgtstpart tau12',myid,tau12(395367)
C       write (6,*)'mgtstpart tau13',myid,tau13(395367)
C       write (6,*)'mgtstpart tau21',myid,tau21(395367)
C       write (6,*)'mgtstpart tau22',myid,tau22(395367)
C       write (6,*)'mgtstpart tau23',myid,tau23(395367)
C       write (6,*)'mgtstpart tau31',myid,tau31(395367)
C       write (6,*)'mgtstpart tau32',myid,tau32(395367)
C       write (6,*)'mgtstpart tau33',myid,tau33(395367)

            t5c = GETSEC(0)
 1000  continue
C------------------------------------------- NEW POSITION

            CALL TSTPOS (NFAM,UPART,
     $                   XPART(1), XPART(NFAM + 1), XPART(2*NFAM + 1),
     $                   DT)
C            write (6,*) 'tstpos finished'
            t6 = GETSEC(0)
            CALL BPART (KK,JJ,II,  X(IP1),   Y(IP1),   Z(IP1),
     $                            DX(IP1),  DY(IP1),  DZ(IP1),
     $                   NFAM,
     $                   XPART(1), XPART(NFAM + 1), XPART(2*NFAM + 1),
     $                   NFRO,NBAC,NRGT,NLFT,NBOT,NTOP)

            t7 = GETSEC(0)

C       write (6,*)'mgtstpart xpart',myid,xpart(100)
C            write (6,*) 'bpart finished'

C------------------------------------------- LOOP OVER CLUSTER FINISHED

C------------------------------------------- OUTPUT?
C        write (99,*) xpart(100),xpart(nfam+100),xpart(nfam*2+100)



      ENDDO

      IF (ITSTEP .EQ. 1) THEN

            tsearch = t2-t1
            tgetvel = t3-t2
            tgetder = t4-t3
            ttstorient = t5-t4
            trr = t5a1-t5
            trrrr = t5a2-t5a1
            tstress = t5a-t5a2
            tset    = t5b2 - t5a
            tparticle_stress = t5b-t5b2
            tsmooth_stress = t5c-t5b
            ttstpos = t6-t5c
            tbpart = t7-t6
            tges   = t7 - t1

        if (myid .eq. 0) then
            write (6,*)'   search index:',tsearch
            write (6,*)' get velocities:',tgetvel
            write (6,*)'get derivatives:',tgetder
            write (6,*)'   orientations:',ttstorient
            write (6,*)'             rr:',trr
            write (6,*)'           rrrr:',trrrr
            write (6,*)'         stress:',tstress
            write (6,*)'     set tau_ij:',tset
            write (6,*)'particle_stress:',tparticle_stress
            write (6,*)'  smooth_stress:',tsmooth_stress
            write (6,*)'  new positions:',ttstpos
            write (6,*)' boundary cond.:',tbpart
            write (6,*)'            all:',tges
        endif

        ENDIF

         RETURN
         END
