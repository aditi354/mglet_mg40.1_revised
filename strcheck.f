










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
      SUBROUTINE STRCHECK
     $ (IDIM3D,IDIM2D,IDIM1D,NBUF,NBND,IDIMA,IDIM1L,IDIM2L,IDIMF,IIDENT)
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C     19. 5.93 (MM) : ORIGINAL
C     27.12.93 (MM) : KOPF GEAENDERT, JETZT AUCH AUFRUF VON SETCOMGRID,
C                     SETCOBOUND UND SETCOLEVEL
C     26. 5.95 (MM) : MPI EINGEFUEHRT
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
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

      PARAMETER ( NBSP  =  8 )
      
      COMMON /COGRDDEF/
     &                 NBX,      NBY,      NBZ,
     &                 NTXB,     NTYB,     NTZB,
     &                 DLX,      DLY,      DLZ,
     &                 NLXB,     NLYB,     NLZB,
     &                 SX,       SY,       SZ,
     &                 NRXB,     NRYB,     NRZB,
     &                 DRX,      DRY,      DRZ

      INTEGER
     &       NBX(MAXGRIDS),       NBY(MAXGRIDS),       NBZ(MAXGRIDS),
     & NTXB(NBSP,MAXGRIDS), NTYB(NBSP,MAXGRIDS), NTZB(NBSP,MAXGRIDS),
     & NLXB(NBSP,MAXGRIDS), NLYB(NBSP,MAXGRIDS), NLZB(NBSP,MAXGRIDS),
     & NRXB(NBSP,MAXGRIDS), NRYB(NBSP,MAXGRIDS), NRZB(NBSP,MAXGRIDS)

      REAL
     &  DLX(NBSP,MAXGRIDS),  DLY(NBSP,MAXGRIDS),  DLZ(NBSP,MAXGRIDS),
     &   SX(NBSP,MAXGRIDS),   SY(NBSP,MAXGRIDS),   SZ(NBSP,MAXGRIDS),
     &  DRX(NBSP,MAXGRIDS),  DRY(NBSP,MAXGRIDS),  DRZ(NBSP,MAXGRIDS)

      
      COMMON /COGRDOLD/
     &                 KMXOLD,  JMXOLD,  IMXOLD,
     &                 LEVELOLD,    LCHILDOLD, IPARENTOLD,
     &                 IPOSOLD, JPOSOLD, KPOSOLD,
     &                 XMINOLD,    YMINOLD,    ZMINOLD,
     &                 XTOTOLD,     YTOTOLD,     ZTOTOLD,   
     &                 NBOCOLD,
     &                 FRONTOLD,    BACKOLD,     RIGHTOLD,     LEFTOLD,
     &                 BOTTOMOLD,   TOPOLD,      CUBEOLD,
     &                 IFRNBROLD, IBANBROLD, IRINBROLD, ILENBROLD,
     &                 IBONBROLD, ITONBROLD,
     &                 IBPOSOLD,   JBPOSOLD,    KBPOSOLD,
     &                 IBANFOLD,  IBENDOLD,
     &                 JBANFOLD,  JBENDOLD,
     &                 KBANFOLD,  KBENDOLD,
     &                 XHOMOGOLD,   YHOMOGOLD,   ZHOMOGOLD,
     &                 LPLEVELOLD

      INTEGER
     &       KMXOLD(MAXGRIDS), JMXOLD(MAXGRIDS), IMXOLD(MAXGRIDS),
     &       LEVELOLD(MAXGRIDS), 
     &       IPARENTOLD(MAXGRIDS),
     & IPOSOLD(MAXGRIDS), JPOSOLD(MAXGRIDS), KPOSOLD(MAXGRIDS),

     &  NBOCOLD(7,MAXGRIDS),
     &  IFRNBROLD(MAXBOCONDS,MAXGRIDS), IBANBROLD(MAXBOCONDS,MAXGRIDS),
     &  IRINBROLD(MAXBOCONDS,MAXGRIDS), ILENBROLD(MAXBOCONDS,MAXGRIDS),
     &  IBONBROLD(MAXBOCONDS,MAXGRIDS), ITONBROLD(MAXBOCONDS,MAXGRIDS),
     &  IBPOSOLD(MAXBOCONDS,7,MAXGRIDS),JBPOSOLD(MAXBOCONDS,7,MAXGRIDS),
     &  KBPOSOLD(MAXBOCONDS,7,MAXGRIDS),
     &  IBANFOLD(MAXBOCONDS,7,MAXGRIDS),IBENDOLD(MAXBOCONDS,7,MAXGRIDS),
     &  JBANFOLD(MAXBOCONDS,7,MAXGRIDS),JBENDOLD(MAXBOCONDS,7,MAXGRIDS),
     &  KBANFOLD(MAXBOCONDS,7,MAXGRIDS),KBENDOLD(MAXBOCONDS,7,MAXGRIDS)


      REAL
     &       XTOTOLD(MAXGRIDS), YTOTOLD(MAXGRIDS), ZTOTOLD(MAXGRIDS),
     &       XMINOLD(MAXGRIDS), YMINOLD(MAXGRIDS), ZMINOLD(MAXGRIDS)

      CHARACTER (LEN=16)
     &   FRONTOLD(MAXBOCONDS,MAXGRIDS),   BACKOLD(MAXBOCONDS,MAXGRIDS),
     &   RIGHTOLD(MAXBOCONDS,MAXGRIDS),   LEFTOLD(MAXBOCONDS,MAXGRIDS),
     &  BOTTOMOLD(MAXBOCONDS,MAXGRIDS),    TOPOLD(MAXBOCONDS,MAXGRIDS),
     &    CUBEOLD(MAXBOCONDS,MAXGRIDS)

      LOGICAL
     &      LCHILDOLD(MAXGRIDS),
     &  XHOMOGOLD(MAXGRIDS), YHOMOGOLD(MAXGRIDS), ZHOMOGOLD(MAXGRIDS),
     &      LPLEVELOLD(MAXGRIDS)

      PARAMETER (NMBODOLD = 100)
      
      COMMON /COBODYOLD/
     &               NBODOLD,  CTYPO,
     &               IB1O,   IB2O,   JB1O,   JB2O,   KB1O,   KB2O,
     &               XB1O,   XB2O,   YB1O,   YB2O,   ZB1O,   ZB2O,
     &               NCOUNO, XMITO,  HEIGHTO, ALPHAO, CDIRO

      INTEGER
     &       NBODOLD,
     &       IB1O(NMBODOLD),   IB2O(NMBODOLD),
     &       JB1O(NMBODOLD),   JB2O(NMBODOLD),
     &       KB1O(NMBODOLD),   KB2O(NMBODOLD),
     &       NCOUNO(NMBODOLD) 

      REAL
     &       XB1O(NMBODOLD),   XB2O(NMBODOLD),
     &       YB1O(NMBODOLD),   YB2O(NMBODOLD),
     &       ZB1O(NMBODOLD),   ZB2O(NMBODOLD),
     &      XMITO(NMBODOLD),HEIGHTO(NMBODOLD),ALPHAO(NMBODOLD)

      CHARACTER (LEN=16) CTYPO(NMBODOLD),CDIRO(NMBODOLD)
 

      INTEGER IIDENT(100)
      INTEGER IFAILALL
C
C
            IFAILALL = 0
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      WRITE (6,*) 
     $'***************************************************************'
      WRITE (6,*) 
     $'       PRUEFUNG DER GITTERDEFINITIONEN AUF KONSISTENZ '
      WRITE (6,*) 
     $'              MIT DEM VORANGEGANGENEN LAUF'
      WRITE (6,*) 
      WRITE (6,*) 
      WRITE (6,*) ' ES WURDEN ',NGRID,'  GITTER DEFINIERT'
      WRITE (6,*) ' ES SIND   ',NGRDOLD,'  GITTER VORHANDEN'

      DO IGRID=1,NGRDOLD

         CALL GRDOUT (IGRID,  KMXOLD,   JMXOLD,    IMXOLD,
     &                      LEVELOLD,LCHILDOLD,IPARENTOLD,
     &                       IPOSOLD,  JPOSOLD,   KPOSOLD,
     &                       XTOTOLD,  YTOTOLD,   ZTOTOLD,
     &             FRONTOLD,BACKOLD,RIGHTOLD,LEFTOLD,
     &             BOTTOMOLD,TOPOLD,CUBEOLD,
     &             XHOMOGOLD,YHOMOGOLD,ZHOMOGOLD,
     &             NXGRAEOLD,NYGRAEOLD,NZGRAEOLD,LPLEVELOLD)

         CALL GRDKON(IGRID,IFAIL,IIDENT)

         IF (IFAIL.GT.0) THEN

            WRITE (6,*) ' BITTE STEUERFILE KORRIGIEREN !!!!!!!!!!!!!!'
            IFAILALL = IFAILALL + 1
C            STOP 
C     $      'STEUERFILE NICHT KONSISTENT MIT ALTEM LAUF, SIEHE OUTPUT'
         ENDIF

         WRITE (6,*) 
         WRITE (6,*) 
         WRITE (6,*) 

      ENDDO


      WRITE (6,*) 
      WRITE (6,*) 

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      WRITE (6,*)
     $'***************************************************************'
      WRITE (6,*)
     $'       PRUEFUNG DER GITTERDEFINITIONEN AUF ZULAESSIGKEIT '

      WRITE (6,*) 
     $ ' ES WERDEN ',NGRDDFD-NGRDOLD,'  GITTER NEU GEWUENSCHT'

      DO IGRID= NGRDOLD+1,NGRDDFD

         CALL GRDCHK (IGRID,IFAIL)

         IF (IFAIL.GT.0) THEN

            WRITE (6,*) ' BITTE STEUERFILE KORRIGIEREN !!!!!!!!!!!!!!'
            WRITE (6,*) ' GITTER NR.',IGRID,'  NICHT ZULAESSIG!!!!!'
            IFAILALL = IFAILALL + 1
C            STOP
C     $      'GEWUENSCHTES GITTER NICHT ZULAESSIG,  SIEHE OUTPUT'
         ELSE
            WRITE (6,*) 'GITTER',IGRID,'  O.K. !'
C           CALL SETGRDPRO (IGRID)
         ENDIF

         WRITE (6,*)

      ENDDO

C
C                         IM FEHLERFALL ABBRUCH
C
      IF ( IFAILALL .GT. 0) CALL ERRR (555,' STRCHECK ')
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

      DO IGRID= 1,NGRDDFD
         GRADPXOLD(IGRID) = GRADPX(IGRID) 
      END DO


      RETURN
C                            SETZEN DER RANDBED.-STEUERUNG

C      DO IGRID= 1,NGRDDFD

C         CALL SETCOBOUND  (IGRID,IGRID)

C         ITST = 0
C         IF (LTST(IGRID)) THEN
C            IF (.NOT.LSLICE(IGRID)) ITST=1
C         ENDIF

C         IVPIT = 0
C         IF (LVP(IGRID)) THEN
C            IF (.NOT.LSLICE(IGRID)) IVPIT=1
C         ENDIF

C         IF (     LSLICE(IGRID)) ISLCED=1
C         IF (.NOT.LSLICE(IGRID)) ISLCED=0

C         IPSDIR=0
C         IF (LPOISSONDIR(IGRID)) THEN
C             IPSDIR=1
C             IVPIT =0
C         ENDIF

C         CALL SETCOLEVEL 
C     $       (IGRID,LEVEL(IGRID),IVPIT,ITST,ISLCED,0,IPSDIR)

C      ENDDO
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                            SETZEN DER DIMENSIONIERUNGEN UND POINTER

C      DO IGRID= 1,NGRDDFD
      
C            CALL SETCOMGRID (IDIM3D,IDIM2D,IDIM1D,NBUF,NBND,
C     $                       IDIMA,IDIM1L,IDIM2L,IDIMF,
C     $                       IGRID,0,0,0,0)

C      ENDDO

C

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                            

      RETURN
      END
