










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
      SUBROUTINE ITSAMPLE (IGRID,IPCGES,DIVGMX,EPSU,EPSV,EPSW,
     $     ESUMG,ESUMS,WSSX,WSSY,WSSZ,WNSX,WNSY,WNSZ,
     $     UBULK,WALLSSX,WALLSSY,WALLSSZ
     $    ,IDIVMAX,JDIVMAX,KDIVMAX,XDIVMAX,YDIVMAX,ZDIVMAX,GRDIVMAX
     $     )


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
 


      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
C                                 FUER DIE KONTROLLE DES RECHENVORGANGES
C
      REAL     EPSU(MAXGRIDS),  EPSV(MAXGRIDS),  EPSW(MAXGRIDS)
      REAL     ESUMG(MAXGRIDS),  ESUMS(MAXGRIDS)
      REAL     WSSX(MAXGRIDS),  WSSY(MAXGRIDS),  WSSZ(MAXGRIDS)
      REAL     WNSX(MAXGRIDS),  WNSY(MAXGRIDS),  WNSZ(MAXGRIDS)
      REAL     WALLSSX(6,MAXGRIDS),
     $         WALLSSY(6,MAXGRIDS),
     $         WALLSSZ(6,MAXGRIDS)
      REAL     UBULK(MAXGRIDS)
      INTEGER IPCGES(MAXGRIDS)
      
      REAL    DIVGMX(MAXGRIDS),XDIVMAX(MAXGRIDS),YDIVMAX(MAXGRIDS),
     $        ZDIVMAX(MAXGRIDS)

      INTEGER 
     $ IDIVMAX(MAXGRIDS),JDIVMAX(MAXGRIDS),KDIVMAX(MAXGRIDS),
     $ GRDIVMAX(MAXGRIDS)


          IF (LSLICE(IGRID)) THEN
                G = 1./FLOAT(NOFSLCHILDS(IGRID))
                GX= 1./FLOAT(NXSLICE(IGRID)) 

                      IPCGES(IGRID) = 0
                      DIVGMX(IGRID) = 0.0
                      EPSU(IGRID) = 0.0
                      EPSV(IGRID) = 0.0
                      EPSW(IGRID) = 0.0
                      ESUMG(IGRID) = 0.0
                      ESUMS(IGRID) = 0.0
                      WSSX(IGRID) = 0.0
                      WSSY(IGRID) = 0.0
                      WSSZ(IGRID) = 0.0
                      WNSX(IGRID) = 0.0
                      WNSY(IGRID) = 0.0
                      WNSZ(IGRID) = 0.0
                      UBULK(IGRID) = 0.0
                      GRADPX(IGRID) = 0.0
                      DO IDIR = 1,6
                         WALLSSX(IDIR,IGRID) = 0.0
                         WALLSSY(IDIR,IGRID) = 0.0
                         WALLSSZ(IDIR,IGRID) = 0.0
                      ENDDO     
C
                      IDIVMAX(IGRID) =  0   
                      JDIVMAX(IGRID) =  0   
                      KDIVMAX(IGRID) =  0   
                      XDIVMAX(IGRID) =  0.0 
                      YDIVMAX(IGRID) =  0.0 
                      ZDIVMAX(IGRID) =  0.0                      
                      GRDIVMAX(IGRID)=  0                       

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
             DO I = 1,NOFSLCHILDS(IGRID)
                ISUB = IGRDOFSLCHILD(I,IGRID)

                      IPCGES(IGRID) = MAX(IPCGES(IGRID),IPCGES(ISUB))
C                      DIVGMX(IGRID) = MAX(DIVGMX(IGRID),DIVGMX(ISUB))
                      EPSU(IGRID) = EPSU(IGRID) + EPSU(ISUB)*G
                      EPSV(IGRID) = EPSV(IGRID) + EPSV(ISUB)*G
                      EPSW(IGRID) = EPSW(IGRID) + EPSW(ISUB)*G
                      ESUMG(IGRID) = ESUMG(IGRID) + ESUMG(ISUB)*G
                      ESUMS(IGRID) = ESUMS(IGRID) + ESUMS(ISUB)*G
                      WSSX(IGRID) = WSSX(IGRID) + WSSX(ISUB)*G
                      WSSY(IGRID) = WSSY(IGRID) + WSSY(ISUB)*G
                      WSSZ(IGRID) = WSSZ(IGRID) + WSSZ(ISUB)*G
                      WNSX(IGRID) = WNSX(IGRID) + WNSX(ISUB)*G
                      WNSY(IGRID) = WNSY(IGRID) + WNSY(ISUB)*G
                      WNSZ(IGRID) = WNSZ(IGRID) + WNSZ(ISUB)*G
                      UBULK(IGRID) = UBULK(IGRID) + UBULK(ISUB)*G
                      GRADPX(IGRID) = GRADPX(IGRID) + GRADPX(ISUB)*G
                     
                      DO IDIR = 1,6
                         WALLSSX(IDIR,IGRID) = WALLSSX(IDIR,IGRID) *G
     $                                       + WALLSSX(IDIR,ISUB)*G
                         WALLSSY(IDIR,IGRID) = WALLSSY(IDIR,IGRID)*G 
     $                                       + WALLSSY(IDIR,ISUB)*G
                         WALLSSZ(IDIR,IGRID) = WALLSSZ(IDIR,IGRID)*G
     $                                       + WALLSSZ(IDIR,ISUB) *G
                      ENDDO

                     IF (ABS(DIVGMX(ISUB)) .GT. ABS(DIVGMX(IGRID))) THEN                     
                       DIVGMX(IGRID)  = DIVGMX(ISUB)
                       IDIVMAX(IGRID) = IDIVMAX(ISUB)
                       JDIVMAX(IGRID) = JDIVMAX(ISUB)
                       KDIVMAX(IGRID) = KDIVMAX(ISUB)
                       XDIVMAX(IGRID) = XDIVMAX(ISUB)
                       YDIVMAX(IGRID) = YDIVMAX(ISUB)
                       ZDIVMAX(IGRID) = ZDIVMAX(ISUB)                       
                       GRDIVMAX(IGRID)= ISUB                                           
                      ENDIF 


                
             ENDDO
          ENDIF        
       RETURN
       END
