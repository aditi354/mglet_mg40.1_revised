










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
      SUBROUTINE DOCHEAD (CIDENT,IIDENT,RIDENT,
     $                    ISELEP,ISAMPP,ISLINP,ISLIDI)
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C          SCHREIBT HEADERINFORMATION AUF ERGEBNISFILE
C          BINAER
C
C
C     16. 5.93 (MM) : ORIGINAL
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

      PARAMETER (NMBODY = 100)
      
      COMMON /COBODY/
     &               NBODY,  CTYP,
     &               IB1,   IB2,   JB1,   JB2,   KB1,   KB2,
     &               XB1,   XB2,   YB1,   YB2,   ZB1,   ZB2,
     &               NCOUN, XMIT,  HEIGHT, ALPHA, CDIR

      INTEGER
     &       NBODY,
     &       IB1(NMBODY),   IB2(NMBODY),
     &       JB1(NMBODY),   JB2(NMBODY),
     &       KB1(NMBODY),   KB2(NMBODY),
     &       NCOUN(NMBODY) 

      REAL
     &       XB1(NMBODY),   XB2(NMBODY),
     &       YB1(NMBODY),   YB2(NMBODY),
     &       ZB1(NMBODY),   ZB2(NMBODY),
     &      XMIT(NMBODY),HEIGHT(NMBODY),ALPHA(NMBODY)

      CHARACTER (LEN=16) CTYP(NMBODY),CDIR(NMBODY)
 

      CHARACTER (LEN=8) CIDENT(10)
      INTEGER     IIDENT(100)
      REAL        RIDENT(100)

      INTEGER     ISELEP(2*752)    ,ISLINP(ISLIDI),
     $            ISAMPP  (752*MAXGRIDS)

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      
      WRITE(4,6020)(CIDENT(N),N=1,10)
      WRITE(4,6025)(IIDENT(N),N=1,100)
      WRITE(4,6026)(RIDENT(N),N=1,100)

      WRITE(4,6025) NGRID,NBODY,NAUFP
      DO I=1,NGRID
         WRITE(4,6025) IMX(I),JMX(I),KMX(I)
         WRITE(4,*) LEVEL(I),LCHILD(I),IPARENT(I)
         WRITE(4,6025) IPOSITION(I),JPOSITION(I),KPOSITION(I)
         WRITE(4,6026) XTOT(I),YTOT(I),ZTOT(I)

            WRITE(4,6024) (NBOCD(J,I),J=1,7)
             
            WRITE(4,6020) ( FRONT(J,I),J=1, NBOCD(1,I))
            WRITE(4,6020) (  BACK(J,I),J=1, NBOCD(2,I))
            WRITE(4,6020) ( RIGHT(J,I),J=1, NBOCD(3,I))
            WRITE(4,6020) (  LEFT(J,I),J=1, NBOCD(4,I))
            WRITE(4,6020) (BOTTOM(J,I),J=1, NBOCD(5,I))
            WRITE(4,6020) (   TOP(J,I),J=1, NBOCD(6,I))
            WRITE(4,6020) (  CUBE(J,I),J=1, NBOCD(7,I))

            WRITE(4,6024) (IFRNBR(J,I),J=1, NBOCD(1,I))
            WRITE(4,6024) (IBANBR(J,I),J=1, NBOCD(2,I))
            WRITE(4,6024) (IRINBR(J,I),J=1, NBOCD(3,I))
            WRITE(4,6024) (ILENBR(J,I),J=1, NBOCD(4,I))
            WRITE(4,6024) (IBONBR(J,I),J=1, NBOCD(5,I))
            WRITE(4,6024) (IBONBR(J,I),J=1, NBOCD(6,I))

            DO ID = 1,7
             WRITE(4,6024) (IBPOS(J,ID,I),J=1,NBOCD(ID,I))
             WRITE(4,6024) (JBPOS(J,ID,I),J=1,NBOCD(ID,I))
             WRITE(4,6024) (KBPOS(J,ID,I),J=1,NBOCD(ID,I))
             WRITE(4,6024) (IBANF(J,ID,I),J=1,NBOCD(ID,I))
             WRITE(4,6024) (JBANF(J,ID,I),J=1,NBOCD(ID,I))
             WRITE(4,6024) (KBANF(J,ID,I),J=1,NBOCD(ID,I))
             WRITE(4,6024) (IBEND(J,ID,I),J=1,NBOCD(ID,I))
             WRITE(4,6024) (JBEND(J,ID,I),J=1,NBOCD(ID,I))
             WRITE(4,6024) (KBEND(J,ID,I),J=1,NBOCD(ID,I))
            ENDDO

         WRITE(4,*) XHOMOG(I),YHOMOG(I),ZHOMOG(I)
         WRITE(4,*) LPLEVEL(I)

      ENDDO

      DO I=1,NBODY
         WRITE(4,6020) CTYP(I)
         WRITE(4,6025)
     $        IB1(I),IB2(I),JB1(I),JB2(I),KB1(I),KB2(I)
         WRITE(4,6026) 
     $        XB1(I),XB2(I),YB1(I),YB2(I),ZB1(I),ZB2(I)
         WRITE(4,6021) 
     $        NCOUN(I),XMIT(I),HEIGHT(I),ALPHA(I),CDIR(I)
      ENDDO


      WRITE(4,6025) (ISELEP(I),I=1,2*752)
      WRITE(4,6025) (ISAMPP(I),I=1,752*NGRID)
      WRITE(4,6015) (ISLINP(N),N=1,NAUFP)

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

      RETURN

 6020 FORMAT(5(A8,2X))
 6021 FORMAT(I9,1X,3(E12.5E3,1X),A8)
 6015 FORMAT(5(I19,1X))
 6024 FORMAT(6(I9,1X))
 6025 FORMAT(5(I9,1X))
 6026 FORMAT(6(E12.5E3,1X))

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                                    FEHLERBEHANDLUNG


      END
