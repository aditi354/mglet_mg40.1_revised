










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
      SUBROUTINE IBFIELD (KK,JJ,II,KMX,JMX,IMX,NBND,
     $                   DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,
     $                   B,BP,BU,BV,BW,IDIMF,NGRP3,HILF	
     $                   ,ZCUB,NCUB,IGRID
     $                   )

C*MGLET*****************************************************************
C   I B F I E L D    BELEGT DAS KOERPERERKENNUNGSFELD (B-FELD)
C                    IN DEN DEFINIERTEN KOERPERN MIT -1 ODER -2
C                    AUSSERHALB WIRD ES MIT EINER SEHR GROSSEN 
C                    ZAHL BELEGT
C*MGLET********************************************* M.MANHART 08.04.92
C                                         GEAENDERT:           18. 3.93
C                                                AUS INILET:    8. 4.93
C
C  PARAMETER
C             KK, JJ, II     - ARRAYGRENZEN
C             KMX, JMX, IMX  - ANZAHL DER GITTERPKTE. IN Z-,Y-,X-DIR.
C             X(II)          - KOORDINATEN IN X-RICHTUNG
C             Y(JJ)          - KOORDINATEN IN Y-RICHTUNG
C             Z(KK)          - KOORDINATEN IN Z-RICHTUNG
C             DDX(II)        - X-KANTENLAENGE DER KONTROLLVOLUMINA
C             DDY(JJ)        - Y-KANTENLAENGE DER KONTROLLVOLUMINA
C             DDZ(KK)        - Z-KANTENLAENGE DER KONTROLLVOLUMINA
C             DX(II)         - X-ABSTAND DER GITTERPUNKTE
C             DY(JJ)         - Y-ABSTAND DER GITTERPUNKTE
C             DZ(KK)         - Z-ABSTAND DER GITTERPUNKTE
C             B(KK,JJ,II)    - WANDERKENNUNGSFELD
C
C           10.02.03 (TB):   BT FIELD INCLUDED, BUT FRED_BODY TREAT-
C                            MENT NOT YET IMPLEMENTED
C
C*STAR******************************************************************
C

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

      COMMON /CONLES/  CONV2S,CMUE,CAPPA,ECONST
      SAVE   /CONLES/
C
C
      INTEGER KK, JJ, II, KMX, JMX, IMX
     $        NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID
      INTEGER I,J,K
C
C
C
C
      REAL     X(II), Y(JJ), Z(KK),
     $        DX(II), DY(JJ), DZ(KK), DDX(II), DDY(JJ), DDZ(KK),
     $        B(KK,JJ,II),BP(KK,JJ,II),BU(KK,JJ,II),
     $        BV(KK,JJ,II),BW(KK,JJ,II) 
C
C
C
C
      REAL        HILF(IDIMF,NGRP3)
C

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C
C                                WERT DES B-FELDES IM KOERPERINNEREN
C                                      SLIP:   -2.0
C                                    NOSLIP:   -1.0
C
      IF( NCUB .EQ. 5) BPAR = -1.0
      IF( NCUB .EQ. 6) BPAR = -2.0


      ZCUB = 0.0

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C
C                                 BELEGUNG DES WANDERKENNUNGSFELDES
C
         CALL SETS   (KK,JJ,II,KMX,JMX,IMX,B,GREAT)
C
      DO 40 IBODY=1,NBODY
      IF (CTYP(IBODY)(1:8).EQ.'COUNINGH') THEN
          ZCUB = HEIGHT(IBODY)
          CALL SBCOUN
     $            (KK,JJ,II,KMX,JMX,IMX,
     $             NCOUN(IBODY),XMIT(IBODY),HEIGHT(IBODY),ALPHA(IBODY),
     $             B,HILF,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $             XTOT(IGRID),YTOT(IGRID),ZTOT(IGRID),BPAR)

      ELSEIF (CTYP(IBODY)(1:8).EQ.'ZYLINDER') THEN

          ZCUB = MIN(XB2(IBODY),YB2(IBODY),ZB2(IBODY))
          CALL SBZYL 
     $            (KK,JJ,II,KMX,JMX,IMX,
     $             XB1(IBODY),XB2(IBODY),YB1(IBODY),YB2(IBODY),
     $             ZB1(IBODY),ZB2(IBODY),CDIR(IBODY),
     $             B,HILF,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,BPAR)

      ELSEIF (CTYP(IBODY)(1:7).EQ.'HEMISPH') THEN

          ZCUB = MIN(XB2(IBODY),YB2(IBODY),ZB2(IBODY))
          CALL SBHEMI
     $            (KK,JJ,II,KMX,JMX,IMX,
     $             XB1(IBODY),XB2(IBODY),YB1(IBODY),YB2(IBODY),
     $             ZB1(IBODY),ZB2(IBODY),
     $             B,HILF,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,BPAR)

      ELSEIF (CTYP(IBODY)(1:5).EQ.'KUBUS') THEN

          IF (IB1(IBODY).EQ.0) THEN
             ZCUB = MIN((XB2(IBODY)-XB1(IBODY)),
     $                  (YB2(IBODY)-YB1(IBODY)),
     $                  (ZB2(IBODY)-ZB1(IBODY)))
          ELSE
             ZCUB = MIN((X(IB2(IBODY))-X(IB1(IBODY))),
     $                  (Y(JB2(IBODY))-Y(JB1(IBODY))),
     $                  (Z(KB2(IBODY))-Z(KB1(IBODY))))
          ENDIF
          CALL SBCUB 
     $            (KK,JJ,II,KMX,JMX,IMX,
     $             IB1(IBODY),IB2(IBODY),JB1(IBODY),JB2(IBODY),
     $             KB1(IBODY),KB2(IBODY),
     $             XB1(IBODY),XB2(IBODY),YB1(IBODY),YB2(IBODY),
     $             ZB1(IBODY),ZB2(IBODY),
     $             B,HILF,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,BPAR)
CTBC  100203: FRED_BODY of Scalar Transport not yet implented
      ELSEIF (CTYP(IBODY)(1:8).EQ.'NEWBLOCK') THEN
         CALL MGBASB(NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
         write(*,*) 'PERFORMING BLOCKING'
C        BU wird als temporaeres FEld verwendet, 
C        da es glich gesetzt wird
         DO I = 1,II
         DO J = 1,JJ
         DO K = 1,KK 
            BU(I + II*(J-1) + II*JJ*(K-1), 1, 1) = BP(K,J,I)
         ENDDO
         ENDDO
         ENDDO
         CALL BLOCKGRID(II,JJ,KK,0,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                 IGRID,NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,
     $        BU,1)
         WRITE(*,*) 'BLOCKING COMPLETED'
         CALL GETQUAD(KK,JJ,II,0,BP,IGRID,4)
         CALL BUBVBW(KK,JJ,II,BU,BV,BW,BP)
         WRITE(*,*) 'FILLING OF THE BLOCKING ARRAYS SUCCESSFUL'
      ELSEIF (CTYP(IBODY)(1:8).EQ.'USEBLOCK') THEN
         CALL GETQUAD(KK,JJ,II,0,BP,IGRID,4)
         CALL BUBVBW(KK,JJ,II,BU,BV,BW,BP)
         WRITE(*,*) 'FILLING OF THE BLOCKING ARRAYS SUCCESSFUL'
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   40 CONTINUE
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
        RETURN
        END
