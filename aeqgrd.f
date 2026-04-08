










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
      SUBROUTINE AEQGRD  (KMX,JMX,IMX,KK,JJ,II,DX,DY,DZ,NBND,
     &                    NXGRAE,NYGRAE,NZGRAE,IGRID)
C*STARLET***************************************************************
C        A E Q G R D      AEQGRD PRUEFT DAS MASCHENGITTER IN ALLEN
C                         DREI KOORDINATENRICHTUNGEN AUF AEQUIDISTANZ.
C*STARLET***************************************************************
C
C PARAM: KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        KK, JJ, II     - ARRAYDIMENSIONEN
C        DX,DY,DZ       - ABSTAND DER GITTERPUNKTE
C        NBND           - ANZAHL DER RANDSCHICHTEN DES BERECHNUNGS-
C                         GEBIETES
C
C UPROG                 : AEQGR1
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        25.09.88 (HW)  : ORIGINAL
C        14.06.89 (HW)  : EPSTOL VON 100.*SMAONE AUF 1000.*SMAONE
C                         WEGEN GRDFMI
C
C*STARLET***************************************************************
C
C
      REAL       DX(II),        DY(JJ),        DZ(KK)
C
C                                 UEBERPRUEFUNG DER DREI KOORDINATEN-
C                                 RICHTUNGEN
C
      WRITE (6,*)
      WRITE (6,*)
      WRITE (6,*) ' ********** INFORMATION AUS SUBR. AEQGRD **********'
      WRITE (6,*)
      WRITE (6,*) ' ********** GITTERNR.:     ',IGRID
      WRITE (6,*)
      CALL AEQGR1  (IMX,II,DX,NBND,NXGRAE,'X')
      CALL AEQGR1  (JMX,JJ,DY,NBND,NYGRAE,'Y')
      CALL AEQGR1  (KMX,KK,DZ,NBND,NZGRAE,'Z')
      WRITE (6,*)
      WRITE (6,*)
      RETURN
      END
      SUBROUTINE AEQGR1  (LMX,LL,DS,NBND,NLGRAE,KOORD)
C*STARLET***************************************************************
C        A E Q G R 1      AEQGR1 PRUEFT DAS MASCHENGITTER IN EINER
C                         KOORDINATENRICHTUNG AUF AEQUIDISTANZ.
C*STARLET***************************************************************
C
C PARAM: LMX            - GRENZE DES BERECHNUNGSGEBIETES (MIT BOUND)
C        LL             - ARRAYDIMENSION
C        DS             - ABSTAND DER GITTERPUNKTE
C        NBND           - ANZAHL DER RANDSCHICHTEN DES BERECHNUNGS-
C                         GEBIETES
C        NLGRAE         + NLGRAE = 1 : MASCHENGITTER IST AEQUIDISTANT
C                         NLGRAE = 0 : MASCHENGITTER IST  N I C H T
C                                      AEQUIDISTANT
C        KOORD          - CHARACTER (LEN=1) VARIABLE (KOORDINATENRICHTUNG)
C
C UPROG                 : ERRR
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        25.09.88 (HW)  : ORIGINAL
C
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
      CHARACTER (LEN=1)  KOORD
      REAL          DS(LL)
C
      NLGRAE = 1
      EPSTOL = SQRT(SMALL)
C
      LBEG   = 1   + NBND
      LEND   = LMX - NBND
C
      DO L = LBEG,LEND
       IF(ABS (DS(L) - DS(L-1)) .GT. EPSTOL) THEN
             NLGRAE = 0
C             WRITE(6,*)'AEQGRID:',DS(L),DS(L-1),EPSTOL
       ENDIF
      ENDDO
C
      IF(NLGRAE .EQ. 1) THEN
         WRITE (6,*) ' DAS GITTER IST IN ',KOORD,'-RICHTUNG ',
     $               'AEQUIDISTANT              -->  N',KOORD,'GRAE = 1'
         GOTO 9999
      END IF
      IF(NLGRAE .EQ. 0) THEN
         WRITE (6,*) ' DAS GITTER IST IN ',KOORD,'-RICHTUNG ',
     $               ' N I C H T  AEQUIDISTANT  -->  N',KOORD,'GRAE = 0'
         GOTO 9999
      END IF
      CALL ERRR (501,' AEQGR1   ')
 9999 RETURN
      END
