










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
       SUBROUTINE ICOGRDCON ()


C*MGLET*****************************************************************
C        I C O G R D C O N INITIALISIEREN DER COMMON-BLOCKES
C                          COGRDCON
C                          ZUM STEUERN DER GITTERVERNETZUNG
C*MGLET*****************************************************************
C
C PARAM: 
C       IVPCHILD(MAXGRIDS)    GITTER, AUS DEM EIN GROBGITTER FUER
C                             DRUCKKORREKTUR HERVORGEHT
C       IPARENT(MAXGRIDS)     PARENT-GITTER (EINE STUFE GROEBER)	
C       ISLPAR(MAXGRIDS)      PARENT-GITTER FUER GEBIETSZERLEGUNG
C                             (GLEICHE VERGROEBERUNGSSTUFE)
C      IPOSITION(MAXGRIDS)    I-POSITION, DES GITTERS IM PARENT-GITTER
C      JPOSITION(MAXGRIDS)    J-POSITION, DES GITTERS IM PARENT-GITTER
C      KPOSITION(MAXGRIDS)    K-POSITION, DES GITTERS IM PARENT-GITTER
C NOFSLCHILDS(MAXGRIDS),IGRDOFSLCHILD(MAXGRIDS,MAXGRIDS),
C    ISLPOS(MAXGRIDS), JSLPOS(MAXGRIDS), KSLPOS(MAXGRIDS),
C       IFRNBR(MAXBOCONDS,MAXGRIDS), IBANBR(MAXBOCONDS,MAXGRIDS),
C       IRINBR(MAXBOCONDS,MAXGRIDS), ILENBR(MAXBOCONDS,MAXGRIDS),
C       IBONBR(MAXBOCONDS,MAXGRIDS), ITONBR(MAXBOCONDS,MAXGRIDS)
C
C
C VERS:  15. 2.94 (MM)  : ORIGINAL
C
C
C*STARLET***************************************************************
C

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )

      
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
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  VORBELEGUNG 
C
      DO IGRID = 1,MAXGRIDS

         IVPCHILD(IGRID) = 0
         IPARENT (IGRID) = 0
         ISLPAR  (IGRID) = 0
         IPOSITION (IGRID) = 0
         JPOSITION (IGRID) = 0
         KPOSITION (IGRID) = 0
         NOFSLCHILDS (IGRID) = 0
         DO I=1,MAXGRIDS
             IGRDOFSLCHILD (I,IGRID) = 0
         ENDDO
         ISLPOS  (IGRID) = 0
         JSLPOS  (IGRID) = 0
         KSLPOS  (IGRID) = 0
         DO I=1,MAXBOCONDS
             IFRNBR  (I,IGRID) = -1
             IBANBR  (I,IGRID) = -1
             IRINBR  (I,IGRID) = -1
             ILENBR  (I,IGRID) = -1
             IBONBR  (I,IGRID) = -1
             ITONBR  (I,IGRID) = -1
         ENDDO

         IPARENT (IGRID) = 0
   
      ENDDO

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
C
      RETURN
      END


