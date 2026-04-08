










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
      SUBROUTINE SETCOLEVEL (IGRID,LEVEL,IVP,ITST,ISCAI,
     $                       ISLCED,ISLICE,IPSDIR)


C*MGLET***************************************************************
C    S E T C O L E V E L     SETZEN DES COMMON-BLOCKES COLEVEL
C
C*MGLET***************************************************************
C
C      MAXLEVEL,   MINLEVEL,     : MAXIMALE UND MINIMALE GITTERLEVEL,
C                                  DIE BELEGT SIND
C      NOFLEVEL,  IGRDOFLEVEL,   : ANZAHL UND NUMMERN DER GITTER
C                                  PRO LEVEL UEBERHAUPT
C      IVP   BELEGT:
C      NOFVPIT ,  IGRDOFVPIT,    : ANZAHL UND NUMMERN DER GITTER
C                                  PRO LEVEL, DIE DRUCKITERATIONEN
C                                  DURCHLAUFEN
C      ITST  BELEGT:
C      NOFTST  ,  IGRDOFTST ,    : ANZAHL UND NUMMERN DER GITTER
C                                  PRO LEVEL, AUF DENEN ZEITSCHRITT
C                                  GEMACHT WIRD
C      ISLCED   BELEGT:
C      NOFSLCED,  IGRDOFSLCED    : ANZAHL UND NUMMERN DER GITTER
C                                  PRO LEVEL, DIE DURCH GEBIETSZERLEGUNG
C                                  ZERSCHNITTEN WERDEN
C      ISLICE   BELEGT:
C      NOFSLICE,  IGRDOFSLICE    : ANZAHL UND NUMMERN DER GITTER
C                                  PRO LEVEL, DIE DURCH GEBIETSZERLEGUNG
C                                  ENTSTANDEN SIND 
C                                  (SUBGITTER VON IGRDOFSLCED)
C
C      IPSDIR  BELEGT:
C      NOFPSDIR,  IGRDOFPSDIR    :ANZAHL UND NUMMERN DER GITTER
C                                  PRO LEVEL, IN DENEN DER DRUCK DURCH
C                                  DIREKTE LOESUNG DER POISSONGL. 
C                                  BERECHNET WERDEN SOLL
C                                  NUR MOEGLICH, FALLS PERIODISCH UND 
C                                  AEQUIDISTANT IN X- UND Y- RICHTUNG

C
C VERS:   6.10.93 (MM)  : ORIGINAL
C        18. 1.94 (MM)  : ERWEITERUNG
C
C
C*STARLET***************************************************************
C

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
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
 
      MAXLEVEL = MAX(MAXLEVEL,LEVEL)
      MINLEVEL = MIN(MINLEVEL,LEVEL)

      IF (MAXLEVEL .GT. MAXLVLMAX) CALL ERRR (501,' SETCOLEVEL ')
      IF (MINLEVEL .LT. MINLVLMIN) CALL ERRR (502,' SETCOLEVEL ')


             NOFLEVEL ( LEVEL ) =  NOFLEVEL ( LEVEL ) + 1

          IGRDOFLEVEL ( NOFLEVEL ( LEVEL ) , LEVEL ) = IGRID

      IF (IVP .EQ. 1) THEN

             NOFVPIT ( LEVEL ) =  NOFVPIT ( LEVEL ) + 1

          IGRDOFVPIT ( NOFVPIT ( LEVEL ) , LEVEL ) = IGRID

      ENDIF

      IF (ITST .EQ. 1) THEN

             NOFTST ( LEVEL ) =  NOFTST ( LEVEL ) + 1

          IGRDOFTST ( NOFTST ( LEVEL ) , LEVEL ) = IGRID

      ENDIF

      IF (ISCAI .EQ. 1) THEN
      
             NOFSCAI ( LEVEL ) =  NOFSCAI ( LEVEL ) + 1

          IGRDOFSCAI ( NOFSCAI ( LEVEL ) , LEVEL ) = IGRID

      ENDIF

      IF (ISLCED .EQ. 1) THEN

             NOFSLCED ( LEVEL ) =  NOFSLICE ( LEVEL ) + 1

          IGRDOFSLCED ( NOFSLCED ( LEVEL ) , LEVEL ) = IGRID

      ENDIF

      IF (ISLICE .EQ. 1) THEN

             NOFSLICE ( LEVEL ) =  NOFSLICE ( LEVEL ) + 1

          IGRDOFSLICE ( NOFSLICE ( LEVEL ) , LEVEL ) = IGRID

      ENDIF

      IF (IPSDIR .EQ. 1) THEN

             NOFPSDIR ( LEVEL ) =  NOFPSDIR ( LEVEL ) + 1

          IGRDOFPSDIR ( NOFPSDIR ( LEVEL ) , LEVEL ) = IGRID

      ENDIF

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C

      RETURN
      END
