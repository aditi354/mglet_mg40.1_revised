










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
      SUBROUTINE WRITE_PASSPART (KK,JJ,II,CHANNEL,RUN,NFAM,
     $     XPART,DISSI,DISTURB,KINDEX,JINDEX,IINDEX,NFAM_MAX,WLOP,
     $     DS,TIMEPH,PARTSCA)
C----------------------------------------------------------
CCREADS POSITION AND ORIENTATION OF PARTICLES
CC1. 8. 1999 (MM): ORIGINAL
CC----------------------------------------------------------
      IMPLICIT NONE
      INTEGER CHANNEL,KK,JJ,II,NFAM,I,WLOP,NFAM_MAX,
     $     RUN(NFAM),KINDEX(NFAM),JINDEX(NFAM),
     $     IINDEX(NFAM)
      REAL     XPART(NFAM_MAX,3)
      REAL     DISSI(NFAM),DS(NFAM),TIMEPH,
     $     DISTURB(NFAM),PARTSCA(NFAM)
      
C      WRITE(6,*)'WRITE_PART:',KK,JJ,II,NFAM,NPART,CHANNEL
C      WRITE(6,*) XPART
      WLOP = WLOP + 1
      DO I=1,NFAM
C         WRITE (CHANNEL,1000,ERR=2000)I, RUN(I),(XPART(I,J),J=1,3),
C     $        DISSI(I)
C     WRITE (CHANNEL,  ERR=2000) IPART(I), (CONF(I,J),J=1,6)
      WRITE(CHANNEL,ERR=2000) I,RUN(I),XPART(I,1),XPART(I,2),XPART(I,3),
     $        DISSI(I),DISTURB(I),DS(I),TIMEPH
C      WRITE(6,*)'DISS TURB: ',I,DISTURB(I),DISSI(I)
      ENDDO
      RETURN
      
 1000 FORMAT (I5,I5, 6(1X,E12.5E3))
      
 2000 CALL ERRR (501,' WRITE_PART')
      
      END
