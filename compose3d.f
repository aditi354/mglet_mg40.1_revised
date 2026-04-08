










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
       SUBROUTINE COMPOSE3D (KKP,JJP,IIP,PHIP,
     &               KKS,JJS,IIS,PHIS,
     &               IPOS,JPOS,KPOS,NBND,
     &               NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB)

C--MGLET----------------------------------------------------------------
C
C                  ZUSAMMENSETZEN EINES FELDES
C
C        PHIP:     PARENT
C        PHIS:     CHILD
C
C        IPOS,JPOS,KPOS:  INDIZES IM PARENT, AUF DENEN SUBGITTERPUNKT
C                         (K,J,I) = (3,3,3) LIEGT
C
C        15. 2.94 (MM)  : AUS SLICE3D ABGELEITET
C
C--MGLET----------------------------------------------------------------
C
C
      REAL
     $        PHIP(KKP,JJP,IIP),PHIS(KKS,JJS,IIS)
C
C-------------------------------------------------- CHECK AUF KONSISTENZ
C
      IF ( IIP .LT. IIS ) CALL ERRR (501,'SLICE1D')
      IF ( JJP .LT. JJS ) CALL ERRR (501,'SLICE1D')
      IF ( KKP .LT. KKS ) CALL ERRR (501,'SLICE1D')
C
C----------------------------------------------- SETZEN DER WERTE
C
C         IF BOUNDARYS OF GRIDS BELONGS TO THE CONNECT-CONDITION
C         THEN THE BOUNDARY VALUES ARE NOT PLACED (A.O. 29.03.1996)

      WRITE(6,*) 'compose3d:NBND',NBND

      ISTAK = 0
      ISTOK = 0
      JSTAK = 0
      JSTOK = 0
      KSTAK = 0
      KSTOK = 0

      IF (NFRO.EQ.7) ISTAK = 1
      IF (NBAC.EQ.7) ISTOK = 1
      IF (NRGT.EQ.7) JSTAK = 1
      IF (NLFT.EQ.7) JSTOK = 1
      IF (NBOT.EQ.7) KSTAK = 1
      IF (NTOP.EQ.7) KSTOK = 1

         ISTART = NBND        + ISTAK
         ISTOP  = IIS-NBND+1  - ISTOK
         JSTART = NBND        + JSTAK
         JSTOP  = JJS-NBND+1  - JSTOK
         KSTART = NBND        + KSTAK
         KSTOP  = KKS-NBND+1  - KSTOK
      write (6,*) 'compose3d, i-Richtung:',istart,istop,ipos
      write (6,*) 'compose3d, j-Richtung:',jstart,jstop,jpos
      write (6,*) 'compose3d, k-Richtung:',kstart,kstop,kpos

      DO IS =   ISTART,ISTOP
         IP = IS+IPOS-3
      DO JS =   JSTART,JSTOP
         JP = JS+JPOS-3
      DO KS =   KSTART,KSTOP
         KP = KS+KPOS-3

         PHIP(KP,JP,IP) = PHIS(KS,JS,IS)

      ENDDO
      ENDDO
      ENDDO

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
         WRITE(6,*)'COMPOSEFIELD, PHIS:',PHIS(3,3,3)
         WRITE(6,*)'COMPOSEFIELD, PHIP:',PHIP(3,3,3)

      RETURN
      END


