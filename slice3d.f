










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
       SUBROUTINE SLICE3D (KKP,JJP,IIP,PHIP,
     &               KKS,JJS,IIS,PHIS,
     &               IPOS,JPOS,KPOS)

C--MGLET----------------------------------------------------------------
C
C                  ZERSCHEIDEN EINES FELDES
C
C        PHIP:     PARENT
C        PHIS:     CHILD
C
C        IPOS,JPOS,KPOS:  INDIZES IM PARENT, AUF DENEN SUBGITTERPUNKT
C                         (K,J,I) = (3,3,3) LIEGT
C
C        15. 2.94 (MM)  : ORIGINAL
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
      DO IS = 1,IIS
         IP = IS+IPOS-3
      DO JS = 1,JJS
         JP = JS+JPOS-3
      DO KS = 1,KKS
         KP = KS+KPOS-3

         PHIS(KS,JS,IS) = PHIP(KP,JP,IP)

      ENDDO
      ENDDO
      ENDDO

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C

      RETURN
      END


