










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
      SUBROUTINE COP3DZERO (KK,JJ,II,UTMP,VTMP,WTMP,U,V,W,P,BP
     $        )

C--MGLET----------------------------------------------------------------
C
C                         KOPIEREN DES FELDES PHI-TMP AUF PHI
C                         UND NULLSETZEN FÜR FREDBODY
C
C        13.03.03 (SE)  : ORIGINAL AUS COP3D UND SETZERO ABGELEITET
C
C--MGLET----------------------------------------------------------------

      IMPLICIT NONE

      INTEGER K,J,I,KK,JJ,II,KM1,JM1,IM1

      REAL  UTMP(KK,JJ,II),VTMP(KK,JJ,II),WTMP(KK,JJ,II),
     &      U   (KK,JJ,II),V   (KK,JJ,II),W   (KK,JJ,II),
     &      P   (KK,JJ,II),BP  (KK,JJ,II)


      IM1 = II - 1
      JM1 = JJ - 1
      KM1 = KK - 1

      DO I = 1,IM1
         DO J = 1,JM1
            DO K = 1,KM1
               P(K,J,I) = P   (K,J,I)*BP(K,J,I)
               U(K,J,I) = UTMP(K,J,I)*BP(K,J,I)*BP(K,J,I+1)
               V(K,J,I) = VTMP(K,J,I)*BP(K,J,I)*BP(K,J+1,I)
               W(K,J,I) = WTMP(K,J,I)*BP(K,J,I)*BP(K+1,J,I)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END


