










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

      SUBROUTINE THOMASI(KK,JJ,II,KL,KU,JL,JU,IL,IU,LCOL,DIAG,
     &                   RCOL,RSGS,FAKTOR,FOUTP
     $                   ,BP,GI,GJ,GK
     $                   )
C*MGLET***************************************************************
C        T H O M A S I   DAS UNTERPROGRAMM THOMAS LOEST EIN TRIDIAGONA-
C                        LES LINEARES GLEICHUNGSSYSTEM NACH DEM
C                        THOMAS - ALGORITHMUS
C*MGLET***************************************************************
C  PARAM:  IMX         - MAXIMALE ANZAHL DER GITTERPUNKTE EINSCHLIESS-
C                        LICH DER RANDSCHICHTEN
C          LCOL        - ENTHAELT KOEFFIZIENTEN DER "LINKEN SPALTE"
C                        PARALLEL ZUR HAUPTDIAGONALEN
C          DIAG        - ENTHAELT KOEFFIZIENTEN DER HAUPTDIAGONALEN
C          RCOL        - ENTHAELT KOEFFIZIENTEN DER "RECHTEN SPALTE"
C                        PARALLEL ZUR HAUPTDIAGONALEN
C          RSGS        - ENTHAELT  WERTE DES KONSTANTEN VEKTORS
C                        AUF DER RECHTEN SEITE DES GL.SYSTEMS
C          FOUTP       - ENTHAELT  DIE LOESUNG DES GL.SYSTEMS   
C
C  DEFINE DIREKTIVEN   : KEINE
C
C  UPROG               : KEINE
C
C     15.10.96 (AM)    : ORIGINAL
C     22.10.03 (NP)    : FRED BODY DEFINE DIREKTIVE, Parallelisierung 
C
C
C Definition
C      GI(K,J,I) = LCOL(I) = COEFFX(I,1)
C      GJ(K,J,I) = DIAG(I) = COEFFX(I,2)
C      GK(K,J,I) = RCOL(I) = COEFFX(I,3)
C
C*MGLET***************************************************************
C
      IMPLICIT NONE

      INTEGER KK,JJ,II,K,J,I,KL,KU,JL,JU,IL,IU,LP,N,IERR

      REAL     LCOL(II), DIAG(II), RCOL(II),RSGS(KK,JJ,II),
     $         FOUTP(KK,JJ,II),
     $         FAKTOR(KK,JJ,II),BP(KK,JJ,II),
     $         GI(KK,JJ,II),GJ(KK,JJ,II),GK(KK,JJ,II)

C
C                                 BILDUNG DER OBEREN DREIECKSMATRIX
C

 
      LP = IL+1


      DO K = KL,KU
         DO J = JL , JU
            DO I = LP,IU
               FAKTOR(K,J,I) = GI(K,J,I) / GJ(K,J,I-1)
               GJ(K,J,I)     = GJ(K,J,I) - (FAKTOR(K,J,I) * GK(K,J,I-1))
            END DO
         END DO
      END DO


      DO J = JL,JU
      DO K = KL,KU
         FOUTP(K,J,IL) = RSGS(K,J,IL)
      END DO
      END DO


       DO K = KL,KU
            DO J = JL , JU
               DO I = LP,IU
         FOUTP(K,J,I) = RSGS(K,J,I) - FAKTOR(K,J,I) * FOUTP(K,J,I-1) 
               END DO
            END DO
       END DO

C
C                                 RUECKSUBSTITUTION
C
      DO 500  J = JL,JU
      DO 400  K = KL,KU
         FOUTP(K,J,IU) = FOUTP(K,J,IU)/GJ(K,J,IU)
  400 CONTINUE
  500 CONTINUE


       DO K = KL,KU
         DO J = JL , JU
            DO N = LP,IU
               I = IU-N+IL
      FOUTP(K,J,I) = (FOUTP(K,J,I) - GK(K,J,I) *
     $                FOUTP(K,J,I+1)) / GJ(K,J,I)   
            END DO
         END DO
      END DO


c      write(*,*) KK,JJ,II,K,J,I,KL,KU,JL,JU,IL,IU,LP,N

CC    WRITE (44,*) 'THOMAS ERFOLGREICH DURCHLAUFEN!!'


      RETURN
      END
