










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
      SUBROUTINE THOMASJ (KK,JJ,II,KL,KU,JL,JU,IL,IU,LCOL,DIAG,
     $                    RCOL,RSGS,FAKTOR,FOUTP
     $                    ,BP,GI,GJ,GK
     $                    )
C*MGLET***************************************************************
C        T H O M A S     DAS UNTERPROGRAMM THOMAS LOEST EIN TRIDIAGONA-
C                        LES LINEARES GLEICHUNGSSYSTEM NACH DEM
C                        THOMAS - ALGORITHMUS
C*MGLET***************************************************************
C  PARAM:  LCOL        - ENTHAELT KOEFFIZIENTEN DER "LINKEN SPALTE"
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
C     30.06.03 (FS)    : 3-DIMENSIONAL
C
C*MGLET***************************************************************

      IMPLICIT NONE

      INTEGER  KK,JJ,II,K,J,I,KL,KU,JL,JU,IL,IU,LP,N

      REAL     LCOL(JJ), DIAG(JJ), RCOL(JJ),RSGS(KK,JJ,II),            
     $         FOUTP(KK,JJ,II),
     $         FAKTOR(KK,JJ,II),BP(KK,JJ,II),
     $         GI(KK,JJ,II),GJ(KK,JJ,II),GK(KK,JJ,II)

C
C                                 BILDUNG DER OBEREN DREIECKSMATRIX
C

      LP = JL+1

      DO I = IL,IU
         DO J = LP,JU
            DO K = KL,KU
               FAKTOR(K,J,I) = GI(K,J,I) / GJ(K,J-1,I)
               GJ(K,J,I)     = GJ(K,J,I) - (FAKTOR(K,J-1,I) * GK(K,J,I))
            END DO
         END DO
      END DO

      DO I = IL,IU
         DO K = KL,KU
            FOUTP(K,JL,I) = RSGS(K,JL,I)
         ENDDO
      ENDDO

      DO I = IL,IU
         DO   J = JL+1,JU
            DO   K = KL,KU
            FOUTP(K,J,I) = RSGS(K,J,I) - FAKTOR(K,J,I) * FOUTP(K,J-1,I)
            ENDDO
         ENDDO
      ENDDO

C                                 RUECKSUBSTITUTION

      DO I = IL,IU
         DO K = KL,KU
            FOUTP(K,JU,I) = FOUTP(K,JU,I)/GJ(K,JU,I)
         ENDDO
      ENDDO

      DO I = IL,IU
         DO  N = LP,JU
            J = JU-N+JL
            DO K = KL,KU
               FOUTP(K,J,I) = (FOUTP(K,J,I) - GK(K,J,I) *
     $                         FOUTP(K,J+1,I)) / GJ(K,J,I)

            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
