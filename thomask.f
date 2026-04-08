










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
      SUBROUTINE THOMASK (KK,JJ,II,KL,KU,JL,JU,IL,IU,
     $                    LCOL,DIAG,RCOL,RSGS,FAKTOR,FOUTP
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
C
C*MGLET***************************************************************

      IMPLICIT NONE

      INTEGER  KK,JJ,II,K,J,I,KL,KU,JL,JU,IL,IU,LP,N

      REAL     LCOL(KK), DIAG(KK), RCOL(KK),     
     $         FOUTP(KK,JJ,II),RSGS(KK,JJ,II),
     $         FAKTOR(KK,JJ,II),BP(KK,JJ,II),
     $         GI(KK,JJ,II),GJ(KK,JJ,II),GK(KK,JJ,II)

C
C                                 BILDUNG DER OBEREN DREIECKSMATRIX
C

      LP = KL+1


      DO I = IL,IU
         DO J = JL,JU
            DO K = LP,KU
               FAKTOR(K,J,I)= GI(K,J,I) / GJ(K-1,J,I)
               GJ(K,J,I)    = GJ(K,J,I) - (FAKTOR(K,J,I) * GK(K-1,J,I))
            END DO
         END DO
      END DO


      DO I = IL,IU
         DO J = JL,JU
            FOUTP(KL,J,I) = RSGS(KL,J,I)
         ENDDO
      ENDDO

      DO I = IL,IU
         DO  J = JL,JU
            DO  K = KL+1,KU
             FOUTP(K,J,I)= RSGS(K,J,I) - FAKTOR(K,J,I) * FOUTP(K-1,J,I)
            ENDDO
         ENDDO
      ENDDO



C
C                                 RUECKSUBSTITUTION
C
      DO I = IL,IU
         DO  J = JL,JU
            FOUTP(KU,J,I) = FOUTP(KU,J,I)/GJ(KU,J,I)
         ENDDO
      ENDDO

      DO I = IL,IU
         DO  J = JL,JU
            DO   N = LP,KU
               K = KU-N+KL
            FOUTP(K,J,I) = (FOUTP(K,J,I) - GK(K,J,I) *
     $                      FOUTP(K+1,J,I)) / GJ(K,J,I)
            ENDDO
         ENDDO
      ENDDO

C     WRITE (44,*) 'THOMAS ERFOLGREICH DURCHLAUFEN!!'

      RETURN
      END
