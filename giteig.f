










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
      subroutine GITEIG (KK,JJ,II, DX,DY,DZ,
     $                    AW,AE,AN,AS,AT,AB,
     $                    BP,
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,AP
     &                        )

**********************************************************************
*     Gittereigenschaft                                              *
*     INPUT : KK, JJ, II     - Gitterpunktzahl in n-Richtung         *
*             dx, dy, dz     - Gitterweite                           *
*             NFRO,NBAC,NRGT,NLFT,NBOT,NTOP  - RANDBEDINGUNGEN       *
*                                              (SIEHE SETCOBOUND)    *
*     RETURN: RAP, AW, AE, AP                                        *
*             AN, AS, AT, AB - Gitterwerte                           *
*                                                                    *
**********************************************************************

C------------------------------- RANDBEDINGUNGEN
C        NEUMANN FUER DRUCK BEI:
C        FIXED (2), NOSLIP (5), SLIP (6), PARENT (8), BLO (9),
C        FLU (11), REC (12)
C
C        DIRICHLET FUER DRUCK BEI:
C        PERIODISCH (1), OPEN WALL (OP1) (3), OP2 (4), CONNECT (7),
C        KONVECTIVE (10), 


      IMPLICIT NONE

      INTEGER KK, JJ, II, K, J, I, NFRO, NBAC, NRGT, NLFT, NBOT, NTOP

      REAL    AP(KK,JJ,II)


      REAL    BP(KK,JJ,II)

      REAL    AW(II), AE(II), AS(JJ), AN(JJ), AB(KK), AT(KK)

      REAL    DX(II), DY(JJ), DZ(KK)


      DO I = 3, II-2
         AE(I) = 2./((DX(I-1)+DX(I))*DX(I))
         AW(I) = 2./((DX(I-1)+DX(I))*DX(I-1  ))
      ENDDO
      DO J = 3, JJ-2
         AN(J) = 2./((DY(J-1)+DY(J))*DY(J))
         AS(J) = 2./((DY(J-1)+DY(J))*DY(J-1  ))
      ENDDO
      DO K = 3, KK-2
         AT(K) = 2./((DZ(K-1)+DZ(K))*DZ(K))
         AB(K) = 2./((DZ(K-1)+DZ(K))*DZ(K-1  ))
      ENDDO

      DO I = 3, II-2
         DO J = 3, JJ-2
            DO K = 3, KK-2
               AP(K,J,I)=-2./(DX(I-1)*DX(I))
     +                   -2./(DY(J-1)*DY(J))
     +                   -2./(DZ(K-1)*DZ(K))

c               AP(K,J,I)    = AP(K,J,I) + AW(I)*(1.0-BP(K,J,I-1))
c               AP(K,J,I)    = AP(K,J,I) + AS(J)*(1.0-BP(K,J-1,I))
c               AP(K,J,I)    = AP(K,J,I) + AB(K)*(1.0-BP(K-1,J,I))
c               AP(K,J,I)    = AP(K,J,I) + AE(I)*(1.0-BP(K,J,I+1))
c               AP(K,J,I)    = AP(K,J,I) + AN(J)*(1.0-BP(K,J+1,I))
c               AP(K,J,I)    = AP(K,J,I) + AT(K)*(1.0-BP(K+1,J,I))

               AP(K,J,I)    = AP(K,J,I) + AW(I)*(1.0-BP(K,J,I-1))
     $                                  + AS(J)*(1.0-BP(K,J-1,I))
     $                                  + AB(K)*(1.0-BP(K-1,J,I))
     $                                  + AE(I)*(1.0-BP(K,J,I+1))
     $                                  + AN(J)*(1.0-BP(K,J+1,I))
     $                                  + AT(K)*(1.0-BP(K+1,J,I))

            ENDDO
         ENDDO
      ENDDO

*********NOSLIP Randbedingung an BOTTOM-FLAECHE

      IF (NBOT .EQ. 2 .OR. NBOT .EQ. 5 .OR. NBOT .EQ. 6 
     $   .OR. NBOT .EQ. 8  .OR. NBOT .EQ. 9 .OR. 
     $   NBOT .EQ. 11 .OR. NBOT .EQ. 12 )         THEN

         DO I = 3, II-2
            DO J =3, JJ-2
               AP(3,J,I)    = AP(3,J,I)    + AB(3) * BP(2,J,I)
            ENDDO
         ENDDO

         AB(3) = 0.0

      ENDIF 

*********NOSLIP Randbedingung an TOP-FLAECHE

      IF (NTOP .EQ. 2 .OR. NTOP .EQ. 5 .OR. NTOP .EQ. 6 
     $   .OR. NTOP .EQ. 8  .OR. NTOP .EQ. 9 .OR. 
     $   NTOP .EQ. 11 .OR. NTOP .EQ. 12 )         THEN

         DO I = 3, II-2
            DO J =3, JJ-2
               AP(KK-2,J,I) = AP(KK-2,J,I) + AT(KK-2) * BP(KK-1,J,I)
            ENDDO
         ENDDO

         AT(KK-2) = 0.0

      ENDIF 

*********NOSLIP Randbedingung an RIGHT-FLAECHE

      IF (NRGT .EQ. 2 .OR. NRGT .EQ. 5 .OR. NRGT .EQ. 6 
     $   .OR. NRGT .EQ. 8  .OR. NRGT .EQ. 9 .OR.
     $     NRGT .EQ. 11 .OR. NRGT .EQ. 12)          THEN

         DO I = 3, II-2
            DO K =3, KK-2

               AP(K,3,I)    = AP(K,3,I)    + AS(3) * BP(K,2,I)

            ENDDO
         ENDDO

         AS(3)    = 0.

      ENDIF
*********NOSLIP Randbedingung an LEFT-FLAECHE
      IF (NLFT .EQ. 2 .OR. NLFT .EQ. 5 .OR. NLFT .EQ. 6 
     $   .OR. NLFT .EQ. 8  .OR. NLFT .EQ. 9 .OR. 
     $   NLFT .EQ. 11 .OR. NLFT .EQ. 12 ) THEN
         DO I = 3, II-2
            DO K =3, KK-2
               AP(K,JJ-2,I) = AP(K,JJ-2,I) + AN(JJ-2) * BP(K,JJ-1,I)
            ENDDO
         ENDDO

         AN(JJ-2) = 0.

      ENDIF
*********NOSLIP Randbedingung an FRONT-FLAECHE
      IF (NFRO .EQ. 2 .OR. NFRO .EQ. 5 .OR. NFRO .EQ. 6 
     $   .OR. NFRO .EQ. 8  .OR. NFRO .EQ. 9 .OR. 
     $   NFRO .EQ. 11 .OR. NFRO .EQ. 12 ) THEN

         DO J = 3, JJ-2
            DO K =3, KK-2
               AP(K,J,3)    = AP(K,J,3)   + AW(3) * BP(K,J,2)
            ENDDO
         ENDDO

         AW(3)    = 0.

      ENDIF
*********NOSLIP Randbedingung an BAC-FLAECHE
      IF (NBAC .EQ. 2 .OR. NBAC .EQ. 5 .OR. NBAC .EQ. 6 
     $   .OR. NBAC .EQ. 8  .OR. NBAC .EQ. 9 .OR. 
     $   NBAC .EQ. 11 .OR. NBAC .EQ. 12 ) THEN

         DO J = 3, JJ-2
            DO K =3, KK-2
               AP(K,J,II-2) = AP(K,J,II-2) + AE(II-2) * BP(K,J,II-1)
            ENDDO
         ENDDO

         AE(II-2) = 0.

      ENDIF

      RETURN
      END
