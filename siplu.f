










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
      SUBROUTINE SIPLU (KK,JJ,II,AW,AE,AN,AS,AT,AB,
     $                 BP,
     $                  AP,LW,LS,LB,LPR,UE,UN,UT)

C****************************************************************
C        S I P L U        SIPLU BELEGT DIE KOEFFIZIENTEN DES SIP
C                         SOLVERS
C****************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C
C UPROG                 : KEINE
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        25.02.03 (GT)  : ORIGINAL
C
C****************************************************************

      IMPLICIT NONE

      INTEGER KK,JJ,II,K,J,I

      REAL P1,P2,P3,ALFA

      REAL BP(KK,JJ,II)

      REAL AW(II),AE(II),AN(JJ),AS(JJ),AT(KK),AB(KK)

      REAL   AP(KK,JJ,II),LB(KK,JJ,II),LW(KK,JJ,II),LS(KK,JJ,II),
     $      LPR(KK,JJ,II),UE(KK,JJ,II),UN(KK,JJ,II),UT(KK,JJ,II)

      PARAMETER ( ALFA=0.92 )

      CALL DPHI0(KK,JJ,II,KK,JJ,II,LB)
      CALL DPHI0(KK,JJ,II,KK,JJ,II,LW)
      CALL DPHI0(KK,JJ,II,KK,JJ,II,LS)
      CALL DPHI0(KK,JJ,II,KK,JJ,II,UE)
      CALL DPHI0(KK,JJ,II,KK,JJ,II,UT)
      CALL DPHI0(KK,JJ,II,KK,JJ,II,LPR)

      DO I=3,II-2
         DO J=3,JJ-2
            DO K=3,KK-2

            LB(K,J,I) = AB(K) * BP(K-1,J,I)
     $                / (1.+ALFA*(UN(K-1,J,I) +UE(K-1,J,I)))

            LW(K,J,I) = AW(I) * BP(K,J,I-1)
     $                / (1.+ALFA*(UN(K,J,I-1) +UT(K,J,I-1)))

            LS(K,J,I) = AS(J) * BP(K,J-1,I)
     $                / (1.+ALFA*(UE(K,J-1,I) +UT(K,J-1,I)))

            P1=ALFA*(LB(K,J,I)*UN(K-1,J,I) + LW(K,J,I)*UN(K,J,I-1))
            P2=ALFA*(LB(K,J,I)*UE(K-1,J,I) + LS(K,J,I)*UE(K,J-1,I))
            P3=ALFA*(LW(K,J,I)*UT(K,J,I-1) + LS(K,J,I)*UT(K,J-1,I))

            LPR(K,J,I) = 1./(AP(K,J,I)+P1+P2+P3-LB(K,J,I)*UT(K-1,J,I)
     $                 - LW(K,J,I)*UE(K,J,I-1)
     $                 - LS(K,J,I)*UN(K,J-1,I)+1.E-20)

            UN(K,J,I)=(AN(J) * BP(K,J+1,I)-P1)*LPR(K,J,I)
            UE(K,J,I)=(AE(I) * BP(K,J,I+1)-P2)*LPR(K,J,I)
            UT(K,J,I)=(AT(K) * BP(K+1,J,I)-P3)*LPR(K,J,I)
            ENDDO
         ENDDO
      ENDDO
      DO I=3,II-2
         DO J=3,JJ-2
          DO K=3,KK-2
               LW(K,J,I)=LW(K,J,I)*LPR(K,J,I)
               LS(K,J,I)=LS(K,J,I)*LPR(K,J,I)
               LB(K,J,I)=LB(K,J,I)*LPR(K,J,I)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
