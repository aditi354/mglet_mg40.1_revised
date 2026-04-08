










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
        SUBROUTINE SIPITER1(KK,JJ,II,PHI,RHS,RES,
     $                 AW,AE,AN,AS,AT,AB,AP,
     $                 BP,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,
     $                 LW,LS,LB,LPR)

C****************************************************************
C        S I P I T E R    Eine Iteration des ILU Algorithmus nach
C                         Stone (SIP)
C****************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C
C UPROG                 : KEINE
C
C DEFINE-DIREKTIVEN     : FRED_BODY
C
C        25.02.03 (GT)  : ORIGINAL
C        31.10.03 (NP)  : FRED BODY define directive
C
C****************************************************************

      IMPLICIT NONE
      
      INTEGER KK,JJ,II,K,J,I,MAXK,MAXJ,MAXI

      REAL PHI(KK,JJ,II),RHS(KK,JJ,II),RES(KK,JJ,II), AP(KK,JJ,II),
     $      LW(KK,JJ,II), LS(KK,JJ,II), LB(KK,JJ,II),LPR(KK,JJ,II)

      REAL  BP(KK,JJ,II)

      REAL AW(II),AE(II),AN(JJ),AS(JJ),AT(KK),AB(KK)

      INTEGER NFRO,NBAC,NRGT,NLFT,NBOT,NTOP
      INTEGER NPROC,ISTART,ISTOP,M8,iproc

C      write(*,*) KK,JJ,II
C      CALL WRITE3DDIAGY (KK,JJ,II,RHS,99,3)
C      stop

      DO I=3,II-2
         DO J=3,JJ-2
            DO K=3,KK-2


               RES(K,J,I)=( RHS(K,J,I)
     $              - AW(I) * PHI(K,J,I-1) * BP(K,J,I-1)
     $              - AE(I) * PHI(K,J,I+1) * BP(K,J,I+1)
     $              - AS(J) * PHI(K,J-1,I) * BP(K,J-1,I)
     $              - AN(J) * PHI(K,J+1,I) * BP(K,J+1,I)
     $              - AB(K) * PHI(K-1,J,I) * BP(K-1,J,I)
     $              - AT(K) * PHI(K+1,J,I) * BP(K+1,J,I)
     $              - AP(K,J,I) * PHI(K,J,I))*  LPR(K,J,I)



               RES(K,J,I) =    RES(K,J,I)
     $              -LW(K,J,I)*RES(K,J,I-1)
     $              -LS(K,J,I)*RES(K,J-1,I)
     $              -LB(K,J,I)*RES(K-1,J,I)
            ENDDO
         ENDDO
      ENDDO

C     PERIODISCHE RANDBED. IN Y-RI.
C     -----------------------------
C

      IF (NRGT.EQ.1) THEN

         DO I = 2,II-1
            DO K = 2,KK-1

            RES(K,2,I)=RES(K,JJ-2,I)
            RES(K,JJ-1,I)=RES(K,3,I)

            ENDDO
         ENDDO

      ENDIF

C     PERIODISCHE RANDBED. IN X-RI.
C     -----------------------------
C

      IF ( NFRO .EQ. 1 ) THEN

         DO J = 2, JJ-1
            DO K = 2, KK-1

               RES(K,J,2)=RES(K,J,II-2)
               RES(K,J,II-1)=RES(K,J,3)

            ENDDO
         ENDDO
      ENDIF

      END


        SUBROUTINE SIPITER2(KK,JJ,II,PHI,RES,
     $                 NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,RESMAX,
     $                 UE,UN,UT)

C****************************************************************
C        S I P I T E R    Eine Iteration des ILU Algorithmus nach
C                         Stone (SIP)
C****************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C
C UPROG                 : KEINE
C
C DEFINE-DIREKTIVEN     : FRED_BODY
C
C        25.02.03 (GT)  : ORIGINAL
C
C****************************************************************

      IMPLICIT NONE
      
      INTEGER KK,JJ,II,K,J,I

      REAL PHI(KK,JJ,II),RES(KK,JJ,II),UE(KK,JJ,II), 
     $     UN (KK,JJ,II),UT (KK,JJ,II)

      REAL RESMAX

      INTEGER NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,RB
      INTEGER NPROC,ISTART,ISTOP,M8,iproc


C
C-----CALCULATE CORRECTION AND UPDATE VARIABLES (BACKWARD SUBST.)
C

      RESMAX = 0.0

C           CALL WRITE3DDIAGY (KK,JJ,II,RES,97,3)

      DO I=II-2,3,-1
         DO J=JJ-2,3,-1
            DO K=KK-2,3,-1

               RES(K,J,I)= RES(K,J,I)
     $              -UN(K,J,I)*RES(K,J+1,I)
     $              -UT(K,J,I)*RES(K+1,J,I)
     $              -UE(K,J,I)*RES(K,J,I+1)


               RESMAX = MAX(RESMAX,ABS(RES(K,J,I)))

               PHI(K,J,I)=PHI(K,J,I)+RES(K,J,I)

            ENDDO
         ENDDO
      ENDDO

C           CALL WRITE3DDIAGY (KK,JJ,II,UE,99,3)
C      write(*,*) RESMAX
C      stop
C                                 PERIODISCHE RANDBED. IN Y-RI.
C                                 -----------------------------
      IF (NRGT .EQ. 1) THEN

         DO I = 2, II-1
            DO k = 2, KK-1

               PHI(K,2,I)=PHI(K,JJ-2,I)
               PHI(K,JJ-1,I)=PHI(K,3,I)

            ENDDO
         ENDDO

      ENDIF



      
C                                 PERIODISCHE RANDBED. IN X-RI.
C                                 -----------------------------
C
      IF (NFRO.EQ.1) THEN

         DO J = 2, JJ-1
            DO K = 2, KK-1

               PHI(K,J,2)=PHI(K,J,II-2)
               PHI(K,J,II-1)=PHI(K,J,3)
            ENDDO
         ENDDO

      ENDIF

      RETURN
      END 
