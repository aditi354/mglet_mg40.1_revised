










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
      SUBROUTINE GSIT(KK,JJ,II,DP,RHS,OMG,RES,IGRID,
     $                 AW,AE,AN,AS,AT,AB,RAP,
     $               NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,DIVGMX)

      IMPLICIT NONE

************************************************************************
*     GAUSS-SEIDEL ALGORITHMUS                                         *
*     INPUT:  KK, JJ, II     - GITTERPUNKTZAHL IN N-RICHTUNG           *
*             RHS            - RECHTE SEITE                            *
*             DX, DY, DZ     - GITTERWEITE                             *
*             P              - DRUCKFELD                               *
*             RAP, AN, AS,                                             *
*             AW, AE, AT, AB - GITTERWERTE                             *
*     RETURN: P              - DRUCKFELD                               *
************************************************************************

      INTEGER KK,JJ,II,K,J,I,KDM,JDM,IDM,NBND,KRB, KSTART,IND,IGRID

      INTEGER NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,ITER

      REAL    DP(KK,JJ,II),RAP(KK,JJ,II),RHS(KK,JJ,II)

      REAL    DIVGMX,OMG,RES

      REAL    AW(KK,JJ,II), AE(KK,JJ,II), AS(KK,JJ,II)
      REAL    AN(KK,JJ,II), AB(KK,JJ,II), AT(KK,JJ,II)

      save iter

***** BERECHNUNG DES RESIDUUMS *****************************************

      iter = iter + 1

      if (iter .eq.6 ) iter = 1

      DIVGMX = 0.
      DO  I = 3, II-2
         DO  J = 3, JJ-2
            DO K = 3, KK-2

C***** RED BLACK ********************************************************
C            DO KRB = 1,2 
C               KSTART = 2 + KRB
C               DO K = KSTART, KK-2,2
C
C***** BERECHNET WIRD RESIDUUM/AP ! *************************************

              RES=
     +            (AW(K,J,I)*DP(K,J,I-1) +
     +             AE(K,J,I)*DP(K,J,I+1) +
     +             AS(K,J,I)*DP(K,J-1,I) +
     +             AN(K,J,I)*DP(K,J+1,I) +
     +             AB(K,J,I)*DP(K-1,J,I) +
     +             AT(K,J,I)*DP(K+1,J,I) 
     +            - RHS(K,J,I))*RAP(K,J,I) + DP(K,J,I)
c                IF (ABS(RES) .GT. DIVGMX) THEN
c                     KDM = K
c                     JDM = J
c                     IDM = I
c                   DIVGMX = ABS(RES)
c                ENDIF

                  DP(K,J,I) = DP(K,J,I) - OMG * RES
                  DIVGMX = MAX(DIVGMX,ABS(RES))
               ENDDO
            ENDDO
         ENDDO

c         write(6,*)'DIVGMX',IGRID,ITER,DIVGMX

****** "EIGENTLICHER L\366SUNGSWEG" **************************************
c      DO I = 3, II-2
c         DO J = 3, JJ-2
c            DO K = 3, KK-2
c               DP(K,J,I) =OMG*((RHS(K,J,I)
c     +              -AW(I)*DP(K,J,I-1)-AE(I)*DP(K,J,I+1)
c     +              -AS(J)*DP(K,J-1,I)-AN(J)*DP(K,J+1,I)
c     +              -AB(K)*DP(K-1,J,I)-AT(K)*DP(K+1,J,I))*RAP(K,J,I))
c     +              +(1.0-OMG)*DP(K,J,I)
c            ENDDO
c         ENDDO
c      ENDDO                     ! RED BLACK



C                                 PERIODISCHE RANDBED. IN Y-RI.
C                                 -----------------------------
      IF (NRGT .EQ. 1) THEN

         DO I = 2, II-1
            DO K = 2, KK-1

               DP(K,2,I)=DP(K,JJ-2,I)
               DP(K,JJ-1,I)=DP(K,3,I)

            ENDDO
         ENDDO

      ENDIF



      
C                                 PERIODISCHE RANDBED. IN X-RI.
C                                 -----------------------------
C
      IF (NFRO.EQ.1) THEN

         DO J = 2, JJ-1
            DO K = 2, KK-1

               DP(K,J,2)=DP(K,J,II-2)
               DP(K,J,II-1)=DP(K,J,3)
            ENDDO
         ENDDO

      ENDIF

      RETURN
      END 
