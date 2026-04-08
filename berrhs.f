










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
      SUBROUTINE BERRHS(KK,JJ,II,HILF,RHS,RHSRES,
     $                   RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                     AW,AE AN,AS,AT,AB,
     $                      AP,LW,LS,
     $                      LB,LPR,UE,UN,UT)


      INTEGER II,JJ,KK
      REAL HILF(KK,JJ,II),RHS(KK,JJ,II),RHSRES(KK,JJ,II)
      REAL AW(KK,JJ,II),AE(KK,JJ,II),AN(KK,JJ,II),
     $     AS(KK,JJ,II),AT(KK,JJ,II),AB(KK,JJ,II)
      REAL RDX(II),RDY(JJ),RDZ(KK),RDDX(II),RDDY(JJ),RDDZ(KK)
      REAL  LW(KK,JJ,II),LS(KK,JJ,II),LB(KK,JJ,II),AP(KK,JJ,II),
     $      LPR(KK,JJ,II),UE(KK,JJ,II),UN(KK,JJ,II),UT(KK,JJ,II) 

      
      DO I=3,II-2                         
         DO J=3,JJ-2                      
            DO K=3,KK-2
c        RHSRES(K,J,I) = RHS(K,J,I) -  LW(K,J,I) * HILF(K,J,I-1)
c     &                             -  UE(K,J,I) * HILF(K,J,I+1)
c     &                             -  LS(K,J,I) * HILF(K,J-1,I)
c     &                             -  UN(K,J,I) * HILF(K,J+1,I)
c     &                             -  LB(K,J,I) * HILF(K-1,J,I)
c     &                             -  UT(K,J,I) * HILF(K+1,J,I)
c     &                             -  HILF(K,J,I)/ LPR(K,J,I)


c          RHSRES(K,J,I) = RHS(K,J,I) -  AW(I) *  HILF(K,J,I-1)
c     &                               -  AE(I) *  HILF(K,J,I+1)
c     &                               -  AS(J) *  HILF(K,J-1,I)
c     &                               -  AN(J) *  HILF(K,J-1,I)
c     &                               -  AB(K) *  HILF(K-1,J,I)
c     &                               -  AT(K) *  HILF(K+1,J,I)
c     &                               -  AP(K,J,I)* HILF(K,J,I)
           RHSRES(K,J,I) = RHS(K,J,I)- 2*HILF(K,J,I)*(
     &                  RDX(i)*RDDX(I)+RDY(J)*RDDY(J)+RDZ(K)*RDDZ(K))



            ENDDO
         ENDDO
      ENDDO




      RETURN
      END
