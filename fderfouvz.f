










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
      SUBROUTINE FDERFOUVZ(KK,JJ,II,KSTART,JSTART,ISTART,
     $     KSTOP,JSTOP,ISTOP,DZ,DDZ,COEFDZ,RSGS,FAKTOR,
     $     F,FIN,FD,NBOT,NTOP,LCOL,DIAG,RCOL
     $     ,BP,GI,GJ,GK
     $     )

C*MGLET***************************************************************
C        F D E R F O U V Z  (First DERive Fourth Order)
C        BERECHNUNG DER ERSTEN ABLEITUNG (DU/DZ,DV/DZ)                    
C*MGLET***************************************************************
C
C PARAM: F(K,J,I)       - BEKANNTE GROESSE AN DEN GITTERKANTEN
C      : FIN(K,J,I)     - INTERPOLIERTE GROESSE
C                         ZWISCHEN DEN BEIDEN KANTEN       
C      : FD (K,J,I)       - ERSTE ABLEITUNG, FIRST DERINATIVE:
C                         (DU/DZ,DV/DZ)
C      : LCOL           - ENTHAELT KOEFFIZIENTEN DER "LINKEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : DIAG           - ENTHAELT KOEFFIZIENTEN DER HAUPTDIAGONALEN
C      : RCOL           - ENTHAELT KOEFFIZIENTEN DER "RECHTEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C
C DEFINE DIREKTIVEN     : KEINE 
C
C UPROG                 : THOMASK     THOMAS-ALGORITHMUS IN Z-RICHTUNG
C
C VERS:  08.01.97(AM)  : ORIGINAL        (KOMPAKT 4. ORDNUNG)
C
C*MGLET***************************************************************

      IMPLICIT NONE

      INTEGER KK,JJ,II,K,J,I,KSTART,JSTART,ISTART,KSTOP,JSTOP,ISTOP,
     $        KA,KE,NBOT,NTOP

      REAL    DZ(KK),DDZ(KK),LCOL(KK),DIAG(KK),RCOL(KK),
     $        F(KK,JJ,II),FIN(KK,JJ,II),FD(KK,JJ,II),RSGS(KK,JJ,II),
     $        COEFDZ(KK,12*3),FAKTOR(KK,JJ,II),A,B,C,D,
     $        GI(KK,JJ,II),GJ(KK,JJ,II),GK(KK,JJ,II),BP(KK,JJ,II)

      IF (NBOT .EQ. 5 .OR. NBOT .EQ. 11) THEN 
         KA = 3
      ELSE
         WRITE(6,*) 'RANDBEDINGUNG IN FDERFOUVZ NOCH NICHT '
         STOP
      ENDIF

      IF (NTOP .EQ. 5) THEN 
         KE = KK-1
      ELSEIF (NTOP .EQ. 3) THEN 
         KE = KK-2
      ELSE
         WRITE(6,*) 'RANDBEDINGUNG IN FDERFOUVZ NOCH NICHT '
         STOP
      ENDIF

C                                      ERSTER PHYSIKALISCHER PUNKT

      K = KA
      DO I = ISTART,ISTOP
         DO J = JSTART, JSTOP
            GI(K,J,I) = COEFDZ(K,7)
            GJ(K,J,I) = COEFDZ(K,8)
            GK(K,J,I) = COEFDZ(K,9)
         ENDDO
      ENDDO

      IF(NBOT.EQ.5) THEN
         DO I = ISTART,ISTOP
            DO J = JSTART, JSTOP

               RSGS(K,J,I) =  COEFDZ(K,10) * FIN(K  ,J,I)
     $                      + COEFDZ(K,11) * F  (K  ,J,I)
     $                      + COEFDZ(K,12) * F  (K+1,J,I)
C               RSGS(K,J,I) = 2*(
C     $                      +  F  (K  ,J,I))*DZ(K)
            ENDDO
         ENDDO
      ENDIF
C                                     LETZTER PHYSIKALISCHER PUNKT
      K = KE
      DO I = ISTART,ISTOP
         DO J = JSTART,JSTOP
            GI(K,J,I) = COEFDZ(K,7)
            GJ(K,J,I) = COEFDZ(K,8)
            GK(K,J,I) = COEFDZ(K,9)
         ENDDO
      ENDDO


      IF(NTOP.EQ.5) THEN
         DO I = ISTART,ISTOP
            DO J = JSTART,JSTOP
               RSGS(K,J,I) =  COEFDZ(K,10) * FIN(K  ,J,I)
     $                      + COEFDZ(K,11) * F  (K-1,J,I)
     $                      + COEFDZ(K,12) * F  (K-2,J,I)
C               RSGS(K,J,I) =  2*(
C     $                      -  F  (K-1,J,I))*DZ(K)
            ENDDO
         ENDDO
      ELSEIF (NTOP.EQ.3) THEN
         DO I = ISTART,ISTOP
            DO J = JSTART,JSTOP
               RSGS(K,J,I) =  COEFDZ(K,10) * F(K  ,J,I)
     $                     +  COEFDZ(K,11) * F(K-1,J,I)
     $                     +  COEFDZ(K,12) * F(K-2,J,I)
            ENDDO
         ENDDO
      ENDIF
C                                      IM BERECHNUNGSGEBIET
C--------------Compact differentiation for inner domain----------


      DO I = ISTART,ISTOP
         DO J = JSTART,JSTOP
            DO K = KA+1 , KE-1
C          A =    BP(K,J,I)
C          B = (1-BP(K,J,I))*    BP(K+1,J,I)
C          C = (1-BP(K,J,I))*    BP(K-1,J,I)
C          D = (1-BP(K,J,I))* (1-BP(K-1,J,I))* (1-BP(K+1,J,I))
          A =    BP(K-1,J,I) *   BP(K,J,I)
          B = (1-BP(K-1,J,I))*   BP(K,J,I)
          C = (1-BP(K,J,I  ))*   BP(K-1,J,I)
          D = (1-BP(K,J,I  ))*(1-BP(K-1,J,I))


          GI(K,J,I) =  COEFDZ(K,7)    * A
     $               + COEFDZ(K,12+7) * B
     $               + COEFDZ(K,24+7) * C
          GJ(K,J,I) =  COEFDZ(K,8)    * A
     $               + COEFDZ(K,12+8) * B
     $               + COEFDZ(K,24+8) * C
     $               + 1.0            * D
          GK(K,J,I) =  COEFDZ(K,9)    * A
     $               + COEFDZ(K,12+9) * B
     $               + COEFDZ(K,24+9) * C

c                RSGS(K,J,I) =
c     $   (COEFDX(I,10)    * F(K,J,I)+
c     $    COEFDX(I,11)    * FIN(K,J,I)+
c     $    COEFDX(I,12)    * F(K,J,I-1))*DDX(I-1) * A
c     $ + (COEFDX(I,12+10) * FIN(K,J,I)+
c     $    COEFDX(I,12+11) * F(K,J,I)+
c     $    COEFDX(I,12+12) * F(K,J,I+1)) *        B
c     $ + (COEFDX(I,24+10) * FIN(K,J,I)+
c     $    COEFDX(I,24+11) * F(K,J,I-1)+
c     $    COEFDX(I,24+12) * F(K,J,I-2)) *        C

c bessere Divergenz
                RSGS(K,J,I) =
     $   (COEFDZ(K,10)    * F(K,J,I)+
     $    COEFDZ(K,11)    * FIN(K,J,I)+
     $    COEFDZ(K,12)    * F(K-1,J,I))*DZ(K-1) * A
     $ + (COEFDZ(K,12+10) * F(K-1,J,I)+
     $    COEFDZ(K,12+11) * F(K,J,I)+
     $    COEFDZ(K,12+12) * F(K+1,J,I)) *        B
     $ + (COEFDZ(K,24+10) * F(K,J,I)+
     $    COEFDZ(K,24+11) * F(K-1,J,I)+
     $    COEFDZ(K,24+12) * F(K-2,J,I)) *        C
            ENDDO
         ENDDO
      ENDDO

       CALL THOMASK(KK,JJ,II,KA,KE,JSTART,JSTOP,ISTART,ISTOP,
     $              LCOL,DIAG,RCOL,RSGS,FAKTOR,FD
     $     ,BP,GI,GJ,GK
     $     )

C*********************  RECHTER RAND   *******************************
C
       IF(NTOP .EQ. 3) THEN
          K = KE+1
          DO I = ISTART,ISTOP
             DO J = JSTART,JSTOP
                FD(K,J,I) = -  COEFDZ(K,7 ) * FD (K-1,J,I)
     $                      +  COEFDZ(K,10) * FIN(K  ,J,I)
     $                      +  COEFDZ(K,11) * F  (K-1,J,I)
     $                      +  COEFDZ(K,12) * F  (K-2,J,I)
             ENDDO
          ENDDO
       ENDIF

       RETURN
       END

