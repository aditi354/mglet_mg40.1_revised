










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
      SUBROUTINE FDERFOUWY(KK,JJ,II,KSTART,JSTART,ISTART,KSTOP,JSTOP,
     $                    ISTOP,DY,DDY,COEFDY,RSGS,FAKTOR,
     $                    F,FIN,FD,NRGT,NLFT,LCOL,DIAG,RCOL
     $     ,BP,GI,GJ,GK
     $     )

C*MGLET***************************************************************
C        F D E R F O U W Y  (First DERive Fourth Order)
C        BERECHNUNG DER ERSTEN ABLEITUNG (DU/DY,DW/DY)                    
C*MGLET***************************************************************
C
C PARAM: F(K,J)         - BEKANNTE GROESSE AN DEN GITTERKANTEN
C      : FIN(K,J)       - INTERPOLIERTE GROESSE IN DER GLEICHEN RICHTUNG
C                         ZWISCHEN DEN BEIDEN KANTEN       
C      : FD (K,J)       - ERSTE ABLEITUNG, FIRST DERINATIVE 
C                         (DU/DY, DW/DY)
C      : LCOL           - ENTHAELT KOEFFIZIENTEN DER "LINKEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : DIAG           - ENTHAELT KOEFFIZIENTEN DER HAUPTDIAGONALEN
C      : RCOL           - ENTHAELT KOEFFIZIENTEN DER "RECHTEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C
C DEFINE DIREKTIVEN     : KEINE 
C
C UPROG                 : THOMASJ     THOMAS-ALGORITHMUS IN Y-RICHTUNG
C
C VERS:  08.01.97(AM)   : ORIGINAL     (KOMPAKT 4. ORDNUNG)
C        20.03.03(FS)   : NUR FUER NRGT = NLEFT = 5 
C        01.07.03(FS)   : 3-Dimensional
C
C*MGLET***************************************************************

      IMPLICIT NONE

      INTEGER KK,JJ,II,K,J,I,KSTART,JSTART,ISTART,JSTOP,KSTOP,ISTOP,
     $        JA,JE,NRGT,NLFT

      REAL    DY(JJ),DDY(JJ),LCOL(JJ),DIAG(JJ),RCOL(JJ),
     $        F(KK,JJ,II),FIN(KK,JJ,II),FD(KK,JJ,II),RSGS(KK,JJ,II),
     $        COEFDY(JJ,12*3),FAKTOR(KK,JJ,II),A,B,C,D,
     $        GI(KK,JJ,II),GJ(KK,JJ,II),GK(KK,JJ,II),BP(KK,JJ,II)


      IF (NRGT.EQ.5)   THEN
         JA = 3
      ELSE
         WRITE(6,*) 'RANDBEDINGUNG IN FDERFOUWY NOCH NICHT '
         STOP
      ENDIF

      IF(NLFT.EQ.5) THEN
         JE = JJ-1
      ELSE
         WRITE(6,*) 'RANDBEDINGUNG IN FDERFOUWY NOCH NICHT '
         STOP
      ENDIF

      J = JA
      DO I = ISTART,ISTOP
         DO K = KSTART,KSTOP
            GI(K,J,I) = COEFDY(J,7)
            GJ(K,J,I) = COEFDY(J,8)
            GK(K,J,I) = COEFDY(J,9)
            RSGS(K,J,I) = COEFDY(J,10)*FIN(K,J,I) 
     $                  + COEFDY(J,11)*F(K,J,I)
     $                  + COEFDY(J,12)*F(K,J+1,I)
         ENDDO
      ENDDO
          
      J = JE
      DO I=ISTART,ISTOP
         DO K = KSTART,KSTOP      
            GI(K,J,I) = COEFDY(J,7)
            GJ(K,J,I) = COEFDY(J,8)
            GK(K,J,I) = COEFDY(J,9)
            RSGS(K,J,I) = COEFDY(J,10)*FIN(K,J,I) 
     $                  + COEFDY(J,11)*F(K,J-1,I)
     $                 +  COEFDY(J,12)*F(K,J-2,I)
         ENDDO
      ENDDO
       
      DO I = ISTART,ISTOP
         DO J = JA+1 , JE-1
            DO K = KSTART,KSTOP
C          A =    BP(K,J,I)
C          B = (1-BP(K,J,I))*    BP(K,J+1,I)
C          C = (1-BP(K,J,I))*    BP(K,J-1,I)
C          D = (1-BP(K,J,I))* (1-BP(K,J-1,I))* (1-BP(K,J+1,I))

          A =    BP(K,J-1,I) *   BP(K,J,I)
          B = (1-BP(K,J-1,I))*   BP(K,J,I)
          C = (1-BP(K,J,I  ))*   BP(K,J-1,I)
          D = (1-BP(K,J,I  ))*(1-BP(K,J-1,I))

          GI(K,J,I) =  COEFDY(J,7)    * A
     $               + COEFDY(J,12+7) * B
     $               + COEFDY(J,24+7) * C
          GJ(K,J,I) =  COEFDY(J,8)    * A
     $               + COEFDY(J,12+8) * B
     $               + COEFDY(J,24+8) * C
     $               + 1.0            * D
          GK(K,J,I) =  COEFDY(J,9)    * A
     $               + COEFDY(J,12+9) * B
     $               + COEFDY(J,24+9) * C

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
     $   (COEFDY(J,10)    * F(K,J,I)+
     $    COEFDY(J,11)    * FIN(K,J,I)+
     $    COEFDY(J,12)    * F(K,J-1,I))*DY(J-1) * A
     $ + (COEFDY(J,12+10) * F(K,J-1,I)+
     $    COEFDY(J,12+11) * F(K,J,I)+
     $    COEFDY(J,12+12) * F(K,J+1,I)) *        B
     $ + (COEFDY(J,24+10) * F(K,J,I)+
     $    COEFDY(J,24+11) * F(K,J-1,I)+
     $    COEFDY(J,24+12) * F(K,J-2,I)) *        C



            ENDDO
         ENDDO
      ENDDO

       CALL THOMASJ(KK,JJ,II,KSTART,KSTOP,JA,JE,ISTART,ISTOP,
     $      LCOL,DIAG,RCOL,RSGS,FAKTOR,FD
     $     ,BP,GI,GJ,GK
     $     )

C*********************   E  N  D       *******************************
       RETURN
       END

