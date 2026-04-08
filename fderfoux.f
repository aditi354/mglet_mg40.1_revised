










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
      SUBROUTINE FDERFOUX(KK,JJ,II,KSTART,JSTART,ISTART,
     $     KSTOP,JSTOP,ISTOP,DX,DDX,COEFDX,RSGS,FAKTOR,
     $     F,FIN,FD,NFRO,NBAC,LCOL,DIAG,RCOL
     $     ,BP,GI,GJ,GK
     $     )
C*MGLET***************************************************************
C        F D E R F O U X  (First DERive Fourth Order)
C        BERECHNUNG DER ERSTEN ABLEITUNG (DU/DX)                    
C*MGLET***************************************************************
C
C PARAM: F(K,J,I)       - BEKANNTE GROESSE AN DEN GITTERKANTEN
C      : FIN(K,J,I)     - INTERPOLIERTE GROESSE
C                         ZWISCHEN DEN BEIDEN KANTEN       
C      : FD (K,J,I)     - ERSTE ABLEITUNG, FIRST DERINATIVE  DU/DX
C
C      : LCOL           - ENTHAELT KOEFFIZIENTEN DER "LINKEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : DIAG           - ENTHAELT KOEFFIZIENTEN DER HAUPTDIAGONALEN
C      : RCOL           - ENTHAELT KOEFFIZIENTEN DER "RECHTEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C
C DEFINE DIREKTIVEN     : KEINE 
C
C UPROG                 :  THOMASI     THOMAS-ALGORITHMUS IN X-RICHTUNG
C
C VERS:  08.01.97(AM)   : ORIGINAL        (KOMPAKT 4. ORDNUNG)
C        26.01.04(NP)   : FRED BODY VERSION
C
C*MGLET***************************************************************

      IMPLICIT NONE

      INTEGER KK,JJ,II,K,J,I,KSTART,JSTART,ISTART,KSTOP,JSTOP,ISTOP,
     $        IA,IE,NFRO,NBAC

      REAL    DX(II),DDX(II),LCOL(II),DIAG(II),RCOL(II),
     $        F(KK,JJ,II),FIN(KK,JJ,II),FD(KK,JJ,II),RSGS(KK,JJ,II)

      REAL    COEFDX(II,12*3),FAKTOR(KK,JJ,II),A,B,C,D,
     $        GI(KK,JJ,II),GJ(KK,JJ,II),GK(KK,JJ,II),BP(KK,JJ,II)


      IF (NFRO .EQ. 5 .OR. NFRO .EQ. 11) THEN
         IA = 3
      ELSEIF(NFRO .EQ.7) THEN
         IA = 3
      ELSE
         WRITE(6,*) 'RB NFRO IN FDERFOUX NICHT IMPLEMENTIERT',NFRO
         STOP
      ENDIF

      IF (NBAC .EQ. 5 .OR. NBAC .EQ. 3) THEN
         IE = II-2
      ELSEIF(NBAC .EQ.7) THEN
         IE = II-2
      ELSE
         WRITE(6,*) 'RANDBEDINGUNG IN FDERFOUX NICHT IMPLEMENTIERT'
         STOP
      ENDIF

C                                    ERSTER PHYSIKALISCHER PUNKT
      I = IA


      DO J = JSTART, JSTOP
         DO K = KSTART, KSTOP

            GI(K,J,I) = COEFDX(I,1)
            GJ(K,J,I) = COEFDX(I,2)
            GK(K,J,I) = COEFDX(I,3)

            RSGS(K,J,I) =  COEFDX(I,4) * F(K,J,I-1)
     $                   + COEFDX(I,5) * F(K,J,I  )
     $                   + COEFDX(I,6) * F(K,J,I+1)
         ENDDO
      ENDDO


C                                    LETZTER PHYSIKALISCHER PUNKT
      I = IE


      DO J = JSTART, JSTOP
         DO K = KSTART, KSTOP

            GI(K,J,I) = COEFDX(I,1)
            GJ(K,J,I) = COEFDX(I,2)
            GK(K,J,I) = COEFDX(I,3)

            RSGS(K,J,I) =  COEFDX(I,4) * F(K,J,I  )
     $                   + COEFDX(I,5) * F(K,J,I-1)
     $                   + COEFDX(I,6) * F(K,J,I-2)
         ENDDO
      ENDDO


C                                      IM BERECHNUNGSGEBIET
C--------------------Compact Diff for inner domain-----------
C       DO I = IA+1 , IE-1
C
C#ifndef 
C          LCOL(I) = COEFDX(I,1)
C          DIAG(I) = COEFDX(I,2)
C          RCOL(I) = COEFDX(I,3)
C#endif
C
C       ENDDO

       DO I = IA+1 , IE-1


          DO J = JSTART,JSTOP
             DO K = KSTART,KSTOP


C          A =    BP(K,J,I)
C          B = (1-BP(K,J,I))*    BP(K,J,I+1)
C          C = (1-BP(K,J,I))*    BP(K,J,I-1)
C          D = (1-BP(K,J,I))* (1-BP(K,J,I-1))* (1-BP(K,J,I+1))
          A =    BP(K,J,I-1) *   BP(K,J,I)
          B = (1-BP(K,J,I-1))*   BP(K,J,I)
          C = (1-BP(K,J,I  ))*   BP(K,J,I-1)
          D = (1-BP(K,J,I  ))*(1-BP(K,J,I-1))

          GI(K,J,I) =  COEFDX(I,1)    * A 
     $               + COEFDX(I,12+1) * B
     $               + COEFDX(I,24+1) * C 
          GJ(K,J,I) =  COEFDX(I,2)    * A
     $               + COEFDX(I,12+2) * B
     $               + COEFDX(I,24+2) * C
     $               + 1.0            * D
          GK(K,J,I) =  COEFDX(I,3)    * A
     $               + COEFDX(I,12+3) * B
     $               + COEFDX(I,24+3) * C

                RSGS(K,J,I) =
     $   (COEFDX(I,4)    * F(K,J,I)+
     $    COEFDX(I,5)    * FIN(K,J,I)+
     $    COEFDX(I,6)    * F(K,J,I-1))*DDX(I) * A
     $ + (COEFDX(I,12+4) * F(K,J,I-1)+
     $    COEFDX(I,12+5) * F(K,J,I)+
     $    COEFDX(I,12+6) * F(K,J,I+1)) *          B
     $ + (COEFDX(I,24+4) * F(K,J,I)+
     $    COEFDX(I,24+5) * F(K,J,I-1)+
     $    COEFDX(I,24+6) * F(K,J,I-2)) *        C

             ENDDO
          ENDDO
       ENDDO



      CALL THOMASI(KK,JJ,II,KSTART,KSTOP,JSTART,JSTOP,IA,IE,
     $     LCOL,DIAG,RCOL,RSGS,FAKTOR,FD
     $     ,BP,GI,GJ,GK
     $     )

********************  EXTRAPOLATION   *******************************

       IF (NBAC.EQ.3) THEN
          I = IE+1
          DO J = JSTART,JSTOP
             DO K = KSTART,KSTOP
                FD(K,J,I) = - COEFDX(I,1) * FD(K,J,I-1)
     $                     +  COEFDX(I,4) * F (K,J,I-1)
     $                     +  COEFDX(I,5) * F (K,J,I-2)
     $                     +  COEFDX(I,6) * F (K,J,I-3)
            ENDDO
         ENDDO
      ENDIF

C*********************   E  N  D       *******************************
       RETURN
       END

