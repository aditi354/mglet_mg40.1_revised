










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
      SUBROUTINE INTERPOLATEVWX(KK,JJ,II,KSTART,JSTART,ISTART,
     $     KSTOP,JSTOP,ISTOP,COEFFX,RSGS,FAKTOR,F,FIN,
     $     NFRO,NBAC,LCOL,DIAG,RCOL
     $     ,BP,GI,GJ,GK
     $     )

C*MGLET***************************************************************
C        I N T E R P O L A T E V W X 
C        INTERPOLIERT V BZW. W IN X-RICHTUNG
C*MGLET***************************************************************
C
C PARAM: F(K,J,I)       - ZU INTERPOLIERENDE GROESSE (V,W)
C      : FIN(K,J,I)     - INTERPOLIERTE GROESSE IN I(X)-RICHTUNG
C
C      : COEFFX(I,1)    - ENTHAELT KOEFFIZIENTEN DER "LINKEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : COEFFX(I,2)    - ENTHAELT KOEFFIZIENTEN DER HAUPTDIAGONALEN
C      : COEFFX(I,3)    - ENTHAELT KOEFFIZIENTEN DER "RECHTEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : COEFFX(I,4),   - KOEFFIZIENTEN AUF DER RECHTEN SEITE DES 
C        COEFFX(I,5),     GLEICHUNGSSYSTEMS
C        COEFFX(I,6)
C      : RSGS           - RECHTE SEITE DES GLEICHUNGSSYSTEMS      
C
C DEFINE DIREKTIVEN     : KEINE 
C
C UPROG                 : THOMASI    : THOMAS-ALGORITHMUS IN I-RICHTUNG
C
C VERS:  10.10.96 (AM)  : ORIGINAL        (KOMPAKT 4. ORDNUNG)
C        22.01.03 (NP)  : FRED_BODY eingefuehrt
C
C*MGLET***************************************************************

      IMPLICIT NONE

      INTEGER KK,JJ,II,K,J,I,KSTART,JSTART,ISTART,KSTOP,JSTOP,ISTOP,
     $        IA,IE,NFRO,NBAC

      REAL    LCOL(II),DIAG(II),RCOL(II),
     $        F(KK,JJ,II),FIN(KK,JJ,II),RSGS(KK,JJ,II)
      REAL    COEFFX(II,12*3),FAKTOR(KK,JJ,II),A,B,C,D,
     $        GI(KK,JJ,II),GJ(KK,JJ,II),GK(KK,JJ,II),BP(KK,JJ,II)

      IF (NFRO .EQ. 5 .OR. NFRO .EQ. 11) THEN
         IA = 3
      ELSEIF(NFRO .EQ.7) THEN
         IA = 3
      ELSE
         WRITE(6,*) 'RANDBEDINGUNG IN INTERPOLATEVWX NOCH NICHT '
         STOP
      ENDIF  
      
      IF (NBAC .EQ. 5 ) THEN
         IE = II-1
      ELSEIF(NBAC .EQ.7) THEN
         IE = II-1
      ELSEIF (NBAC .EQ. 3) THEN
         IE = II-2
c#ifdef 
c         WRITE(6,*) 'RB IN INTERPOLATEVWX NOCH NICHT FUER FRED_BODY'
c         STOP
c#endif
      ELSE
         WRITE(6,*) 'RANDBEDINGUNG IN INTERPOLATEVWX NOCH NICHT '
         STOP
      ENDIF


CC                                    ERSTER PHYSIKALISCHER PUNKT
      I = IA


      DO J = JSTART, JSTOP
         DO K = KSTART, KSTOP

            GI(K,J,I) = COEFFX(I,7)
            GJ(K,J,I) = COEFFX(I,8)
            GK(K,J,I) = COEFFX(I,9)

            RSGS(K,J,I) = COEFFX(I,10) * F(K,J,I-1)
     $                  + COEFFX(I,11) * F(K,J,I  )
     $                  + COEFFX(I,12) * F(K,J,I+1)

         ENDDO
      ENDDO

CC                                    LETZTER PHYSIKALISCHER PUNKT
       I = IE


       DO  J = JSTART, JSTOP
          DO  K = KSTART, KSTOP

             GI(K,J,I) = COEFFX(I,7)
             GJ(K,J,I) = COEFFX(I,8)
             GK(K,J,I) = COEFFX(I,9)

             RSGS(K,J,I) = COEFFX(I,10) * F(K,J,I  )
     $                   + COEFFX(I,11) * F(K,J,I-1)
     $                   + COEFFX(I,12) * F(K,J,I-2)

          ENDDO
       ENDDO
C                                     KOMPAKTER ANSATZ IM GEBIET


       DO I = IA+1, IE-1


          DO J = JSTART, JSTOP
             DO K = KSTART, KSTOP


C          A =    BP(K,J,I)
C          B = (1-BP(K,J,I))*    BP(K,J,I+1)
C          C = (1-BP(K,J,I))*    BP(K,J,I-1)
C          D = (1-BP(K,J,I))* (1-BP(K,J,I-1))* (1-BP(K,J,I+1))
          A =    BP(K,J,I-1) *   BP(K,J,I)
          B = (1-BP(K,J,I-1))*   BP(K,J,I)
          C = (1-BP(K,J,I  ))*   BP(K,J,I-1)
          D = (1-BP(K,J,I  ))*(1-BP(K,J,I-1))

             GI(K,J,I) =  COEFFX(I,7)    * A
     $                  + COEFFX(I,12+7) * B
     $                  + COEFFX(I,24+7) * C
             GJ(K,J,I) =  COEFFX(I,8)    * A
     $                  + COEFFX(I,12+8) * B
     $                  + COEFFX(I,24+8) * C
     $                  + 1.0            * D
             GK(K,J,I) =  COEFFX(I,9)    * A
     $                  + COEFFX(I,12+9) * B
     $                  + COEFFX(I,24+9) * C

                RSGS(K,J,I) =
     $   (COEFFX(I,10)    * F(K,J,I)+
     $    COEFFX(I,11)    * F(K,J,I-1)) *  A
     $ + (COEFFX(I,12+10) * F(K,J,I-1)+
     $    COEFFX(I,12+11) * F(K,J,I)+
     $    COEFFX(I,12+12) * F(K,J,I+1)) *  B
     $ + (COEFFX(I,24+10) * F(K,J,I)+
     $    COEFFX(I,24+11) * F(K,J,I-1)+
     $    COEFFX(I,24+12) * F(K,J,I-2)) *  C



             ENDDO
          ENDDO
       ENDDO

       CALL THOMASI(KK,JJ,II,KSTART,KSTOP,JSTART,JSTOP,IA,IE,
     $              LCOL,DIAG,RCOL,RSGS,FAKTOR,FIN
     $     ,BP,GI,GJ,GK
     $     )


C*********************  EXTRAPOLATION   *******************************

       IF (NBAC.EQ.3) THEN
          I = IE + 1
          DO J = JSTART,JSTOP
             DO K = KSTART,KSTOP
                FIN(K,J,I) = - COEFFX(I,7 ) * FIN(K,J,I-1)
     $                       + COEFFX(I,10) * F  (K,J,I-1)
     $                       + COEFFX(I,11) * F  (K,J,I-2)
     $                       + COEFFX(I,12) * F  (K,J,I-3)
             ENDDO
          ENDDO
       ENDIF
C
C*********************   E  N  D       *******************************

       RETURN
       END
