










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
      SUBROUTINE INTERPOLATEUX(KK,JJ,II,KSTART,JSTART,ISTART,
     $     KSTOP,JSTOP,ISTOP,COEFFX,RSGS,FAKTOR,U,UINI,
     $     NFRO,NBAC,LCOL,DIAG,RCOL
     $     ,BP,GI,GJ,GK
     $     )
C*MGLET***************************************************************
C        I N T E R P O L A T E U X      
C        INTERPOLIERT U IN X-RICHTUNG (UINI)                            
C*MGLET***************************************************************
C
C PARAM: U(K,J,I)       - ZU INTERPOLIERENDE GROESSE(U)
C      : UINI(K,J,I)    - INTERPOLIERTE U IN I(X)-RICHTUNG
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
C UPROG                 : THOMASI: THOMAS-ALGORITHMUS IN I-RICHTUNG
C
C VERS:  10.10.96 (AM)  : ORIGINAL        (KOMPAKT 4. ORDNUNG)
C
C*MGLET***************************************************************

      IMPLICIT NONE

      INTEGER KK,JJ,II,K,J,I,KSTART,JSTART,ISTART,KSTOP,JSTOP,ISTOP,
     $        IA,IE,NFRO,NBAC

      REAL    LCOL(II),DIAG(II),RCOL(II),
     $        U(KK,JJ,II),UINI(KK,JJ,II),RSGS(KK,JJ,II)

      REAL    FAKTOR(KK,JJ,II),COEFFX(II,12*3),A,B,C,D,
     $        GI(KK,JJ,II),GJ(KK,JJ,II),GK(KK,JJ,II),BP(KK,JJ,II)


      IF (NFRO .EQ. 5 .OR. NFRO .EQ. 11) THEN
         IA = 3
      ELSEIF(NFRO .EQ.7) THEN
         IA = 3
      ELSE
         WRITE(6,*) 'RANDBEDINGUNG IN INTERPOLATEUX NOCH NICHT '
         STOP
      ENDIF  
      
      IF (NBAC .EQ. 5 .OR. NBAC .EQ. 3) THEN
         IE = II-2
      ELSEIF(NBAC .EQ.7) THEN
         IE = II-2
      ELSE
         WRITE(6,*) 'RANDBEDINGUNG IN INTERPOLATEUX NOCH NICHT '
         STOP
      ENDIF

C                                    ERSTER PHYSIKALISCHER PUNKT

      I = IA


      DO J = JSTART, JSTOP
         DO K = KSTART, KSTOP

            GI(K,J,I) = COEFFX(I,1)
            GJ(K,J,I) = COEFFX(I,2)
            GK(K,J,I) = COEFFX(I,3)

            RSGS(K,J,I) = COEFFX(I,4) * U(K,J,I-1)
     $                  + COEFFX(I,5) * U(K,J,I  )
     $                  + COEFFX(I,6) * U(K,J,I+1)

         ENDDO
      ENDDO



CC                                    LETZTER PHYSIKALISCHER PUNKT
       I = IE


       DO J = JSTART, JSTOP
          DO K = KSTART, KSTOP

             GI(K,J,I) = COEFFX(I,1)
             GJ(K,J,I) = COEFFX(I,2)
             GK(K,J,I) = COEFFX(I,3)

             RSGS(K,J,I) = COEFFX(I,4) * U(K,J,I  )
     $                   + COEFFX(I,5) * U(K,J,I-1)
     $                   + COEFFX(I,6) * U(K,J,I-2)

          ENDDO
       ENDDO


C                                     KOMPAKTER ANSATZ IM GEBIET
C-----------Compact scheme for inner doamin---------------------
       DO I = IA+1, IE-1


          DO J = JSTART, JSTOP
             DO K = KSTART, KSTOP



c          A =    BP(K,J,I)
c          B =    BP(K,J,I+1) * (1-BP(K,J,I))
c          C = (1-BP(K,J,I))*    BP(K,J,I-1)
c          D = (1-BP(K,J,I))* (1-BP(K,J,I-1))*(1-BP(K,J,I+1))

          A =    BP(K,J,I-1) *   BP(K,J,I)
          B = (1-BP(K,J,I-1))*   BP(K,J,I)
          C = (1-BP(K,J,I  ))*   BP(K,J,I-1)
          D = (1-BP(K,J,I  ))*(1-BP(K,J,I-1))


             GI(K,J,I) =  COEFFX(I,1)    * A 
     $                  + COEFFX(I,12+1) * B
     $                  + COEFFX(I,24+1) * C
             GJ(K,J,I) =  COEFFX(I,2)    * A
     $                  + COEFFX(I,12+2) * B
     $                  + COEFFX(I,24+2) * C
     $                  + 1.0            * D
             GK(K,J,I) =  COEFFX(I,3)    * A
     $                  + COEFFX(I,12+3) * B
     $                  + COEFFX(I,24+3) * C

                RSGS(K,J,I) = 
     $   (COEFFX(I,4)    * U(K,J,I)+
     $    COEFFX(I,5)    * U(K,J,I-1)) *   A
     $ + (COEFFX(I,12+4) * U(K,J,I-1)+
     $    COEFFX(I,12+5) * U(K,J,I)+
     $    COEFFX(I,12+6) * U(K,J,I+1)) *   B
     $ + (COEFFX(I,24+4) * U(K,J,I)+
     $    COEFFX(I,24+5) * U(K,J,I-1)+
     $    COEFFX(I,24+6) * U(K,J,I-2)) *   C

             ENDDO
          ENDDO
       ENDDO

       CALL THOMASI(KK,JJ,II,KSTART,KSTOP,JSTART,JSTOP,IA,IE,
     $              LCOL,DIAG,RCOL,RSGS,FAKTOR,UINI
     $     ,BP,GI,GJ,GK
     $     )


C*********************  EXTRAPOLATION   *******************************


       IF (NBAC.EQ.3) THEN
          I = IE + 1
          DO J = JSTART,JSTOP
             DO K = KSTART,KSTOP
                UINI(K,J,I) = - COEFFX(I,1) * UINI(K,J,I-1)
     $                        + COEFFX(I,4) * U   (K,J,I-1)
     $                        + COEFFX(I,5) * U   (K,J,I-2)
     $                        + COEFFX(I,6) * U   (K,J,I-3)
             ENDDO
          ENDDO
       ENDIF


C*********************   E  N  D       *******************************


       RETURN
       END
