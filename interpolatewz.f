










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
      SUBROUTINE INTERPOLATEWZ(KK,JJ,II,KSTART,JSTART,ISTART,
     $     KSTOP,JSTOP,ISTOP,COEFFZ,RSGS,FAKTOR,W,WINK,
     $     NBOT,NTOP,LCOL,DIAG,RCOL
     $     ,BP,GI,GJ,GK
     $     )

C*MGLET***************************************************************
C        I N T E R P O L A T E W Z      
C        INTERPOLIERT W IN Z-RICHTUNG (WINK)                            
C*MGLET***************************************************************
C
C PARAM: W(K,J,I)         - ZU INTERPOLIERENDE GROESSE(W)
C      : WINK(K,J,I)      - INTERPOLIERTE W IN K(Z)-RICHTUNG
C
C      : COEFFZ(I,1)    - ENTHAELT KOEFFIZIENTEN DER "LINKEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : COEFFZ(I,2)    - ENTHAELT KOEFFIZIENTEN DER HAUPTDIAGONALEN
C      : COEFFZ(I,3)    - ENTHAELT KOEFFIZIENTEN DER "RECHTEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : COEFFZ(I,4),   - KOEFFIZIENTEN AUF DER RECHTEN SEITE DES 
C        COEFFZ(I,5),     GLEICHUNGSSYSTEMS
C        COEFFZ(I,6)
C      : RSGS           - RECHTE SEITE DES GLEICHUNGSSYSTEMS      
C
C DEFINE DIREKTIVEN     : KEINE 
C
C UPROG                 : THOMASK    : THOMAS-ALGORITHMUS IN K-RICHTUNG
C
C VERS:  07.10.97 (AM)  : ORIGINAL        (KOMPAKT 4. ORDNUNG)
C        24.04.03  (FS) : Komplett ueberarbeitet
C
C*MGLET***************************************************************

      IMPLICIT NONE

      INTEGER KK,JJ,II,K,J,I,KSTART,JSTART,ISTART,KSTOP,JSTOP,ISTOP,
     $        KA,KE,NBOT,NTOP

      REAL    LCOL(KK),DIAG(KK),RCOL(KK),
     $        W(KK,JJ,II),WINK(KK,JJ,II),RSGS(KK,JJ,II),
     $        FAKTOR(KK,JJ,II),COEFFZ(KK,12*3),A,B,C,D,
     $        GI(KK,JJ,II),GJ(KK,JJ,II),GK(KK,JJ,II),BP(KK,JJ,II)



C -------------debug for testing interpolation of W in Z
      
C      DO K = 1, KK
C      WRITE(6,*)  COEFFZ(K,1),COEFFZ(K,2),COEFFZ(K,3),COEFFZ(K,4),
C     $            COEFFZ(K,5) ,COEFFZ(K,6) 
C      ENDDO





      IF (NBOT .EQ. 5 .OR. NBOT .EQ. 11) THEN
         KA = 3
      ELSE
         WRITE(6,*) 'RANDBEDINGUNG IN INTERPOLATEWZ NOCH NICHT '
         STOP
      ENDIF

      IF (NTOP .EQ. 5 .OR. NTOP .EQ. 3) THEN
         KE = KK-2
      ELSE
         WRITE(6,*) 'RANDBEDINGUNG IN INTERPOLATEWZ NOCH NICHT '
         STOP
      ENDIF

C                                    ERSTER PHYSIKALISCHER PUNKT

      K = KA


      DO I = ISTART,ISTOP
         DO J = JSTART, JSTOP

            GI(K,J,I) = COEFFZ(K,1)
            GJ(K,J,I) = COEFFZ(K,2)
            GK(K,J,I) = COEFFZ(K,3)

            RSGS(K,J,I) = COEFFZ(K,4) * W(K-1,J,I)
     $                  + COEFFZ(K,5) * W(K  ,J,I)
     $                  + COEFFZ(K,6) * W(K+1,J,I)
C            RSGS(K,J,I) = 0.5 * W(K,J,I)

         ENDDO
      ENDDO

CC                                    LETZTER PHYSIKALISCHER PUNKT
      K = KE

      DO I = ISTART,ISTOP
         DO J = JSTART, JSTOP
             GI(K,J,I) = COEFFZ(K,1)
             GJ(K,J,I) = COEFFZ(K,2)
             GK(K,J,I) = COEFFZ(K,3)
            RSGS(K,J,I) = COEFFZ(K,4) * W(K  ,J,I)
     $                  + COEFFZ(K,5) * W(K-1,J,I)
     $                  + COEFFZ(K,6) * W(K-2,J,I)
C            RSGS(K,J,I) = 0.5 * W(K-1  ,J,I)
         ENDDO
      ENDDO


C                                     KOMPAKTER ANSATZ IM GEBIET
C------------------Compact scheme for innerdomain--------------
      DO I = ISTART,ISTOP
         DO J = JSTART,JSTOP
            DO K = KA+1, KE-1

C          A =    BP(K,J,I)
C          B = (1-BP(K,J,I))*    BP(K+1,J,I)
C          C = (1-BP(K,J,I))*    BP(K-1,J,I)
C          D = (1-BP(K,J,I))* (1-BP(K-1,J,I))*(1-BP(K+1,J,I))
          A =    BP(K-1,J,I) *   BP(K,J,I)
          B = (1-BP(K-1,J,I))*   BP(K,J,I)
          C = (1-BP(K,J,I  ))*   BP(K-1,J,I)
          D = (1-BP(K,J,I  ))*(1-BP(K-1,J,I))

             GI(K,J,I) =  COEFFZ(K,1)    * A
     $                  + COEFFZ(K,12+1) * B
     $                  + COEFFZ(K,24+1) * C
             GJ(K,J,I) =  COEFFZ(K,2)    * A
     $                  + COEFFZ(K,12+2) * B
     $                  + COEFFZ(K,24+2) * C
     $                  + 1.0            * D
             GK(K,J,I) =  COEFFZ(K,3)    * A
     $                  + COEFFZ(K,12+3) * B
     $                  + COEFFZ(K,24+3) * C

                RSGS(K,J,I) =
     $   (COEFFZ(K,4)    * W(K,J,I)+
     $    COEFFZ(K,5)    * W(K-1,J,I)) *   A
     $ + (COEFFZ(K,12+4) * W(K-1,J,I)+
     $    COEFFZ(K,12+5) * W(K,J,I)+
     $    COEFFZ(K,12+6) * W(K+1,J,I)) *   B
     $ + (COEFFZ(K,24+4) * W(K,J,I)+
     $    COEFFZ(K,24+5) * W(K-1,J,I)+
     $    COEFFZ(K,24+6) * W(K-2,J,I)) *   C


            ENDDO
         ENDDO
      ENDDO

      CALL THOMASK(KK,JJ,II,KA,KE,JSTART,JSTOP,ISTART,ISTOP,
     $             LCOL,DIAG,RCOL,RSGS,FAKTOR,WINK
     $     ,BP,GI,GJ,GK
     $     )


C*********************  EXTRAPOLATION   *******************************

      IF (NTOP.EQ. 3) THEN
         K = KE +1
         DO I = ISTART,ISTOP
            DO J = JSTART,JSTOP
               WINK(K,J,I) = - COEFFZ(K,1) * WINK(K-1,J,I)
     $                       + COEFFZ(K,4) * W   (K-1,J,I)
     $                       + COEFFZ(K,5) * W   (K-2,J,I)
     $                       + COEFFZ(K,6) * W   (K-3,J,I)
            ENDDO
         ENDDO
      ENDIF


C*********************   E  N  D       *******************************

       RETURN
       END
