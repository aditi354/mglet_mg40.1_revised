










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
      SUBROUTINE INTERPOLATEVY(KK,JJ,II,KSTART,JSTART,ISTART,KSTOP,
     $                        JSTOP,ISTOP,COEFFY,RSGS,FAKTOR,
     $                        V,VINJ,NRGT,NLFT,LCOL,DIAG,RCOL
     $     ,BP,GI,GJ,GK
     $     )

C*MGLET***************************************************************
C        I N T E R P O L A T E V Y      
C        INTERPOLIERT V IN Y-RICHTUNG (VINJ)                            
C*MGLET***************************************************************
C
C PARAM: V(K,J)         - ZU INTERPOLIERENDE GROESSE(V)
C      : VINJ(K,J)      - INTERPOLIERTE V IN J(Y)-RICHTUNG
C
C      : COEFFY(I,1)    - ENTHAELT KOEFFIZIENTEN DER "LINKEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : COEFFY(I,2)    - ENTHAELT KOEFFIZIENTEN DER HAUPTDIAGONALEN
C      : COEFFY(I,3)    - ENTHAELT KOEFFIZIENTEN DER "RECHTEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : COEFFY(I,4),   - KOEFFIZIENTEN AUF DER RECHTEN SEITE DES 
C        COEFFY(I,5),     GLEICHUNGSSYSTEMS
C        COEFFY(I,6)
C      : RSGS           - RECHTE SEITE DES GLEICHUNGSSYSTEMS      
C
C DEFINE DIREKTIVEN     : KEINE 
C
C UPROG                 : THOMASJ: THOMAS-ALGORITHMUS IN J(Y)-RICHTUNG
C
C VERS:  30.06.03 (FS)  : ORGINAL
C
C*MGLET***************************************************************
C
      IMPLICIT NONE

      INTEGER KK,JJ,II,K,J,I,KSTART,JSTART,ISTART,
     $     ISTOP,KSTOP,JSTOP,NRGT,NLFT,JA,JE

      REAL    LCOL(JJ),DIAG(JJ),RCOL(JJ),
     $        V(KK,JJ,II),VINJ(KK,JJ,II),RSGS(KK,JJ,II),
     $        FAKTOR(KK,JJ,II),COEFFY(JJ,12*3),A,B,C,D,
     $        GI(KK,JJ,II),GJ(KK,JJ,II),GK(KK,JJ,II),BP(KK,JJ,II)

      IF(NRGT.EQ.5) THEN
         JA = 3
      ELSE
         WRITE(6,*) 'RANDBEDINGUNG IN INTERPOLATEVY NOCH NICHT '
         STOP
      ENDIF

      IF(NLFT.EQ.5) THEN
         JE = JJ - 2
      ELSE
         WRITE(6,*) 'RANDBEDINGUNG IN INTERPOLATEVY NOCH NICHT '
         STOP
      ENDIF

C-------------------------ERSTER PHYSIKALISCHER PUNKT ----------


      J = JA
      DO I = ISTART,ISTOP
         DO K = KSTART, KSTOP
            GI(K,J,I) = COEFFY(J,1)
            GJ(K,J,I) = COEFFY(J,2)
            GK(K,J,I) = COEFFY(J,3)
            RSGS(K,J,I) = COEFFY(J,4) * V(K,J-1,I) +
     $                    COEFFY(J,5) * V(K,J,I)   + 
     $                    COEFFY(J,6) * V(K,J+1,I)

         ENDDO
      ENDDO
C-------------------------LETZTER PHYSIKALISCHER PUNKT ---------

      J = JE
      DO I = ISTART,ISTOP
         DO K = KSTART,KSTOP
            GI(K,J,I) = COEFFY(J,1)
            GJ(K,J,I) = COEFFY(J,2)
            GK(K,J,I) = COEFFY(J,3)
            RSGS(K,J,I) = COEFFY(J,4) * V(K,J,I) + 
     $                    COEFFY(J,5) * V(K,J-1,I) +  
     $                    COEFFY(J,6) * V(K,J-2,I)
         ENDDO
      ENDDO
C-------------------------KOMPAKTER ANSATZ IM GEBIET ----------
C--------------------Compact scheme for inner domain-----------
      DO I = ISTART,ISTOP
         DO J = JA+1, JE-1
            DO K = KSTART,KSTOP


C          A =    BP(K,J,I)
C          B = (1-BP(K,J,I))*    BP(K,J+1,I)
C          C = (1-BP(K,J,I))*    BP(K,J-1,I)
C          D = (1-BP(K,J,I))* (1-BP(K,J-1,I))*(1-BP(K,J+1,I))
          A =    BP(K,J-1,I) *   BP(K,J,I)
          B = (1-BP(K,J-1,I))*   BP(K,J,I)
          C = (1-BP(K,J,I  ))*   BP(K,J-1,I)
          D = (1-BP(K,J,I  ))*(1-BP(K,J-1,I))

             GI(K,J,I) =  COEFFY(J,1)    * A
     $                  + COEFFY(J,12+1) * B
     $                  + COEFFY(J,24+1) * C
             GJ(K,J,I) =  COEFFY(J,2)    * A
     $                  + COEFFY(J,12+2) * B
     $                  + COEFFY(J,24+2) * C
     $                  + 1.0            * D
             GK(K,J,I) =  COEFFY(J,3)    * A
     $                  + COEFFY(J,12+3) * B
     $                  + COEFFY(J,24+3) * C

                RSGS(K,J,I) =
     $   (COEFFY(J,4)    * V(K,J,I)+
     $    COEFFY(J,5)    * V(K,J-1,I)) *   A
     $ + (COEFFY(J,12+4) * V(K,J-1,I)+
     $    COEFFY(J,12+5) * V(K,J,I)+
     $    COEFFY(J,12+6) * V(K,J+1,I)) *   B
     $ + (COEFFY(J,24+4) * V(K,J,I)+
     $    COEFFY(J,24+5) * V(K,J-1,I)+
     $    COEFFY(J,24+6) * V(K,J-2,I)) *   C


            ENDDO
         ENDDO
      ENDDO

      CALL THOMASJ(KK,JJ,II,KSTART,KSTOP,JA,JE,ISTART,ISTOP,
     $             LCOL,DIAG,RCOL,RSGS,FAKTOR,VINJ
     $     ,BP,GI,GJ,GK
     $     )

C********************* EXTRAPOLATION *******************************
C
      IF (NLFT .EQ. 3) THEN
         J = JE + 1
         DO I = ISTART,ISTOP
            DO K = KSTART,KSTOP
               VINJ(K,J,I) = - COEFFY(J,1)*VINJ(K,J-1,I) 
     $                       + COEFFY(J,4)*V   (K,J-1,I) 
     $                       + COEFFY(J,5)*V   (K,J-2,I)
     $                       + COEFFY(J,6)*V   (K,J-3,I)

            ENDDO
         ENDDO    
      ENDIF

C*********************   E  N  D       *******************************

      RETURN
      END

