










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
      SUBROUTINE FDERFOVY(KK,JJ,II,KSTART,JSTART,ISTART,KSTOP,JSTOP,
     $                    ISTOP,DY,DDY,COEFDY,RSGS,FAKTOR,
     $                    F,FIN,FD,NRGT,NLFT,LCOL,DIAG,RCOL
     $     ,BP,GI,GJ,GK
     $     )

C*MGLET***************************************************************
C        F D E R F O V Y  (First DERive Fourth Order)
C        BERECHNUNG DER ERSTEN ABLEITUNG (DV/DY)                    
C*MGLET***************************************************************
C
C PARAM: F(K,J)         - BEKANNTE GROESSE AN DEN GITTERKANTEN
C      : FIN(K,J)       - INTERPOLIERTE GROESSE IN DER GLEICHEN RICHTUNG
C                         ZWISCHEN DEN BEIDEN KANTEN       
C      : FD (K,J)       - ERSTE ABLEITUNG, FIRST DERINATIVE  DV/DY
C
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
C VERS:  08.01.97(AM)  : ORIGINAL     (KOMPAKT 4. ORDNUNG)
C        01.07.03 (FS) : ORGINAL AUS AM-ROUNTINE
C
C*MGLET***************************************************************

      IMPLICIT NONE
      INTEGER KK,JJ,II,K,J,I,KSTART,JSTART,ISTART,KSTOP,JSTOP,ISTOP,
     $        JA,JE,NRGT,NLFT

      REAL    DY(JJ),DDY(JJ),LCOL(JJ),DIAG(JJ),RCOL(JJ),
     $        F(KK,JJ,II),FIN(KK,JJ,II),FD(KK,JJ,II),RSGS(KK,JJ,II),
     $        COEFDY(JJ,12*3),FAKTOR(KK,JJ,II),A,B,C,D,
     $        GI(KK,JJ,II),GJ(KK,JJ,II),GK(KK,JJ,II),BP(KK,JJ,II)

      IF(NRGT.EQ.5) THEN
         JA = 3
      ELSE
         WRITE(6,*) 'RANDBEDINGUNG IN FDERFOVY NOCH NICHT '
         STOP
      ENDIF

      IF(NLFT.EQ.5) THEN
         JE = JJ-2
      ELSE
         WRITE(6,*) 'RANDBEDINGUNG IN FDERFOVY NOCH NICHT '
         STOP
      ENDIF
C                                    ERSTER PHYSIKALISCHER PUNKT      
      J = JA
      DO I=ISTART,ISTOP
         DO K = KSTART, KSTOP
            GI(K,J,I) = COEFDY(J,1)
            GJ(K,J,I) = COEFDY(J,2)
            GK(K,J,I) = COEFDY(J,3)
            RSGS(K,J,I) = COEFDY(J,4)*F(K,J-1,I) 
     $                + COEFDY(J,5)*F(K,J,I) 
     $                + COEFDY(J,6)*F(K,J+1,I)
         ENDDO
      ENDDO
C                                    LETZTER PHYSIKALISCHER PUNKT
      J = JE
      DO I=ISTART,ISTOP
         DO K = KSTART, KSTOP
            GI(K,J,I) = COEFDY(J,1)
            GJ(K,J,I) = COEFDY(J,2)
            GK(K,J,I) = COEFDY(J,3)
            RSGS(K,J,I) = COEFDY(J,4)*F(K,J,I)
     $                + COEFDY(J,5)*F(K,J-1,I)
     $                + COEFDY(J,6)*F(K,J-2,I)
         ENDDO
      ENDDO
C                                                       IM GEBIET
C-------------------------Compact differentiation for inner domain---
      DO I = ISTART,ISTOP
         DO J = JA+1 , JE-1
            DO K = KSTART,KSTOP

C          A =     BP(K,J,I)
C          B = (1.-BP(K,J,I))*      BP(K,J+1,I)
C          C = (1.-BP(K,J,I))*      BP(K,J-1,I)
C          D = (1.-BP(K,J,I))* (1.-BP(K,J-1,I))* (1.-BP(K,J+1,I))
          A =    BP(K,J-1,I) *   BP(K,J,I)
          B = (1-BP(K,J-1,I))*   BP(K,J,I)
          C = (1-BP(K,J,I  ))*   BP(K,J-1,I)
          D = (1-BP(K,J,I  ))*(1-BP(K,J-1,I))

          GI(K,J,I) =  COEFDY(J,1)    * A
     $               + COEFDY(J,12+1) * B
     $               + COEFDY(J,24+1) * C
          GJ(K,J,I) =  COEFDY(J,2)    * A
     $               + COEFDY(J,12+2) * B
     $               + COEFDY(J,24+2) * C
     $               + 1.0            * D
          GK(K,J,I) =  COEFDY(J,3)    * A
     $               + COEFDY(J,12+3) * B
     $               + COEFDY(J,24+3) * C

                RSGS(K,J,I) =
     $   (COEFDY(J,4)    * F(K,J,I)+
     $    COEFDY(J,5)    * FIN(K,J,I)+
     $    COEFDY(J,6)    * F(K,J-1,I))*DDY(J) * A
     $ + (COEFDY(J,12+4) * F(K,J-1,I)+
     $    COEFDY(J,12+5) * F(K,J,I)+
     $    COEFDY(J,12+6) * F(K,J+1,I)) *          B
     $ + (COEFDY(J,24+4) * F(K,J,I)+
     $    COEFDY(J,24+5) * F(K,J-1,I)+
     $    COEFDY(J,24+6) * F(K,J-2,I)) *        C
            ENDDO                   
         ENDDO          
      ENDDO         

       CALL THOMASJ(KK,JJ,II,KSTART,KSTOP,JA,JE,ISTART,ISTOP,
     $              LCOL,DIAG,RCOL,RSGS,FAKTOR,FD
     $     ,BP,GI,GJ,GK
     $     )

C*********************  EXTRAPOLATION *******************************

       IF(NLFT .EQ. 3) THEN
          J = JE+1
          DO I = ISTART,ISTOP
             DO K = KSTART,KSTOP
                FD(K,J,I) = - COEFDY(J,1)*FD(K,J-1,I) + 
     $                        COEFDY(J,4)*F(K,J-1,I) + 
     $                        COEFDY(J,5)*F(K,J-2,I) + 
     $                        COEFDY(J,6)*F(K,J-3,I)
             ENDDO  
          ENDDO            
      ENDIF

C*********************   E  N  D       *******************************
       RETURN
       END

