










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
      SUBROUTINE INTERPOLATEUX3 (KK,JJ,II,KSTART,JSTART,ISTART,KSTOP,
     $                         JSTOP,ISTOP,COEFFX,RSGS,FAKTOR,
     $                         U,UINI,NFRO,NBAC,LCL,DIL,RCL,
     $                         DX,DDX,B,LCR,DIR,RCR)
C*MGLET***************************************************************
C        I N T E R P O L A T E U X 3      
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
C UPROG                 :  THOMASI: THOMAS-ALGORITHMUS IN I(X)-RICHTUNG
C
C
C VERS:  10.10.99 (AM)  : ORIGINAL    (KOMPAKT-UPWIND 3. ORDNUNG)
C
C*MGLET***************************************************************
C

      REAL       DX(II),DDX(II)
C
      REAL       FAKTOR(KK,JJ,II)
      REAL       LCL(KK,JJ,II),DIL(KK,JJ,II),RCL(KK,JJ,II)
      REAL       LCR(KK,JJ,II),DIR(KK,JJ,II),RCR(KK,JJ,II)
      REAL       LCOL,DIAG,RCOL,RSIDE1
C
      REAL       U(KK,JJ,II),UINI(KK,JJ,II),RSGS(KK,JJ,II),
     $           COEFFX(II,12), B(KK,JJ,II)
C
C
CC                                    ERSTER PHYSIKALISCHER PUNKT
        I = ISTART
        DO  J = JSTART, JSTOP
        DO  K = KSTART, KSTOP
C
         LCL(K,J,I) = 0.0
         DIL(K,J,I) = 1.0
         RCL(K,J,I) = 0.0
C        RSGS(K,J,I) = 0.5 * (U(K,J,I) + U(K,J,I-1))
         RSGS(K,J,I) = 1.875 * U(K,J,I) -
     $          1.25 * U(K,J,I+1)  +0.375*U(K,J,I+2)

        ENDDO                
        ENDDO                
C
C
        I = ISTOP+1
        DO  J = JSTART, JSTOP
        DO  K = KSTART, KSTOP
C
         LCL(K,J,I) = 0.0
         DIL(K,J,I) = 1.0
         RCL(K,J,I) = 0.0
         RSGS(K,J,I) = 0.5 * (U(K,J,I) + U(K,J,I-1))
C        RSGS(K,J,I) = 1.875 * U(K,J,I) -
C    $          1.25 * U(K,J,I+1)  +0.375*U(K,J,I+2)

        ENDDO
        ENDDO
C                                     KOMPAKTER ANSATZ IM GEBIET (FLUID)
C
      DO  I = ISTART+1, ISTOP
C
C                                     PROBLEM: if 'COUNINGH',
C                                     calculate 2nd order in this segment
CKOMPCOUN       IKB = 3
       IKB = 20
       IF (I .LT. IKB)  THEN
C
        DO  J = JSTART, JSTOP
        DO  K = KSTART, KSTOP
C
             LCL(K,J,I) = 0.0
             DIL(K,J,I) = 1.0
             RCL(K,J,I) = 0.0
        RSGS (K,J,I)  = (U(K,J,I-1)+(U(K,J,I)-U(K,J,I-1))
     $                  *0.5*DX(I-1)/DDX(I))
        ENDDO
        ENDDO
C
       ELSE

        DO  J = JSTART, JSTOP
        DO  K = KSTART, KSTOP
CCCCC                                        Kompakt-Upwind- 3th
C
         IF (U(K,J,I) .GE. 0.) THEN
         CALL INTERCOEFKU3(I,II,DX,LCOL,DIAG,RCOL,RSIDE1,1.0)
C
 	     LCL(K,J,I) = LCOL
 	     DIL(K,J,I) = 1.0
 	     RCL(K,J,I) = RCOL
             RSGS(K,J,I) = RSIDE1*U(K,J,I-1)
C
         ELSEIF (U(K,J,I) .LT. 0.) THEN
         CALL INTERCOEFKU3(I,II,DX,LCOL,DIAG,RCOL,RSIDE1,2.0)
C
 	     LCL(K,J,I) = LCOL
 	     DIL(K,J,I) = 1.0
 	     RCL(K,J,I) = RCOL
             RSGS(K,J,I) = RSIDE1*U(K,J,I)
C
         ENDIF
C
        ENDDO                
        ENDDO                
C
       ENDIF
C
      ENDDO                
C                                      KOERPER IM GEBIET ?????????
CCC
C
        DO  I = ISTART+1, ISTOP-1
        DO  J = JSTART, JSTOP
        DO  K = KSTART, KSTOP
C
C                           Fluidzelle und rechts davon Koerperzelle
C
          IF (B(K,J,I) .GT. 0.0 .AND. B(K,J,I+1) .LT. 0.0)  THEN
         LCL(K,J,I) = 0.5
         DIL(K,J,I) = 1.0
         RCL(K,J,I) = -(1./6.)
         RSGS(K,J,I) = (4./3.)*U(K,J,I-1)
C
C        RSGS(K,J,I) = 2.1875 * (U(K,J,I-1) - U(K,J,I-2)) +
C    $          1.3125 * U(K,J,I-3)  + 0.3125*U(K,J,I-4)
C
C
C                           Fluidzelle und links davon Koerperzelle
C
          ELSEIF (B(K,J,I) .GT. 0.0 .AND. B(K,J,I-1) .LT. 0.0)  THEN
             LCL(K,J,I) = -(1./6.)
             DIL(K,J,I) = 1.0
             RCL(K,J,I) = 0.5
             RSGS(K,J,I) = (4./3.)*U(K,J,I)
C
C        LCL(K,J,I) = 0.0
C        DIL(K,J,I) = 1.0
C        RCL(K,J,I) = 0.0
C        RSGS(K,J,I) = 2.1875 * (U(K,J,I) - U(K,J,I+1)) +
C    $          1.3125 * U(K,J,I+2)  + 0.3125*U(K,J,I+3)

C        RSGS(K,J,I) = 0.5 * (U(K,J,I) + U(K,J,I-1))

          ELSEIF (B(K,J,I) .LT. 0.0)  THEN
            LCL(K,J,I) = 0.0
            DIL(K,J,I) = 1.0
            RCL(K,J,I) = 0.0
            RSGS(K,J,I) = 0.0
CCCC
C
          ENDIF
C
        ENDDO                
        ENDDO                
        ENDDO                
C
        CALL THOMASI(KK,JJ,II,KSTART,KSTOP,JSTART,JSTOP,ISTART,
     $            ISTOP+1,LCL,DIL,RCL,RSGS,FAKTOR,UINI)
C
C

C

C*********************  RECHTER RAND   *******************************
C    wird im Thomas-Algorithmus mitberechnet
C
C     I = ISTOP +1
C     DO J = JSTART,JSTOP
C     DO K = KSTART,KSTOP
C
C      IF (NBAC.EQ.7) THEN
C
C       UINI(K,J,I) = 0.5 * ( U(K,J,I-1) + U(K,J,I) )
C
C      ELSE
C
C       UINI(K,J,I) = 0.5 * ( U(K,J,I-1) + U(K,J,I) )
C	UINI(K,J,I) = - COEFFX(I,1)*UINI(K,J,I-1) + COEFFX(I,4)*
C    $  U(K,J,I-1) + COEFFX(I,5)*U(K,J,I-2) + COEFFX(I,6)*U(K,J,I-3)
C
C      ENDIF
C
C       ENDDO              
C     ENDDO              
C
C
C*********************   E  N  D       *******************************
       RETURN
       END

