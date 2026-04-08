










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
      SUBROUTINE TRIZYK3DJ (KK,JJ,II,KL,KU,JL,JU,IL,IU,
     $                           UD,HD,OD,UZ,RSP,RS,X)
C***********************************************************************
C
C LOESUNG EINES LINEAREN GLEICHUNGSSYSTEMS MIT ZYKLISCH
C TRIDIAGONALER MATRIX (3-D FALL)
C
C  10. 05. 97 (FS):   - FUER J-RICHTUNG                     
C***********************************************************************
C
      IMPLICIT NONE

      INTEGER KK,JJ,II,K,J,I,KL,KU,JL,JU,IL,IU,IFEHL

      REAL    UD(JJ), HD(JJ), OD(JJ), RS(KK,JJ,II),  X(KK,JJ,II)
      REAL    UZ(JJ), RSP(JJ)
C
C                                 LOESUNG DES GLEICHUNGSSYSTEMS
C
      CALL FZYKTR3DJ (KK,JJ,II,KL,KU,JL,JU,IL,IU,
     $     UD,HD,OD,UZ,RSP,RS,X,IFEHL)
C
      IF(IFEHL .NE. 0) THEN
         WRITE (6,*) 'FEHLER AUFGETRETEN, trizyk3dj.src'
         STOP ' TRIZYK ERROR'
      END IF

       RETURN
       END

C-----------------------------------------------------------------------
      SUBROUTINE FZYKTR3DJ (KK,JJ,II,KL,KU,JL,JU,IL,IU,
     $     UD,HD,OD,UZ,RSP,RS,X,IFEHL)

      IMPLICIT NONE

      INTEGER KK,JJ,II,K,J,I,KL,KU,JL,JU,IL,IU,IFEHL,N

      REAL    UD(JJ),   HD(JJ),   OD(JJ),   UZ(JJ),   RSP(JJ),
     $     RS(KK,JJ,II),   X (KK,JJ,II)

      N = JU - JL +1
      IFEHL  = -1
      IF(N .LT. 3) RETURN

      CALL FZYKTZ3DJ (JJ,JL,JU,UD,HD,OD,UZ,RSP,IFEHL)

      IF(IFEHL .EQ. 0) THEN

      CALL FZYKTL3D_SJ (KK,JJ,II,KL,KU,JL,JU,IL,IU,
     $	             UD,HD,OD,UZ,RSP,RS,X)

      END IF

      RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE FZYKTZ3DJ (JJ,JL,JU,UD,HD,OD,UZ,RSP,IFEHL)

      IMPLICIT NONE

      INTEGER JJ,J,JL,JU,IFEHL,N

      REAL    UD(JJ),   HD(JJ),   OD(JJ),   UZ(JJ),   RSP(JJ)

      REAL HILF,S

      N = JU - JL +1
      IFEHL  = -1
      IF(N .LT. 3) RETURN

      UD (JL)   = 0.0
      OD (JU)   = 0.0
      UZ (JU-1) = 0.0
      UZ (JU)   = 0.0
      RSP(JU-1) = 0.0
      RSP(JU)   = 0.0
      IF(HD(JL) .EQ. 0.0) THEN
         IFEHL = 1
         RETURN
      END IF

      HILF   = 1.0 / HD(JL)
      OD(JL)  = OD(JL) * HILF
      RSP(JL) = RSP(JL) * HILF

      DO J = JL+1,JU-2
         HD(J) = HD(J)-UD(J)*OD(J-1)
         IF(HD(J) .EQ. 0.0) THEN
            IFEHL = J
            RETURN
         END IF
         HILF  = 1.0 / HD(J)
         OD(J)  = OD(J) * HILF
         RSP(J) = -UD(J)*RSP(J-1)*HILF
      ENDDO

      HD(JU-1) = HD(JU-1) - UD(JU-1)*OD(JU-2)
      IF(HD(JU-1) .EQ. 0.0) THEN
         IFEHL  = JU-1
         RETURN
      END IF
C
      DO J = JL+1,JU-2
         UZ(J) = -UZ(J-1)*OD(J-1)
      ENDDO

      UD(JU)    = UD(JU) - UZ(JU-2)*OD(JU-2)
      OD(JU-1)  = (OD(JU-1) - UD(JU-1)*RSP(JU-2)) / HD(JU-1)
      S        = 0.0
C
      DO J = JL,JU-2
         S = S - UZ(J)*RSP(J)
      ENDDO
C
      HD(JU)   = HD(JU) + S - UD(JU)*OD(JU-1)
      IF(HD(JU) .EQ. 0.0) THEN
         IFEHL  = N
         RETURN
      END IF
      IFEHL = 0
      RETURN
      END

C----------------------------------------------------------------
      SUBROUTINE FZYKTL3DJ (KK,JJ,II,KL,KU,JL,JU,IL,IU,
     $                     UD,HD,OD,UZ,RSP,RS,X)

      IMPLICIT NONE

      INTEGER KK,JJ,II,K,J,I,KL,KU,JL,JU,IL,IU

      REAL    UD(JJ),   HD(JJ),   OD(JJ),   UZ(JJ),   RSP(JJ),
     $        RS(KK,JJ,II),   X (KK,JJ,II)

      REAL S


      DO I = IL, IU
         DO K = KL, KU

            RS(K,JL,I)   = RS(K,JL,I) / HD(JL)

         ENDDO
      ENDDO

      DO  I = IL,IU
         DO J = JL+1, JU-1
            DO K = KL, KU
               RS(K,J,I) = (RS(K,J,I) - RS(K,J-1,I) * UD(J)) / HD(J)
            ENDDO
         ENDDO
      ENDDO

      DO I = IL, IU
         DO K = KL, KU
            S       = 0.0
            DO J = JL,JU-2
               S     = S - UZ(J) * RS(K,J,I)
            ENDDO
C
            RS(K,JU,I)  = 
     $           (RS(K,JU,I) + S - UD(JU) * RS(K,JU-1,I)) / HD(JU)
         ENDDO
      ENDDO
C
C                                 BERECHNUNG DER LOESUNGEN DURCH
C                                 RUECKWAERTSELIMINATION
C
      DO I = IL, IU
         DO K = KL, KU
            X(K,JU,I)    = RS(K,JU,I)
            X(K,JU-1,I)  = RS(K,JU-1,I) - X(K,JU,I) * OD(JU-1)
         ENDDO
      ENDDO

      DO I = IL,IU
         DO J = JU-1, JL,-1
            DO K = KL, KU
               X(K,J,I)  = RS(K,J,I) - OD(J) * X(K,J+1,I) - 
     $              RSP(J) * X(K,JU,I)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
c-----------------------------------------------------------------
      SUBROUTINE FZYKTL3D_SJ (KK,JJ,II,KL,KU,JL,JU,IL,IU,
     $                     UD,HD,OD,UZ,RSP,RS,X)

      IMPLICIT NONE

      INTEGER KK,JJ,II,K,J,I,KL,KU,JL,JU,IL,IU

      REAL    UD(JJ),   HD(JJ),   OD(JJ),   UZ(JJ),   RSP(JJ),
     $        RS(KK,JJ,II),   X (KK,JJ,II)

      REAL S

      REAL SUM(1030)

      IF (KK .GT. 1030) CALL ERRR(501,'FZYKTL3D_S')
C
C                                 VORWAERTSELIMINATION
C
C-------------------------------------------------------
      DO I = IL, IU
         DO K = KL, KU
            RS(K,JL,I)   = RS(K,JL,I) / HD(JL)
         ENDDO
      ENDDO

      DO  I = IL,IU
         DO J = JL+1, JU-1
            DO K = KL, KU
               RS(K,J,I) = (RS(K,J,I) - RS(K,J-1,I) * UD(J)) / HD(J)
            ENDDO
         ENDDO
      ENDDO

      DO I = IL, IU
         DO K = KL, KU
            SUM(K)       = 0.0
         ENDDO
         DO J = JL,JU-2
            DO K = KL, KU
               SUM(K)     = SUM(K) - UZ(J) * RS(K,J,I)
            ENDDO
         ENDDO

         DO K = KL, KU
            S = SUM(K)
            RS(K,JU,I)  = 
     $           (RS(K,JU,I) + S - UD(JU) * RS(K,JU-1,I)) / HD(JU)
         ENDDO
      ENDDO
C
C                                 BERECHNUNG DER LOESUNGEN DURCH
C                                 RUECKWAERTSELIMINATION
C

      DO I = IL, IU
         DO K = KL, KU
            X(K,JU,I)    = RS(K,JU,I)
            X(K,JU-1,I)  = RS(K,JU-1,I) - X(K,JU,I) * OD(JU-1)
         ENDDO
      ENDDO

      DO I = IL,IU
         DO J = JU-1, JL,-1
            DO K = KL, KU
               X(K,J,I)  = RS(K,J,I) - OD(J) * X(K,J+1,I) - 
     $              RSP(J) * X(K,JU,I)
            ENDDO
         ENDDO
      ENDDO
C
C-------------------------------------------------------
C
      RETURN
      END

