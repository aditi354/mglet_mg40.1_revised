










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
      SUBROUTINE TRIZYK3D (KK,JJ,II,KL,KU,JL,JU,IL,IU,
     $                           UD,HD,OD,UZ,RSP,RS,X)
C***********************************************************************
C
C LOESUNG EINES LINEAREN GLEICHUNGSSYSTEMS MIT ZYKLISCH
C TRIDIAGONALER MATRIX (3-D FALL)
C
C  10. 05. 97 (AM):   - ORIGINAL                      
C***********************************************************************
C
      REAL    UD(II), HD(II), OD(II), RS(KK,JJ,II),  X(KK,JJ,II)
      REAL    UZ(II), RSP(II)
C
C                                 LOESUNG DES GLEICHUNGSSYSTEMS
C
      CALL FZYKTR3D (KK,JJ,II,KL,KU,JL,JU,IL,IU,
     $               UD,HD,OD,UZ,RSP,RS,X,IFEHL)
C
      IF(IFEHL .NE. 0) THEN
         WRITE (6,*) 'FEHLER AUFGETRETEN, trizyk3d.src '
         STOP ' TRIZYK ERROR'
      END IF

       RETURN
       END

C-----------------------------------------------------------------------
      SUBROUTINE FZYKTR3D (KK,JJ,II,KL,KU,JL,JU,IL,IU,
     $                     UD,HD,OD,UZ,RSP,RS,X,IFEHL)

      REAL    UD(II),   HD(II),   OD(II),   UZ(II),   RSP(II),
     $        RS(KK,JJ,II),   X (KK,JJ,II)

      N = IU - IL +1
      IFEHL  = -1
      IF(N .LT. 3) RETURN
C
C                                 ZERLEGUNG DER MATRIX  A
C
      CALL FZYKTZ3D (II,IL,IU,UD,HD,OD,UZ,RSP,IFEHL)
C
C                                 FALLS IFEHL=0: VORWAERTS- UND
C                                 RUECKWAERTSELIMINATION
C
      IF(IFEHL .EQ. 0) THEN

      CALL FZYKTL3D_S (KK,JJ,II,KL,KU,JL,JU,IL,IU,
     $	             UD,HD,OD,UZ,RSP,RS,X)

      END IF

      RETURN
      END
C$
      SUBROUTINE FZYKTZ3D (II,IL,IU,UD,HD,OD,UZ,RSP,IFEHL)
C***********************************************************************
C
C***********************************************************************
C
      REAL    UD(II),   HD(II),   OD(II),   UZ(II),   RSP(II)
C
      N = IU - IL +1
      IFEHL  = -1
      IF(N .LT. 3) RETURN

      UD (IL)   = 0.0
      OD (IU)   = 0.0
      UZ (IU-1) = 0.0
      UZ (IU)   = 0.0
      RSP(IU-1) = 0.0
      RSP(IU)   = 0.0
      IF(HD(IL) .EQ. 0.0) THEN
         IFEHL = 1
         RETURN
      END IF
      HILF   = 1.0 / HD(IL)
      OD(IL)  = OD(IL) * HILF
      RSP(IL) = RSP(IL) * HILF
C
      DO 10 I = IL+1,IU-2
         HD(I) = HD(I)-UD(I)*OD(I-1)
         IF(HD(I) .EQ. 0.0) THEN
            IFEHL = I
            RETURN
         END IF
         HILF  = 1.0 / HD(I)
         OD(I)  = OD(I) * HILF
         RSP(I) = -UD(I)*RSP(I-1)*HILF
   10   CONTINUE
      HD(IU-1) = HD(IU-1) - UD(IU-1)*OD(IU-2)
      IF(HD(IU-1) .EQ. 0.0) THEN
         IFEHL  = IU-1
         RETURN
      END IF
C
      DO 20 K = IL+1,IU-2
         UZ(K) = -UZ(K-1)*OD(K-1)
   20 CONTINUE
C
      UD(IU)    = UD(IU) - UZ(IU-2)*OD(IU-2)
      OD(IU-1)  = (OD(IU-1) - UD(IU-1)*RSP(IU-2)) / HD(IU-1)
      S        = 0.0
C
      DO 30 J = IL,IU-2
         S = S - UZ(J)*RSP(J)
   30 CONTINUE
C
      HD(IU)   = HD(IU) + S - UD(IU)*OD(IU-1)
      IF(HD(IU) .EQ. 0.0) THEN
         IFEHL  = N
         RETURN
      END IF
      IFEHL = 0
      RETURN
      END
C$ 
      SUBROUTINE FZYKTL3D (KK,JJ,II,KL,KU,JL,JU,IL,IU,
     $                     UD,HD,OD,UZ,RSP,RS,X)
C***********************************************************************
C
C***********************************************************************
C
      REAL    UD(II),   HD(II),   OD(II),   UZ(II),   RSP(II),
     $        RS(KK,JJ,II),   X (KK,JJ,II)
C
C                                  J J J J J J J J J J J J J J J J J J
C                                  K K K K K K K K K K K K K K K K K K
C                                 VORWAERTSELIMINATION
      DO J = JL, JU
         DO K = KL, KU
C
            RS(K,J,IL)   = RS(K,J,IL) / HD(IL)

         ENDDO
      ENDDO

      DO  I = IL+1,IU-1
         DO J = JL, JU
            DO K = KL, KU
               RS(K,J,I) = (RS(K,J,I) - RS(K,J,I-1) * UD(I)) / HD(I)
            ENDDO
         ENDDO
      ENDDO
C
      DO J = JL, JU
         DO K = KL, KU
            S       = 0.0
            DO I = IL,IU-2
               S     = S - UZ(I) * RS(K,J,I)
            ENDDO
C
            RS(K,J,IU)  = 
     $           (RS(K,J,IU) + S - UD(IU) * RS(K,J,IU-1)) / HD(IU)
         ENDDO
      ENDDO
C
C                                 BERECHNUNG DER LOESUNGEN DURCH
C                                 RUECKWAERTSELIMINATION
C
      DO J = JL, JU
         DO K = KL, KU
            X(K,J,IU)    = RS(K,J,IU)
            X(K,J,IU-1)  = RS(K,J,IU-1) - X(K,J,IU) * OD(IU-1)
         ENDDO
      ENDDO
C     
      DO I = IU-1,IL,-1
         DO J = JL, JU
            DO K = KL, KU
               X(K,J,I)  = RS(K,J,I) - OD(I) * X(K,J,I+1) - 
     $              RSP(I) * X(K,J,IU)
            ENDDO
         ENDDO
      ENDDO
C
C                                  K K K K K K K K K K K K K K K K K K
C                                  J J J J J J J J J J J J J J J J J J
      RETURN
      END
      SUBROUTINE FZYKTL3D_S (KK,JJ,II,KL,KU,JL,JU,IL,IU,
     $                     UD,HD,OD,UZ,RSP,RS,X)
C***********************************************************************
C
C         SKALARE VERSION
C
C***********************************************************************
C
      REAL    UD(II),   HD(II),   OD(II),   UZ(II),   RSP(II),
     $        RS(KK,JJ,II),   X (KK,JJ,II)
C
      REAL SUM(1030)
      IF (KK .GT. 1030) CALL ERRR(501,'FZYKTL3D_S')
C
C                                 VORWAERTSELIMINATION
C
C-------------------------------------------------------
C
      DO J = JL, JU
         DO K = KL, KU
C
            RS(K,J,IL)   = RS(K,J,IL) / HD(IL)

         ENDDO
      ENDDO

      DO  I = IL+1,IU-1
         DO J = JL, JU
            DO K = KL, KU
               RS(K,J,I) = (RS(K,J,I) - RS(K,J,I-1) * UD(I)) / HD(I)
            ENDDO
         ENDDO
      ENDDO
C
      DO J = JL, JU
         DO K = KL, KU
            SUM(K)       = 0.0
         ENDDO
         DO I = IL,IU-2
            DO K = KL, KU
               SUM(K)     = SUM(K) - UZ(I) * RS(K,J,I)
            ENDDO
         ENDDO

         DO K = KL, KU
            S = SUM(K)
            RS(K,J,IU)  = 
     $           (RS(K,J,IU) + S - UD(IU) * RS(K,J,IU-1)) / HD(IU)
         ENDDO
      ENDDO
C
C                                 BERECHNUNG DER LOESUNGEN DURCH
C                                 RUECKWAERTSELIMINATION
C

      DO J = JL, JU
         DO K = KL, KU
            X(K,J,IU)    = RS(K,J,IU)
            X(K,J,IU-1)  = RS(K,J,IU-1) - X(K,J,IU) * OD(IU-1)
         ENDDO
      ENDDO
C     
      DO I = IU-1,IL,-1
         DO J = JL, JU
            DO K = KL, KU
               X(K,J,I)  = RS(K,J,I) - OD(I) * X(K,J,I+1) - 
     $              RSP(I) * X(K,J,IU)
            ENDDO
         ENDDO
      ENDDO
C
C-------------------------------------------------------
C
      RETURN
      END
