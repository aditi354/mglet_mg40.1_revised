










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
      SUBROUTINE WRITE3DMPI (KK,JJ,II,X,Y,Z,PHI1,CHANNEL,
     $                         SLICE,APPEND,ITTOT)

C--MGLET----------------------------------------------------------------
C
C                  SCHREIBT SCHEIBE AUS 3D-FELD AUF KANAL
C                  ZU DEBUGGZWECKEN
C
C        18. 9.03 (FS)  : ORIGINAL
C
C--MGLET----------------------------------------------------------------
C
C
      INTEGER KK,JJ,II,SLICE,CHANNEL,APPEND,ITTOT
      REAL
     $        PHI1(KK,JJ,II),X(II),Y(JJ),Z(KK)
      CHARACTER (LEN=6) NUM,PREFIX
      CHARACTER (LEN=12) FILEOUT
C

      PREFIX="ittot_"
      WRITE(NUM(1:6),'(i6)') ITTOT
      DO I = 1,6
         FILEOUT(I:I) = PREFIX(I:I)
      END DO
      DO I = 1,6
         IF (NUM(I:I).NE.' ') THEN
         FILEOUT(6+I:6+I) = NUM(I:I)
         ELSE
         FILEOUT(6+I:6+I) = "0"
         END IF
      END DO 

      IF(APPEND.EQ.1) THEN
                 OPEN(UNIT=CHANNEL,FILE=FILEOUT,FORM="FORMATTED",
     $                        POSITION="REWIND",STATUS="REPLACE")
      ELSE
                 OPEN(UNIT=CHANNEL,FILE=FILEOUT,FORM="FORMATTED",
     $                                          ACCESS="APPEND")
      ENDIF

c      OPEN(CHANNEL)
c      write(*,*) ITTOT
c      DO K = 1,KK
       K = int(KK/2)
         DO I = 1,II
c            WRITE(CHANNEL,*) K,I,PHI1(K,SLICE,I)
            WRITE(CHANNEL,*) X(I),PHI1(K,SLICE,I)
         ENDDO
         WRITE(CHANNEL,*)
c      ENDDO
      CLOSE(CHANNEL)

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C

      RETURN
      END

      SUBROUTINE WRITE3DDIAGX (KK,JJ,II,PHI1,CHANNEL,SLICE)

C--MGLET----------------------------------------------------------------
C
C                  SCHREIBT SCHEIBE AUS 3D-FELD AUF KANAL
C                  ZU DEBUGGZWECKEN
C
C        18. 9.03 (FS)  : ORIGINAL
C
C--MGLET----------------------------------------------------------------
C
C
      INTEGER KK,JJ,II,SLICE,CHANNEL
      REAL
     $        PHI1(KK,JJ,II)
C
      DO K = 1,KK
       DO J = 1,JJ 
          WRITE(CHANNEL,*)K,J,PHI1(K,J,SLICE)
         ENDDO
         WRITE(CHANNEL,*)
      ENDDO

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C

      RETURN
      END
      SUBROUTINE WRITE3DDIAGZ (KK,JJ,II,PHI1,CHANNEL,SLICE)

C--MGLET----------------------------------------------------------------
C
C                  SCHREIBT SCHEIBE AUS 3D-FELD AUF KANAL
C                  ZU DEBUGGZWECKEN
C
C        18. 9.03 (FS)  : ORIGINAL
C
C--MGLET----------------------------------------------------------------
C
C
      INTEGER KK,JJ,II,SLICE,CHANNEL
      REAL
     $        PHI1(KK,JJ,II)
C
      DO I = 1,II
       DO J = 1,JJ
          WRITE(CHANNEL,*)I,J,PHI1(SLICE,J,I)
         ENDDO
         WRITE(CHANNEL,*)
      ENDDO

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C

      RETURN
      END


      SUBROUTINE WRITE3DDIAGY (KK,JJ,II,PHI1,CHANNEL,SLICE)

      INTEGER KK,JJ,II,SLICE,CHANNEL
      REAL
     $        PHI1(KK,JJ,II)
C
       OPEN(CHANNEL)
      DO K = 1,KK
       DO I = 1,II
            WRITE(CHANNEL,*)K,I,PHI1(K,SLICE,I)
         ENDDO
         WRITE(CHANNEL,*)
      ENDDO
       CLOSE(CHANNEL)

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      RETURN
      END


      SUBROUTINE WRITE1D_TMIX(KK,JJ,II,T,CH,
     $                         TIMEPH)
      INTEGER I,CH
      INTEGER KK,JJ,II
      REAL    T(KK,JJ,II),X(II),TIMEPH

       
      WRITE(CH,ERR=2000)TIMEPH,T(50,76,10),T(50,76,19)
     $   ,T(50,76,28),T(50,76,37)
     $   ,T(50,76,46),T(50,76,55)
     $   ,T(50,76,73)
     $   ,T(50,76,91)
C     $   ,T(100,150,180),T(100,150,229)

      RETURN

 1000 FORMAT (9(1X,E12.5E3))
 2000 CALL ERRR (501,' WRITE1D_TMIX')
      END
CCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72


      SUBROUTINE WRITEKOORD(KK,JJ,II,X,CH)

      INTEGER I,KK,JJ,II
      REAL X(II)
      INTEGER CH

      WRITE(CH,ERR=2000)X(10),X(19),X(28),X(37),X(46),X(55),X(73)
C     $      X(91),X(137),X(180),X(229)


      RETURN

 1000 FORMAT (I5,I5, 6(1X,E12.5E3))
 2000 CALL ERRR (501,' WRITEKOORD')

      END
c---------------------------------------------------------------------72

      SUBROUTINE WRITEORRSOMER(KK,JJ,II,X,DX,Z,DZ,U,SLICE,CHANNEL)

      IMPLICIT NONE
      INTEGER I,J,K,KK,JJ,II,SLICE,CHANNEL
      REAL U(KK,JJ,II),Z(KK),DZ(KK),X(II),DX(II)
      REAL x_p,x_m,z_p,z_m,u_fluk

      DO I=3,II-2
      DO K=3,KK-2
       x_p = X(I) + 0.5*DX(I)
       x_m = X(I) - 0.5*DX(I-1)
       z_p = Z(K) + 0.5*DZ(K)
       z_m = Z(K) - 0.5*DZ(K-1)
       U(K,SLICE,I) = U(K,SLICE,I) +
     $     ((1.0/3.0*z_p**3-z_p**2) - (1.0/3.0*z_m**3-z_m**2))/
     $                (z_p - z_m )
      ENDDO
      ENDDO
      DO I=3,II-2
      DO K=3,KK-2
       write(CHANNEL,*)I,Z(K),U(K,SLICE,I)
      ENDDO
      write(CHANNEL,*)
      ENDDO
      RETURN

 1000 FORMAT (I5,I5, 6(1X,E12.5E3))
 2000 CALL ERRR (501,' WRITEKOORD')

      END

