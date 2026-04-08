










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
      SUBROUTINE SMOOTH_STRESS  (KK,JJ,II,X,Y,Z,
     $                    DX,DY,DZ,DDX,DDY,DDZ,
     $                    RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                    TAU,
     $                    SAMPLES,HELP)
C***MGLET***************************************************************
C        S M O O T H _ S T R E S S
C                         GLAETTUNG DES SPANNUNGSTENSORS NACH 
C                         INTERPOLATION VON PARTIKELPOSITIONEN
C***MGLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        X, Y, Z        - KOORDINATEN DER ZELLDEFINITIONSPUNKTE
C        DX,DY,DZ       - ABSTAND DER GITTERPUNKTE
C        DDX,DDY,DDZ    - KANTENLAENGEN DER KONTROLLVOLUMINA
C        TAU11..TAU33   - ELEMENTE DES SPANNUNGSTENSORS
C        SAMPLES         - ANZAHL DER TREFFER PRO ZELLE (PARTIKEL)
C
C VERS:  02.05.99 (MM)  : ORIGINAL
C
C DEFINE-DIREKTIVEN     : 
C
C*MGLET*****************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
C
      REAL        X(II),         Y(JJ),         Z(KK),
     $           DX(II),        DY(JJ),        DZ(KK),
     $          DDX(II),       DDY(JJ),       DDZ(KK),
     $          RDX(II),       RDY(JJ),       RDZ(KK),
     $         RDDX(II),      RDDY(JJ),      RDDZ(KK)
C
      REAL 	TAU(KK,JJ,II), HELP(KK,JJ,II)
C      REAL 	TAU11(KK,JJ,II), TAU12(KK,JJ,II), TAU13(KK,JJ,II),
C     $          TAU21(KK,JJ,II), TAU22(KK,JJ,II), TAU23(KK,JJ,II),
C     $          TAU31(KK,JJ,II), TAU32(KK,JJ,II), TAU33(KK,JJ,II)
C
      REAL SAMPLES(KK,JJ,II)

      integer KK,JJ,II
C
C      real nverteilung(0:99)

C
C      write (6,*) 'smooth_stress,',
C     $ 'myid, tau...',myid,help(5,5,5),samples(5,5,5)

      NBND = 2

      KM2    = KK-NBND
      JM2    = JJ-NBND
      IM2    = II-NBND

C     
C
      DO I=3,II-2
         DO J=3,JJ-2
            DO K=3,KK-2
               
               RN = 1.0/(SAMPLES(K,J,I) + 0.0000001)

               TAU(K,J,I)= HELP(K,J,I) * RN

            ENDDO
         ENDDO
      ENDDO

C------------------------------------------------------------
C                                    BOUNDARY CONDITIONS
C                                    WE ASSUME PERIODIC IN X- AND Y
C                                    AND NOSLIP IN Z

      
C------------------------------------------------------------
C                       Z-DIRECTION, NOSLIP-WALL --> NEUMANN FOR STRESSES

      DO I=3,II-2
        DO J=3,JJ-2

          TAU(2,J,I) = TAU(3,J,I)

          TAU(KK-1,J,I) = TAU(KK-2,J,I)

        ENDDO
      ENDDO
C------------------------------------------------------------
C                      Y-DIRECTION, PERIODIC

      DO I=3,II-2
        DO K=2,KK-1

          TAU(K,   2,I) = TAU(K,JJ-2,I)
          TAU(K,JJ-1,I) = TAU(K,   3,I)

        ENDDO
      ENDDO

C------------------------------------------------------------
C                     X-DIRECTION, PERIODIC

      DO J=3,JJ-2
        DO K=2,KK-1

          TAU(K,J,   2) = TAU(K,J,II-2)
          TAU(K,J,II-1) = TAU(K,J,   3)

        ENDDO
      ENDDO

C------------------------------------------------------------
C      write (6,*) 'smooth_stress,',
C     $ 'myid, tau...',myid,tau(5,5,5),samples(5,5,5)
C
c      if (myid .eq. 0) then
c      write (6,*) 'smooth_stress'
c      do i=1,ii
c      write (6,1000) 'tau',i,tau(5,5,i),samples(5,5,i)
c      enddo
c
c      endif
 1000 FORMAT (A5,I6,2(1X,E12.5E3,1X))


C------------------------------------------------------------
C      N= 0

C      do i=0,99
C        nverteilung(i) = 0
C      enddo

C      DO I=3,II-2
C       DO J=3,JJ-2
C        DO K=3,KK-2
C
C           N = N + SAMPLES(K,J,I)
C
C           nverteilung(SAMPLES(K,J,I)) = nverteilung(SAMPLES(K,J,I)) + 1
C
C        ENDDO
C       ENDDO
C      ENDDO

C      if (myid .eq. 0) then
C
C      write (6,*)'gesamtanzahl in smooth_stress',n
C      do i=0,25
C         write (6,*) 'nvert',i,nverteilung(i)
C      enddo
C
C      endif

C      DO I=2,II-1
C         DO J=2,JJ-1
C            DO K=2,KK-1
C
C               TAU(K,J,I)= 0.0
C
C            ENDDO
C         ENDDO
C      ENDDO

C------------------------------------------------------------

      RETURN
      END

