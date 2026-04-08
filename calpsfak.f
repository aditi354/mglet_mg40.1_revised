










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
      SUBROUTINE CALPSFAK (KK,JJ,II,NBND,
     $                      DX, DY, DZ,
     $                     DDX,DDY,DDZ,
     $                     FPSFAK,FCOSMY,FCOSNY,
     &                     NXGRAE,NYGRAE,NZGRAE,
     $                     NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,IGRID)

C**MGLET*********************************************************
C
C      C A L P S F A K        BERECHNET DIE FAKTOREN, DIE
C                             IN DER DISKRETEN POISSONGLEICHUNGEN
C                             STEHEN
C                             UND DURCH VOLU/AX DIVIDIERT
C
C      DIV/DT  *  VOLU/AX*DDX  =
C
C               P(I+1)     -    2*P(I)      +     P(I-1)
C
C      FAK(2)*( P(J+1)     -    2*P(J)      +     P(J-1) )
C
C      FAK(3)*  P(K+1)     - FAK(4)*P(K)    + FAK(5)*P(K-1)
C
C                             FOURIERTRANSFORMIERT IN X- UND Y-RICHTUNG
C
C      FT(DIV)/DT  * VOLU/AX*DDX =
C
C                                             FAK(5)*FT(P(K-1))
C
C         (FAK(4) + FAK(2))FCOSNY + FCOSMY) *        FT(P(K))
C
C                                             FAK(3)*FT(P(K+1))
C
C      
C     FPSFAK(K,1):    VOLU/AX*DDX
C     FPSFAK(K,2):    FAK(2) = VOLU/AX*DDX * AY/VOLV/DDY
C
C     FPSFAK(K,3):    FAK(3) = VOLU/AX*DDX *   AZ/VOLZ(K)/DDZ(K)
C     FPSFAK(K,4):    FAK(4) = VOLU/AX*DDX * 
C                              -(AZ/VOLZ(K)/DDZ(K)+AZ/VOLZ(K-1)/DDZ(K-1))
C     FPSFAK(K,5):    FAK(5) = VOLU/AX*DDX *   AZ/VOLZ(K-1)/DDZ(K-1)
C

C
C     09.03.1994 (MM):   ORIGINAL 
C                        UNTER KRAEFTIGER HILFE VON FRIEDEMANN UNGER
C
C*MGLET************************************************************
C
C
      REAL
     $      DX(II), DY(JJ), DZ(KK),
     $     DDX(II),DDY(JJ),DDZ(KK)

      REAL
     $     FPSFAK(KK,5),FCOSMY(II),FCOSNY(JJ)
C
      DATA  PI /3.141592653589/
C
C---- ----------------------------  IST FOURIERTANSF. UEBERHAUPT
C                                   MOEGLICH
C
      IF (NXGRAE .NE. 1) CALL ERRR (501,' CALPSFAK')
      IF (NYGRAE .NE. 1) CALL ERRR (502,' CALPSFAK')
C
C---- ----------------------------  PERIODIC BOUNDARY-CONDITIONS?
C
      IF (NFRO .NE. 1 ) CALL ERRR (503,' CALPSFAK')
      IF (NBAC .NE. 1 ) CALL ERRR (504,' CALPSFAK')
      IF (NRGT .NE. 1 ) CALL ERRR (505,' CALPSFAK')
      IF (NLFT .NE. 1 ) CALL ERRR (506,' CALPSFAK')
C
C---- ----------------------------  VORBELEGUNG
C

      DO J = 1,5
      DO K = 1,KK

         FPSFAK(K,J) = 0.0

      ENDDO
      ENDDO
C
      DO I=1,II

         FCOSMY(I) = 0.0

      ENDDO

      DO J=1,JJ

         FCOSNY(J) = 0.0

      ENDDO


C
C---- -----------------------------------------   VOLU/AX*DDX
C

      DO K=1,KK

         FPSFAK(K,1) = DX(3) * DDX(3)

      ENDDO

C
C---- ------------------------------------   VOLU/AX*DDX * AY/VOLV/DDY
C

      DO K=1,KK

         FPSFAK(K,2) = DX(3)*DDX(3)/(DY(3)*DDY(3))
      
      ENDDO

C
C---- ------------------------------ VOLU/AX*DDX *   AZ/VOLZ(K)/DDZ(K)
C

      DO K=1,KK

         FPSFAK(K,3) = DX(3)*DDX(3) / (DZ(K)*DDZ(K))

      ENDDO

C
C---- ----------------------------- VOLU/AX*DDX * 
C                          -(AZ/VOLZ(K)/DDZ(K) + AZ/VOLZ(K-1)/DDZ(K-1))
C
      DO K=2,KK

         FPSFAK(K,4) = DX(3)*DDX(3) * (-1.0) * 
     $                 ( 1.0/(DZ(K)*DDZ(K)) + 1.0/(DZ(K-1)*DDZ(K)))

      ENDDO

C
C---- ----------------------------  VOLU/AX*DDX *   AZ/VOLZ(K-1)/DDZ(K)
C

      DO K=2,KK

         FPSFAK(K,5) = DX(3)*DDX(3) / (DZ(K-1)*DDZ(K))

      ENDDO

C
C---- -----------------------------------  BOUNDARY CONDITIONS
C
      IF ((NBOT .EQ. 5) .OR. (NBOT .EQ. 6)) THEN
C
C---- -----------------------------------  NEUMANN FOR PRESSURE  ON BOTTOM
C                                          DRUECKE AUSSERHALB WERDEN
C                                          GLEICH DEN DRUECKEN INNERHALB
C                                          GESETZT
C
      FPSFAK(   NBND+1,4) = FPSFAK(   NBND+1,4) + FPSFAK(   NBND+1,5)
      FPSFAK(   NBND+1,5) = 0.0

      ELSE 
C
C---- -----------------------------------  DIRICHLET FOR PRESSURE 
C                                          PRESSURE IN FIRST GRID POINT
C                                          OUT OF THE AREA IS SET TO ZERO
      FPSFAK(   NBND+1,5) = 0.0
      ENDIF

      IF ((NTOP .EQ. 5) .OR. (NTOP .EQ. 6)) THEN
C
C---- -----------------------------------  NEUMANN FOR PRESSURE  ON TOP

      FPSFAK(KK-NBND  ,4) = FPSFAK(KK-NBND  ,4) + FPSFAK(KK-NBND  ,3)
      FPSFAK(KK-NBND  ,3) = 0.0

      ELSE 
C
C---- -----------------------------------  DIRICHLET FOR PRESSURE 
      FPSFAK(KK-NBND  ,3) = 0.0
      ENDIF


C
C---- -----------------------------------   KONTROLLAUSDRUCK
C

      WRITE (6,*) '***********************************************'
      WRITE (6,*) '*'
      WRITE (6,*) '* CALPSFAK:'
      WRITE (6,*) '*'
      WRITE (6,*) '*'
      WRITE (6,*) 'FPSFAK(1) FPSFAK(2) FPSFAK(3) FPSFAK(4) FPSFAK(5)'
      WRITE (6,*) '*'
      WRITE (6,*) '*'

      DO K=1,KK

         WRITE (6,1000) FPSFAK(K,1),FPSFAK(K,2),FPSFAK(K,3),
     $                  FPSFAK(K,4),FPSFAK(K,5)

      ENDDO

      WRITE (6,*) ' '
      WRITE (6,*) ' '
 1000 FORMAT (6(F10.5))

C
C---- ------------  FAKTOREN FUER FOURIER-TRANSFORIERTE POISSONGLEICHUNG
C
C
C---- ------------- FCOSMY:   TRANSF. IN X-RICHTUNG
C

      TERM1=2.*PI
      TERM2=TERM1/FLOAT(II-2*NBND)

      DO I=1,II-2*NBND

         TERM3 = FLOAT(I-1)

         FCOSMY(I) = 2.0 * (COS(TERM2*TERM3)-1.0)

      ENDDO
C
C---- ------------- FCOSNY:   TRANSF. IN Y-RICHTUNG
C

      TERM1=2.*PI
      TERM2=TERM1/FLOAT(JJ-2*NBND)

      DO J=1,JJ-2*NBND

         TERM3 = FLOAT(J-1)

         FCOSNY(J) = 2.0 * (COS(TERM2*TERM3)-1.0)

      ENDDO

C
C---- ------------- END OF CALPSFAK
C
      RETURN
      END
