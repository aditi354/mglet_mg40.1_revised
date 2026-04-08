










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
      SUBROUTINE SETIDX(KK,JJ,II,KMX,JMX,IMX,B,JUMP,NBND,MAXCBP,
     1                  ISTA,IEND,JSTA,JEND,KSTA,KEND,
     2                  LALFAX,LGAMAX,LALFAY,LGAMAY,LALFAZ,LGAMAZ,
     3                  IDXAKT,ICUBPT,JCUBPT,KCUBPT,IDXIGD,LCUB)

      INTEGER II,JJ,KK,NBND,I,J,K,JUMP,IMX,JMX,KMX,IDXAKT
      INTEGER ISTA,IEND,JSTA,JEND,KSTA,KEND
      INTEGER ICUBPT(MAXCBP),JCUBPT(MAXCBP),KCUBPT(MAXCBP)
      INTEGER IDXIGD
      REAL    LCUB(8,MAXCBP)
      REAL B(KK,JJ,II)
      REAL LALFAX(II),LGAMAX(II),LALFAY(JJ),LGAMAY(JJ),
     1     LALFAZ(KK),LGAMAZ(KK)

C      WRITE (*,*) IDXAKT,ISTA,IEND,JSTA,JEND,KSTA,KEND
      DO I = ISTA,IEND,2*JUMP
        DO J = JSTA,JEND,2*JUMP
	  DO K = KSTA,KEND,2*JUMP

            FIP =  FSENSI(KK,JJ,II,IMX,K,J,I,NBND, 1,JUMP,B)
            FIM =  FSENSI(KK,JJ,II,IMX,K,J,I,NBND,-1,JUMP,B)
            FJP =  FSENSJ(KK,JJ,II,JMX,K,J,I,NBND, 1,JUMP,B)
            FJM =  FSENSJ(KK,JJ,II,JMX,K,J,I,NBND,-1,JUMP,B)
            FKP =  FSENSK(KK,JJ,II,KMX,K,J,I,NBND, 1,JUMP,B)
            FKM =  FSENSK(KK,JJ,II,KMX,K,J,I,NBND,-1,JUMP,B)

            IF ((FIP+FIM+FJP+FJM+FKP+FKM) .LT. 5.5) THEN
              IDXAKT = IDXAKT + 1
	      ICUBPT(IDXAKT) = I
	      JCUBPT(IDXAKT) = J
	      KCUBPT(IDXAKT) = K
	      IF ((FIP + FIM) .LT. 0.5) FIP = 1.0
	      IF ((FJP + FJM) .LT. 0.5) FJP = 1.0
	      IF ((FKP + FKM) .LT. 0.5) FKP = 1.0
	      LCUB(1,IDXAKT) = FIM * LALFAX(I)
	      LCUB(2,IDXAKT) = FIP * LGAMAX(I)
	      LCUB(3,IDXAKT) = FJM * LALFAY(J)
	      LCUB(4,IDXAKT) = FJP * LGAMAY(J)
	      LCUB(5,IDXAKT) = FKM * LALFAZ(K)
	      LCUB(6,IDXAKT) = FKP * LGAMAZ(K)
	    ENDIF
	  ENDDO
	ENDDO
      ENDDO

      IDXIGD = IDXAKT

      RETURN
      END



