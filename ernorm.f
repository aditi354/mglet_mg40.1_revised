










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
	SUBROUTINE ERNORM(KK,JJ,II,KMX,JMX,IMX,NBND,PHI1,PHI2,EPSPHI)
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C               LIEFERT DIE FEHLERNORM
C        VERSION:  22. 4.92 (MM): ORIGINAL
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

	REAL 
     $       PHI1(KK,JJ,II)
     $      ,PHI2(KK,JJ,II)

        RNPT = 1./FLOAT((IMX-2*NBND)*(JMX-2*NBND)*(KMX-2*NBND))

        AVPHI2 = 0.0
        SIGPHI = 0.0
        EPSPHI = 0.0

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                      MITTELWERT VON PHI2

       DO 14 I=NBND+1,IMX-NBND
       DO 14 J=NBND+1,JMX-NBND
       DO 14 K=NBND+1,KMX-NBND

		AVPHI2 = AVPHI2 + PHI2(K,J,I)

   14    CONTINUE

		AVPHI2 = AVPHI2*RNPT

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                       STANDARDABWEICHUNG VON PHI2

       DO 29 I=NBND+1,IMX-NBND
       DO 29 J=NBND+1,JMX-NBND
       DO 29 K=NBND+1,KMX-NBND

		SIGPHI = SIGPHI + (PHI2(K,J,I) - AVPHI2)**2

   29    CONTINUE

      IF(SIGPHI.EQ.0.0) SIGPHI = 1.0

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                       NORMIERTE ABWEICHUNG PHI1 VON PHI2

        DO 44 I=NBND+1,IMX-NBND
        DO 44 J=NBND+1,JMX-NBND
        DO 44 K=NBND+1,KMX-NBND

		EPSPHI = EPSPHI + (PHI1(K,J,I)-PHI2(K,J,I))**2

   44     CONTINUE

		EPSPHI = EPSPHI/SIGPHI

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

       RETURN
       END

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

