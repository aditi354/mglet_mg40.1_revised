










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
      SUBROUTINE INTERPOLATE3(KK,JJ,II,KSTART,JSTART,ISTART,
     $                        KSTOP,JSTOP,ISTOP,DX,DY,DZ,DDX,
     $                        DDY,DDZ,F,FINI,FINJ,FINK)
C*MGLET***************************************************************
C        I N T E R P O L A T E 3                               
C*MGLET***************************************************************
C
C PARAM: F(K,J,I)       - ZU INTERPOLIERENDE GROESSE
C      : FINI(K,J,I)    - INTERPOLIERTE GROESSE IN I-RICHTUNG
C      : FINJ(K,J,I)    - INTERPOLIERTE GROESSE IN J-RICHTUNG
C      : FINK(K,J,I)    - INTERPOLIERTE GROESSE IN K-RICHTUNG
C
C DEFINE DIREKTIVEN     : 
C
C UPROG                 :
C
C VERS:  10.10.96 (AM)  : ORIGINAL (ZENTRAL 2. ORDNUNG
C                                  +KOMPAKT 4. ORDNUNG)
C
C*MGLET***************************************************************
C

      REAL       DX(II),        DY(JJ),        DZ(KK),
     $          DDX(II),       DDY(JJ),       DDZ(KK)

      REAL F(KK,JJ,II),FINJ(KK,JJ,II),FINK(KK,JJ,II),
     $     FINI(KK,JJ,II)

CCC                                    HIER ZENTRALE INTERPOLATION !!!

       DO I = ISTART, ISTOP
	 DO J = JSTART, JSTOP
	    DO K = KSTART, KSTOP


       FINJ(K,J,I) = 0.5*(F(K,J,I)   + F(K,J,I+1))
       FINK(K,J,I) = 0.5*(F(K,J,I)   + F(K,J+1,I))
       FINI(K,J,I) = 0.5*(F(K,J,I)*(DZ(K+1)/DDZ(K+1)) 
     $                  + F(K+1,J,I)*(DZ(K)/DDZ(K+1)))


	    ENDDO                   
	 ENDDO                   
       ENDDO                   

CCC                                    HIER KOMPAKTE INTERPOLATION !!!

       RETURN
       END

