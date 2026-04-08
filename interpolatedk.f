










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
      SUBROUTINE INTERPOLATEDK(KK,JJ,II,KMX,JMX,IMX,KSTART,JSTART,
     $			       ISTART,
     $                        KSTOP,JSTOP,ISTOP,DX,DY,DZ,DDX,DDY,DDZ,
     $                        LCOL,DIAG,RCOL,RSIDE1,RSIDE2,
     $                        RSGS,FAKTOR,F,FIN,FD)
C*MGLET***************************************************************
C        I N T E R P O L A T E D K     
C        BERECHNUNG DER ERSTEN ABLEITUNG IN K BZW. Z-RICHTUNG
C*MGLET***************************************************************
C
C PARAM: F(K,J,I)       - BEKANNTE GROESSE AN DEN GITTERKANTEN
C      : FIN(K,J,I)     - INTERPOLIERTE GROESSE IN I(X)-RICHTUNG
C                         ZWISCHEN DEN BEIDEN KANTEN aus interpolatek)      
C      : FD (K,J,I)     - ERSTE ABLEITUNG, First Dirivative
C
C      : LCOL           - ENTHAELT KOEFFIZIENTEN DER "LINKEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : DIAG           - ENTHAELT KOEFFIZIENTEN DER HAUPTDIAGONALEN
C      : RCOL           - ENTHAELT KOEFFIZIENTEN DER "RECHTEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : RSIDE1,RSIDE2  - KOEFFIZIENTEN AUF DER RECHTEN SEITE DES 
C                         GLEICHUNGSSYSTEMS
C
C DEFINE DIREKTIVEN     : KEINE 
C
C UPROG                 :  THOMASK    : THOMAS-ALGORITHMUS IN I-RICHTUNG
C                          INTERCOEF3 : BERECHNUNG DER KOEFFIZIENTEN
C                                       VIERTER ORDNUNG
C                          INTERCOEF2 : BERECHNUNG DER KOEFFIZIENTEN
C                                       DRITTER ORDNUNG
C
C VERS: (18.10.96 AM)  : ORIGINAL        (KOMPAKT 4. ORDNUNG)
C
C*MGLET***************************************************************
C

      REAL       DX(II),        DY(JJ),        DZ(KK),
     $          DDX(II),       DDY(JJ),       DDZ(KK)

      REAL F(KK,JJ,II),FIN(KK,JJ,II),RSGS(KK,JJ,II),
     $     FD(KK,JJ,II)
C
      REAL LCOL(KK),RCOL(KK),DIAG(KK),FAKTOR(KK),
     $     RSIDE1(KK),RSIDE2(KK)


       DO I = ISTART, ISTOP
	 DO J = JSTART, JSTOP
	    DO K = KSTART, KSTOP

            IF (K.EQ.KSTART) THEN
            CALL RANDWERTAND (KK,K,DZ,DDZ,LCOL(K),DIAG(K),
     $	    RCOL(K),RSIDE1(K),RSIDE2(K))
 
            RSGS(K,J,I) = RSIDE1(K)*F(K,J,I) + RSIDE2(K)*
     $			F(K+1,J,I) + F(K+2,J,I)
C           WRITE (44,*) K,LCOL(K),DIAG(K),RCOL(K),
C    $                  RSIDE1(K),RSIDE2(K)
C
            ELSEIF (K.EQ.KSTOP)THEN
 
 	    CALL RANDWERTD (KK,K,DZ,DDZ,LCOL(K),DIAG(K),
     $      RCOL(K),RSIDE1(K),RSIDE2(K))
C
C
            RSGS(K,J,I) = RSIDE1(K)*F(K,J,I) + RSIDE2(K)*
     $			F(K-1,J,I) + F(K-2,J,I)
C           WRITE (44,*) K,LCOL(K),DIAG(K),RCOL(K),
C    $                  RSIDE1(K),RSIDE2(K)
            ELSE
C
C                                BERECHNUNG VIERTER ORDNUNG !!
C
           CALL INTERCOEF3(KK,K,DZ,DDZ,LCOL(K),DIAG(K),
     $     RCOL(K),RSIDE1(K),RSIDE2(K))


           RSGS(K,J,I) = RSIDE1(K)*F(K+1,J,I) + RSIDE2(K)*
     $			 FIN(K,J,I) + F(K,J,I)
C          WRITE (44,*) K,LCOL(K),DIAG(K),RCOL(K),
C    $                  RSIDE1(K),RSIDE2(K)


            ENDIF
C
	    ENDDO                   
	 ENDDO                   
       ENDDO                   


CC     WRITE (44,*) 'THOMAS WIRD IN INTERPOLATED  AUFGERUFEN'        
       CALL THOMASK (KK,JJ,II,KSTART,KSTOP,JSTART,JSTOP,ISTART,ISTOP,
     $               LCOL,DIAG,RCOL,RSGS,FAKTOR,FD)
CC     WRITE (44,*) 'THOMAS WURDE IN INTERPOLATED ERFOLGREICH BEENDET '     

       RETURN
       END

