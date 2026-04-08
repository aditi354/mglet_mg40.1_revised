










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
      SUBROUTINE INTERPOLATEKK(KK,JJ,II,KMX,JMX,IMX,KSTART,JSTART,
     $                        ISTART,
     $                        KSTOP,JSTOP,ISTOP,DX,DY,DZ,DDX,DDY,DDZ,
     $                        LCOL,DIAG,RCOL,RSIDE1,RSIDE2,RSGS,FAKTOR,
     $                        F,FIN)
C*MGLET***************************************************************
C        I N T E R P O L A T E K K   
C        INTERPOLIERT EINE BELIEBIGE GROESSE IN K BZW. Z-RICHTUNG
C*MGLET***************************************************************
C
C PARAM: F(K,J,I)       - ZU INTERPOLIERENDE GROESSE
C      : FIN(K,J,I)     - INTERPOLIERTE GROESSE IN K(Z)-RICHTUNG
C
C      : LCOL           - ENTHAELT KOEFFIZIENTEN DER "LINKEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : DIAG           - ENTHAELT KOEFFIZIENTEN DER HAUPTDIAGONALEN
C      : RCOL           - ENTHAELT KOEFFIZIENTEN DER "RECHTEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : RSIDE1,RSIDE2  - KOEFFIZIENTEN AUF DER RECHTEN SEITE DES 
C                         GLEICHUNGSSYSTEMS
C      : RSGS           - RECHTE SEITE DES GLEICHUNGSSYSTEMS      
C
C DEFINE DIREKTIVEN     : KEINE 
C
C UPROG                 :  THOMASK    : THOMAS-ALGORITHMUS IN K-RICHTUNG
C                          INTERCOEF1 : BERECHNUNG DER KOEFFIZIENTEN
C                          RANDWERT   : BERECHNUNG DER KOEFFIZIENTEN
C                                       FUER DEN LETZTEN PUNKT      
C                          RANDWERTAN : BERECHNUNG DER KOEFFIZIENTEN
C                                       FUER DEN ERSTEN PUNKT      
C
C VERS:  10.10.96 (AM)  : ORIGINAL        (KOMPAKT 4. ORDNUNG)
C
C*MGLET***************************************************************
C

      REAL       DX(II),        DY(JJ),        DZ(KK),
     $          DDX(II),       DDY(JJ),       DDZ(KK)

      REAL F(KK,JJ,II),FIN(KK,JJ,II),LCOL(KK),FAKTOR(KK),
     $     RCOL(KK),DIAG(KK),RSIDE1(KK),RSIDE2(KK),
     $     RSGS(KK,JJ,II)

C     WRITE (44,*) 'K,LCOL(K),DIAG(K),RCOL(K),RSIDE1(K),RSIDE2(K)'

       DO 30 I = ISTART, ISTOP
 	 DO 20 J = JSTART, JSTOP
	    DO 10 K = KSTART, KSTOP

         IF (K.EQ.KSTART) THEN
CC                                    KOMPAKTER ANSATZ FUER DEN 
CC                                    ERSTEN PHYSIKALISCHEN PUNKT

         CALL RANDWERTAN(KK,K,DZ,DDZ,LCOL(K),DIAG(K),
     $                 RCOL(K),RSIDE1(K),RSIDE2(K))
           RSGS(K,J,I) = RSIDE1(K) * F(K,J,I) + 
     $ 	   RSIDE2(K) * F(K+1,J,I) + 1.0 * F(K+2,J,I)
C          WRITE (44,*) K,LCOL(K),DIAG(K),RCOL(K),
C    $                  RSIDE1(K),RSIDE2(K)
C          WRITE (44,*) 'RSGS(',K,',',J,',',I,')',RSGS(K,J,I)
C
         ELSE IF (K.EQ.KSTOP) THEN
CC                                    KOMPAKTER ANSATZ FUER DEN 
CC                                    LETZTEN PHYSIKALISCHEN PUNKT

         CALL RANDWERT(KK,K,DZ,LCOL(K),DIAG(K),
     $                 RCOL(K),RSIDE1(K),RSIDE2(K))
           RSGS(K,J,I) = RSIDE1(K) * F(K,J,I) + 
     $ 	   RSIDE2(K) * F(K-1,J,I) + 1.0 * F(K-2,J,I)
C          WRITE (44,*) K,LCOL(K),DIAG(K),RCOL(K),
C    $                  RSIDE1(K),RSIDE2(K)
C          WRITE (44,*) 'RSGS(',K,',',J,',',I,')',RSGS(K,J,I)

         ELSE
CC                                    KOMPAKTER ANSATZ IM GEBIET

         CALL INTERCOEF1(KK,K,DZ,DDZ,LCOL(K),DIAG(K),
     $                 RCOL(K),RSIDE1(K),RSIDE2(K))
           RSGS(K,J,I) = RSIDE1(K)*F(K+1,J,I) + F(K,J,I)
C          WRITE (44,*) K,LCOL(K),DIAG(K),RCOL(K),
C    $                  RSIDE1(K),RSIDE2(K)
C          WRITE (44,*) 'RSGS(',K,',',J,',',I,')',RSGS(K,J,I)

         ENDIF

   10      CONTINUE                
   20    CONTINUE                
   30  CONTINUE                

C      WRITE (44,*) 'THOMAS WIRD JETZT AUFGERUFEN'        
       CALL THOMASK(KK,JJ,II,KSTART,KSTOP,JSTART,JSTOP,ISTART,ISTOP,
     $               LCOL,DIAG,RCOL,RSGS,FAKTOR,FIN)
C      WRITE (44,*) 'THOMAS WURDE ERFOLGREICH ABGESCHLOSSEN !!!'

       RETURN
       END

