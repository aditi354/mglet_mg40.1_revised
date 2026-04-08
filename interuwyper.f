










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
      SUBROUTINE INTERUWY_PER(KK,JJ,KSTART,JSTART,KSTOP,JSTOP,
     $                        COEFFY,RSGS,F,FIN,LCOL,DIAG,RCOL,UZ,RSP)
C*MGLET***************************************************************
C        I N T E R U W Y P E R
C                  = = =
C        INTERPOLIERT U, W IN Y-RICHTUNG (VINJ, WINJ) PERIODISCHE RB
C*MGLET***************************************************************
C
C PARAM: F(K,J)         - ZU INTERPOLIERENDE GROESSE(U,W)
C      : FIN (K,J)      - INTERPOLIERTE U,W IN J(Y)-RICHTUNG
C
C      : COEFFY(I,7)    - ENTHAELT KOEFFIZIENTEN DER "LINKEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : COEFFY(I,8)    - ENTHAELT KOEFFIZIENTEN DER HAUPTDIAGONALEN
C      : COEFFY(I,9)    - ENTHAELT KOEFFIZIENTEN DER "RECHTEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : COEFFY(I,10),   - KOEFFIZIENTEN AUF DER RECHTEN SEITE DES 
C        COEFFY(I,11),     GLEICHUNGSSYSTEMS
C        COEFFY(I,12) = 0.0
C
C DEFINE DIREKTIVEN     : KEINE 
C
C UPROG                 : TRIZYK
C
C VERS:  07.10.97 (AM)  : ORIGINAL        (KOMPAKT 4. ORDNUNG)
C
C*MGLET***************************************************************
C

      REAL       LCOL(JJ),DIAG(JJ),RCOL(JJ),UZ(JJ),RSP(JJ),
     $           F(KK,JJ),FIN(KK,JJ),COEFFY(JSTART:JSTOP+1,12),
     $           RSGS(KK,JJ)
C
       JJSTOP  = JSTOP + 1
C
       DO 20 J = JSTART, JJSTOP
C
C                                     BILDUNG DER MATRIX  A              
	 LCOL(J) = 1./6.       
	 DIAG(J) = 1.           
	 RCOL(J) = 1./6.        
C
       DO 30 K = KSTART, KSTOP
C                                     RECHTE SEITE DES GLEICHUNGSSYSTEMS
C
         RSGS(K,J) = (2./3.) * (F(K,J) + F(K,J-1))
C
   30  CONTINUE                
   20  CONTINUE                

       UZ (1) = 1./6.
       RSP(1) = 1./6.

C                                     LOESUNG DES GLEICHUNGSSYSTEMS      
       DO 50 K = KSTART, KSTOP
       CALL TRIZYK (K,KK,JJ,JSTART,JJSTOP,
     $              LCOL,DIAG,RCOL,UZ,RSP,RSGS,FIN)
   50  CONTINUE

C
C*MGLET***************************************************************
C
       RETURN
       END

