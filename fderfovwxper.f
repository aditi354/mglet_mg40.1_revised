










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
      SUBROUTINE FDERFOVWX_PER(KK,JJ,II,KSTART,JSTART,ISTART,KSTOP,
     $                        JSTOP,
     $                        ISTOP,COEFDX,RSGS,DX,DDX,F,FINI,FD,
     $                        LCOL,DIAG,RCOL,UZ,RSP)
C*MGLET***************************************************************
C        F D E R F O V W X P E R 
C                    = = = 
C        BERECHNUNG DER ERSTEN ABLEITUNG FDVI,FDWI / PERIODISCHE RB
C*MGLET***************************************************************
C
C PARAM: F(K,J,I)         - GESCH.-KOMPONENTE AN DEN KANTEN
C      : FINI(K,J,I)      - INTERPOLIERTE F IN I(X)-RICHTUNG
C      : FD  (K,J,I)      - ERSTE ABLEITUNG IN X-RICHTUNG    
C
C      : COEFDX(I,7)    - ENTHAELT KOEFFIZIENTEN DER "LINKEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : COEFDX(I,8)    - ENTHAELT KOEFFIZIENTEN DER HAUPTDIAGONALEN
C      : COEFDX(I,9)    - ENTHAELT KOEFFIZIENTEN DER "RECHTEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : COEFDX(I,10),  - KOEFFIZIENTEN AUF DER RECHTEN SEITE DES 
C        COEFDX(I,11),    GLEICHUNGSSYSTEMS
C        COEFDX(I,12)
C
C DEFINE DIREKTIVEN     : KEINE 
C
C UPROG                 : TRIZYK
C
C VERS:  14.02.97 (AM)  : ORIGINAL        (KOMPAKT 4. ORDNUNG)
C
C*MGLET***************************************************************
C

      REAL       LCOL(II),DIAG(II),RCOL(II),UZ(II),RSP(II),
     $           DX(II),DDX(II)
      REAL       F(KK,JJ,II),FINI(KK,JJ,II),COEFDX(ISTART:ISTOP+1,12),
     $           RSGS(KK,JJ,II),FD(KK,JJ,II)
C


       DO 30 J = JSTART, JSTOP
       DO 40 K = KSTART, KSTOP
C
C                                     KOMPAKTER ANSATZ IM GEBIET
 	 DO 20 I = ISTART, ISTOP+1
C
         RSGS(K,J,I) = (24./22.) *( F(K,J,I)  - F(K,J,I-1))*
     $                 DX(I-1)
C
	 LCOL(I) = 1./22.      
	 DIAG(I) = 1.          
	 RCOL(I) = 1./22.       
C
   20    CONTINUE                

       UZ (1) = 1./22.
       RSP(1) = 1./22.


C
       CALL TRIZYK (II,ISTART,ISTOP+1,LCOL,DIAG,RCOL,
     $              RSGS(K,J,1),FD(K,J,1),UZ,RSP)

   40  CONTINUE                
   30  CONTINUE                
C
C*MGLET***************************************************************
C
       RETURN
       END

