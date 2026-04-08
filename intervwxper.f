










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
      SUBROUTINE INTERVWX_PER(KK,JJ,II,KSTART,JSTART,ISTART,KSTOP,JSTOP,
     $                   ISTOP,COEFFX,RSGS,F,FINI,LCOL,DIAG,RCOL,UZ,RSP)
C*MGLET***************************************************************
C        I N T E R V W X P E R
C                  = = =
C        INTERPOLIERT F IN X-RICHTUNG (VINI,WINI) PERIODISCHE RB
C*MGLET***************************************************************
C
C PARAM: F(K,J,I)       - ZU INTERPOLIERENDE GROESSE(V,W)
C      : FINI (K,J,I)   - INTERPOLIERTE V und W IN I(X)-RICHTUNG
C
C      : COEFFX(I,1)    - ENTHAELT KOEFFIZIENTEN DER "LINKEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : COEFFX(I,2)    - ENTHAELT KOEFFIZIENTEN DER HAUPTDIAGONALEN
C      : COEFFX(I,3)    - ENTHAELT KOEFFIZIENTEN DER "RECHTEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : COEFFX(I,4),   - KOEFFIZIENTEN AUF DER RECHTEN SEITE DES 
C        COEFFX(I,5),     GLEICHUNGSSYSTEMS
C        COEFFX(I,6)
C
C DEFINIE DIREKTIVEN     : KEINE 
C
C UPROG                 : TRIZYK
C
C VERS:  07.10.97 (AM)  : ORIGINAL        (KOMPAKT 4. ORDNUNG)
C
C*MGLET***************************************************************
C

      REAL       LCOL(II),DIAG(II),RCOL(II),UZ(II),RSP(II),
     $           F(KK,JJ,II),FINI(KK,JJ,II),COEFFX(ISTART:ISTOP+1,12),
     $           RSGS(KK,JJ,II)
C
       DO 30 J = JSTART, JSTOP
       DO 40 K = KSTART, KSTOP
C
C                                     KOMPAKTER ANSATZ IM GEBIET
 	 DO 20 I = ISTART, ISTOP+1
C
         RSGS(K,J,I) = (2./3.) * (F(K,J,I) + F(K,J,I-1))
C
	 LCOL(I) = 1./6.       
	 DIAG(I) = 1.          
	 RCOL(I) = 1./6.       
C
   20    CONTINUE                

       UZ (1) = 1./6.
       RSP(1) = 1./6.

C
       CALL TRIZYK (II,ISTART,ISTOP+1,LCOL,DIAG,RCOL,
     $              RSGS(K,J,1),FINI(K,J,1),UZ,RSP)

   40  CONTINUE                
   30  CONTINUE                
C
C*MGLET***************************************************************
C
       RETURN
       END

