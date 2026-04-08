










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
      SUBROUTINE FDERFOUWY_PER (KK,JJ,KSTART,JSTART,KSTOP,JSTOP,
     $           COEFDY,RSGS,DY,DDY,F,FD,LCOL,DIAG,RCOL,UZ,RSP)
C*MGLET***************************************************************
C        F D E R F O U W Y P E R   (First DERive Fourth Order)
C        BERECHNUNG DER ERSTEN ABLEITUNG (DU/DY, DW/DY) BEI 
C        PERIODISCHEN RANDBEDINGUNGEN
C*MGLET***************************************************************
C
C PARAM: F(K,J)         - BEKANNTE GROESSE AN DEN GITTERKANTEN
C      : FIN(K,J)       - INTERPOLIERTE GROESSE IN DER GLEICHEN RICHTUNG
C                         ZWISCHEN DEN BEIDEN KANTEN       
C      : FD (K,J)       - ERSTE ABLEITUNG, FIRST DERINATIVE  DV/DY
C
C      : LCOL           - ENTHAELT KOEFFIZIENTEN DER "LINKEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : DIAG           - ENTHAELT KOEFFIZIENTEN DER HAUPTDIAGONALEN
C      : RCOL           - ENTHAELT KOEFFIZIENTEN DER "RECHTEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C
C DEFINE DIREKTIVEN     : KEINE 
C
C UPROG                 : TRIZYK
C
C VERS:  20.01.97(AM)   : ORIGINAL     (KOMPAKT 4. ORDNUNG)
C
C*MGLET***************************************************************
C

      REAL    DY(JJ),       DDY(JJ),
     $        LCOL(JJ),DIAG(JJ),RCOL(JJ),
     $        UZ(JJ),       RSP(JJ)
C
      REAL    F(KK,JJ),FD(KK,JJ),COEFDY(JSTART:JSTOP+1,12),
     $        RSGS(KK,JJ)            
C
       JJSTOP = JSTOP + 1

       DO 20 J = JSTART, JJSTOP
C
	    LCOL(J) = 1./22.      
	    DIAG(J) = 1.
	    RCOL(J) = 1./22.         
C
       DO 30 K = KSTART, KSTOP
C
            RSGS(K,J) = (24./22.)*(F(K,J)-F(K,J-1))*(DY(J-1))
C
   30  CONTINUE                   
   20  CONTINUE                   
C
           UZ (1) = 1./22.
           RSP(1) = 1./22.
	   
C
       DO 40 K = KSTART, KSTOP
       CALL TRIZYK (K,KK,JJ,JSTART,JJSTOP,
     $              LCOL,DIAG,RCOL,UZ,RSP,RSGS,FD)

   40  CONTINUE
C
C*MGLET***************************************************************
       RETURN
       END

