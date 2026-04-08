










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
      SUBROUTINE TAUINC  (KK,JJ,II,
     $                    DDX,DDY,DDZ,
     $                    U,V,W,FINC,
     $                    TAU11,TAU12,TAU13,
     $                    TAU22,TAU23,TAU33)
C***MGLET***************************************************************
C        T A U I N C       BERECHNUNG DES SGS-LES
C                         SPANNUNGSTENSORS
C***MGLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        X, Y, Z        - KOORDINATEN DER ZELLDEFINITIONSPUNKTE
C        DDX,DDY,DDZ    - KANTENLAENGEN DER KONTROLLVOLUMINA
C        U(KK,JJ,II)    - GESCHWINDIGKEITSFELD
C        V(KK,JJ,II)    - GESCHWINDIGKEITSFELD
C        W(KK,JJ,II)    - GESCHWINDIGKEITSFELD
C        TAU11..TAU33   - ELEMENTE DES SPANNUNGSTENSORS
C
C UPROG                 : ERRR,
C
C VERS:  13.04.99 (MM)  : ORIGINAL
C        23.04.99 (CB)  : MODIF
C        14.03.00 (CB)  : MODIF
C
C
C*STARLET***************************************************************
C
C
C
      REAL        DDX(II),       DDY(JJ),       DDZ(KK),
     $            U(KK,JJ,II),   V(KK,JJ,II),   W(KK,JJ,II)
C
      REAL 	TAU11(KK,JJ,II), TAU12(KK,JJ,II), TAU13(KK,JJ,II),
     $                           TAU22(KK,JJ,II), TAU23(KK,JJ,II),
     $                           FINC(KK,JJ,II),  TAU33(KK,JJ,II)
C
C
      integer KK,JJ,II, k,j,i
C
C    


!!!!!!   INITIALISATION    
      do i = 1,ii
      do j = 1,jj
      do k = 1,kk
      TAU11(k,j,i) = 0.0
      TAU12(k,j,i) = 0.0
      TAU13(k,j,i) = 0.0
      TAU22(k,j,i) = 0.0
      TAU23(k,j,i) = 0.0
      TAU33(k,j,i) = 0.0
      FINC (k,j,i) = 0.0
      end do
      end do
      end do
c     return      !!!!!    POUR  "NO MODEL"  ENLEVER LE COMMENTAIRE
      CALL LEON (KK,JJ,II,U,V,W,DDX,DDY,DDZ,FINC)
      CALL INCH (KK,JJ,II,U,V,W,TAU11,TAU22,TAU33)
      CALL DYN (KK,JJ,II,TAU11,TAU22,TAU33,FINC)
      CALL INC (KK,JJ,II,U,V,W,TAU11,TAU22,TAU33)
      CALL MODEL (KK,JJ,II,FINC,TAU11,TAU22,TAU12)
      CALL MODEL (KK,JJ,II,FINC,TAU11,TAU33,TAU13)
      CALL MODEL (KK,JJ,II,FINC,TAU22,TAU33,TAU23)
      CALL MODEL (KK,JJ,II,FINC,TAU11,TAU11,TAU11)
      CALL MODEL (KK,JJ,II,FINC,TAU22,TAU22,TAU22)
      CALL MODEL (KK,JJ,II,FINC,TAU33,TAU33,TAU33)
C
      RETURN
      END


