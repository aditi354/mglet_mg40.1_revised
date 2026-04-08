










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
      SUBROUTINE EXCHAN  (KK,JJ,II,KMX,JMX,IMX,U,V,W,UO,VO,WO,IEX)
C*STARLET***************************************************************
C        E X C H A N      AUSTAUSCH DER SPEICHERINHALTE (IEX = 1)
C                         UO = U, VO = V, WO = W        (IEX = 2)
C*STARLET***************************************************************
C
C VERS:  23.08.85 (HW)  : ORIGINAL
C
C UPROG                 : ERRR
C
C*STARLET***************************************************************
C
      REAL   U(KK,JJ,II),   V(KK,JJ,II),   W(KK,JJ,II),
     $       UO(KK,JJ,II),  VO(KK,JJ,II),  WO(KK,JJ,II)
C
      IF(IEX .NE. 1) GOTO 2100
C
      DO 10 I=1,IMX
         DO 20 J=1,JMX
            DO 30 K=1,KMX
               U1        = U(K,J,I)
               V1        = V(K,J,I)
               W1        = W(K,J,I)
C
               U(K,J,I)  = UO(K,J,I)
               V(K,J,I)  = VO(K,J,I)
               W(K,J,I)  = WO(K,J,I)
               UO(K,J,I) = U1
               VO(K,J,I) = V1
   30          WO(K,J,I) = W1
   20    CONTINUE
   10 CONTINUE
C
      U(3    ,1,1) = UO(3    ,1,1)
      U(KMX-2,1,1) = UO(KMX-2,1,1)
C
      RETURN
C
 2100 IF(IEX .NE. 2) CALL ERRR(501,' EXCHAN   ')
C
      DO 100 I=1,IMX
         DO 110 J=1,JMX
            DO 120 K=1,KMX
               UO(K,J,I) = U(K,J,I)
               VO(K,J,I) = V(K,J,I)
  120          WO(K,J,I) = W(K,J,I)
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
