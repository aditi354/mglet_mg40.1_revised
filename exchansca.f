










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
C#ifdef _TSCAL_
      SUBROUTINE EXCHANSCA  (KK,JJ,II,KMX,JMX,IMX,T,T0,IEX)

C*STARLET***************************************************************
C        E X C H A N S C A     AUSTAUSCH DER SPEICHERINHALTE (IEX = 1)
C                              T0 = T                        (IEX = 2)
C*STARLET***************************************************************
C
C VERS:  27.11.02 (TB)  : ORIGINAL FOR SCALAR TRANSPORT
C
C UPROG                 : ERRR
C
C*STARLET***************************************************************
C
      REAL   T(KK,JJ,II), T0(KK,JJ,II)
C
      IF(IEX .NE. 1) GOTO 2100
C
      DO 10 I=1,IMX
         DO 20 J=1,JMX
            DO 30 K=1,KMX
               T1        = T(K,J,I)
               T(K,J,I)  = T0(K,J,I)
               T0(K,J,I) = T1
   30       CONTINUE   
   20    CONTINUE
   10 CONTINUE
C
C      T(3    ,1,1) = T0(3    ,1,1)
C      T(KMX-2,1,1) = T0(KMX-2,1,1)
C
      RETURN
C
 2100 IF(IEX .NE. 2) CALL ERRR(501,' EXCHANSCA')
C
      DO 100 I=1,IMX
         DO 110 J=1,JMX
            DO 120 K=1,KMX
               T0(K,J,I) = T(K,J,I)
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C#endif
