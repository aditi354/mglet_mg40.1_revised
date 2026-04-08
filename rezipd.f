










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
      SUBROUTINE REZIPD(KK,JJ,II,KMX,JMX,IMX,
     $                  DX,DY,DZ,DDX,DDY,DDZ,
     $                  RDX,RDY,RDZ,RDDX,RDDY,RDDZ)
C*STARLET***************************************************************
C  R E Z I P D    BERECHNET REZIPROKWERTE DER GEOMETRIEFAKTOREN
C*STARLET***************************************************************
C
C  PARAMETER KK, JJ, II           - ARRAYGRENZEN
C            KMX,JMX,IMX          - GRENZE D. BER.-GEB.(MIT BOUND)
C            DX(I),DY(J),DZ(K)    - ABSTAND DER GITTERPUNKTE
C            DDX(I),DDY(J),DDZ(K) - KANTENLAENGE DER KONTROLLVOLUMINA
C            RDX(I),RDY(J),RDZ(K)    - REZIPROKWERTE
C            RDDX(I),RDDY(J),RDDZ(K) - REZIPROKWERTE
C
C  VERS:  10. 6.92 (MM)  : ORIGINAL
C
C
C*STAR******************************************************************
C
C
C
      IMPLICIT NONE
      INTEGER II,JJ,KK,KMX,JMX,IMX,I,J,K
      REAL    DX(II),DY(JJ),DZ(KK),DDX(II),DDY(JJ),DDZ(KK),
     $        RDX(II),RDY(JJ),RDZ(KK),RDDX(II),RDDY(JJ),RDDZ(KK)


       DO I=1,IMX
          RDX(I) = 1.0/DX(I)
          RDDX(I) = 1.0/DDX(I)
       ENDDO

       DO 200 J=1,JMX

          RDY(J) = 1./DY(J)
          RDDY(J) = 1./DDY(J)

  200   CONTINUE

       DO 300 K=1,KMX

          RDZ(K) = 1./DZ(K)
          RDDZ(K) = 1./DDZ(K)

  300   CONTINUE

       RETURN
       END
