










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
      SUBROUTINE FLUCTUATIONS (KK,JJ,II,
     $                   U,V,W,H3D1,H3D2,H3D3,
     $                   UOLD,VOLD,WOLD,DREAD,DZ,Z)
      
C
C
      INTEGER      I,J,K,KK,JJ,II
      LOGICAL      DREAD
      REAL         U(KK,JJ,II),   V(KK,JJ,II),   W(KK,JJ,II)
      REAL         UOLD(KK,JJ,II),   VOLD(KK,JJ,II),   WOLD(KK,JJ,II)
      REAL         H3D1(KK,JJ,II),   H3D2(KK,JJ,II),   H3D3(KK,JJ,II)
      REAL         DZ(KK),Z(KK)
C
C
C
      IF (DREAD) THEN
      DO I=1,II
       DO J=1,JJ
        DO K=1,KK
         H3D1(K,J,I) = -UOLD(K,J,I) + U(K,J,I) 
C-UOLD(K,J,I) + U(K,J,I)
         H3D2(K,J,I) = 0.0
C-VOLD(K,J,I) + V(K,J,I)
         H3D3(K,J,I) = W(K,J,I)
C-WOLD(K,J,I) + W(K,J,I)
        ENDDO
       ENDDO
      ENDDO
      ELSE
      DO I=2,II-2
       DO J=2,JJ-2
        DO K=2,KK-2
              H3D1(K,J,I) = U(K,J,I) +
     $      ((1.0/3.0*(Z(K)+0.5*DZ(K))**3-(Z(K)+0.5*DZ(K))**2) -
     $     (1.0/3.0*(Z(K)-0.5*DZ(K-1))**3-(Z(K)-0.5*DZ(K-1))**2))/
     $                (0.5*(DZ(K)+DZ(K-1)))
         H3D2(K,J,I) =  V(K,J,I)
         H3D3(K,J,I) =  W(K,J,I)
        ENDDO
       ENDDO
      ENDDO
      ENDIF

      RETURN
      END
