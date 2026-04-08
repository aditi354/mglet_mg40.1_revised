










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
C$    LIST(ALL=0)      
      SUBROUTINE COEFFOP32 (II,K,DX,ALFA,BETA,GAMA)              
C*MGLET***************************************************************
C
C     UNTERPROGRAMM ZU BERECHNUNG DER KOEFFIZIENTEN 
C     FUER DIE EXTRAPOLATION FUER DEN AUSSTROEMRAND
C     BEIM ZENTRALDIFFERENZENVERFAHREN  !!!!!
C     GESCHWINDIGKEITSKOMPONENTE V UND W
C
C*MGLET***************************************************************
C
C
C DEFINE DIREKTIVEN     : KEINE
C
C UPROG                 : KEINE                                      
C                                                            
C 24.09.1997            : ADNAN MERI (ORIGINAL)
C
C*MGLET***************************************************************
C

      REAL ALFA,BETA,GAMA,a,b,c,DX(II)

      a = DX(K-1)
      b = DX(K-1) + DX(K-2)
      c = DX(K-1) + DX(K-2) + DX(K-3)



      GAMA = (a*b)/(b*(a-c) + c*(c-a))
      BETA = (a+ GAMA*(c-a)) / (a-b)
      alfa = 1.0 - GAMA - BETA
  
c
c*******************************************************************
c
      RETURN
      END



