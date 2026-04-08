










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
      SUBROUTINE SKONLE   (GMOL,RHO,ITCON)
C*STARLET***************************************************************
C        S K O N L E      FESTLEGUNG DER KONSTANTEN FUER DIE
C                         LES-SIMULATION
C*STARLET***************************************************************
C
C PARAM: GMOL           - MOLEKULARE DYNAMISCHE VISKOSITAET
C        RHO            - DICHTE DES FLUIDS (= CONST)
C        ITCON          - INTEGERKONSTANTE ZUR AUSWAHL DES ENTSPRECHEN-
C                         DEN SATZES VON TURBULENZMODELLKONSTANTEN
C                         ITCON = 1 : KONSTANTEN NACH DEARDORFF
C
C VERS:  10.09.85 (HW)  : BERECHNUNG WICHTIGER KONSTANTEN FUER DAS
C                         1/7 POTENZGESETZ
C        14.01.86 (HW)  : ERWEITERUNG DES KONSTANTENSATZES (CPO11,CPO12)
C
C UPROG                 : ERRR
C
C*STARLET***************************************************************
C

      COMMON /CONLES/  CONV2S,CMUE,CAPPA,ECONST
      SAVE   /CONLES/

      COMMON /CONPOT/  CWA,CWB,CPO1,CPO2,CPO3,CPO4,CPO5,CPO6,CPO7,CPO8,
     $                 CPO9,CPO10,CPO11,CPO12
      SAVE   /CONPOT/
C
      IF(ITCON .NE. 1) GOTO 2100
C
      CONV1S = 0.1
      CONV2S = 0.094
      CMUE   = 0.09
C
C                                 KONSTANTEN DER WANDFUNKTION:
C
      CAPPA  = 0.4
      ECONST = 9.0
C
C                                 KONSTANTEN FUER DIE BERECHNUNG VON
C                                 TAUW, UTAU UND DUDZ MITTELS DES 1/7
C                                 POTENZGESETZES
C
C                                 KONSTANTEN A  UND  B  DES 1/7 POTENZ-
C                                 GESETZES
C
      CWA    = 8.3
      CWB    = 1.0/7.0
C
C                                 ABGELEITETE KONSTANTEN
C
      CPO1   = 1.0 - CWB
      CPO2   = 1.0 + CWB
      CPO3   = 1.0 / CPO2
      CPO4   = CPO2/CWA*(GMOL/RHO)**CWB
      CPO5   = 0.5*CPO1*(CWA**(CPO2/CPO1))*(GMOL/RHO)**CPO2
      CPO6   = 0.5*GMOL/RHO*CWA**(2.0/CPO1)
      CPO7   = CWA*CWB/((GMOL/RHO)**CWB)
      CPO8   = 2.0 / CPO2
      CPO9   = -CPO1*CPO6/CPO2
      CPO10  = GMOL/RHO * CWA**(1.0/CPO1)
      CPO11  = 2.0*CPO6*CWB*CPO1 / (2.0** CWB)
      CPO12  = CWB*CPO2          / (2.0**(CWB-1.0))
C
      RETURN
C
 2100 CALL ERRR(510,' SKONLE   ')
      RETURN
      END
