










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
      SUBROUTINE ENSTRO  (KK,JJ,II,KMX,JMX,IMX,OX,OY,OZ,O2)
C*STARLET***************************************************************
C        E N S T R O      BERECHNUNG DER ENSTROPHY:
C                         ENSTROPHY = 1/2 * <OMEGA_I * OMEGA_I>;
C                         (SUMMATIONSKONVENTION GILT !)
C                         KUBUSFOERMIGE EINBAUTEN WERDEN BERUECKSICHTIGT
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        OX(KK,JJ,II)   - X-KOMPONENTE DER VORTICITY
C        OY(KK,JJ,II)   - Y-KOMPONENTE DER VORTICITY
C        OZ(KK,JJ,II)   - Z-KOMPONENTE DER VORTICITY
C        O2(KK,JJ,II)   + ENSTROPHY (O2 -> OMEGA **2)
C
C UPROG                 : KEINE
C
C DEFINE-DIREKTIVEN     : KEINE
C
C VERS:  22.08.88 (HW)  : ORIGINAL
C
C
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/

      COMMON /CONLES/  CONV2S,CMUE,CAPPA,ECONST
      SAVE   /CONLES/
C
      CHARACTER (LEN=1)  ITYP
C
      REAL          OX(KK,JJ,II),  OY(KK,JJ,II),  OZ(KK,JJ,II),
     $              O2(KK,JJ,II)
C
C
      KM2    = KMX-2
      JM2    = JMX-2
      IM2    = IMX-2
C
      DO 100 I = 2,IM2
         IKST   =      MAX0(I      ,3)
         IKSTPS = MIN0(MAX0(I+1    ,3),IM2)
         DO 100 J = 2,JM2
            JKST   =      MAX0(J      ,3)
            JKSTPS = MIN0(MAX0(J+1    ,3),JM2)
C
            KSTART = 3
C
            DO 100 K = KSTART,KM2
               O2 (K,J,I) = 0.5 * (
     $                      (0.5*(OX(K  ,J  ,I  ) + OX(K  ,J  ,I+1)))**2
     $                    + (0.5*(OY(K  ,J  ,I  ) + OY(K  ,J+1,I  )))**2
     $                    + (0.5*(OZ(K  ,J  ,I  ) + OZ(K+1,J  ,I  )))**2
     $                    )
  100 CONTINUE
C
      RETURN
      END
