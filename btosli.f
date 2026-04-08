










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
      SUBROUTINE BTOSLI (KMX,JMX,IMX,
     $                    U,V,W,P,G,
     $                   ISTA,ISTO,JSTA,JSTO,
     $                   ITYP,IRB,JRB,KRB
     $                   )     
C*STARLET***************************************************************
C        B T O S L I   SETZEN DER RANDBEDINGUNGEN FUER DIE
C                         GESCHWINDIGKEITSFELDER AND SCALAR T.
C                         (LARGE-EDDY-SIMULATION)
C*STARLET***************************************************************
C
C PARAM: 
C        KMX, JMX, IMX  - DIMENSIONEN DES GITTERS 
C        U(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        V(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        W(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        T(KMX,JMX,IMX) + SCALAR T FIELD
C        P(KMX,JMX,IMX) - DRUCKFELD
C        G(KMX,JMX,IMX) - EFFEKTIVE DYNAMISCHE VISKOSITAET (= MUE)
C        ITYP           - HOLLERITH-KONSTANTE:
C                         'P'  VOR DER DRUCKKORREKTUR
C                         'T'  VOR DEM ZEITSCHRITT
C
C VERS:   8. 3.93 (MM)  : ORIGINAL
C                         AUS BTOPLE (BTOPL) ABGELEITET
C        28.02.03 (TB)  : SCALAR BOUNDARY TREATMENT IMPLEMENTED
C
C DEFINE-DIREKTIVEN     : QUD 
C
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
      CHARACTER (LEN=1)  ITYP
C
      REAL
     $            U(KMX,JMX,IMX),   V(KMX,JMX,IMX),   W(KMX,JMX,IMX),
     $            P(KMX,JMX,IMX),   G(KMX,JMX,IMX)
C
C
      KM1 = KMX-1
      KM2 = KMX-2
      KM3 = KMX-3
      KM4 = KMX-4
C
      NBND = 2
      IF (IRB .EQ. 0) THEN
         ISTART = 1
         ISTOP  = IMX
      ELSEIF (IRB .EQ. 1) THEN
         ISTART = NBND + 1
         ISTOP  = IMX/2
      ELSEIF (IRB .EQ. 2) THEN
         ISTART = IMX/2 + 1
         ISTOP  = IMX - NBND
      ELSE
         CALL ERRR (504,'IRB BTOSLI')
      ENDIF   
C
      IF (JRB .EQ. 0) THEN
         JSTART = 1
         JSTOP  = JMX
      ELSEIF (JRB .EQ. 1) THEN
         JSTART =   2  + 1
         JSTOP  = JMX/2
      ELSEIF (JRB .EQ. 2) THEN
         JSTART = JMX/2 + 1
         JSTOP  = JMX -   2 
      ELSE
         CALL ERRR (504,'JRB BBACON')
      ENDIF   
C
C                                 **************************************
C                                 REIBUNGSLOSE, UNDURCHLAESSIGE WAND
C                                 **************************************
C
C
C
C                                 FUER ZEITSCHRITT
      IF (ITYP.EQ.'T') THEN

         DO I=ISTA,ISTO
         DO J=JSTA,JSTO

               U(KM1,J,I) = 0.0
               V(KM1,J,I) = 0.0
               W(KM2,J,I) = 0.0
               P(KM1,J,I) = P(KM2,J,I)
               G(KM1,J,I) = 2.0*SMALL

         ENDDO
         ENDDO


C                                 FUER DRUCKKORREKTUR
C                                 NUR NORMALKOMPONENTE MUSS 
C                                 GESETZT WERDEN
       ELSEIF (ITYP.EQ.'P') THEN

         DO I=ISTART,ISTOP
         DO J=JSTART,JSTOP

               W(KM2,J,I) = 0.0

         ENDDO
         ENDDO

 
       ELSE
                CALL ERRR (501,' BTOSLI')
       ENDIF

      RETURN
      END

