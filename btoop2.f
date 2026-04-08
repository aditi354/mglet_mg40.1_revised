










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
      SUBROUTINE BTOOP2 (KMX,JMX,IMX,
     $                    DZ,DDZ,
     $                    U,V,W,P,G,
     $                    ISTART,ISTOP,JSTART,JSTOP,ITYP)
C*STARLET***************************************************************
C        B T O O P 2   SETZEN DER RANDBEDINGUNGEN FUER DIE
C                         GESCHWINDIGKEITSFELDER.
C                         (LARGE-EDDY-SIMULATION)
C*STARLET***************************************************************
C
C PARAM: 
C        KMX, JMX, IMX  - DIMENSIONEN DES GITTERS 
C        U(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        V(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        W(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        P(KMX,JMX,IMX) - DRUCKFELD
C        G(KMX,JMX,IMX) - EFFEKTIVE DYNAMISCHE VISKOSITAET (= MUE)
C        ITYP           - HOLLERITH-KONSTANTE:
C                         'P'  VOR DER DRUCKKORREKTUR
C                         'T'  VOR DEM ZEITSCHRITT
C
C VERS:   8. 3.93 (MM)  : ORIGINAL
C                         AUS BTOPLE (BTOPL) ABGELEITET
C
C DEFINE-DIREKTIVEN     : QUD 
C
C*STARLET***************************************************************
C
      CHARACTER (LEN=1)  ITYP
C
      REAL
     $              DZ(  KMX   ),    DDZ(  KMX   ),
     $            U(KMX,JMX,IMX),   V(KMX,JMX,IMX),   W(KMX,JMX,IMX),
     $            P(KMX,JMX,IMX),   G(KMX,JMX,IMX)
C
      KM1 = KMX-1
      KM2 = KMX-2
      KM3 = KMX-3
      KM4 = KMX-4
C
C
C                                 **************************************
C                                 OPEN-WALL, 2. ABLEITUNG = 0
C                                 **************************************
C
C
C
         DO I=ISTART,ISTOP
         DO J=JSTART,JSTOP

               U(KM1,J,I) = U(KM2,J,I) + (U(KM2,J,I)-U(KM3,J,I))
     $                    *   DZ(KM2)          / DZ(KM3)
               V(KM1,J,I) = V(KM2,J,I) + (V(KM2,J,I)-V(KM3,J,I))
     $                    *   DZ(KM2)          / DZ(KM3)
               W(KM2,J,I) = W(KM3,J,I) + (W(KM3,J,I)-W(KM4,J,I))
     $                    *  DDZ(KM2)          /DDZ(KM3)
               G(KM1,J,I) = G(KM2,J,I) + (G(KM2,J,I)-G(KM3,J,I))
     $                    *   DZ(KM2)          / DZ(KM3)


         ENDDO
         ENDDO


      RETURN
      END

