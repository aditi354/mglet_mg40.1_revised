










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
      SUBROUTINE BBAPER (KMX,JMX,IMX,     
     $                    U,V,W,P,G,
     $                    JSTART,JSTOP,KSTART,KSTOP,ITYP
     $                   )       
C*STARLET***************************************************************
C        B B A P E R    SETZEN DER RANDBEDINGUNGEN FUER DIE
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
C                         AUS BFROLE (BFROL) ABGELEITET
C        28.02.03 (TB)  : SCALAR BOUNDARY TREATMENT IMPLEMENTED
C
C DEFINE-DIREKTIVEN     : QUD 
C
C*STARLET***************************************************************
C
      CHARACTER (LEN=1)  ITYP
C
      REAL
     $            U(KMX,JMX,IMX),   V(KMX,JMX,IMX),   W(KMX,JMX,IMX),
     $            P(KMX,JMX,IMX),   G(KMX,JMX,IMX)
C
C
      IM1 = IMX-1
      IM2 = IMX-2
      IM3 = IMX-3
C
C                                 **************************************
C                                 PERIODISCHE RANDBEDINGUNGEN IN
C                                 X-RICHTUNG
C                                 **************************************
C
C                                 FUER ZEITSCHRITT
C                                 NORMALKOMPONENTE WIRD EINEN INDEX
C                                 AUSSERHALB  GESETZT
      IF (ITYP.EQ.'T') THEN

         DO J=JSTART,JSTOP
         DO K=KSTART,KSTOP

               U(K,J,IM1) = U(K,J,  3)
               V(K,J,IM1) = V(K,J,  3)
               W(K,J,IM1) = W(K,J,  3)
               P(K,J,IM1) = P(K,J,  3)
               G(K,J,IM1) = G(K,J,  3)

         ENDDO
         ENDDO
C

C                                 FUER DRUCKKORREKTUR
C                                 NORMALKOMPONENTE MUSS NICHT
C                                 GESETZT WERDEN
       ELSEIF (ITYP.EQ.'P') THEN

         RETURN
 
       ELSE
                CALL ERRR (501,' BBAPER')
       ENDIF


      RETURN
      END

