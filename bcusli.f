










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
      SUBROUTINE BCUSLI (KMX,JMX,IMX,
     $                    U,V,W,P,G,B,HILF,
     $                    IC1,IC2,JC1,JC2,KC1,KC2,ITYP
     $                   )      
C*STARLET***************************************************************
C        B C U S L I   SETZEN DER RANDBEDINGUNGEN FUER DIE
C                         GESCHWINDIGKEITSFELDER AND SCALAR T.
C                         (LARGE-EDDY-SIMULATION)
C*STARLET***************************************************************
C
C PARAM: 
C        KMX, JMX, IMX  - DIMENSIONEN DES GITTERS 
C        U(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        V(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        W(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        T(KK,JJ,II)    + SCALAR T FIELD
C        P(KMX,JMX,IMX) - DRUCKFELD
C        G(KMX,JMX,IMX) - EFFEKTIVE DYNAMISCHE VISKOSITAET (= MUE)
C        ITYP           - HOLLERITH-KONSTANTE:
C                         'P'  VOR DER DRUCKKORREKTUR
C                         'T'  VOR DEM ZEITSCHRITT
C
C VERS:   8. 3.93 (MM)  : ORIGINAL
C                         AUS BCUBLE (BCUBL) ABGELEITET
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
     $            P(KMX,JMX,IMX),   G(KMX,JMX,IMX),   B(KMX,JMX,IMX),
     $         HILF(KMX,JMX)
C
C
C
C                                 **************************************
C                                 KUBUS-WAENDE REIBUNGSBEHAFTET UND
C                                 UNDURCHLAESSIG (NOSLIP)
C                                 **************************************
C
      GPAR   = 2.0*SMALL

C
C                                 BELEGUNG DES KOERPERINNEREN
C
      DO 10 I= IC1,IC2
         DO 50 J= JC1,JC2
            DO 40 K= KC1,KC2
               FACTOR    = 1.0 + MIN( 0.0 , MAX( -1.0 , B(K,J,I) ) )
               U(K  ,J  ,I  ) = U(K  ,J  ,I  )*FACTOR   
               U(K  ,J  ,I-1) = U(K  ,J  ,I-1)*FACTOR   
               V(K  ,J  ,I  ) = V(K  ,J  ,I  )*FACTOR   
               V(K  ,J-1,I  ) = V(K  ,J-1,I  )*FACTOR   
               W(K  ,J  ,I  ) = W(K  ,J  ,I  )*FACTOR   
               W(K-1,J  ,I  ) = W(K-1,J  ,I  )*FACTOR   
               G(K  ,J  ,I  ) = G(K  ,J  ,I  )*FACTOR   
     $                        + GPAR          *(1.0-FACTOR   )
   40       CONTINUE
   50    CONTINUE
   10 CONTINUE
      RETURN
      END

