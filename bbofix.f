










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
      SUBROUTINE BBOFIX (KMX,JMX,IMX,
     $                    U,V,W,P,G,UBO,VBO,WBO,
     $                    ISTART,ISTOP,JSTART,JSTOP,ITYP
     $                   )     
C*STARLET***************************************************************
C        B B O F I X   SETZEN DER RANDBEDINGUNGEN FUER DIE
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
c      WBO(KMX,JMX, 2 ) - EINSTROEMFELD
C        ITYP           - HOLLERITH-KONSTANTE:
C                         'P'  VOR DER DRUCKKORREKTUR
C                         'T'  VOR DEM ZEITSCHRITT
C
C VERS:   6. 6.97 (MM)  : ORIGINAL
C                         AUS BBONOS und BFRFIX ABGELEITET
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
     $            P(KMX,JMX,IMX),   G(KMX,JMX,IMX), 
     $            UBO(JMX,IMX,2),VBO(JMX,IMX,2),WBO(JMX,IMX,2)
C
C
C
C                                 **************************************
C                                 FIXED-CONDITIONS AN DER BOTTOM FLAECHE
C                                 **************************************
C
C                                 DIE W-KOMPONENTEN WERDEN BEIM EINLESEN
C                                 VOM DATENTRAEGER IM FELD "WBO" GESPEI-
C                                 CHERT. VOR DEM ZEITSCHRITT WERDEN (WIE
C                                 BEI DER DRUCKKORREKTUR) NUR DIE W-KOMP
C                                 NOCHMALS KORRIGIERT. DIE U- UND V-KOMP
C                                 SIND BEIM EINLESEN V. DATENTRAEGER BE-
C                                 REITS RICHTIG GESETZT.
C
C
C
      IF ( ITYP .EQ. 'T' ) THEN

         DO I=ISTART,ISTOP
         DO J=JSTART,JSTOP

             U(2,J,I) = UBO(J,I,2)
             V(2,J,I) = VBO(J,I,2)
             W(2,J,I) = WBO(J,I,2)
             P(2,J,I) =   P(3,J,I)

         ENDDO
         ENDDO

      ELSEIF ( ITYP .EQ. 'P') THEN

         DO I=ISTART,ISTOP
         DO J=JSTART,JSTOP

             W(2,J,I) = WBO(J,I,2)

         ENDDO
         ENDDO

      ENDIF

      RETURN
      END

