










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
      SUBROUTINE BBACONBP (
     $     KMX,JMX,IMX,BP,
     $     KMXN,JMXN,IMXN,BPN,
     $     JSTA,JSTO,KSTA,KSTO,
     $     IPOS,JPOS,KPOS,
     $     IPROC,IPROCNBR)
C     *MGLET*****************************************************************
C     B B A C O N B P   SETZEN DER RANDBEDINGUNGEN FUER DIE
C     GESCHWINDIGKEITSFELDER. AUS NACHBARN
C     *MGLET*****************************************************************
C     
C  KREUZINGER UND MANHART TURBULENZ GmbH
C  JK      11. 8.2005
C  aus BBACON
C
C     PARAM: 
C     KMX, JMX, IMX  - DIMENSIONEN DES GITTERS 
C     U(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C     JSTART,JSTOP
C     KSTART,KSTOP   - BEREICH IN DEM RANDBEDINGUNG GESETZT WIRD
C     IPOS,JPOS,KPOS - POSITION IM NACHBARGITTER, AUF DER DER PUNKT
C     MIT DEN INDIZES (3,3,3) ZU LIEGEN KOMMT
C     
C     *STARLET***************************************************************
C     
C     
      REAL
     $     BP(KMX,JMX,IMX)
      REAL
     $     BPN(KMXN,JMXN,IMXN)
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
      IM1 = IMX-1
      IM2 = IMX-2
      IM3 = IMX-3
C     
C     IN = IM1 + IPOS - IM2
      IN = IPOS +1


C     
         IF ((JSTA.NE.1).OR.(JSTO.NE.JMX)) STOP 'bbaconbp A'
         IF ((KSTA.NE.1).OR.(KSTO.NE.KMX)) STOP 'bbaconbp A'
C     
C     FUER ZEITSCHRITT
C     NORMALKOMPONENTE WIRD EINEN INDEX
C     AUSSERHALB  GESETZT

         DO J=JSTA,JSTO
            DO K=KSTA,KSTO

               JN = J + JPOS - 3
               KN = K + KPOS - 3

               BP(K,J,IM1) = BPN(KN,JN,IN)
               BP(K,J,IMX) = BPN(KN,JN,IN+1)

            ENDDO
         ENDDO
C     


      RETURN
      END

