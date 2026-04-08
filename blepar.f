










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
      SUBROUTINE BLEPAR (KMX,JMX,IMX,U,V,W,P,G,UFR,
     $                   KMXN,JMXN,IMXN,UN,VN,WN,PN,GN,
     $                   ISTA,ISTO,KSTA,KSTO,
     $                   IPOS,JPOS,KPOS,ITYP,IRB,JRB,KRB,
     $                   IPROC,IPROCNBR
     $                  )     
C*MGLET*****************************************************************
C        B L E P A R    SETZEN DER RANDBEDINGUNGEN FUER DIE
C                         GESCHWINDIGKEITSFELDER. AUS NACHBARGITTER
C*MGLET*****************************************************************
C
C PARAM: 
C        KMX, JMX, IMX  - DIMENSIONEN DES GITTERS 
C        U(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        V(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        W(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        T(KK,JJ,II)    + SCALAR T FIELD
C        P(KMX,JMX,IMX) - DRUCKFELD
C        G(KMX,JMX,IMX) - EFFEKTIVE DYNAMISCHE VISKOSITAET (= MUE)
C        ...N           - ENTSPRECHENDE GROESSE AUS DEM NACHBARGITTER
C        JSTART,JSTOP
C        KSTART,KSTOP   - BEREICH IN DEM RANDBEDINGUNG GESETZT WIRD
C        IPOS,JPOS,KPOS - POSITION IM NACHBARGITTER, AUF DER DER PUNKT
C                         MIT DEN INDIZES (3,3,3) ZU LIEGEN KOMMT
C        ITYP           - HOLLERITH-KONSTANTE:
C                         'P'  VOR DER DRUCKKORREKTUR
C                         'T'  VOR DEM ZEITSCHRITT
C        IRB            - INTEGER-KENNZAHL I-RICHTUNG RED-BLACK
C        JRB            - INTEGER-KENNZAHL J-RICHTUNG RED-BLACK
C        KRB            - INTEGER-KENNZAHL K-RICHTUNG RED-BLACK
C
C VERS:  27.07.95 (MM)  : ORIGINAL AUS BLECON ABGELEITET
C        28.02.03 (TB)  : SCALAR BOUNDARY TREATMENT IMPLEMENTED
C        
C
C*MGLET*****************************************************************
C
C
      CHARACTER (LEN=1)  ITYP
C
      REAL
     $     U(KMX,JMX,IMX),  V(KMX,JMX,IMX),  W(KMX,JMX,IMX),
     $     P(KMX,JMX,IMX),  G(KMX,JMX,IMX)
      REAL
     $    UN(KMXN,JMXN,IMXN), VN(KMXN,JMXN,IMXN),  WN(KMXN,JMXN,IMXN),
     $    PN(KMXN,JMXN,IMXN),  GN(KMXN,JMXN,IMXN)
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CRED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-
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
         CALL ERRR (504,'IRB BLEPAR')
      ENDIF   
C
C
      IF (JRB .EQ. 0) THEN
         JSTART = 1
         JSTOP  = JMX
      ELSEIF (JRB .EQ. 1) THEN
          IF (ITYP.EQ.'P') RETURN
      ELSEIF (JRB .EQ. 2) THEN
         JSTART = JMX/2 + 1
         JSTOP  = JMX - NBND
      ELSE
         CALL ERRR (504,'JRB BLEPAR')
      ENDIF   
C
C
      IF (KRB .EQ. 0) THEN
         KSTRIDE = 1
         KP      = 0
      ELSEIF (KRB .EQ. 1) THEN
         KSTRIDE = 2
         KP      = 0
      ELSEIF (KRB .EQ. 2) THEN
         KSTRIDE = 2
         KP      = 1
      ELSE
         CALL ERRR (504,'KRB BLEPAR')
      ENDIF   
C

C
CRED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-
C
C                         FINE-GRID: J-INDEX=1 AND 2 (FRONT)
               JF = JMX - 1
C                         J-INDEX OF COARSE-GRID
C                                   (=PARENT, CALLED NEIGHBOUR)
            JC = JPOS - 1 + (JF-1)/2

C
C
C
C                                 FUER ZEITSCHRITT
      IF (ITYP.EQ.'T') THEN

         DO IF=ISTA,ISTO,2
               IC = IPOS - 1 + (IF-1)/2
         DO KF=KSTA,KSTO,2
            KC = KPOS - 1 + (KF-1)/2

               U(KF  ,JF  ,IF  ) = 0.5*( UN(KC,JC,IC) + UN(KC,JC,IC-1))
               U(KF+1,JF  ,IF  ) = 0.5*( UN(KC,JC,IC) + UN(KC,JC,IC-1))
               U(KF  ,JF  ,IF+1) = UN(KC,JC,IC)
               U(KF+1,JF  ,IF+1) = UN(KC,JC,IC)

               V(KF  ,JF  ,IF  ) = 0.5*( VN(KC,JC,IC) + VN(KC,JC-1,IC))
               V(KF+1,JF  ,IF  ) = 0.5*( VN(KC,JC,IC) + VN(KC,JC-1,IC))
               V(KF  ,JF  ,IF+1) = 0.5*( VN(KC,JC,IC) + VN(KC,JC-1,IC))
               V(KF+1,JF  ,IF+1) = 0.5*( VN(KC,JC,IC) + VN(KC,JC-1,IC))

               W(KF  ,JF  ,IF  ) = 0.5*( WN(KC,JC,IC)  + WN(KC-1,JC,IC))
               W(KF  ,JF  ,IF+1) = 0.5*( WN(KC,JC,IC)  + WN(KC-1,JC,IC))
               W(KF+1,JF  ,IF  ) = WN(KC,JC,IC)
               W(KF+1,JF  ,IF+1) = WN(KC,JC,IC)

               P(KF  ,JF  ,IF  ) = PN(KC,JC,IC)
               P(KF+1,JF  ,IF  ) = PN(KC,JC,IC)
               P(KF  ,JF  ,IF+1) = PN(KC,JC,IC)
               P(KF+1,JF  ,IF+1) = PN(KC,JC,IC)

               G(KF  ,JF  ,IF  ) = GN(KC,JC,IC)
               G(KF+1,JF  ,IF  ) = GN(KC,JC,IC)
               G(KF  ,JF  ,IF+1) = GN(KC,JC,IC)
               G(KF+1,JF  ,IF+1) = GN(KC,JC,IC)



         ENDDO
         ENDDO

C                                 FUER DRUCKKORREKTUR
C                                 NUR NORMALKOMPONENTE MUSS GESETZT
C                                 WERDEN
       ELSEIF (ITYP.EQ.'P') THEN

C                   NORMAL-KOMPONENT MUST NOT BE SET
C                   BECAUSE PRESSURE IN THE BOUNDARY-CELLS
C                   IS NOT MODIFIED, SO NORMAL-KOMPONENTS MUST
C                   BE ABLE TO DEVELOP FREELY
               RETURN
C                         FINE-GRID: I-INDEX=IMX-2 +1(FRONT)
               JF = JMX - 2
C                         J-INDEX OF COARSE-GRID
C                                   (=PARENT, CALLED NEIGHBOUR)
            JC = JPOS - 1 + (JF-1)/2

         DO IF=ISTART,ISTOP,2
               IC = IPOS - 1 + (IF-1)/2
         DO KF=KSTA+KP,KSTO,2
            KC = KPOS - 1 + (KF-1)/2

               V(KF  ,JF  ,IF  ) = VN(KC,JC,IC)
               V(KF  ,JF  ,IF+1) = VN(KC,JC,IC)

         ENDDO
         ENDDO

       ELSE
                CALL ERRR (501,' BLEPAR')
       ENDIF


      RETURN
      END

