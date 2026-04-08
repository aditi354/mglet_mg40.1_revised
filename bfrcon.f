










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
      SUBROUTINE BFRCON (
     $                   KMX,JMX,IMX,U,V,W,P,G,UFR,
     $                   KMXN,JMXN,IMXN,UN,VN,WN,PN,GN,
     $                   JSTA,JSTO,KSTA,KSTO,
     $                   IPOS,JPOS,KPOS,ITYP,IRB,JRB,KRB,
     $                   IPROC,IPROCNBR)
C*MGLET*****************************************************************
C        B F R C O N    SETZEN DER RANDBEDINGUNGEN FUER DIE
C                         GESCHWINDIGKEITSFELDER. AUS NACHBARGITTER
C*MGLET*****************************************************************
C
C PARAM: 
C        KMX, JMX, IMX  - DIMENSIONEN DES GITTERS 
C        U(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        V(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        W(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        T(KMX,JMX,IMX) + SCALAR T FIELD
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
C                         'S' BEFORE SIP/GAUS-SEIDEL PRESSURE CORRECTION
C            IRB                  - INTEGER-KENNZAHL I-RICHTUNG RED-BLACK
C            JRB                  - INTEGER-KENNZAHL J-RICHTUNG RED-BLACK
C            KRB                  - INTEGER-KENNZAHL K-RICHTUNG RED-BLACK
C
C VERS:  19.12.93 (MM)  : ORIGINAL AUS BFRPER ABGELEITET
C        28.02.03 (TB)  : SCALAR TRANSPORT IMPLEMENTED
C        11.08.03 (TB)  : SIP/GAUS-SEIDEL SUPPORT IMPLEMENTED
C
C DEFINE-DIREKTIVEN     : QUD 
C
C*STARLET***************************************************************

C#ifdef _IRIS_
      CHARACTER (LEN=1)  ITYP
C#else
C      CHARACTER (LEN=1)  ITYP
C#endif

      REAL
     $     U(KMX,JMX,IMX),  V(KMX,JMX,IMX),  W(KMX,JMX,IMX),
     $     P(KMX,JMX,IMX),  G(KMX,JMX,IMX)
      REAL
     $    UN(KMXN,JMXN,IMXN), VN(KMXN,JMXN,IMXN),  WN(KMXN,JMXN,IMXN),
     $    PN(KMXN,JMXN,IMXN),  GN(KMXN,JMXN,IMXN)
C




CRED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-
C
      IF (IRB .EQ. 0) THEN
         ISTART = 1
         ISTOP  = IMX
      ELSEIF (IRB .EQ. 1) THEN
         ISTART =   2  + 2
         ISTOP  = IMX/2
      ELSEIF (IRB .EQ. 2) THEN
          IF (ITYP.EQ.'P') RETURN
      ELSE
         CALL ERRR (504,'IRB BFRCON')
      ENDIF   
C
C
      IF (JRB .EQ. 0) THEN
         JSTRIDE = 1
         JSTART  = 1
         JSTOP   = JMX
      ELSEIF (JRB .EQ. 1) THEN
         JSTART = MAX(  2  + 1,JSTA)
         JSTOP  = MIN(JMX/2   ,JSTO)
      ELSEIF (JRB .EQ. 2) THEN
         JSTART = MAX(JMX/2 + 1 ,JSTA)
         JSTOP  = MIN(JMX -   2 ,JSTO)
      ELSE
         CALL ERRR (504,'JRB BFRCON')
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
         CALL ERRR (504,'KRB BFRCON')
      ENDIF   
C
CRED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-
C
C               IN = 2 + IPOS - 3
      IN = IPOS -1
C
C
C
C                                 FUER ZEITSCHRITT
      IF (ITYP.EQ.'T') THEN

         DO J=JSTA,JSTO
         DO K=KSTA,KSTO

               JN = J + JPOS - 3
               KN = K + KPOS - 3

               U(K,J,  2) = UN(KN,JN,IN)
               V(K,J,  2) = VN(KN,JN,IN)
               W(K,J,  2) = WN(KN,JN,IN)
               P(K,J,  2) = PN(KN,JN,IN)
               G(K,J,  2) = GN(KN,JN,IN)
               

C jk 18.7.2005
               U(K,J,  1) = UN(KN,JN,IN-1)
               V(K,J,  1) = VN(KN,JN,IN-1)
               W(K,J,  1) = WN(KN,JN,IN-1)
               P(K,J,  1) = PN(KN,JN,IN-1)
               G(K,J,  1) = GN(KN,JN,IN-1)

         ENDDO
         ENDDO

C                                 FUER DRUCKKORREKTUR
C                                 NUR NORMALKOMPONENTE MUSS GESETZT
C                                 WERDEN
       ELSEIF (ITYP.EQ.'P') THEN

         DO J=JSTART,JSTOP

         DO K=KSTA+KP,KSTO,KSTRIDE

               JN = J + JPOS - 3
               KN = K + KPOS - 3

               U(K,J,  2) = UN(KN,JN,IN)

         ENDDO
         ENDDO

C                           FUER DRUCKKORREKTUR MIT GAUSS-SEIDEL ODER SIP

       ELSEIF (ITYP.EQ.'S') THEN

         DO J=JSTA,JSTO
         DO K=KSTA,KSTO

               JN = J + JPOS - 3
               KN = K + KPOS - 3

               P(K,J,  2) = PN(KN,JN,IN)

         ENDDO
         ENDDO
        
       ELSE
                CALL ERRR (501,' BFRCON')
       ENDIF
C

      RETURN
      END
