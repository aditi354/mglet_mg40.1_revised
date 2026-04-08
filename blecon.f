










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
      SUBROUTINE BLECON (
     $                   KMX,JMX,IMX,U,V,W,P,G,UFR,
     $                   KMXN,JMXN,IMXN,UN,VN,WN,PN,GN,
     $                   ISTA,ISTO,KSTA,KSTO,
     $                   IPOS,JPOS,KPOS,ITYP,IRB,JRB,KRB,
     $                   IPROC,IPROCNBR)
C*MGLET*****************************************************************
C        B L E C O N    SETZEN DER RANDBEDINGUNGEN FUER DIE
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
C        IRB            - INTEGER-KENNZAHL I-RICHTUNG RED-BLACK
C        JRB            - INTEGER-KENNZAHL J-RICHTUNG RED-BLACK
C        KRB            - INTEGER-KENNZAHL K-RICHTUNG RED-BLACK
C
C VERS:   6. 1.94 (MM)  : ORIGINAL AUS BRICON ABGELEITET
C        11.05.95 (MM,AO): EINFUEHRUNG DER KOMMUNIKATION UNTER MPI
C        28.02.03 (TB)  : SCALAR TRANSPORT IMPLEMENTED
C        11.08.03 (TB)  : SIP/GAUS-SEIDEL SUPPORT IMPLEMENTED
C        
C
C DEFINE-DIREKTIVEN     : QUD 
C
C*STARLET***************************************************************
C
C
C#ifdef _IRIS_
      CHARACTER (LEN=1)  ITYP
C#else
C      CHARACTER (LEN=1)  ITYP
C#endif
C
      REAL
     $     U(KMX,JMX,IMX),  V(KMX,JMX,IMX),  W(KMX,JMX,IMX),
     $     P(KMX,JMX,IMX),  G(KMX,JMX,IMX)
      REAL
     $    UN(KMXN,JMXN,IMXN), VN(KMXN,JMXN,IMXN),  WN(KMXN,JMXN,IMXN),
     $    PN(KMXN,JMXN,IMXN),  GN(KMXN,JMXN,IMXN)
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
         CALL ERRR (504,'IRB BLECON')
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
         CALL ERRR (504,'JRB BLECON')
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
         CALL ERRR (504,'KRB BLECON')
      ENDIF   
C

C
CRED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-
C
      JM1 = JMX-1
      JM2 = JMX-2
      JM3 = JMX-3

C       JN = JPOS + 1 naechste Zeile ausgeschrieben
C       JN = JM1 + JPOS - JM2
       JN = JPOS +1
C
C
C
C                                 FUER ZEITSCHRITT
      IF (ITYP.EQ.'T') THEN


         DO I=ISTA,ISTO
         DO K=KSTA,KSTO

               IN = I + IPOS - 3
               KN = K + KPOS - 3


               U(K,JM1,  I) = UN(KN,JN,IN)
               V(K,JM1,  I) = VN(KN,JN,IN)
               W(K,JM1,  I) = WN(KN,JN,IN)
               P(K,JM1,  I) = PN(KN,JN,IN)
               G(K,JM1,  I) = GN(KN,JN,IN)
C jk 20.7.2005
               U(K,JMX,  I) = UN(KN,JN+1,IN)
               V(K,JMX,  I) = VN(KN,JN+1,IN)
               W(K,JMX,  I) = WN(KN,JN+1,IN)
               P(K,JMX,  I) = PN(KN,JN+1,IN)
               G(K,JMX,  I) = GN(KN,JN+1,IN)
c$$$C MERI 180199
c$$$               U(K,JMX,  I) = U(KN,JN+1,I)
c$$$               V(K,JMX,  I) = V(KN,JN+1,I)
c$$$               W(K,JMX,  I) = W(KN,JN+1,I)
c$$$               P(K,JMX,  I) = P(KN,JN+1,I)
c$$$               G(K,JMX,  I) = G(KN,JN+1,I)
c$$$#ifdef _TSCAL_
c$$$               T(K,JMX,  I) = T(KN,JN+1,I)
c$$$#endif                

         ENDDO
         ENDDO

C                                 FUER DRUCKKORREKTUR
C                                 NUR NORMALKOMPONENTE MUSS GESETZT
C                                 WERDEN
       ELSEIF (ITYP.EQ.'P') THEN

         DO I=ISTART,ISTOP

         DO K=KSTA+KP,KSTO,KSTRIDE

               IN = I + IPOS - 3
               KN = K + KPOS - 3

               V(K,JM2,  I) = VN(KN,JN-1,IN)
C               V(K,JM2,  I) = VN(KN,JN,IN)

         ENDDO
         ENDDO
C                    FUER DRUCKKORREKTUR GAUS_SEIDEL ODER SIP

      ELSEIF (ITYP.EQ.'S') THEN

         DO I=ISTA,ISTO
         DO K=KSTA,KSTO

               IN = I + IPOS - 3
               KN = K + KPOS - 3

               P(K,JM1,  I) = PN(KN,JN,IN)

            ENDDO
         ENDDO

      ELSE
               CALL ERRR (501,' BLECON')
      ENDIF


      RETURN
      END

