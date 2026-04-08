










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
      SUBROUTINE VPLESC(KK,JJ,II,KMX,JMX,IMX,
     $                  DX,DY,DZ,DDX,DDY,DDZ,
     $                  RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                  IPCORR,OMBETA,DT,U,V,W,P,B,
     $                  BP,BU,BV,BW,SDIV,HILF,
     $                  RHO,DIVGMX,DIVG,WSOR,NBND,
     $                  NFRPER,NBAPER,NRIPER,NLEPER,GEOVP,
     $                  IALGO)
C*STAR******************************************************************
C  V P L E S C    DRUCK UND GESCHWINDIGKEITSKORREKTUR (EINE ITERATION)
C                 (PERIODISCHE RANDBEDINGUNGEN IN X- UND Y-RICHTUNG
C                 SIND MOEGLICH)
C*STAR******************************************************************
C
C  PARAMETER KK, JJ, II           - ARRAYGRENZEN
C            KMX,JMX,IMX          - GRENZE D. BER.-GEB.(MIT BOUND)
C            DX(I),DY(J),DZ(K)    - ABSTAND DER GITTERPUNKTE
C            DDX(I),DDY(J),DDZ(K) - KANTENLAENGE DER KONTROLLVOLUMINA
C            OMBETA               - FAKTOR FUER DRUCKKORREKTUR
C            DT                   - ZEITINKREMENT
C            U(KK,JJ,II)          + GESCHWINDIGKEITSFELD
C            V(KK,JJ,II)          + GESCHWINDIGKEITSFELD
C            W(KK,JJ,II)          + GESCHWINDIGKEITSFELD
C            P(KK,JJ,II)          + DRUCKFELD
C            DIVGMX               + MAXIMALWERT DER DIVERGENZ DIVG IM
C                                   FELD
C            DIVG(KK,JJ)          + DIVERGENZFELD
C            WSOR                 - WICHTUNGSFAKTOR FUER DEN QUELLTERM
C            NBND                 - ANZAHL DER RANDSCHICHTEN
C
C  VERS:  18.03.86 (HW)  : ORIGINAL
C         02.09.86 (HW)  : ERSATZ FUER ISAMAX VON (FB) AUS ABOX86
C                          UEBERNOMMEN
C         24.04.89 (HW)  : SAVE-STATEMENT FUER NBAPER, NLEPER
C         30.10.90 (HW)  : UEBERGABE ERWEITERT
C         01.04.92 (MM)  : AENDERUNG DER SCHLEIFENSTEUERUNG, JETZT
C                          WIRD DURCH KOERPER DURCHITERIERT
C         25. 3.93 (MM)  : PERIODISCHE RANDBEDINGUNGEN WERDEN
C                          DURCH NFRPER,NBAPER,NLEPER UND NRIPER
C                          ERKANNT, UND IM KOPF UEBERGEBEN
C
C UPROG                  : ERRR
C
C  DEFINE-DIREKTIVEN     : KEINE
C
C*STAR******************************************************************
C
C
      INTEGER KK, JJ, II
C
      REAL    OMBETA,
     $        DX(II),DY(JJ),DZ(KK),DDX(II),DDY(JJ),DDZ(KK),
     $        RDX(II),RDY(JJ),RDZ(KK),RDDX(II),RDDY(JJ),RDDZ(KK),
     $        U(KK,JJ,II), V(KK,JJ,II), W(KK,JJ,II), P(KK,JJ,II),
     $        B(KK,JJ,II),HILF(KK,JJ,II),
     $        DT,DIVGMX,DIVG(KK,JJ),BP(KK,JJ,II),BU(KK,JJ,II),
     $        BV(KK,JJ,II),BW(KK,JJ,II),SDIV(KK,JJ,II)

C
C
      IM2  = IMX-NBND
      JM2  = JMX-NBND
      KM2  = KMX-NBND
C
      BETA = OMBETA/WSOR
      RFAK = DT/RHO*WSOR
C
CRED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-
C                                  GERADE IN X ?????
C
      IF (2*(IMX/2) .EQ. IMX ) THEN 
C
C                                 X-RICHTUNG KANN IN RED-BLACK ORDNUNG 
C                                 DURCHLAUFEN WERDEN
           NOUTER = 2
           ISTRIDE= 2

      ELSE

          NOUTER = 1
          ISTRIDE= 1

      ENDIF
CRED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-
C
C
C                                  VORBELEGUNG VON DIVG
C        
      DO 5 J=1,JJ
      DO 5 K=1,KK
    5      DIVG(K,J) = 0.0
C
C
C                                 ERMITTLUNG DER DIVERGENZ DIVG UND DES
C                                 KORREKTURDRCKES DP. KORREKTUR DES DR
C                                 UND GESCHWINDIGKEITSFELDES.ERMITTLUNG
C                                 DER MAXIMALDIVERGENZ IM FELD
      DIVGMX = 0.0
C
CRED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-
C                                 I-SCHLEIFE IN RED-BLACK ORDER
C
C      DO IOUTER = 1,NOUTER
C
C      DO 10 I=IOUTER+NBND,IM2,ISTRIDE
CRED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-
      IF (IALGO.NE.1) THEN
      DO 10 I=1+NBND,IM2

         DO 20 J=1+NBND,JM2
            KSTART = 1+NBND
CDIR$ IVDEP

            DO 30 K=KSTART,KM2,2
C
               FACTOR = BP(K,J,I)
                DIVG(K,J) = ((U(K,J,I) - U(K,J,I-1))*RDDX(I)
     $               + (V(K,J,I) - V(K,J-1,I))*RDDY(J)
     $               + (W(K,J,I) - W(K-1,J,I))*RDDZ(K))
     $               *FACTOR
                DP = BETA * DIVG(K,J)
                DP = DP / ((RDX(I)*RDX(I-1))+(RDY(J)*RDY(J-1))+
     $                     (RDZ(K)*RDZ(K-1)))
C
C                 WITH FACTOR MULTIPLIED (A.O.) 14.08.1996
               P(K,J,I) = P(K,J,I) + DP 
C              P(K,J,I) = P(K,J,I) + DP
               DP = DP * RFAK
C
               U(K,J,I)    =
     $         U(K,J,I  )  + DP*RDX(I  )

               U(K,J,I-1)  =
     $         U(K,J,I-1)  - DP*RDX(I-1)

               V(K,J  ,I)  =
     $         V(K,J  ,I)  + DP*RDY(J  )

               V(K,J-1,I)  =
     $         V(K,J-1,I)  - DP*RDY(J-1)

               W(K  ,J,I)  =
     $         W(K  ,J,I)  + DP*RDZ(K  )

               W(K-1,J,I)  =
     $         W(K-1,J,I)  - DP*RDZ(K-1)
C
   30       CONTINUE
C
            KSTART = KSTART + 1
CDIR$ IVDEP

            DO 35 K=KSTART,KM2,2
C
               FACTOR = BP(K,J,I)
                DIVG(K,J) = ((U(K,J,I) - U(K,J,I-1))*RDDX(I)
     $               + (V(K,J,I) - V(K,J-1,I))*RDDY(J)
     $               + (W(K,J,I) - W(K-1,J,I))*RDDZ(K))
     $               * FACTOR
C
                DP = BETA * DIVG(K,J) 
                DP = DP / ((RDX(I)*RDX(I-1))+(RDY(J)*RDY(J-1))+
     $                     (RDZ(K)*RDZ(K-1)))
C
               P(K,J,I) = P(K,J,I) + DP
C
               DP = DP * RFAK
C
               U(K,J,I)   = U(K,J,I)   + DP*RDX(I)
               U(K,J,I-1)  = U(K,J,I-1)  - DP*RDX(I-1)
               V(K,J,I)   = V(K,J,I)   + DP*RDY(J)
               V(K,J-1,I)  = V(K,J-1,I)  - DP*RDY(J-1)
               W(K,J,I)   = W(K,J,I)   + DP*RDZ(K)
               W(K-1,J,I)  = W(K-1,J,I)  - DP*RDZ(K-1)
C
   35       CONTINUE

         IF(NRIPER.EQ.1) THEN
C                                 PERIODISCHE RANDBED. IN Y-RI.
C                                 -----------------------------
C
            IF(J .EQ. (1+NBND)) THEN
C
               IF(NLEPER .NE. 1) CALL ERRR(501,' VPLESC   ')
C
            KSTART = 1+NBND
C
               DO 100 K=KSTART,KM2
  100             V(K,JM2,I) = V(K,2,I)
C
            ENDIF
C
         ENDIF
C
         KSTART = 1+NBND
         DO 50 K = KSTART,KM2
   50       DIVG(K,1) = AMAX1(DIVG(K,1),ABS(DIVG(K,J)))
C
   20    CONTINUE
C
         IF(NFRPER.EQ.1) THEN
C
C                                 PERIODISCHE RANDBED. IN X-RI.
C                                 -----------------------------
C
      IF(I .EQ. (1+NBND)) THEN
         IF(NBAPER .NE. 1) CALL ERRR(502,' VPLESC   ')
C
         DO 200 J=1+NBND,JM2
            KSTART = 1+NBND
C
            DO 205 K=KSTART,KM2
  205          U(K,J,IM2) = U(K,J,2)
  200    CONTINUE
C
      ENDIF
         ENDIF
C
   10 CONTINUE

CRED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-
C     ENDDO
CRED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-
C
      DO 60 K = 1+NBND,KM2
   60    DIVGMX = AMAX1(DIVGMX,DIVG(K,1))
C
      ELSE
         DO 15 I=1+NBND,IM2

         DO 25 J=1+NBND,JM2
            KSTART = 1+NBND
CDIR$ IVDEP
            DO 40 K=KSTART,KM2,2
C
               FACTOR = BP(K,J,I)
                DIVG(K,J) = ((U(K,J,I)*BU(K,J,I)
     $                - U(K,J,I-1)*BU(K,J,I-1))*RDDX(I)
     $               + (V(K,J,I)*BV(K,J,I)
     $               - V(K,J-1,I)*BV(K,J-1,I))*RDDY(J)
     $               + (W(K,J,I)*BW(K,J,I)
     $              - W(K-1,J,I)*BW(K-1,J,I))*RDDZ(K)+
     $                SDIV(K,J,I))
     $               *FACTOR
                DP = BETA * DIVG(K,J)
                DP = DP / ((RDX(I)*RDX(I-1))+(RDY(J)*RDY(J-1))+
     $                     (RDZ(K)*RDZ(K-1)))
C
C                 WITH FACTOR MULTIPLIED (A.O.) 14.08.1996
               P(K,J,I) = P(K,J,I) + DP
C              P(K,J,I) = P(K,J,I) + DP
               DP = DP * RFAK
C
               U(K,J,I)    =
     $         U(K,J,I  )  + DP*RDX(I  )

               U(K,J,I-1)  =
     $         U(K,J,I-1)  - DP*RDX(I-1)

               V(K,J  ,I)  =
     $         V(K,J  ,I)  + DP*RDY(J  )

               V(K,J-1,I)  =
     $         V(K,J-1,I)  - DP*RDY(J-1)

               W(K  ,J,I)  =
     $         W(K  ,J,I)  + DP*RDZ(K  )

               W(K-1,J,I)  =
     $         W(K-1,J,I)  - DP*RDZ(K-1)
C
   40       CONTINUE
C
            KSTART = KSTART + 1
CDIR$ IVDEP
            DO 45 K=KSTART,KM2,2
C
               FACTOR = BP(K,J,I)
                DIVG(K,J) = ((U(K,J,I)*BU(K,J,I)
     $                - U(K,J,I-1)*BU(K,J,I-1))*RDDX(I)
     $               + (V(K,J,I)*BV(K,J,I)
     $               - V(K,J-1,I)*BV(K,J-1,I))*RDDY(J)
     $               + (W(K,J,I)*BW(K,J,I)
     $              - W(K-1,J,I)*BW(K-1,J,I))*RDDZ(K)+
     $                SDIV(K,J,I))

     $               * FACTOR
C
                DP = BETA * DIVG(K,J)
                DP = DP / ((RDX(I)*RDX(I-1))+(RDY(J)*RDY(J-1))+
     $                     (RDZ(K)*RDZ(K-1)))
C
               P(K,J,I) = P(K,J,I) + DP
C
               DP = DP * RFAK
C
               U(K,J,I)   = U(K,J,I)   + DP*RDX(I)
               U(K,J,I-1)  = U(K,J,I-1)  - DP*RDX(I-1)
               V(K,J,I)   = V(K,J,I)   + DP*RDY(J)
               V(K,J-1,I)  = V(K,J-1,I)  - DP*RDY(J-1)
               W(K,J,I)   = W(K,J,I)   + DP*RDZ(K)
               W(K-1,J,I)  = W(K-1,J,I)  - DP*RDZ(K-1)
C
   45       CONTINUE

         IF(NRIPER.EQ.1) THEN
C                                 PERIODISCHE RANDBED. IN Y-RI.
C                                 -----------------------------
C
            IF(J .EQ. (1+NBND)) THEN
C
               IF(NLEPER .NE. 1) CALL ERRR(501,' VPLESC   ')
C
            KSTART = 1+NBND
C
               DO 105 K=KSTART,KM2
  105             V(K,JM2,I) = V(K,2,I)
C
            ENDIF
C
         ENDIF
C
         KSTART = 1+NBND
         DO 55 K = KSTART,KM2
   55       DIVG(K,1) = AMAX1(DIVG(K,1),ABS(DIVG(K,J)))
C
   25    CONTINUE
C
         IF(NFRPER.EQ.1) THEN
C
C                                 PERIODISCHE RANDBED. IN X-RI.
C                                 -----------------------------
C
      IF(I .EQ. (1+NBND)) THEN
         IF(NBAPER .NE. 1) CALL ERRR(502,' VPLESC   ')
C
         DO 210 J=1+NBND,JM2
            KSTART = 1+NBND
C
            DO 215 K=KSTART,KM2
  215          U(K,J,IM2) = U(K,J,2)
  210    CONTINUE
C
      ENDIF
         ENDIF
C
   15 CONTINUE

CRED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-
C     ENDDO
CRED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-
C
      DO 65 K = 1+NBND,KM2
   65    DIVGMX = AMAX1(DIVGMX,DIVG(K,1))
C
      END IF
 
      RETURN
      END
