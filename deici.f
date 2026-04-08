










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
      SUBROUTINE DEICI   (KK,JJ,II,KMX,JMX,IMX,XTOT,YTOT,ZTOT,
     $                    CIDEND,IIDEND,RIDEND,KMXOLD,JMXOLD,IMXOLD,
     $                    UFR,VFR,WFR,PFR,GFR,
     $                    UI1,VI1,WI1,UI2,VI2,WI2,GI1,GI2,
     $                    UGRID,ITSTEP,
     $                    ITTOT,TIMEPH,DT)
C*STARLET***************************************************************
C        D E I C I        EINLESEN DER GESCHWINDIGKEITSKOMPONENTEN U,V,W
C                         UND DER EFFEKTIVEN DYN. VISKOSITAET  G
C                         AM EINTRITTSRAND "FRONT" (FORMATIERT AUF
C                         KANAL 13)
C*STARLET***************************************************************
C
C PARAM: KK,JJ,II       - ARRAYGRENZEN
C
C VERS:  14.01.87 (HW)  : ORIGINAL AUS DEIBI ABGELEITET
C                         NEUE VERSION KOMMT OHNE 'BACKSPACE' AUS.
C        30.12.88 (HW)  : DTEPS  WIRD MASCHINEN"UN"ABHAENGIG BELEGT,
C                         UEBERPRUEFUNGEN WERDEN DURCHGEFUEHRT
C
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
      CHARACTER (LEN=8)     CIDEND(10)
      INTEGER    IIDEND(100)
      REAL       RIDEND(100)
C
CSR35      REAL         U(KK,JJ,II),    V(KK,JJ,II),    W(KK,JJ,II),
CSR35     $             P(KK,JJ,II),    G(KK,JJ,II),
           REAL
     $           UI1(KK,JJ, 2),  VI1(KK,JJ, 2),  WI1(KK,JJ, 2),
     $           UI2(KK,JJ, 2),  VI2(KK,JJ, 2),  WI2(KK,JJ, 2),
     $           GI1(KK,JJ),        GI2(KK,JJ),
     $           UFR(KK,JJ, 2),  VFR(KK,JJ, 2),  WFR(KK,JJ, 2),
     $           PFR(KK,JJ, 2),  GFR(KK,JJ, 2)

C
      SAVE       NCOUNT
      DATA       NCOUNT                     /0/
C
      SAVE       TP1,TP2
C
      JM2    = JMX - 2
      KM2    = KMX - 2
C
      DTEPS  = AMAX1 ((DT / 10000.0) , (100.0 * SMAONE))
      DTVH   = DT / DTEPS
      IF((DTVH .LT. 100.0)  .AND. NCOUNT .EQ. 0) THEN
         WRITE (6,*)
         WRITE (6,*)
         WRITE (6,*) ' ! ! ! ! !   W A R N U N G  AUS SUBR. DEICI ',
     $               ' ! ! ! ! !'
         WRITE (6,*)
         WRITE (6,*) ' DAS VERHAELTNIS  DT / DT_EPS  BETRAEGT : ',
     $               DTVH
         WRITE (6,*) ' DT = ',DT,'   DTEPS = ',DTEPS,'   SMAONE = ',
     $               SMAONE
         WRITE (6,*) ' DAS RICHTIGE EINLESEN DER GESCHWINDIGKEITS',
     $               'PROFILE AM EINTRITT SOLLTE UEBERPRUEFT WERDEN'
         NCOUNT = NCOUNT + 1
      END IF
C
      IF(ITSTEP .GT. 0) GOTO 2100
C
C                                  HIER: DER DATENSATZ WIRD AN DEN
C                                  ANFANG GESETZT, WICHTIGE KENNDATEN
C                                  WERDEN EINGELESEN.
C
 2000    REWIND 13
         READ(13,6010) (CIDEND(N),N=1,10)
         READ(13,6020) (IIDEND(N),N=1,100)
         READ(13,6030) (RIDEND(N),N=1,100)
         READ(13,6040)  KKOLD,JJOLD,IIOLD,KMXOLD,JMXOLD,IMXOLD
C
C
C                                  PLAUSIBILITAETSTEST
C
         IF(KMXOLD     .NE. KMX ) CALL ERRR (501,' DEICI    ')
         IF(JMXOLD     .NE. JMX ) CALL ERRR (502,' DEICI    ')
         IF(RIDEND(10) .NE. ZTOT) CALL ERRR (503,' DEICI    ')
         IF(RIDEND(11) .NE. YTOT) CALL ERRR (504,' DEICI    ')
C
C                                  BEREITSTELLUNG DER ZYKLUSZEIT UND
C                                  DES ZEITSCHRITTES, DER BEI DER
C                                  DATENERZEUGUNG VERWENDET WURDE
C
         TCYCLE  = RIDEND(31)
         DTD     = RIDEND( 2)
C
         IF(INT((TCYCLE+DTEPS)/DTD) .LT. 2) CALL ERRR (510,' DEICI    ')
C
C                                  EINLESEN DER ERSTEN BEIDEN DATEN-
C                                  SAETZE IN DIE PUFFERFELDER P1 UND P2
C
         READ(13,6030) UI1,VI1,WI1,GI1
         TP1    = 0.0

         READ(13,6030)  UI2,VI2,WI2,GI2
         TP2    = DTD

C
C                                  HIER: DIE BENOETIGTE ZEITEBENE LIEGT
C                                  ZWISCHEN DEN BEIDEN ZEITEBENEN TP1
C                                  UND TP2. ES GILT:
C                                  TP1 .LE. TIMEPD .LT. TP2
C
C                                  DER ZUM ZEITPUNKT TP1 GEHOERIGE
C                                  DATENSATZ IST IM PUFFERFELD  P1
C                                  ABGELEGT, DER ZUM ZEITPUNKT TP2
C                                  GEHOERIGE DATENSATZ IST IN  P2.
C
 2100 TCYCLE  = RIDEND(31)
      DTD     = RIDEND( 2)
C
      TIMEPD  = MOD(TIMEPH+DTEPS,TCYCLE)
C
C                                  FALLS TIMEPD ZWISCHEN TP1 UND TP2
C                                  LIEGT, KANN SOFORT INTERPOLIERT
C                                  WERDEN. ES MUESSEN KEINE DATEN
C                                  VON TAPE 13 GELESEN WERDEN.
C
      IF((TIMEPD .GE. TP1) .AND. (TIMEPD .LT. TP2)) GOTO 2999
C
      NREAD  = INT((TIMEPD-TP2)/DTD) + 1
C
C                                  IST DIE ANZAHL DER ZU LESENDEN
C                                  DATENSAETZE 'NREAD' KLEINER OD.
C                                  GLEICH NULL, DANN WIRD VON VORNE
C                                  GELESEN
C
      IF(NREAD .LE. 0) GOTO 2000
C
         DO 110 J = 1,JMXOLD
            DO 120 K = 1,KMXOLD
               UI1(K,J,1) = UI2(K,J,I)
               UI1(K,J,2) = UI2(K,J,I)
               VI1(K,J,1) = VI2(K,J,I)
               VI1(K,J,2) = VI2(K,J,I)
               WI1(K,J,1) = WI2(K,J,I)
               WI1(K,J,2) = WI2(K,J,I)
               GI1(K,J  ) = GI2(K,J  )
  120       CONTINUE
  110    CONTINUE
C
      TP1    = TP2
C
      IF(NREAD .GT. 1) THEN
         DO 150 N = 1,NREAD-1
            READ(13,6030) UI1,VI1,WI1,GI1
  150       TP1   = TP1 + DTD
      ENDIF
C
      READ(13,6030,END=9999) UI2,VI2,WI2,GI2
      TP2    = TP1 + DTD
C
 2999 CONTINUE
C
      IF(TIMEPD .LT. TP1 .OR. TIMEPD .GE. TP2) THEN
C
C                                  IM FEHLERFALL:
C
         WRITE (6,6000) TP1, TIMEPD, TP2, ITSTEP, ITTOT, TIMEPH,
     $                  TCYCLE, DT
         CALL ERRR (515,' DEICI    ')
      ENDIF
C
C                                  DIE PUFFER SIND MIT DEN RICHTIGEN
C                                  WERTEN GEFUELLT, ES WIRD INTERPOLIERT
C
      AMULT   = (TIMEPD-DTEPS-TP1)/DTD
C
      DO 300 I=1,2
         DO 310 J=1,JMX
            DO 320 K=1,KMX
               UFR(K,J,I) = UI1(K,J,I)+(UI2(K,J,I)-UI1(K,J,I))*AMULT
               VFR(K,J,I) = VI1(K,J,I)+(VI2(K,J,I)-VI1(K,J,I))*AMULT
  320          WFR(K,J,I) = WI1(K,J,I)+(WI2(K,J,I)-WI1(K,J,I))*AMULT
  310    CONTINUE
  300 CONTINUE
C
      DO 350 J=1,JMX
         DO 360 K=1,KMX
            GFR(K,J,  2) = GI1(K,J)+(GI2(K,J)-GI1(K,J))*AMULT
  360       GFR(K,J,  1) = GFR(K,J,  2)
  350 CONTINUE
C
C                                  GALILEI-TRANSFORMATION UND
C                                  BELEGUNG DES  UFR-FELDES
C                                  (NOTWENDIG Z. SETZEN DER RANDBED.)
C
      DO 600 I=1,2
         DO 610 J=1,JMX
            DO 620 K=1,KMX
  620          UFR(K,J,I) = UFR(K,J,I) - UGRID
  610    CONTINUE
  600 CONTINUE
C
CSR35      IF(ITSTEP .GT. 0) THEN
CSR35         DO 650 I = 1,2
CSR35            DO 660 J = 3,JM2
CSR35               DO 670 K = 3,KM2
CSR35  670             U(K,J,I) = UFR(K,J,I)
CSR35  660       CONTINUE
CSR35  650    CONTINUE
CSR35      ENDIF
C
      RETURN
C
 9999 WRITE (6,6000) TP1, TIMEPD, TP2, ITSTEP, ITTOT, TIMEPH,
     $               TCYCLE, DT
      STOP ' DEICI, EOF ERREICHT!'
C
 6000 FORMAT(//,1X,'  FEHLER IN SUBR. DEICI !',/,1X,
     $       '  TP1    = ',1PE15.8,/,1X,'  TIMEPD = ',1PE15.8,/,1X,
     $       '  TP2    = ',1PE15.8,/,1X,'  ITSTEP = ',I7,/,1X,
     $       '  ITTOT  = ',I7,/,1X,'  TIMEPH = ',1PE15.8,/,1X,
     $       '  TCYCLE = ',1PE15.8,/,1X,'  DT     = ',1PE15.8,/)
 6010 FORMAT (5(A8,2X))
 6020 FORMAT (5(I9,1X))
 6030 FORMAT (4(E19.12E3,1X))
 6040 FORMAT (6(I9,1X))
      END
