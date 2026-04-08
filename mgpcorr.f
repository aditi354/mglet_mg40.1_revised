










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
      SUBROUTINE MGPCORR(KK,JJ,II,KMX,JMX,IMX,
     $                   DX,DY,DZ,DDX,DDY,DDZ,
     $                   RDX,RDY,RDZ,RDDX,RDDY,RDDZ,
     $                   U,V,W,P,B,DP,
     $                   BP,BU,BV,BW,
     $                   RHO,DT,WSOR,NBND,
     $                   NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,IDP)
C*MGLET*****************************************************************
C  M G P C O R R  DRUCKKORREKTUR BEI GEGEBENEM DP
C                 (PERIODISCHE RANDBEDINGUNGEN IN X- UND Y-RICHTUNG
C                 SIND MOEGLICH)
C*MGLET*****************************************************************
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
C            WSOR                 - WICHTUNGSFAKTOR FUER DEN QUELLTERM
C            NBND                 - ANZAHL DER RANDSCHICHTEN
C            IDP                  - IF IDP = 1, THEN DP IS NEW PRESSURE
C
C  VERS:  28.06.93 (MM)  : ORIGINAL AUS VPLESC ABGELEITET
C          2. 8.95 (MM)  : DP CAN CONTAIN NEW PRESSURE, SO THAT
C                          DELTA_P = DP - P
C         14.08.03 (TB)  : CHANGED FOR SIP AND OLD BODY TREATMENT
C
C UPROG                  : ERRR
C
C  DEFINE-DIREKTIVEN     : KEINE
C
C*STAR******************************************************************
c
      IMPLICIT NONE

      INTEGER KK, JJ, II, K, J, I, KMX, JMX, IMX, KM2, JM2, IM2,
     &     NBND, NFRO, NBAC, NRGT, NLFT, NBOT, NTOP, IDP,
     &     NFU,  NRV,  NBW,  NTW,  NBU

      REAL    DX(II),DY(JJ),DZ(KK),DDX(II),DDY(JJ),DDZ(KK),
     $        RDX(II),RDY(JJ),RDZ(KK),RDDX(II),RDDY(JJ),RDDZ(KK)

      REAL    U (KK,JJ,II), V(KK,JJ,II), W(KK,JJ,II), P(KK,JJ,II),
     $        B (KK,JJ,II),DP(KK,JJ,II)

      REAL    BP(KK,JJ,II),BU(KK,JJ,II),BV(KK,JJ,II),BW(KK,JJ,II)

      REAL RFAK, RHO, WSOR, DT

      IM2  = IMX-NBND
      JM2  = JMX-NBND
      KM2  = KMX-NBND

      RFAK = DT/RHO*WSOR

C                                 BEI LOKAL VERFEINERTEM GITTER
C                                 WIRD DIE NORMALKOMPONENTE DER ERSTEN
C                                 ZELLE BERECHNET, ANSONSTEN WIRD SIE IN
C                                 BOUNDMG GESETZT
      NFU = 0
      IF ( NFRO .EQ. 8 .OR. NFRO .EQ. 3 ) NFU = 1
      NRV = 0
      IF ( NRGT .EQ. 8 .OR. NRGT .EQ. 3 ) NRV = 1
      NBW = 0
      IF ( NBOT .EQ. 8 .OR. NBOT .EQ. 3 ) NBW = 1
      NTW = 0
      NBU = 0

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (IDP .EQ. 1) THEN
C
C                                CALCULATION OF DELTA_P
C
      DO I=2,IMX-1
      DO J=2,JMX-1
      DO K=2,KMX-1
  
               DP(K,J,I)    =   DP(K,J,I) - P(K,J,I)
CTBC  200803: BP NOW USED FOR OLD BODY, TOO!
C#ifdef 
               DP(K,J,I)    =   DP(K,J,I)*BP(K,J,I)
C#endif
      ENDDO
      ENDDO
      ENDDO
C
      ENDIF
C
C                                 KORREKTUR DES DRUCKES
C
C                                 PPPPPPPPPPPPPPPPPPPP
C

      DO I=2,IMX-1
      DO J=2,JMX-1
      DO K=2,KMX-1
  
               P(K,J,I)    =   P(K,J,I) + DP(K,J,I)
     $                                  * BP(K,J,I) 
      ENDDO
      ENDDO
      ENDDO
C
C                                 KORREKTUR DES GESCHWINDIGKEITSFELDES
C
C                                 UUUUUUUUUUUUUUUUUUUU
C
      DO I=3-NFU,IMX-2-NBU
      DO J=3,JMX-2
      DO K=3,KMX-2

           U(K,J,I)    =  ( U(K,J,I) + (DP(K,J,I) - DP(K,J,I+1))
     $                     *BU(K,J,I)
     $                     *RDX(I)*RFAK)
      ENDDO
      ENDDO
      ENDDO
C
C
C                                 VVVVVVVVVVVVVVVVVVVV
C
      DO I=3,IMX-2
      DO J=3-NRV,JMX-2
      DO K=3,KMX-2

            V(K,J,I)    =   (V(K,J,I) + (DP(K,J,I) - DP(K,J+1,I))
     $                       *BV(K,J,I)
     $                       *RDY(J)*RFAK)
      ENDDO
      ENDDO
      ENDDO
C
C                                 WWWWWWWWWWWWWWWWWWWW
C
      DO I=3,IMX-2
      DO J=3,JMX-2
      DO K=3-NBW,KMX-2-NTW

            W(K,J,I)    =   (W(K,J,I) + (DP(K,J,I) - DP(K+1,J,I))
     $                       *BW(K,J,I)
     $                       *RDZ(K)*RFAK)
      ENDDO
      ENDDO
      ENDDO
      RETURN
      END
