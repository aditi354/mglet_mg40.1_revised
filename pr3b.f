










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
      SUBROUTINE PR3B(KK,JJ,II,V,IV,KSTAG,JSTAG,ISTAG,IDUMMY,
     $                IC1,ICS,IC2,M1,MS,M2,N1,NS,N2,DX,DY,DZ,X,Y,Z,
     $                KANAL,ILPRLE,NPR,ISTPR,JSTPR,KSTPR,OPFORM)
C***********************************************************************
C   P R 2 M    DRUCKT BLOECKE AUS 3D-MATRIX
C*************************************************** J.FERSTL 19.11.82 *
C                                                    F.BAETKE 05.10.84 *
C                                                    M.M      01.03.92 *
C                                                    M.M      10.03.92 *
C
C  PARAMETER  KK,JJ,II  - FELDGRENZEN
C             V         - VARIABLE DER FORM V(KK,JJ,II)
C             IV        - CHARACTER-KONSTANTE (CHARACTER (LEN=16)) ZUR IDEN-
C                         TIFIKATION DES AUSZUGEBENDEN FELDES
C             IC        - CHARACTER-KONSTANTE DER FORM 'K','J' ODER 'I'
C
C             NPR      - ANZAHL DER AUSZUDRUCKENDEN BLOECKE MIT 4*4*4
C                        PUNKTEN
C             ISTPR    - I-INDEX DER STARTPUNKTE (BIS ZU 8 INDIZES)
C             JSTPR    - J-INDEX DER STARTPUNKTE (BIS ZU 8 INDIZES)
C             KSTPR    - K-INDEX DER STARTPUNKTE (BIS ZU 8 INDIZES)
C
C             CSP1     - CHARACTER-FELD ZUM PUFFERN EINES OUTPUT-BLOCKES
C             CSP2     - CHARACTER-FELD ZUM PUFFERN EINES OUTPUT-BLOCKES
C
C
C  VERS: 02.09.86 (HW)  : KSTAG,JSTAG,ISTAG EINGEFUEHRT
C  VERS: 25.02.92 (MM)  : KANAL  EINGEFUEHRT
C  VERS: 01.03.92 (MM)  : VOELLIG NEUER OUTPUT
C
C  UPROG                : ERRR, PR3
C
C***********************************************************************
      CHARACTER (LEN=1)   IC,  IK
      CHARACTER (LEN=16)  IV
      CHARACTER (LEN=16)  OPFORM
      CHARACTER (LEN=40)  CSP1(10),CSP2(10)
      INTEGER ISTPR(8),JSTPR(8),KSTPR(8)
      INTEGER KK,JJ,II,IC1,ICS,IC2,M1,MS,M2,N1,NS,N2
      REAL    V(KK,JJ,II), DX(II), DY(JJ), DZ(KK), X(II), Y(JJ), Z(KK)

C     DATA  NPR      /  8  /
C     DATA  ISTPR/   1,   13,    1,   13,    1,   13,    1,   13/
C     DATA  JSTPR/   1,    1,   13,   13,    1,    1,   13,   13/
C     DATA  KSTPR/   1,    1,    1,    1,    9,    9,    9,    9/
C     DATA  NPR      /  4  /
C     DATA  ISTPR/   3,   11,    3,   11,    3,   11,    3,   11/
C     DATA  JSTPR/   3,    3,   11,   11,    3,    3,   11,   11/
C     DATA  KSTPR/   3,    3,    3,    3,    9,    9,    9,    9/
C     DATA  NPR      /  4  /
C     DATA  ISTPR/   5,    9,    5,    9,    1,   13,    1,   13/
C     DATA  JSTPR/   5,    5,    9,    9,    1,    1,   13,   13/
C     DATA  KSTPR/   4,    4,    4,    4,    9,    9,    9,    9/

      IC='J'

      DO 10 LPR=1,NPR

       IF (OPFORM(1:5).EQ.'SHORT') LPAGE=56
       IF (OPFORM(1:4).EQ.'LONG') LPAGE=48
	 IF (ILPRLE.GE.LPAGE) THEN
	       DO 11 I=1,66-ILPRLE
	       WRITE(KANAL,*)
   11            CONTINUE
	       ILPRLE = 0
         ENDIF

	 IP1 = ISTPR(LPR)

C                               AB HIER OUTPUT-FORMAT 'SHORT'

      IF(OPFORM(1:5).EQ.'SHORT') THEN

      DO 20 JPR = 1,4,2
	 JP1 = JPR-1+JSTPR(LPR)
	 JP2 = JPR  +JSTPR(LPR)
         WRITE(KANAL,6010) IV,IC,JP1,IV,IC,JP2

       DO 20 K=4,1,-1
	      KV=K-1+KSTPR(LPR)
	      WRITE(KANAL,6000) KV
     $                         ,(V(KV,JP1,I),I=IP1,IP1+3)
     $                         ,(V(KV,JP2,I),I=IP1,IP1+3)

   20 CONTINUE

	      WRITE(KANAL,6020) (I,I=IP1,IP1+3),(I,I=IP1,IP1+3)
      ILPRLE=ILPRLE+14

C                               AB HIER OUTPUT-FORMAT 'LONG'

      ELSEIF(OPFORM(1:4).EQ.'LONG') THEN

      DO 21 JPR = 1,4
	 JP1 = JPR-1+JSTPR(LPR)
	 JP2 = JPR  +JSTPR(LPR)
         WRITE(KANAL,6010) IV,IC,JP1

       DO 21 K=4,1,-1
	      KV=K-1+KSTPR(LPR)
	      WRITE(KANAL,7000) KV
     $                         ,(V(KV,JP1,I),I=IP1,IP1+3)

   21 CONTINUE

	      WRITE(KANAL,7020) (I,I=IP1,IP1+3)
      ILPRLE=ILPRLE+24

      ENDIF

   10 CONTINUE
C
C
 6000 FORMAT(1X,I3,2X,4(F8.5),4X,4(F8.5))
 7000 FORMAT(1X,I3,2X,4(E16.9E2))
 6010    FORMAT(1X,5(1H-),1X,A10,'-FELD ',4H -- ,A1,3H = ,I3,
     $          1X,5(1H-),1X,A10,'-FELD ',4H -- ,A1,3H = ,I3)
 6020 FORMAT(/,1X,'  K/ I:',I3,3I8,4X,4I8,/)
 7020 FORMAT(/,1X,'  K/ I:',I3,3I16,/)
 6100 FORMAT(/,3X,'Y(J)',3X,'DY(J)',6X,16(2X,I3,2X))
 6200 FORMAT(/,3X,'Z(K)',3X,'DZ(K)',6X,16(2X,I3,2X))
 6400 FORMAT(/,15X,'DX(I) ',16F7.5)
 6500 FORMAT(16X,'X(I) ',16F7.4)
 6600 FORMAT(/,15X,'DY(J) ',16F7.5)
 6700 FORMAT(16X,'Y(J) ',16F7.4)
 6800 FORMAT(/)
C
      RETURN
      END
