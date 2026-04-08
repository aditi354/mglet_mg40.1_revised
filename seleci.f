










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
      SUBROUTINE SELECI  (KK,JJ,II,KKA,JJA,IIA,ISELEP,KKXL,KKYL,KKZL,
     $                    KSXL,KSYL,KSZL,KH0L,JJ1L,JJ2L,KLASS,ILIMXP,
     $                    ISLIDI,IDIMF,NDGL,NMAX,NRZUL)
C*STARLET***************************************************************
C        S E L E C I      NACH DEM DURCHLAUFEN VON SUBR. SELECA WIRD
C                         DEM BENUTZER EINE INFORMATION UEBER DIE FELDER
C                         AUSGEGEBEN, DIE KLEINER DIMENSIONIERT WERDEN
C                         KOENNEN.
C                         AUSSERDEM WIRD DER BENOETIGTE SPEICHERPLATZ
C                         ERMITTELT (UNGEFAEHRE ANGABE), DER SICH BEI
C                         OPTIMALER FELDDIMENSIONIERUNG ERGIBT.
C  A C H T U N G:         DIESE ROUTINE MUSS ERWEITERT WERDEN, FALLS
C                         DIE STATIST. AUSWERTUNG ERWEITERT WIRD !
C                         DIE ZUSAETZLICHEN PROGRAMMSTATEMENTS SIND
C                         ENTSPRECHEND DEM SPEICHERPLATZ DER BE-
C                         TREFFENDEN GROESSE IM ISELEP-FELD EINZU-
C                         ORDNEN
C*STARLET***************************************************************
C
C PARAM: KK ,JJ ,II     - ARRAYDIMENSIONEN
C        KKA,JJA,IIA    - ARRAYDIMENSIONEN DER AUSWERTEFELDER
C        ISELEP(2,752)  - STEUERFELD FUER DIE AUSWERTUNG
C
C UPROG                 : KEINE
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        05.08.86 (HW)  : ORIGINAL
C        29.08.86 (HW)  : ERWEITERUNG: BEI DEN KOMPONENTEN DES SCHUB-
C                         SPANNUNGSTENSORS KANN DIE SUMME AUS GROB-
C                         UND FEINSTRUKTURANTEIL (ENSEMBLE-MITTELWERT
C                         DER REYNOLDSSPANNUNG) AUSGEGEBEN WERDEN.
C                         DIE SUBR. SELECI IST DAVON NICHT BETROFFEN.
C        30.06.88 (HW)  : ERWEITERUNG: AUSWERTUNG DER DREI KOMPONENTEN
C                         DER VORTICITY, DER ENSTROPHIE UND DER
C                         HELIZITAET
C        02.12.88 (HW)  : DIMENSIONIERUNG DES ISELEA- UND
C                         ISELEP-FELDES VON 250 --> 752
C        06.12.88 (HW)  : ERWEITERUNG: AUSWERTUNG VON KORRELATIONS-
C                         FUNKTIONEN, -KOEFFIZIENTEN, LEISTUNGS-
C                         DICHTESPEKTREN UND HAEUFIGKEITSVERTEILUNGEN
C                         A C H T U N G:  UEBERGABE GEAENDERT
C        24.02.89 (HW)  : ISLINI(ILIMXP) --> ISLINI(ISLIDI)
C        21.04.89 (HW)  : ERWEITERUNGEN WEGEN NEUEM GITTERGENERATOR
C                         GRDFMI
C        12.02.03 (TB)  : NO CALL OF SELECI IN >= MG38
C
C*STARLET***************************************************************
C
      CHARACTER *7  FNAM
C
      INTEGER       ISELEP(2,752)
C
C                                 NBSP WIE IM PARAMETER-STATEMENT DES
C                                 HAUPTPROGRAMMES DEFINIERT
C
      NBSP    = 6
C
      MEMORY  = KK    * JJ    * II
      MEMA    = KKA   * JJA   * IIA
C
C                                 KORRELATIONSFUNKTIONEN UND -KOEFF.
C
      MEMKXA  = KKXL  * JJ2L  * ILIMXP
      MEMKXS  = KKXL  * JJ1L  * ILIMXP
      MEMKYA  = KKYL  * JJ2L  * ILIMXP
      MEMKYS  = KKYL  * JJ1L  * ILIMXP
      MEMKZA  = KKZL  * JJ2L  * ILIMXP
      MEMKZS  = KKZL  * JJ1L  * ILIMXP
C
C                                 LEISTUNGSDICHTESPEKTREN
C
      MEMSX   = KSXL  * JJ1L  * ILIMXP
      MEMSY   = KSYL  * JJ1L  * ILIMXP
      MEMSZ   = KSZL  * JJ1L  * ILIMXP
C
C                                 HAEUFIGKEITSVERTEILUNGEN
C
      MEMH   = KH0L  * JJ2L  * ILIMXP
C
C                                 MAXIMALE DIMENSION DER LINIEN-FELDER
C
      KKML   = MAX0 (KKXL, KKYL, KKZL, KSXL, KSYL, KSZL, KH0L)
      MEMML  = KKML  * JJ2L  * ILIMXP
C
      MEMSUC  = 0
      MEMSUI  = 0
      MEMSUR  = 0
C
      WRITE (6,6000)
C
C                                 AUSGABE DES ISELEP(1,...)-FELDES
C
      WRITE (6,6001)
C
      DO 10 I = 1,752,8
         WRITE (6,6005) I,ISELEP(1,I),I+1,ISELEP(1,I+1),I+2,
     $                  ISELEP(1,I+2),I+3,ISELEP(1,I+3),I+4,
     $                  ISELEP(1,I+4),I+5,ISELEP(1,I+5),I+6,
     $                  ISELEP(1,I+6),I+7,ISELEP(1,I+7)
   10 CONTINUE
C
      WRITE (6,6006) KK, JJ, II, KKA, JJA, IIA,KKXL, KKYL, KKZL,
     $               KSXL, KSYL, KSZL, KH0L, JJ1L, JJ2L, KLASS,
     $               ILIMXP, NMAX
C
C                                 GROESSEN DER U-KOMPONENTE
C                                 -------------------------
C
      IF(ISELEP(1,  1) .EQ. 0) THEN
         FNAM ='U      '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
      IF(ISELEP(1,  2) .EQ. 0) THEN
         FNAM ='UO     '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
*     IF(ISELEP(1,  3) .EQ. 0) THEN
*        FNAM ='UOO    '
*        WRITE (6,6010) FNAM, FNAM
*     ELSE
*        MEMSUR = MEMSUR + MEMORY
*     ENDIF
      IF(ISELEP(1,  4) .EQ. 0) THEN
         FNAM ='AU     '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,  5) .EQ. 0) THEN
         FNAM ='SU     '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,  6) .EQ. 0) THEN
         FNAM ='UFG    '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
      IF(ISELEP(1,  7) .EQ. 0) THEN
         FNAM ='AURG   '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,  8) .EQ. 0) THEN
         FNAM ='SURG   '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 10) .EQ. 0) THEN
         FNAM ='AUSKG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 11) .EQ. 0) THEN
         FNAM ='SUSKG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 12) .EQ. 0) THEN
         FNAM ='AUFLG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 13) .EQ. 0) THEN
         FNAM ='SUFLG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 14) .EQ. 0) THEN
         FNAM ='ARXUU  '
         WRITE (6,6042) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXA
      ENDIF
      IF(ISELEP(1, 15) .EQ. 0) THEN
         FNAM ='SRXUU  '
         WRITE (6,6041) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXS
      ENDIF
      IF(ISELEP(1, 16) .EQ. 0) THEN
         FNAM ='ARYUU  '
         WRITE (6,6044) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYA
      ENDIF
      IF(ISELEP(1, 17) .EQ. 0) THEN
         FNAM ='SRYUU  '
         WRITE (6,6043) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYS
      ENDIF
      IF(ISELEP(1, 18) .EQ. 0) THEN
         FNAM ='ARZUU  '
         WRITE (6,6046) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZA
      ENDIF
      IF(ISELEP(1, 19) .EQ. 0) THEN
         FNAM ='SRZUU  '
         WRITE (6,6045) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZS
      ENDIF
      IF(ISELEP(1, 20) .EQ. 0) THEN
         FNAM ='ASXUU  '
         WRITE (6,6050) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSX
      ENDIF
      IF(ISELEP(1, 21) .EQ. 0) THEN
         FNAM ='SSXUU  '
         WRITE (6,6050) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSX
      ENDIF
      IF(ISELEP(1, 22) .EQ. 0) THEN
         FNAM ='ASYUU  '
         WRITE (6,6052) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSY
      ENDIF
      IF(ISELEP(1, 23) .EQ. 0) THEN
         FNAM ='SSYUU  '
         WRITE (6,6052) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSY
      ENDIF
      IF(ISELEP(1, 24) .EQ. 0) THEN
         FNAM ='ASZUU  '
         WRITE (6,6054) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSZ
      ENDIF
      IF(ISELEP(1, 25) .EQ. 0) THEN
         FNAM ='SSZUU  '
         WRITE (6,6054) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSZ
      ENDIF
C
C                                 GROESSEN DER V-KOMPONENTE
C                                 -------------------------
C
      IF(ISELEP(1, 26) .EQ. 0) THEN
         FNAM ='V      '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
      IF(ISELEP(1, 27) .EQ. 0) THEN
         FNAM ='VO     '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
*     IF(ISELEP(1, 28) .EQ. 0) THEN
*        FNAM ='VOO    '
*        WRITE (6,6010) FNAM, FNAM
*     ELSE
*        MEMSUR = MEMSUR + MEMORY
*     ENDIF
      IF(ISELEP(1, 29) .EQ. 0) THEN
         FNAM ='AV     '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 30) .EQ. 0) THEN
         FNAM ='SV     '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 31) .EQ. 0) THEN
         FNAM ='VFG    '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
      IF(ISELEP(1, 32) .EQ. 0) THEN
         FNAM ='AVRG   '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 33) .EQ. 0) THEN
         FNAM ='SVRG   '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 35) .EQ. 0) THEN
         FNAM ='AVSKG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 36) .EQ. 0) THEN
         FNAM ='SVSKG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 37) .EQ. 0) THEN
         FNAM ='AVFLG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 38) .EQ. 0) THEN
         FNAM ='SVFLG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 39) .EQ. 0) THEN
         FNAM ='ARXVV  '
         WRITE (6,6042) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXA
      ENDIF
      IF(ISELEP(1, 40) .EQ. 0) THEN
         FNAM ='SRXVV  '
         WRITE (6,6041) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXS
      ENDIF
      IF(ISELEP(1, 41) .EQ. 0) THEN
         FNAM ='ARYVV  '
         WRITE (6,6044) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYA
      ENDIF
      IF(ISELEP(1, 42) .EQ. 0) THEN
         FNAM ='SRYVV  '
         WRITE (6,6043) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYS
      ENDIF
      IF(ISELEP(1, 43) .EQ. 0) THEN
         FNAM ='ARZVV  '
         WRITE (6,6046) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZA
      ENDIF
      IF(ISELEP(1, 44) .EQ. 0) THEN
         FNAM ='SRZVV  '
         WRITE (6,6045) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZS
      ENDIF
      IF(ISELEP(1, 45) .EQ. 0) THEN
         FNAM ='ASXVV  '
         WRITE (6,6050) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSX
      ENDIF
      IF(ISELEP(1, 46) .EQ. 0) THEN
         FNAM ='SSXVV  '
         WRITE (6,6050) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSX
      ENDIF
      IF(ISELEP(1, 47) .EQ. 0) THEN
         FNAM ='ASYVV  '
         WRITE (6,6052) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSY
      ENDIF
      IF(ISELEP(1, 48) .EQ. 0) THEN
         FNAM ='SSYVV  '
         WRITE (6,6052) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSY
      ENDIF
      IF(ISELEP(1, 49) .EQ. 0) THEN
         FNAM ='ASZVV  '
         WRITE (6,6054) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSZ
      ENDIF
      IF(ISELEP(1, 50) .EQ. 0) THEN
         FNAM ='SSZVV  '
         WRITE (6,6054) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSZ
      ENDIF
C
C                                 GROESSEN DER W-KOMPONENTE
C                                 -------------------------
C
      IF(ISELEP(1, 51) .EQ. 0) THEN
         FNAM ='W      '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
      IF(ISELEP(1, 52) .EQ. 0) THEN
         FNAM ='WO     '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
*     IF(ISELEP(1, 53) .EQ. 0) THEN
*        FNAM ='WOO    '
*        WRITE (6,6010) FNAM, FNAM
*     ELSE
*        MEMSUR = MEMSUR + MEMORY
*     ENDIF
      IF(ISELEP(1, 54) .EQ. 0) THEN
         FNAM ='AW     '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 55) .EQ. 0) THEN
         FNAM ='SW     '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 56) .EQ. 0) THEN
         FNAM ='WFG    '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
      IF(ISELEP(1, 57) .EQ. 0) THEN
         FNAM ='AWRG   '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 58) .EQ. 0) THEN
         FNAM ='SWRG   '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 60) .EQ. 0) THEN
         FNAM ='AWSKG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 61) .EQ. 0) THEN
         FNAM ='SWSKG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 62) .EQ. 0) THEN
         FNAM ='AWFLG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 63) .EQ. 0) THEN
         FNAM ='SWFLG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 64) .EQ. 0) THEN
         FNAM ='ARXWW  '
         WRITE (6,6042) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXA
      ENDIF
      IF(ISELEP(1, 65) .EQ. 0) THEN
         FNAM ='SRXWW  '
         WRITE (6,6041) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXS
      ENDIF
      IF(ISELEP(1, 66) .EQ. 0) THEN
         FNAM ='ARYWW  '
         WRITE (6,6044) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYA
      ENDIF
      IF(ISELEP(1, 67) .EQ. 0) THEN
         FNAM ='SRYWW  '
         WRITE (6,6043) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYS
      ENDIF
      IF(ISELEP(1, 68) .EQ. 0) THEN
         FNAM ='ARZWW  '
         WRITE (6,6046) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZA
      ENDIF
      IF(ISELEP(1, 69) .EQ. 0) THEN
         FNAM ='SRZWW  '
         WRITE (6,6045) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZS
      ENDIF
      IF(ISELEP(1, 70) .EQ. 0) THEN
         FNAM ='ASXWW  '
         WRITE (6,6050) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSX
      ENDIF
      IF(ISELEP(1, 71) .EQ. 0) THEN
         FNAM ='SSXWW  '
         WRITE (6,6050) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSX
      ENDIF
      IF(ISELEP(1, 72) .EQ. 0) THEN
         FNAM ='ASYWW  '
         WRITE (6,6052) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSY
      ENDIF
      IF(ISELEP(1, 73) .EQ. 0) THEN
         FNAM ='SSYWW  '
         WRITE (6,6052) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSY
      ENDIF
      IF(ISELEP(1, 74) .EQ. 0) THEN
         FNAM ='ASZWW  '
         WRITE (6,6054) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSZ
      ENDIF
      IF(ISELEP(1, 75) .EQ. 0) THEN
         FNAM ='SSZWW  '
         WRITE (6,6054) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSZ
      ENDIF
C
C                                 GROESSEN DER P-KOMPONENTE
C                                 -------------------------
C
      IF(ISELEP(1, 76) .EQ. 0) THEN
         FNAM ='P      '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
*     IF(ISELEP(1, 77) .EQ. 0) THEN
*        FNAM ='PO     '
*        WRITE (6,6010) FNAM, FNAM
*     ELSE
*        MEMSUR = MEMSUR + MEMORY
*     ENDIF
*     IF(ISELEP(1, 78) .EQ. 0) THEN
*        FNAM ='POO    '
*        WRITE (6,6010) FNAM, FNAM
*     ELSE
*        MEMSUR = MEMSUR + MEMORY
*     ENDIF
      IF(ISELEP(1, 79) .EQ. 0) THEN
         FNAM ='AP     '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 80) .EQ. 0) THEN
         FNAM ='SP     '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 81) .EQ. 0) THEN
         FNAM ='PFG    '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
      IF(ISELEP(1, 82) .EQ. 0) THEN
         FNAM ='APRG   '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 83) .EQ. 0) THEN
         FNAM ='SPRG   '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 85) .EQ. 0) THEN
         FNAM ='APSKG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 86) .EQ. 0) THEN
         FNAM ='SPSKG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 87) .EQ. 0) THEN
         FNAM ='APFLG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1, 88) .EQ. 0) THEN
         FNAM ='SPFLG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
C
C                                 GROESSEN DER G-KOMPONENTE
C                                 -------------------------
C
      IF(ISELEP(1,101) .EQ. 0) THEN
         FNAM ='G      '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
C
C                                 KINETISCHE ENERGIE DER SCHWANKUNGS-
C                                 GESCHWINDIGKEITEN
C                                 -----------------------------------
C
      IF(ISELEP(1,106) .EQ. 0) THEN
         FNAM ='AEFG   '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,107) .EQ. 0) THEN
         FNAM ='SEFG   '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,109) .EQ. 0) THEN
         FNAM ='AEFS   '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,110) .EQ. 0) THEN
         FNAM ='SEFS   '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
C
C                                 KOMPONENTEN DES SPANNUNGSTENSORS
C                                 --------------------------------
C
C                                 HIER: U - W
C                                 -----------
C
      IF(ISELEP(1,146) .EQ. 0) THEN
         FNAM ='AUFWFG '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,147) .EQ. 0) THEN
         FNAM ='SUFWFG '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,149) .EQ. 0) THEN
         FNAM ='AUFWFS '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,150) .EQ. 0) THEN
         FNAM ='SUFWFS '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,151) .EQ. 0) THEN
         FNAM ='AUFWFM '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,152) .EQ. 0) THEN
         FNAM ='SUFWFM '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
C
C                                 HIER: V - W
C                                 -----------
C
      IF(ISELEP(1,156) .EQ. 0) THEN
         FNAM ='AVFWFG '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,157) .EQ. 0) THEN
         FNAM ='SVFWFG '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,159) .EQ. 0) THEN
         FNAM ='AVFWFS '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,160) .EQ. 0) THEN
         FNAM ='SVFWFS '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,161) .EQ. 0) THEN
         FNAM ='AVFWFM '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,162) .EQ. 0) THEN
         FNAM ='SVFWFM '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
C
C                                 HIER: U - V
C                                 -----------
C
      IF(ISELEP(1,166) .EQ. 0) THEN
         FNAM ='AUFVFG '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,167) .EQ. 0) THEN
         FNAM ='SUFVFG '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,169) .EQ. 0) THEN
         FNAM ='AUFVFS '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,170) .EQ. 0) THEN
         FNAM ='SUFVFS '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,171) .EQ. 0) THEN
         FNAM ='AUFVFM '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,172) .EQ. 0) THEN
         FNAM ='SUFVFM '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
C
C                                 HIER: OMEGA-X
C                                 -------------
C
      IF(ISELEP(1,176) .EQ. 0) THEN
         FNAM ='OX     '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
      IF(ISELEP(1,177) .EQ. 0) THEN
         FNAM ='AOX    '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,178) .EQ. 0) THEN
         FNAM ='SOX    '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,179) .EQ. 0) THEN
         FNAM ='OXFG   '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
      IF(ISELEP(1,180) .EQ. 0) THEN
         FNAM ='AOXRG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,181) .EQ. 0) THEN
         FNAM ='SOXRG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
C
C                                 HIER: OMEGA-Y
C                                 -------------
C
      IF(ISELEP(1,182) .EQ. 0) THEN
         FNAM ='OY     '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
      IF(ISELEP(1,183) .EQ. 0) THEN
         FNAM ='AOY    '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,184) .EQ. 0) THEN
         FNAM ='SOY    '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,185) .EQ. 0) THEN
         FNAM ='OYFG   '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
      IF(ISELEP(1,186) .EQ. 0) THEN
         FNAM ='AOYRG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,187) .EQ. 0) THEN
         FNAM ='SOYRG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
C
C                                 HIER: OMEGA-Z
C                                 -------------
C
      IF(ISELEP(1,188) .EQ. 0) THEN
         FNAM ='OZ     '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
      IF(ISELEP(1,189) .EQ. 0) THEN
         FNAM ='AOZ    '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,190) .EQ. 0) THEN
         FNAM ='SOZ    '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,191) .EQ. 0) THEN
         FNAM ='OZFG   '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
      IF(ISELEP(1,192) .EQ. 0) THEN
         FNAM ='AOZRG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,193) .EQ. 0) THEN
         FNAM ='SOZRG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
C
C                                 HIER: ENSTROPHIE
C                                 ----------------
C
      IF(ISELEP(1,195) .EQ. 0) THEN
         FNAM ='AO2    '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,196) .EQ. 0) THEN
         FNAM ='SO2    '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,197) .EQ. 0) THEN
         FNAM ='O2FG   '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
      IF(ISELEP(1,198) .EQ. 0) THEN
         FNAM ='AO2RG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,199) .EQ. 0) THEN
         FNAM ='SO2RG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
C
C                                 HIER: HELIZITAET
C                                 ----------------
C
      IF(ISELEP(1,200) .EQ. 0) THEN
         FNAM ='AHE    '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,201) .EQ. 0) THEN
         FNAM ='SHE    '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,202) .EQ. 0) THEN
         FNAM ='HEFG   '
         WRITE (6,6010) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMORY
      ENDIF
      IF(ISELEP(1,203) .EQ. 0) THEN
         FNAM ='AHERG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
      IF(ISELEP(1,204) .EQ. 0) THEN
         FNAM ='SHERG  '
         WRITE (6,6020) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMA
      ENDIF
C
C                                  HIER: KREUZKORRELATIONSFUNKTIONEN
C                                  ---------------------------------
C
      IF(ISELEP(1,210) .EQ. 0) THEN
         FNAM ='ARXUV  '
         WRITE (6,6042) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXA
      ENDIF
      IF(ISELEP(1,211) .EQ. 0) THEN
         FNAM ='SRXUV  '
         WRITE (6,6041) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXS
      ENDIF
      IF(ISELEP(1,212) .EQ. 0) THEN
         FNAM ='ARYUV  '
         WRITE (6,6044) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYA
      ENDIF
      IF(ISELEP(1,213) .EQ. 0) THEN
         FNAM ='SRYUV  '
         WRITE (6,6043) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYS
      ENDIF
      IF(ISELEP(1,216) .EQ. 0) THEN
         FNAM ='ARXUW  '
         WRITE (6,6042) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXA
      ENDIF
      IF(ISELEP(1,217) .EQ. 0) THEN
         FNAM ='SRXUW  '
         WRITE (6,6041) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXS
      ENDIF
      IF(ISELEP(1,218) .EQ. 0) THEN
         FNAM ='ARYUW  '
         WRITE (6,6044) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYA
      ENDIF
      IF(ISELEP(1,219) .EQ. 0) THEN
         FNAM ='SRYUW  '
         WRITE (6,6043) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYS
      ENDIF
      IF(ISELEP(1,224) .EQ. 0) THEN
         FNAM ='ARYVW  '
         WRITE (6,6044) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYA
      ENDIF
      IF(ISELEP(1,225) .EQ. 0) THEN
         FNAM ='SRYVW  '
         WRITE (6,6043) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYS
      ENDIF
C
C                                  HIER: KORRELATIONSKOEFFIZIENTEN
C                                  -------------------------------
C
      IF(ISELEP(1,234) .EQ. 0) THEN
         FNAM ='ACXUW  '
         WRITE (6,6042) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXA
      ENDIF
      IF(ISELEP(1,235) .EQ. 0) THEN
         FNAM ='SCXUW  '
         WRITE (6,6041) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXS
      ENDIF
      IF(ISELEP(1,238) .EQ. 0) THEN
         FNAM ='ACZUW  '
         WRITE (6,6046) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZA
      ENDIF
      IF(ISELEP(1,239) .EQ. 0) THEN
         FNAM ='SCZUW  '
         WRITE (6,6045) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZS
      ENDIF
C
C                                  OBERE DREIECKSMATRIX DES TENSORS
C                                  DER KORRELATIONSFUNKTIONEN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  --------------------------------
C
      IF(ISELEP(1,246) .EQ. 0) THEN
         FNAM ='ARXOXX '
         WRITE (6,6042) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXA
      ENDIF
      IF(ISELEP(1,247) .EQ. 0) THEN
         FNAM ='SRXOXX '
         WRITE (6,6041) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXS
      ENDIF
      IF(ISELEP(1,248) .EQ. 0) THEN
         FNAM ='ARXOXY '
         WRITE (6,6042) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXA
      ENDIF
      IF(ISELEP(1,249) .EQ. 0) THEN
         FNAM ='SRXOXY '
         WRITE (6,6041) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXS
      ENDIF
      IF(ISELEP(1,250) .EQ. 0) THEN
         FNAM ='ARXOXZ '
         WRITE (6,6042) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXA
      ENDIF
      IF(ISELEP(1,251) .EQ. 0) THEN
         FNAM ='SRXOXZ '
         WRITE (6,6041) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXS
      ENDIF
      IF(ISELEP(1,252) .EQ. 0) THEN
         FNAM ='ARXOYY '
         WRITE (6,6042) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXA
      ENDIF
      IF(ISELEP(1,253) .EQ. 0) THEN
         FNAM ='SRXOYY '
         WRITE (6,6041) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXS
      ENDIF
      IF(ISELEP(1,254) .EQ. 0) THEN
         FNAM ='ARXOYZ '
         WRITE (6,6042) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXA
      ENDIF
      IF(ISELEP(1,255) .EQ. 0) THEN
         FNAM ='SRXOYZ '
         WRITE (6,6041) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXS
      ENDIF
      IF(ISELEP(1,256) .EQ. 0) THEN
         FNAM ='ARXOZZ '
         WRITE (6,6042) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXA
      ENDIF
      IF(ISELEP(1,257) .EQ. 0) THEN
         FNAM ='SRXOZZ '
         WRITE (6,6041) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKXS
      ENDIF
      IF(ISELEP(1,258) .EQ. 0) THEN
         FNAM ='ARYOXX '
         WRITE (6,6044) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYA
      ENDIF
      IF(ISELEP(1,259) .EQ. 0) THEN
         FNAM ='SRYOXX '
         WRITE (6,6043) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYS
      ENDIF
      IF(ISELEP(1,260) .EQ. 0) THEN
         FNAM ='ARYOXY '
         WRITE (6,6044) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYA
      ENDIF
      IF(ISELEP(1,261) .EQ. 0) THEN
         FNAM ='SRYOXY '
         WRITE (6,6043) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYS
      ENDIF
      IF(ISELEP(1,262) .EQ. 0) THEN
         FNAM ='ARYOXZ '
         WRITE (6,6044) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYA
      ENDIF
      IF(ISELEP(1,263) .EQ. 0) THEN
         FNAM ='SRYOXZ '
         WRITE (6,6043) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYS
      ENDIF
      IF(ISELEP(1,264) .EQ. 0) THEN
         FNAM ='ARYOYY '
         WRITE (6,6044) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYA
      ENDIF
      IF(ISELEP(1,265) .EQ. 0) THEN
         FNAM ='SRYOYY '
         WRITE (6,6043) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYS
      ENDIF
      IF(ISELEP(1,266) .EQ. 0) THEN
         FNAM ='ARYOYZ '
         WRITE (6,6044) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYA
      ENDIF
      IF(ISELEP(1,267) .EQ. 0) THEN
         FNAM ='SRYOYZ '
         WRITE (6,6043) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYS
      ENDIF
      IF(ISELEP(1,268) .EQ. 0) THEN
         FNAM ='ARYOZZ '
         WRITE (6,6044) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYA
      ENDIF
      IF(ISELEP(1,269) .EQ. 0) THEN
         FNAM ='SRYOZZ '
         WRITE (6,6043) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKYS
      ENDIF
      IF(ISELEP(1,270) .EQ. 0) THEN
         FNAM ='ARZOXX '
         WRITE (6,6046) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZA
      ENDIF
      IF(ISELEP(1,271) .EQ. 0) THEN
         FNAM ='SRZOXX '
         WRITE (6,6045) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZS
      ENDIF
      IF(ISELEP(1,272) .EQ. 0) THEN
         FNAM ='ARZOXY '
         WRITE (6,6046) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZA
      ENDIF
      IF(ISELEP(1,273) .EQ. 0) THEN
         FNAM ='SRZOXY '
         WRITE (6,6045) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZS
      ENDIF
      IF(ISELEP(1,274) .EQ. 0) THEN
         FNAM ='ARZOXZ '
         WRITE (6,6046) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZA
      ENDIF
      IF(ISELEP(1,275) .EQ. 0) THEN
         FNAM ='SRZOXZ '
         WRITE (6,6045) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZS
      ENDIF
      IF(ISELEP(1,276) .EQ. 0) THEN
         FNAM ='ARZOYY '
         WRITE (6,6046) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZA
      ENDIF
      IF(ISELEP(1,277) .EQ. 0) THEN
         FNAM ='SRZOYY '
         WRITE (6,6045) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZS
      ENDIF
      IF(ISELEP(1,278) .EQ. 0) THEN
         FNAM ='ARZOYZ '
         WRITE (6,6046) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZA
      ENDIF
      IF(ISELEP(1,279) .EQ. 0) THEN
         FNAM ='SRZOYZ '
         WRITE (6,6045) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZS
      ENDIF
      IF(ISELEP(1,280) .EQ. 0) THEN
         FNAM ='ARZOZZ '
         WRITE (6,6046) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZA
      ENDIF
      IF(ISELEP(1,281) .EQ. 0) THEN
         FNAM ='SRZOZZ '
         WRITE (6,6045) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMKZS
      ENDIF
C
C                                  LEISTUNGSDICHTESPEKTREN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  ----------------------------
C
      IF(ISELEP(1,282) .EQ. 0) THEN
         FNAM ='ASXOXX '
         WRITE (6,6050) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSX
      ENDIF
      IF(ISELEP(1,283) .EQ. 0) THEN
         FNAM ='SSXOXX '
         WRITE (6,6050) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSX
      ENDIF
      IF(ISELEP(1,284) .EQ. 0) THEN
         FNAM ='ASXOYY '
         WRITE (6,6050) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSX
      ENDIF
      IF(ISELEP(1,285) .EQ. 0) THEN
         FNAM ='SSXOYY '
         WRITE (6,6050) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSX
      ENDIF
      IF(ISELEP(1,286) .EQ. 0) THEN
         FNAM ='ASXOZZ '
         WRITE (6,6050) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSX
      ENDIF
      IF(ISELEP(1,287) .EQ. 0) THEN
         FNAM ='SSXOZZ '
         WRITE (6,6050) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSX
      ENDIF
      IF(ISELEP(1,288) .EQ. 0) THEN
         FNAM ='ASYOXX '
         WRITE (6,6052) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSY
      ENDIF
      IF(ISELEP(1,289) .EQ. 0) THEN
         FNAM ='SSYOXX '
         WRITE (6,6052) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSY
      ENDIF
      IF(ISELEP(1,290) .EQ. 0) THEN
         FNAM ='ASYOYY '
         WRITE (6,6052) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSY
      ENDIF
      IF(ISELEP(1,291) .EQ. 0) THEN
         FNAM ='SSYOYY '
         WRITE (6,6052) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSY
      ENDIF
      IF(ISELEP(1,292) .EQ. 0) THEN
         FNAM ='ASYOZZ '
         WRITE (6,6052) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSY
      ENDIF
      IF(ISELEP(1,293) .EQ. 0) THEN
         FNAM ='SSYOZZ '
         WRITE (6,6052) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMSY
      ENDIF
C
C                                  HAEUFIGKEITSVERTEILUNG DER
C                                  INKLINATIONSWINKEL DER VORTICITY-
C                                  VEKTOREN
C                                  ---------------------------------
C
      IF(ISELEP(1,306) .EQ. 0) THEN
         FNAM ='AHOZOY '
         WRITE (6,6060) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMH
      ENDIF
      IF(ISELEP(1,307) .EQ. 0) THEN
         FNAM ='SHOZOY '
         WRITE (6,6060) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMH
      ENDIF
      IF(ISELEP(1,308) .EQ. 0) THEN
         FNAM ='AHOZOX '
         WRITE (6,6060) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMH
      ENDIF
      IF(ISELEP(1,309) .EQ. 0) THEN
         FNAM ='SHOZOX '
         WRITE (6,6060) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMH
      ENDIF
      IF(ISELEP(1,310) .EQ. 0) THEN
         FNAM ='AHOYOX '
         WRITE (6,6060) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMH
      ENDIF
      IF(ISELEP(1,311) .EQ. 0) THEN
         FNAM ='SHOYOX '
         WRITE (6,6060) FNAM, FNAM
      ELSE
         MEMSUR = MEMSUR + MEMH
      ENDIF
C
C                                 ERMITTLUNG DES SPEICHERPLATZES
C                                 DER FELDER 'HILF', 'HILFL' UND
C                                 'GRHF'.
C                                 MAX0 - ANWEISUNG WEGEN EQUIVALENCE-
C                                 STATEMENT IM HAUPTPROGRAMM
C
      MEMSUR = MEMSUR + MAX0 (MEMORY, MEMML, (IDIMF*(IDIMF+3)))
C
C                                 X,DX,DDX, Y,DY,DDY, Z,DZ,DDZ:
C
      MEMSUR = MEMSUR + 3*(II + JJ + KK)
C
C                                 WCU, WCV, WCW, DIVG:
C
      MEMSUR = MEMSUR + 4*(KK * JJ)
C
C                                 RIDENT/D
C
      MEMSUR = MEMSUR + 200
C
C                                 IIDENT/D
C
      MEMSUI = MEMSUI + 200
C
C                                 CIDENT/D
C
      MEMSUC = MEMSUC +  20
C
C                                 ISELEP, ISELEA:
C
      MEMSUI = MEMSUI + 3008
C
C                                 UFR, GI1 UND GI2:
C
      MEMSUR = MEMSUR + 4*(KK * JJ)
C
C
C                                 IGRHF
C
      MEMSUI = MEMSUI + IDIMF*(IDIMF+1)
C
C                                 ISLINP, ISLINA, ISLINI
C
      MEMSUI = MEMSUI + ISLIDI + ISLIDI + ISLIDI
C
C                                 FVT-FELD,  HVX-FELD
C
      MEMSUR = MEMSUR + (IDIMF+1) + (KLASS * JJ2L)
C
C                                 REAL-BEDARF FUER DEN GITTERGENERATOR
C                                 GRDFMI
C
      MEMSUR = MEMSUR + 2*IDIMF + NBSP*18
C
C                                 GESAMTBEDARF FUER DAS PARTIKEL-
C                                 PROGRAMM VON G. EDER
C
      MEMSPC = 0
      MEMSPI = NMAX + 2*NRZUL
      MEMSPR = NMAX*NDGL + 7*NDGL
C
      WRITE (6,6030) MEMSUC,        MEMSUI,        MEMSUR,
     $               MEMSUC+MEMSPC, MEMSUI+MEMSPI, MEMSUR+MEMSPR
C
      RETURN
 6000 FORMAT (/,4X,10(1H*),'  INFORMATION UEBER DIE FELD',
     $        'DIMENSIONIERUNG (SUBR. SELECI)  ',10(1H*),/)
 6001 FORMAT (4X,'INHALT DES ISELEP(1,...) - FELDES : ',/)
 6005 FORMAT (4X,8('(1,',I3,') = ',I1,4X))
 6006 FORMAT (/,4X,'KK    = ',I4,4X,'JJ    = ',I4,4X,'II    = ',I4,/,
     $           4X,'KKA   = ',I4,4X,'JJA   = ',I4,4X,'IIA   = ',I4,/,
     $           4X,'KKXL  = ',I4,4X,'KKYL  = ',I4,4X,'KKZL  = ',I4,/,
     $           4X,'KSXL  = ',I4,4X,'KSYL  = ',I4,4X,'KSZL  = ',I4,/,
     $           4X,'KH0L  = ',I4,4X,'JJ1L  = ',I4,4X,'JJ2L  = ',I4,/,
     $           4X,'KLASS = ',I4,4X,'ILIMXP= ',I4,4X,'NMAX  = ',I6,/)
 6010 FORMAT (4X,'STATT ',A7,'(KK   , JJ  , II    ) GENUEGT FOLGENDE',
     $        ' DIMENSIONIERUNG : ',A7,'( 1, 1, 1)')
 6020 FORMAT (4X,'STATT ',A7,'(KKA  , JJA , IIA   ) GENUEGT FOLGENDE',
     $        ' DIMENSIONIERUNG : ',A7,'( 1, 1, 1)')
 6030 FORMAT (/,4X,'BENOETIGTER SPEICHERPLATZ (IN WORTEN)',
     $        ' BEI OPTIMALER FELDDIMENSIONIERUNG FUER:',/,
     $        4X,'CHARACTER-FELDER',4X,'INTEGER-FELDER',4X,
     $        'REAL-FELDER',/,6X,I12,7X,I12,5X,I12,3X,
     $        '(OHNE PARTIKELPROGRAMM)',/,
     $        6X,I12,7X,I12,5X,I12,3X,'(MIT  PARTIKELPROGRAMM)',/)
 6041 FORMAT (4X,'STATT ',A7,'(KKXL , JJ1L, ILIMXP) GENUEGT FOLGENDE',
     $        ' DIMENSIONIERUNG : ',A7,'( 1, 1, 1)')
 6042 FORMAT (4X,'STATT ',A7,'(KKXL , JJ2L, ILIMXP) GENUEGT FOLGENDE',
     $        ' DIMENSIONIERUNG : ',A7,'( 1, 1, 1)')
 6043 FORMAT (4X,'STATT ',A7,'(KKYL , JJ1L, ILIMXP) GENUEGT FOLGENDE',
     $        ' DIMENSIONIERUNG : ',A7,'( 1, 1, 1)')
 6044 FORMAT (4X,'STATT ',A7,'(KKYL , JJ2L, ILIMXP) GENUEGT FOLGENDE',
     $        ' DIMENSIONIERUNG : ',A7,'( 1, 1, 1)')
 6045 FORMAT (4X,'STATT ',A7,'(KKZL , JJ1L, ILIMXP) GENUEGT FOLGENDE',
     $        ' DIMENSIONIERUNG : ',A7,'( 1, 1, 1)')
 6046 FORMAT (4X,'STATT ',A7,'(KKZL , JJ2L, ILIMXP) GENUEGT FOLGENDE',
     $        ' DIMENSIONIERUNG : ',A7,'( 1, 1, 1)')
 6050 FORMAT (4X,'STATT ',A7,'(KSXL , JJ1L, ILIMXP) GENUEGT FOLGENDE',
     $        ' DIMENSIONIERUNG : ',A7,'( 1, 1, 1)')
 6052 FORMAT (4X,'STATT ',A7,'(KSYL , JJ1L, ILIMXP) GENUEGT FOLGENDE',
     $        ' DIMENSIONIERUNG : ',A7,'( 1, 1, 1)')
 6054 FORMAT (4X,'STATT ',A7,'(KSZL , JJ1L, ILIMXP) GENUEGT FOLGENDE',
     $        ' DIMENSIONIERUNG : ',A7,'( 1, 1, 1)')
 6060 FORMAT (4X,'STATT ',A7,'(KH0L , JJ2L, ILIMXP) GENUEGT FOLGENDE',
     $        ' DIMENSIONIERUNG : ',A7,'( 1, 1, 1)')
      END
