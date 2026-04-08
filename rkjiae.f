










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
      SUBROUTINE RKJIAE  (KK,JJ,II,KMX,JMX,IMX,X,Y,Z,KAUFP,JAUFP,IAUFP,
     $                    KDIR,JDIR,IDIR,NBND,RKOMXP,
     $                    ILI,KJIBEG,KJIEND,KJILAG)
C*STARLET***************************************************************
C        R K J I A E      IN RKJIAE WERDEN DIE ANFANGS- UND ENDINDIZES
C                         -> KJIBEG -> KJIEND FUER DIE BESTIMMUNG DER
C                         KORRELATION ZWEIER GROESSEN BESTIMMT.
C                         DIE ROUTINE FUNKTIONIERT FUER ALLE DREI
C                         KOORDINATENRICHTUNGEN. DIE INFORMATION, IN
C                         WELCHER DER DREI RICHTUNGEN DIE KORRELATION
C                         BESTIMMT WERDEN SOLL, WIRD MITTELS 'KDIR',
C                         'JDIR' UND 'IDIR' UEBERGEBEN.
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        X,Y,Z          - KOORDINATEN DER ZELLMITTELPUNKTE
C        K-, J-, IAUFP  - INDIZES DES AUFPUNKTES FUER DIE KORRELATIONEN
C        K-, J-, IDIR   - NDIR = 1: IN DIESER KOORDINATENRICHTUNG
C                                   WIRD DIE KORRELATION GEBILDET.
C                         NDIR = 0: IN DIESER KOORDINATENRICHTUNG
C                                   WIRD DIE KORRELATION  N I C H T
C                                   GEBILDET. BEACHTE: DIE SUMME
C                                   KDIR+JDIR+IDIR MUSS STETS = 1 SEIN !
C        NBND           - ANZAHL DER RANDSCHICHTEN DES BERECHNUNGS-
C                         GEBIETES
C        RKOMXP         - MAXIMAL ZULAESSIGE KORRELATIONSLAENGE
C        ILI            - KENNZEICHNET DIE GERADE ABGEARBEITETE I-LINIE
C        KJIBEG         + INDEX DES STARTPUNKTES FUER DIE KORRELATION
C        KJIEND         + INDEX DES ENDPUNKTES FUER DIE KORRELATION
C        KJILAG         + GIBT BEI HOMOGENER RICHTUNG  U N D  AEQUI-
C                         DISTANTEM GITTER DEN MAXIMALEN LAG (ANZAHL
C                         DER GITTERPUNKTE ZWISCHEN AUFPUNKT UND DEM
C                         AM WEITESTEN ENTFERNTEN KORRELATIONSPUNKT )
C                         AN. INTERESSANT IST DIESE GROESSE, FALLS
C                         KJILAG < (NMX/2) ; N = K, J, I.
C
C UPROG                 : ERRR
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        11.09.88 (HW)  : ORIGINAL
C
C*STARLET***************************************************************
C
C
      REAL     X(II),  Y(JJ),  Z(KK)
C
C                                 UEBERPRUEFUNG, OB DER AUFPUNKT
C                                 AUSSERHALB DES STROEMUNGSFELDES LIEGT
C
      NBP1   = NBND + 1
      NBM1   = NBND - 1
C
      IF(IAUFP .LE. NBND  .OR.  IAUFP .GE. (IMX-NBM1)) THEN
         WRITE (6,6010)
         WRITE (6,6020) 'I', ILI-1
         CALL ERRR (501,' RKJIAE   ')
      END IF
      IF(JAUFP .LE. NBND  .OR.  JAUFP .GE. (JMX-NBM1)) THEN
         WRITE (6,6010)
         WRITE (6,6020) 'J', ILI-1
         CALL ERRR (502,' RKJIAE   ')
      END IF
      IF(KAUFP .LE. NBND  .OR. KAUFP .GE. (KMX-NBM1)) THEN
         WRITE (6,6010)
         WRITE (6,6020) 'K', ILI-1
         CALL ERRR (503,' RKJIAE   ')
      END IF
C
C                                 UEBERPRUEFUNG, OB DIE KORRELATION
C                                 TATSAECHLICH NUR IN EINER KOORD.-RI.
C                                 GEBILDET WIRD.
C
      IF(IDIR .NE. 0  .AND.  IDIR .NE. 1) CALL ERRR (511,' RKJIAE   ')
      IF(JDIR .NE. 0  .AND.  JDIR .NE. 1) CALL ERRR (512,' RKJIAE   ')
      IF(KDIR .NE. 0  .AND.  KDIR .NE. 1) CALL ERRR (513,' RKJIAE   ')
      IF((KDIR + JDIR + IDIR) .NE. 1)     CALL ERRR (514,' RKJIAE   ')
C
      KJIAUF = KAUFP*KDIR + JAUFP*JDIR + IAUFP*IDIR
      KJIM   = (KMX-NBM1)*KDIR + (JMX-NBM1)*JDIR + (IMX-NBM1)*IDIR
      NEND   = MAX0 (KMX,JMX,IMX)
C
C                                 BESTIMMUNG VON KJIBEG
C
      DO 100 N = 1,NEND
         KJIBEG = KJIAUF - N
         JKS    = JAUFP*(1-JDIR) + KJIBEG*JDIR
         JKS    = MIN0 (MAX0 (JKS,NBP1), JMX-NBND)
         IKS    = IAUFP*(1-IDIR) + KJIBEG*IDIR
         IKS    = MIN0 (MAX0 (IKS,NBP1), IMX-NBND)
         KSTBEG = 3
         KSTAUF = KAUFP * MAX0 (IDIR,JDIR)  +  KJIBEG*KDIR
C
C                                 UEBERPRUEFUNG DER MAXIMALEN
C                                 KORRELATIONSLAENGE
C
         IXA    = MIN0 (MAX0 (KJIAUF,NBND), IMX-NBM1)
         IXB    = MIN0 (MAX0 (KJIBEG,NBND), IMX-NBM1)
         IYA    = MIN0 (MAX0 (KJIAUF,NBND), JMX-NBM1)
         IYB    = MIN0 (MAX0 (KJIBEG,NBND), JMX-NBM1)
         IZA    = MIN0 (MAX0 (KJIAUF,NBND), KMX-NBM1)
         IZB    = MIN0 (MAX0 (KJIBEG,NBND), KMX-NBM1)
         RKOR   = ABS( (X(IXA) - X(IXB)) * FLOAT(IDIR)
     $          +      (Y(IYA) - Y(IYB)) * FLOAT(JDIR)
     $          +      (Z(IZA) - Z(IZB)) * FLOAT(KDIR) )
         IF((KSTBEG .GT. KSTAUF)  .OR.  (KJIBEG .LE. NBND)  .OR.
     $      (RKOR   .GT. RKOMXP)) THEN
            KJIBEG = KJIBEG + 1
            GOTO 2100
         END IF
  100 CONTINUE
C
C                                 BESTIMMUNG VON KJIEND
C
 2100 DO 110 N = 1,NEND
         KJIEND = KJIAUF + N
         JKS    = JAUFP*(1-JDIR) + KJIEND*JDIR
         JKS    = MIN0 (MAX0 (JKS,NBP1), JMX-NBND)
         IKS    = IAUFP*(1-IDIR) + KJIEND*IDIR
         IKS    = MIN0 (MAX0 (IKS,NBP1), IMX-NBND)
         KSTEND = 3
         KSTAUF = KAUFP * MAX0 (IDIR,JDIR)  +  KJIEND*KDIR
C
C                                 UEBERPRUEFUNG DER MAXIMALEN
C                                 KORRELATIONSLAENGE
C
         IXA    = MIN0 (MAX0 (KJIAUF,NBND), IMX-NBM1)
         IXB    = MIN0 (MAX0 (KJIEND,NBND), IMX-NBM1)
         IYA    = MIN0 (MAX0 (KJIAUF,NBND), JMX-NBM1)
         IYB    = MIN0 (MAX0 (KJIEND,NBND), JMX-NBM1)
         IZA    = MIN0 (MAX0 (KJIAUF,NBND), KMX-NBM1)
         IZB    = MIN0 (MAX0 (KJIEND,NBND), KMX-NBM1)
         RKOR   = ABS( (X(IXA) - X(IXB)) * FLOAT(IDIR)
     $          +      (Y(IYA) - Y(IYB)) * FLOAT(JDIR)
     $          +      (Z(IZA) - Z(IZB)) * FLOAT(KDIR) )
         IF((KSTEND .GT. KSTAUF)  .OR.  (KJIEND .GE. KJIM)  .OR.
     $      (RKOR   .GT. RKOMXP)) THEN
            KJIEND = KJIEND - 1
            GOTO 2200
         END IF
  110 CONTINUE
C
C                                 BESTIMMUNG VON KJILAG
C
 2200 KJILAG = MAX0 (IABS(KJIAUF-KJIBEG), IABS(KJIAUF-KJIEND))
C
      RETURN
 6010 FORMAT (' ********** FEHLERMELDUNG AUS SUBR. RKJIAE',
     $        ' **********')
 6020 FORMAT (' DER ',A1,'-AUFPUNKT DES ',I4,'. VOM BENUTZER',
     $        ' IN SUBR. SELAUF VORGEGEBENEN AUFPUNKTES',/,
     $        ' LIEGT AUSSERHALB DES PHYSIKALISCHEN GEBIETES !')
      END
