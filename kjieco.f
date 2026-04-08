










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
      SUBROUTINE KJIECO  (ISLINP,ISLIDI,ISLINI,ILIMXP,ILI,IL0,ILX,ILY,
     $                    ILZ,CD,KAUP,JAUP,IAUP,IB1,IB2,JB1,JB2,KB,
     $                    NBND,IE1,IE2,IE3,IE4,IE5)
C*STARLET***************************************************************
C        K J I E C O      IN KJIECO WIRD DIE INFORMATION UEBER DEN
C                         AUFPUNKT UND EINE EVENTUELL BENOETIGTE KOOR-
C                         DINATENRICHTUNG IN EINE EINZIGE 10-STELLIGE
C                         (3*NSTELL + 1) INTEGER-ZAHL KODIERT.
C*STARLET***************************************************************
C
C PARAM: ISLINP(ISLIDI) + ENTHAELT DIE AUFPUNKTE IN KODIERTER FORM
C        ISLIDI         - ARRAYDIMENSION
C        ISLINI(ISLIDI) + REINES INFORMATIONSFELD UM DIE ZUORDNUNG DER
C                         AUFPUNKTE WIE SIE IN  SELAUF  ANGEGEBEN WERDEN
C                         UND WIE SIE IM ISLINP-FELD ABGELEGT WERDEN
C                         HERZUSTELLEN
C        ILIMXP         - DIE MAXIMAL ZULAESSIGE ANZAHL VON EIN-
C                         DIMENSIONALEN "LINIEN" IST ILIMXP-1
C                         (ILI = 2 ... ILIMXP)
C        ILI            + ZAEHLT DEN IM ISLINP-FELD BENOETIGTEN
C                         SPEICHERPLATZ MIT (ZU KONTROLLZWECKEN !)
C        IL0            + ZAEHLT DIE ANZAHL DER AUFPUNKTE OHNE
C                         EINE DEFINIERTE RICHTUNGSINFORMATION (CD='0')
C        ILX            + ZAEHLT DIE ANZAHL DER AUFPUNKTE MIT
C                         DER RICHTUNGSINFORMATION  X  (CD='X')
C        ILY            + ZAEHLT DIE ANZAHL DER AUFPUNKTE MIT
C                         DER RICHTUNGSINFORMATION  Y  (CD='Y')
C        ILZ            + ZAEHLT DIE ANZAHL DER AUFPUNKTE MIT
C                         DER RICHTUNGSINFORMATION  Z  (CD='Z')
C        CD             - CHARACTER (LEN=1) VARIABLE (ENTHAELT DIE
C                         RICHTUNGSINFORMATION
C        K-, J-, IAUP   - INDIZES DES AUFPUNKTES FUER DIE KORRELATIONEN
C        IB1            - I-INDEX DER ERSTEN  MASCHENZELLE  I M  KUBUS
C        IB2            - I-INDEX DER LETZTEN MASCHENZELLE  I M  KUBUS
C        JB1            - J-INDEX DER ERSTEN  MASCHENZELLE  I M  KUBUS
C        JB2            - J-INDEX DER LETZTEN MASCHENZELLE  I M  KUBUS
C        KB             - K-INDEX DER LETZTEN MASCHENZELLE  I M  KUBUS
C        NBND           - ANZAHL DER RANDSCHICHTEN
C        IE.            + INTEGERZAHL ZUR FEHLERERKENNUNG
C                         IE. =   0 : ALLES IN ORDNUNG
C                         IE1 =   1 : K-, J-, ODER IAUFP IST < 0 !
C                         IE2 =   1 : K-, J-, ODER IAUFP BENOETIGT
C                                     MEHR ALS 'NSTELL' STELLEN
C                                     (NORMALERWEISE 3-STELLIG)
C                         IE3 =   1 : ISLIDI IST ZU KLEIN, UM ALLE
C                                     VOM BENUTZER GEWUENSCHTEN AUF-
C                                     PUNKTE SPEICHERN ZU KOENNEN.
C                         IE4 =   1 : ILIMXP IST ZU KLEIN, UM ALLE
C                                     VOM BENUTZER GEWUENSCHTEN AUF-
C                                     PUNKTE SPEICHERN ZU KOENNEN.
C                         IE5 =   1 : DIE RICHTUNGSINFORMATION IST
C                                     UNZULAESSIG
C
C UPROG                 : ERRR,  KJIDCO
C
C DEFINE-DIREKTIVEN     : XHOMOG, YHOMOG, ZHOMOG
C
C        10.10.88 (HW)  : ORIGINAL
C        24.11.88 (HW)  : VERBESSERUNG DER UEBERPRUEFUNGEN
C        24.02.89 (HW)  : ISLINI(ILIMXP) --> ISLINI(ISLIDI)
C
C*STARLET***************************************************************
C
      CHARACTER (LEN=1)  CD,  CDT,  CDV
      INTEGER        ISLINP (ISLIDI),  ISLINI (ISLIDI)
      SAVE           NCA
      DATA           NCA           /0/
C
C                                 DIE X-RICHTUNG IST  N I C H T  HOMOGEN
      NXRIHO = 0
C
C                                 DIE Y-RICHTUNG IST  N I C H T  HOMOGEN
      NYRIHO = 0
C
C                                 DIE Z-RICHTUNG IST  N I C H T  HOMOGEN
      NZRIHO = 0
C
C                                 NCA ZAEHLT DIE AUFRUFE VON KJIECO MIT
C
      NCA     = NCA + 1
C
      KAUFP   = KAUP
      JAUFP   = JAUP
      IAUFP   = IAUP
C
C                                 ANZAHL DER FUER K-, J-, ODER IAUFP
C                                 ZULAESSIGEN STELLEN
C
      NSTELL = 3
      IL0M   = 0
      ILXM   = 0
      ILYM   = 0
      ILZM   = 0
      IDIRCO = 9 * 10**(3*NSTELL)
C
C                                 MASCHENINDIZES IN KUBUSMITTE
C
      IF(KB .GT. NBND) THEN
         IBMIT  = (IB1    + IB2) / 2
         JBMIT  = (JB1    + JB2) / 2
         KBMIT  = (NBND+1 + KB ) / 2
      END IF
C
C                                 IN DEN HOMOGENEN RICHTUNGEN KOENNEN
C                                 DIE AUFPUNKTE BELIEBIG VERSCHOBEN
C                                 WERDEN. UM UNSINNIGE AUFPUNKTE IN
C                                 SUBR. RKJIAE ENTDECKEN ZU KOENNEN
C                                 WIRD DER AUFPUNKT IN KUBUSMITTE
C                                 VERSCHOBEN.
C
      IF(KB .GT. NBND  .AND.  NXRIHO .EQ. 1) IAUFP = IBMIT
      IF(KB .GT. NBND  .AND.  NYRIHO .EQ. 1) JAUFP = JBMIT
      IF(KB .GT. NBND  .AND.  NZRIHO .EQ. 1) KAUFP = KBMIT
C
C                                 UEBERPRUEFUNGEN
C
      MINKJI = MIN0 (KAUFP, JAUFP, IAUFP)
      MAXKJI = MAX0 (KAUFP, JAUFP, IAUFP)
C
      IF(MINKJI .LT.   0) THEN
         WRITE (6,6010)
         WRITE (6,6020) ILI
         IE1 = 1
      END IF
C
C                                 KODIERUNG DER RICHTUNGSINFORMATION
C
      IF(CD .EQ. '0') THEN
         IDIRCO = 0 * 10**(3*NSTELL)
         IL0M   = 1
         GOTO 2000
      END IF
      IF(CD .EQ. 'X') THEN
         IDIRCO = 1 * 10**(3*NSTELL)
         ILXM   = 1
         GOTO 2000
      END IF
      IF(CD .EQ. 'Y') THEN
         IDIRCO = 2 * 10**(3*NSTELL)
         ILYM   = 1
         GOTO 2000
      END IF
      IF(CD .EQ. 'Z') THEN
         IDIRCO = 3 * 10**(3*NSTELL)
         ILZM   = 1
         GOTO 2000
      END IF
C
C                                 HIER: DIE RICHTUNGSINFORMATION
C                                 IST UNZULAESSIG !
C
         WRITE (6,6010)
         WRITE (6,6025) ILI, CD
         IE5 = 1
C
 2000 KJIZUL = 0
      DO 100 N = 1,NSTELL
  100    KJIZUL = KJIZUL + 9 * 10**(N-1)
C
      IF(MAXKJI .GT. KJIZUL) THEN
         WRITE (6,6010)
         WRITE (6,6030) ILI, NSTELL, NSTELL+1
         IE2 = 1
      END IF
C
C                                 BESTIMMUNG DER INTEGER-KODE-ZAHL
C
      ICODE  = KAUFP * (10**(2*NSTELL))
     $       + JAUFP * (10**(1*NSTELL))
     $       + IAUFP * (10**(0*NSTELL))
     $       + IDIRCO
C
C                                 DIE EBEN CODIERTE INTEGER-ZAHL
C                                 WIRD PROBEWEISE DECODIERT. ES MUSS
C                                 UEBEREINSTIMMUNG MIT DEN
C                                 ....T - VARIABLEN VORHANDEN SEIN
C
      CALL KJIDCO (ICODE, CDT, IDIRCT, KAUFT, JAUFT, IAUFT)
C
      IF(ABS(KAUFT) .NE. KAUFP  .OR.  ABS(JAUFT) .NE. JAUFP .OR.
     $   ABS(IAUFT) .NE. IAUFP  .OR.      CD     .NE. CDT       ) THEN
         WRITE (6,6010)
         WRITE (6,6035) 3*NSTELL + 1
         CALL ERRR (501,' KJIECO   ')
      END IF
C
C                                 UEBERPRUEFUNG, OB DER GERADE BETRACH-
C                                 TETE AUFPUNKT NICHT DAS GLEICHE ERGEB-
C                                 NIS LIEFERN WUERDE WIE EIN BEREITS
C                                 DEFINIERTER AUFPUNKT. IST DIES DER
C                                 FALL, SO WIRD DER GERADE BETRACHTETE
C                                 AUFPUNKT IGNORIERT.
C
      IF(ILI .GE. 2) THEN
         DO 200 IL = 2,ILI
            CALL KJIDCO (ISLINP(IL), CDV, IDIRCV, KAUFV, JAUFV, IAUFV)
            IF(CDV .EQ. CD) THEN
               IF(NXRIHO .EQ. 1  .AND.  NYRIHO .EQ. 1  .AND.
     $            NZRIHO .EQ. 1)                              THEN
                  WRITE (6,6050)           NCA, ISLINI (IL), NCA
                  GOTO 9999
               END IF
               IF(NYRIHO .EQ. 1  .AND.  NZRIHO .EQ. 1)        THEN
                  IF(IAUFV .EQ. IAUFP) THEN
                     WRITE (6,6060) 'Y', 'Z', NCA, ISLINI (IL), NCA
                     GOTO 9999
                  END IF
                  GOTO 199
               END IF
               IF(NXRIHO .EQ. 1  .AND.  NZRIHO .EQ. 1)        THEN
                  IF(JAUFV .EQ. JAUFP) THEN
                     WRITE (6,6060) 'X', 'Z', NCA, ISLINI (IL), NCA
                     GOTO 9999
                  END IF
                  GOTO 199
               END IF
               IF(NXRIHO .EQ. 1  .AND.  NYRIHO .EQ. 1)        THEN
                  IF(KAUFV .EQ. KAUFP) THEN
                     WRITE (6,6060) 'X', 'Y', NCA, ISLINI (IL), NCA
                     GOTO 9999
                  END IF
                  GOTO 199
               END IF
               IF(NZRIHO .EQ. 1)                              THEN
                  IF(IAUFV .EQ. IAUFP  .AND.  JAUFV .EQ. JAUFP) THEN
                     WRITE (6,6070)      'Z', NCA, ISLINI (IL), NCA
                     GOTO 9999
                  END IF
                  GOTO 199
               END IF
               IF(NYRIHO .EQ. 1)                              THEN
                  IF(IAUFV .EQ. IAUFP  .AND.  KAUFV .EQ. KAUFP) THEN
                     WRITE (6,6070)      'Y', NCA, ISLINI (IL), NCA
                     GOTO 9999
                  END IF
                  GOTO 199
               END IF
               IF(NXRIHO .EQ. 1)                              THEN
                  IF(JAUFV .EQ. JAUFP  .AND.  KAUFV .EQ. KAUFP) THEN
                     WRITE (6,6070)      'X', NCA, ISLINI (IL), NCA
                     GOTO 9999
                  END IF
                  GOTO 199
               END IF
               IF(IAUFV .EQ. IAUFP  .AND.  JAUFV .EQ. JAUFP  .AND.
     $            KAUFV .EQ. KAUFP)                           THEN
                  WRITE (6,6040)                ISLINI (IL), NCA
                  GOTO 9999
               END IF
            END IF
  199       CONTINUE
  200    CONTINUE
      END IF
C
      ILI    = ILI + 1
      IL0    = IL0 + IL0M
      ILX    = ILX + ILXM
      ILY    = ILY + ILYM
      ILZ    = ILZ + ILZM
      IL0XYZ = IL0*IL0M + ILX*ILXM + ILY*ILYM + ILZ*ILZM
      IF(ILI .LE. ISLIDI) THEN
         IF(IL0XYZ .LE. ILIMXP) THEN
            ISLINP (ILI) = ICODE
            ISLINI (ILI) = NCA
         ELSE
            IE4 = 1
         END IF
      ELSE
         IE3 = 1
      END IF
C
 9999 RETURN
C
 6010 FORMAT (' ********** FEHLERMELDUNG AUS SUBR. KJIECO',
     $        ' **********')
 6020 FORMAT (' DER ',I4,'. VOM BENUTZER ANGEGEBENE AUFPUNKT',
     $        ' ENTHAELT EINEN K-, J-, ODER I-INDEX',/,
     $        ' DER KLEINER ALS NULL IST --> ABBRUCH')
 6025 FORMAT (' DER ',I4,'. VOM BENUTZER ANGEGEBENE AUFPUNKT',
     $        ' ENTHAELT EINE UNZULAESSIGE',/,
     $        ' RICHTUNGSINFORMATION (CD = ',A1,')  --> ABBRUCH')
 6030 FORMAT (' DER FEHLER TRAT BEIM ',I4,'. VOM BENUTZER',
     $        ' ANGEGEBENEN AUFPUNKT EIN !',/,
     $        ' DER K-, J-, ODER I-AUFPUNKT DARF MAXIMAL ',
     $        I2,'-STELLIG SEIN !',/,
     $        ' FALLS AUF ',I2,' ODER MEHR STELLEN',
     $        ' ERWEITERT WIRD, SO IST ZU BERUECKSICHTIGEN,',/,
     $        ' DASS DIE RESULTIERENDE INTEGER-ZAHL DREIMAL',
     $        ' SOVIELE STELLEN HABEN WIRD !',/,
     $        ' DIESE AENDERUNG BETRIFFT AUCH DIE SUBR.',
     $        ' KJIDCO !!!')
 6035 FORMAT (' BEI DER PROBEWEISEN DEKODIERUNG TRAT EIN FEHLER AUF !',
     $        /,' FOLGENDE PUNKTE SOLLTEN ABGECHECKT WERDEN :',/,
     $        ' - DARF AUF DER VERWENDETEN MASCHINE EINE INTEGER-ZAHL',
     $        ' AUS ',I3,' DEZIMALZIFFERN',/,'   BESTEHEN ?',/,
     $        ' - STIMMT DER WERT DER VARIABLEN  NSTELL  IN SUBR.',
     $        ' KJIECO UND KJIDCO',/,'   UEBEREIN ?')
 6040 FORMAT (' DER ',I4,'. VOM BENUTZER (IN SUBR. SELAUF)',
     $        ' ANGEGEBENEN AUFPUNKT IST',/,
     $        ' IDENTISCH ZUM ',I4,'. AUFPUNKT. LETZTERER',
     $        ' WIRD IGNORIERT !')
 6050 FORMAT (' DA DIE X-, Y- UND Z-RICHTUNG HOMOGEN IST, WUERDE',
     $        ' DER ',I4,'. IN SELAUF ',/,' ANGEGEBENE AUFPUNKT DAS',
     $        ' GLEICHE ERGEBNIS LIEFERN WIE DER ',I4,'.',/,' DER ',
     $        I4,'. AUFPUNKT WIRD DAHER IGNORIERT !')
 6060 FORMAT (' DA DIE ',A1,'- UND ',A1,'-RICHTUNG HOMOGEN IST, WUERDE',
     $        ' DER ',I4,'. IN SELAUF ',/,' ANGEGEBENE AUFPUNKT DAS',
     $        ' GLEICHE ERGEBNIS LIEFERN WIE DER ',I4,'.',/,' DER ',
     $        I4,'. AUFPUNKT WIRD DAHER IGNORIERT !')
 6070 FORMAT (' DA DIE ',A1,'-RICHTUNG HOMOGEN IST, WUERDE',
     $        ' DER ',I4,'. IN SELAUF ',/,' ANGEGEBENE AUFPUNKT DAS',
     $        ' GLEICHE ERGEBNIS LIEFERN WIE DER ',I4,'.',/,' DER ',
     $        I4,'. AUFPUNKT WIRD DAHER IGNORIERT !')
      END
