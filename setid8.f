










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
      SUBROUTINE SETID8  (CIDENT,IIDENT,RIDENT,RNAME,VERS,MTSTEP,
     $                    NRRUN,MPCORR,KB,JB1,JB2,IB1,IB2,RE,DT,EPCORR,
     $                    ZTOT,YTOT,XTOT,BETA,VREF,EXPON,TU,RHO,GMOL,
     $                    DREAD,ZCUB
     $		         )
C*STAR******************************************************************
C*STAR*  S E T I D        ABSPEICHERN WICHTIGER STEUERPARAMETER UND
C*STAR*                   KENNGROESSEN
C*STAR******************************************************************
C CIDENT(10) ENTHAELT KENNUNG
C
C             VARIABLE          WIRD BELEGT IN ROUTINE       ANMERKUNGEN
C
C                                K-EPS           LES
C
C  CIDENT( 1) = "FB/HW"          SETID          SETID
C  CIDENT( 2) = PROGRAMMNAME       "              "
C  CIDENT( 4) = DATE   ()          "              "
C  CIDENT( 5) = CLOCK  ()          "              "
C  CIDENT( 6) = DATEC  ()          ?              ?
C  CIDENT( 7) = CLOCKC ()          ?              ?
C
C  IIDENT(100) ENTHAELT INTEGER-KONSTANTEN
C
C  IIDENT( 1) = MTSTEP           SETID          SETID
C  IIDENT( 2) = MPCORR             "              "
C  IIDENT( 3) = NRRUN              "              "
C  IIDENT( 4) = NRRUNC             ?              ?
C  IIDENT( 5) = VERSIONSNUMMER     ?            SETID      FUER EINLESEN
C  IIDENT( 6) = VERSIONSNUMMER     ?            SETID      FUER SCHREIBEN
C  IIDENT( 7) = KK                             SETREF
C  IIDENT( 8) = JJ                               ""
C  IIDENT( 9) = II                               ""
C  IIDENT(10) = KB               SETID          SETID
C  IIDENT(11) = JB1                "              "
C  IIDENT(12) = JB2                "              "
C  IIDENT(13) = IB1                "              "
C  IIDENT(14) = IB2                "              "
C  IIDENT(15) = KKA              -----          SETREF       DIM. D.
C                                                            AUSWERTEF.
C  IIDENT(16) = JJA              -----          SETREF             ""
C  IIDENT(17) = IIA              -----          SETREF             ""
C  IIDENT(20) = KPP              SETCPP         SETREF
C  IIDENT(21) = JPP              SETCPP         SETREF
C  IIDENT(22) = IPP              SETCPP         SETREF
C  IIDENT(23) = KMXA             -----            ""
C  IIDENT(24) = JMXA             -----            ""
C  IIDENT(25) = IMXA             -----            ""
C  IIDENT(30) = ITTOT            -----          HAUPTPROG.
C  IIDENT(31) = NPRTOT           -----          HAUPTPROG.
C  IIDENT(75) = 1                -----          HAUPTPROG.  ISELEP-FELD
C                                                           ERWEITERT,
C                                                           ISLINP-FELD
C                                                           NEU EINGEF.
C  IIDENT(50) = SCALAR FIELD 
C               IDENTIFIER       SETID          DIC/DIB
C
C  RIDENT(100) ENTHAELT REAL-KONSTANTEN
C
C  RIDENT( 1) = RE               SETID          SETREF
C  RIDENT( 2) = DT                 ""            SETID
C  RIDENT( 3) = EPCORR             ""              ""
C  RIDENT(10) = ZTOT               "              "
C  RIDENT(11) = YTOT               "              "
C  RIDENT(12) = XTOT               "              "
C  RIDENT(13) = BETA               "              "
C  RIDENT(14) = VREF               "              "
C  RIDENT(15) = EXPON              "              "
C  RIDENT(16) = TU                 ""           ANDERE BED.
C  RIDENT(16) = GRADPX           ANDERE BED.    SETID
C  RIDENT(17) = RHO              SETID          SETID
C  RIDENT(18) = GMOL               "              "
C  RIDENT(20) = RHO              SETCPP         SETREF
C  RIDENT(21) = VREF             SETCPP         ANDERE BED.
C  RIDENT(21) = UREF             ANDERE BED.    SETREF
C  RIDENT(22) = CPP              SETCPP         -----
C  RIDENT(23) = ALREF            -----          SETREF
C  RIDENT(24) = TRF              -----          SETREF
C  RIDENT(25) = ZCUB             -----          SETREF
C  RIDENT(30) = PHYS.
C               GESAMTZEIT       -----          HAUPTPROG.
C  RIDENT(31) = TCYCLE           -----          VORBELEGUNG IN DEOBI
C                                               ENDGUELTIG IN PROG.
C                                               "FILOPI"
C  RIDENT(40) = ?????            SET??          -----
C  RIDENT(41) = ?????            SET??          -----
C  RIDENT(42) = ?????            SET??          -----
C
C  RIDENT(50) = LAMDA            SETID          SETID
C  RIDENT(51) = PRMOL              "              "
C  RIDENT(52) = PRTURB             "              "
C  RIDENT(53) = EXPONT             "              "
C  RIDENT(53) = TREF               "              "
C
C
C  RIDENT(75) = UREF             -----          SETREF       GESCHW.
C  RIDENT(76) = ALREF            -----          SETREF       LAENGE
C  RIDENT(77) = TRF              -----          SETREF       ZEIT
C  RIDENT(78) = EREF             -----          SETREF       ENERGIE
C  RIDENT(79) = GREF             -----          SETREF       VISKOSIT.
C  RIDENT(80) = PREF             -----          SETREF       DRUCK
C  RIDENT(81) = TAURE            -----          SETREF       SCHUBSPG.
C
C DEFINE-DIREKTIVEN     : KEINE
C
C VERS:  07.12.84 (HW)  : ORIGINAL
C        29.03.85 (HW)  : DATE UND CLOCK DIREKT IN CIDENT-FELD
C                         GESCHRIEBEN
C        14.11.85 (HW)  : IIDENT UND RIDENT MIT NULL VORBELEGT
C        09.07.86 (HW)  : ERLAEUTERUNGEN AUF NEUESTEN STAND GEBRACHT
C        30.09.86 (HW)  : SETIDL AUS SETID FUER LES ABGELEITET.
C                         DREAD WIRD ZUSAETZLICH UEBERGEBEN, SO DASS IM
C                         RIDENT-FELD REFERENZGROESSEN DES VORAN-
C                         GEGANGENEN LAUFES GERETTET WERDEN KOENNEN,
C                         SOFERN SIE IM IIDENT- UND RIDENT-FELD AUF DEN
C                         SPEICHERPLAETZEN 75 ... 100 STEHEN !!!
C        13.12.88 (HW)  : IIDENT (75) WIRD MIT 1 BELEGT
C        21.02.89 (HW)  : GETAET EINGEFUEHRT (ACTUAL ELAPSED TIME)
C        27.02.89 (HW)  : INTEGER GETAET, GETDAT FUER CRAY !
C         6. 4.92 (MM)  : VERSIONSNUMMER IIDENT(5) EINGEFUEHRT
C        10. 4.92 (MM)  : ES WERDEN, FALLS DREAD NUR NOCH BESTIMMTE
C                         GROESSEN BELEGT
C        29.01.03 (TB)  : _KSR_ REMOVED
C        29.01.03 (TB)  : SCALAR PARAMETERS LAMDA, PRMOL,PRTURB,EXPONT 
C                         TO RIDENT(50-54)
C
C*STAR******************************************************************
C
      CHARACTER (LEN=8) RNAME, VERS, CIDENT(10)
      CHARACTER (LEN=8) GETAET, GETDAT
      INTEGER      IIDENT(100)
      REAL         RIDENT(100)
      LOGICAL      DREAD
C
      IF(.NOT.DREAD) THEN
C
C                                     VORBELEGUNG
C
      DO 90  K = 1,10
   90    CIDENT(K) = '        '
C
      DO 100 K=1,74
         IIDENT(K) = 0
  100    RIDENT(K) = 0.0
C
         DO 110 K = 75,100
            IIDENT(K) = 0
  110       RIDENT(K) = 0.0
C
C
      IIDENT(5)  = 40
      IIDENT(10) = KB
      IIDENT(11) = JB1
      IIDENT(12) = JB2
      IIDENT(13) = IB1
      IIDENT(14) = IB2
CTBA1 070203 IIDENT(5) SET TO 1 -> T field will be written to result
C
      RIDENT(10) = ZTOT
      RIDENT(11) = YTOT
      RIDENT(12) = XTOT
      RIDENT(25) = ZCUB
C
      ELSE
C
      ZTOT = RIDENT(10)
      YTOT = RIDENT(11)	
      XTOT = RIDENT(12)	
      ZCUB = RIDENT(25)
      KB = IIDENT(10)
      JB1 = IIDENT(11)
      JB2 = IIDENT(12)
      IB1 = IIDENT(13)
      IB2 = IIDENT(14)
      
C
      ENDIF
C
      CIDENT(1)  = RNAME
      CIDENT(2)  = VERS
      WRITE(CIDENT(4),'(A8)') GETDAT(0)
      WRITE(CIDENT(5),'(A8)') GETAET(0)
C
      IIDENT(1)  = MTSTEP
      IIDENT(2)  = MPCORR
      IIDENT(3)  = NRRUN
C                                           VERSIONSNUMMER
C                                           MUSS VOR RAUSSCHREIBEN AUF
C                                           IIDENT(5) UEBERTRAGEN WERDEN
C                                           (IN SUBR. SETREF)
      IIDENT(6)  = 40
C
C     RIDENT(1)  = RE
      RIDENT(2)  = DT
      RIDENT(3)  = EPCORR
      RIDENT(13) = BETA
      RIDENT(14) = VREF
      RIDENT(15) = EXPON
      RIDENT(16) = TU
      RIDENT(17) = RHO
      RIDENT(18) = GMOL
C
      RETURN
      END
