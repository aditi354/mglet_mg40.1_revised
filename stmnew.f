










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
      SUBROUTINE STMNEW  (KKA,JJA,IIA,KMXA,JMXA,IMXA,NPROLD,
     $                    NPRRUN,DREAD,DCONT,NPRNEU,FPRNEU,APHI,PHIS)
C*STARLET***************************************************************
C        S T M N E W      IN STMNEW WERDEN DIE ALTEN ENSEMBLE-
C                         MITTELWERTE VERBESSERT. DIE GEWICHTUNG DER
C                         ALTEN (WAEHREND DER VORANGEGANGENEN LAEUFE
C                         ERZEUGT) UND DER NEUEN ENSEMBLE-MITTELWERTE
C                         (WAEHREND DES MOMENTANEN LAUFES ERZEUGT)
C                         ERGIBT SICH AUS DEM VERHAELTNIS DER ANZAHL
C                         DER STICHPROBEN, DIE ZU DEN JEWEILIGEN
C                         MITTELWERTEN FUEHRTEN.
C*STARLET***************************************************************
C
C PARAM: KKA,JJA,IIA    - ARRAYDIMENSIONEN DER AUSWERTEFELDER
C                         BEISPIEL: FALLS J-RICHTUNG NICHT HOMOGEN IST,
C                         GILT JJA = JJ. ANDERNFALLS IST JJA = 1.
C        KMXA,JMXA,IMXA - GRENZEN DER AUSWERTEFELDER (MIT BOUND)
C        NPROLD        ( )ANZAHL DER STATIST. STICHPROBEN, DIE
C                         WAEHREND ALLER VORANGEGANGENEN LAEUFE
C                         GENOMMEN WURDEN
C                         NPROLD WIRD UNTER UMSTAENDEN VERAENDERT !!
C                         (SIEHE  DCONT)
C        NPRRUN         - ANZAHL DER STATIST. STICHPROBEN, DIE
C                         WAEHREND DES MOMENTANEN LAUFES GE-
C                         NOMMEN WURDEN
C        DREAD          - DREAD = .T. : DIE WAEHREND DES VORANGEGANGENEN
C                         LAUFES ERZEUGTEN ENSEMBLE-MITTELWERTE WERDEN
C                         EINGELESEN.
C                         DREAD = .F. : ES FINDET EINE VORBELEGUNG MIT
C                         "BESTMOEGLICHEN" ENSEMBLE-MITTELWERTEN STATT.
C        DCONT          - DCONT = .T. .AND. DREAD = .T. : DIE GEWICH-
C                         TUNG DER "ALTEN" ENSEMBLE-MITTELWERTE ERFOLGT
C                         ENTSPRECHEND DER ANZAHL DER "ALTEN" STICH-
C                         PROBEN. NPROLD BLEIBT U N V E R A E N D E R T
C
C                         DCONT = .F. .AND. DREAD = .T. : DIE "ALTEN"
C                         ENSEMBLE-MITTELWERTE WERDEN MIT NPRNEU
C                         STICHPROBEN GEWICHTET. NPROLD WIRD DAHER
C                         MIT NPRNEU BELEGT UND  V E R A E N D E R T
C                         IN DAS RUFENDE UNTERPROG. STMIMP ZURUECK-
C                         GEGEBEN.
C        NPRNEU         - SIEHE BESCHREIBUNG  "DCONT"
C        APHI(KKA,JJA,IIA)  + ENSEMBLE-MITTELWERT DER ALLGEMEINEN
C                             VARIABLEN PHI
C        PHIS(KKA,JJA,IIA)  - FELD, DAS DIE SUMME ALLER STICHPROBEN,
C                             DIE WAEHREND DES MOMENTANEN LAUFES ER-
C                             ZEUGT WURDEN, ENTHAELT.
C
C UPROG                 : KEINE
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        29.08.86 (HW)  : ORIGINAL
C        25.03.98 (MM)  : Falls DCONT=.T. werden alte Stichproben
C                         mit FPRNEU*NPROLD gewichtet
C
C*STARLET***************************************************************
C
      REAL     APHI(KKA,JJA,IIA),      PHIS(KKA,JJA,IIA)
C
      LOGICAL  DREAD,  DCONT
C
      IF(DREAD) THEN
         IF( .NOT. DCONT) THEN
            NPROLD = NPRNEU
         ELSE
            NPROLD = INT( FLOAT(NPROLD) * FPRNEU )
         ENDIF
      ENDIF
C
C                                 WICHTUNGSFAKTOREN
C
      WFOLD  = FLOAT(NPROLD) / FLOAT(NPROLD + NPRRUN)
      RPRTOT = 1.0           / FLOAT(NPROLD + NPRRUN)
C
      DO 100 I = 1,IMXA
         DO 110 J = 1,JMXA
            DO 120 K = 1,KMXA
  120          APHI(K,J,I) = APHI(K,J,I)*WFOLD + PHIS(K,J,I)*RPRTOT
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
