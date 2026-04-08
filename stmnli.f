










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
      SUBROUTINE STMNLI  (APHI,JJAL,SPHI,JJSL,KKNL,ILIMXP,
     $                    IBEG,IEND,JBEG,JEND,KBEG,KEND,NPROLD,
     $                    NPRRUS,DREAD,DCONT,NPRNEU)
C*STARLET***************************************************************
C        S T M N L I      IN STMNLI WERDEN DIE ALTEN ENSEMBLE-
C                         MITTELWERTE VON KORRELATIONSFUNKTIONEN,
C                         -KOEFFIZIENTEN, LEISTUNGSDICHTESPEKTREN UND
C                         HAEUFIGKEITSVERTEILUNGEN (ALLGEMEIN: ALLE
C                         "LINIEN"-FELDER) VERBESSERT. DIE GEWICHTUNG
C                         DER ALTEN (WAEHREND DER VORANGEGANGENEN LAEUFE
C                         ERZEUGT) UND DER NEUEN ENSEMBLE-MITTELWERTE
C                         (WAEHREND DES MOMENTANEN LAUFES ERZEUGT)
C                         ERGIBT SICH AUS DEM VERHAELTNIS DER ANZAHL
C                         DER STICHPROBEN, DIE ZU DEN JEWEILIGEN
C                         MITTELWERTEN FUEHRTEN.
C*STARLET***************************************************************
C
C PARAM: APHI(KKNL,     + ENSEMBLE-MITTELWERT DER ALLGEMEINEN
C        JJAL,ILIMXP)     VARIABLEN PHI
C        JJAL           - ARRAYDIM. IN J-RICHTUNG DES A_VERAGE-FELDES
C        SPHI(KKNL,     - FELD, DAS DIE SUMME ALLER STICHPROBEN,
C        JJSL,ILIMXP)     DIE WAEHREND DES MOMENTANEN LAUFES ER-
C                         ZEUGT WURDEN, ENTHAELT.
C        JJSL           - ARRAYDIM. IN J-RICHTUNG DES S_UMMATIONS-FELDES
C        KKNL,ILIMXP    - ARRAYDIMENSIONEN IN K- UND I-RICHTUNG
C        IBEG,IEND      - START- UND ENDINDEX FUER DIE "I" DO-SCHLEIFE
C        JBEG,JEND      - START- UND ENDINDEX FUER DIE "J" DO-SCHLEIFE
C        KBEG,KEND      - START- UND ENDINDEX FUER DIE "K" DO-SCHLEIFE
C        NPROLD        (+)ANZAHL DER STATIST. STICHPROBEN, DIE
C                         WAEHREND ALLER VORANGEGANGENEN LAEUFE
C                         GENOMMEN WURDEN
C                         NPROLD WIRD UNTER UMSTAENDEN VERAENDERT !!
C                         (SIEHE  DCONT)
C        NPRRUS         - ANZAHL DER STATIST. STICHPROBEN, DIE
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
C                         DER ZUSATZ  .AND. (NPROLD .GT. 0)  VERMEIDET,
C                         DASS MIT 0.0 VORBELEGTE WERTE EINGEHEN.
C        NPRNEU         - SIEHE BESCHREIBUNG  "DCONT"
C
C UPROG                 : ERRR
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        14.12.88 (HW)  : ORIGINAL AUS STMNEW ABGELEITET
C
C*STARLET***************************************************************
C
      LOGICAL  DREAD,  DCONT
      REAL     APHI(KKNL,JJAL,ILIMXP),      SPHI(KKNL,JJSL,ILIMXP)
C
C
      DO 100 I = IBEG,IEND
         NPRRUN = IFIX (SPHI (KKNL-13, 1, I))
         IF(NPRRUN .NE. NPRRUS) THEN
C
C                                 IM FEHLERFALL WIRD DEM BENUTZER
C                                 EINE MOEGLICHST UMFANGREICHE INFOR-
C                                 MATION GEGEBEN
C
            WRITE (6,*)
            WRITE (6,*) ' ********** FEHLERMELDUNG AUS SUBR.: STMNLI',
     $                  ' **********'
            WRITE (6,*) ' SOLLWERT VON NPRRUN = ',NPRRUS,
     $                  '  ISTWERT VON NPRRUN = ',NPRRUN
            WRITE (6,*)
            WRITE (6,*) '  I = ',I,'  NPRRUN = ',NPRRUN
            WRITE (6,*) ' KKNL = ',KKNL,' JJSL = ',JJSL,' ILIMXP = ',
     $                  ILIMXP,'    JJAL = ',JJAL
            WRITE (6,*)
            WRITE (6,*) '   J     K              APHI             SPHI'
            DO 101 J = JBEG,JEND
               DO 101 K = KBEG,KKNL
  101             WRITE (6,6010) J, K, APHI(K,J,I), SPHI(K,J,I)
 6010       FORMAT (4X,I2,3X,I3,4X,2(2X,F14.7))
            CALL ERRR (501,' STMNLI   ')
         END IF
C
         NPROLD = IFIX (APHI (KKNL-13, 1, I))
         IF(DREAD .AND. (.NOT. DCONT) .AND. (NPROLD .GT. 0)) THEN
            NPROLD = NPRNEU
         END IF
C
C                                 WICHTUNGSFAKTOREN
C
         WFOLD  = FLOAT(NPROLD) / FLOAT(NPROLD + NPRRUN)
         RPRTOT = 1.0           / FLOAT(NPROLD + NPRRUN)
C
         DO 110 J = JBEG,JEND
            DO 120 K = KBEG,KEND
  120          APHI(K,J,I) = APHI(K,J,I)*WFOLD + SPHI(K,J,I)*RPRTOT
  110    CONTINUE
C
C                                 DIE NEUE GESAMTZAHL VON STICHPROBEN
C                                 WIRD VERMERKT
C
         APHI(KKNL-13, 1, I) = FLOAT (NPRRUN + NPROLD)
  100 CONTINUE
C
      RETURN
      END
