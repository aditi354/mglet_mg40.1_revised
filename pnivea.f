










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
      SUBROUTINE PNIVEA  (KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,
     $                    KMXA,JMXA,IMXA,ISELEP,P, AP, KPP, JPP, IPP)
C*STARLET***************************************************************
C        P N I V E A      IN PNIVEA WIRD DAS DRUCKNIVEAU SO ANGEPASST,
C                         DASS AN DER STELLE  KPP, JPP, IPP  DER
C                         ENSEMBLE-MITTELWERT DES DRUCKES ZU NULL WIRD.
C                         DIESER REFERENZDRUCK WIRD AUCH VOM MOMENTANEN
C                         DRUCKFELD SUBTRAHIERT.
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KKA,JJA,IIA    - ARRAYDIMENSIONEN DER AUSWERTEFELDER
C                         BEISPIEL: FALLS J-RICHTUNG NICHT HOMOGEN IST,
C                         GILT JJA = JJ. ANDERNFALLS IST JJA = 1.
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        KMXA,JMXA,IMXA - GRENZEN DER AUSWERTEFELDER (MIT BOUND)
C        ISELEP(2,752)  - STEUERFELD F. DIE STAT. AUSWERTUNG
C                         (ENTHAELT DIE STEUERDATEN DES MOMEN-
C                         TANEN LAUFES)
C        P  (KK,JJ,II)  + MOMENTANES DRUCKFELD
C        AP (KKA,JJA,   + ENSEMBLE-MITTELWERT DES DRUCKFELDES
C            IIA     )
C
C UPROG                 : KEINE
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        20.12.89 (HW)  : ORIGINAL
C
C*STARLET***************************************************************
C
      INTEGER  KK,JJ,II,KKA,JJA,IIA,KMX,JMX,IMX,KMXA,JMXA,IMXA
      INTEGER  ISELEP (2,752)
      REAL     P  (KK,JJ,II),       AP  (KKA,JJA,IIA)
C
C                                 INDIZES FUER DEN ORT DES BEZUGSDRUCKES
C
      KPPR   = KPP
      KPPRP  = KPP + 1
      JPPR   = JPP
      IPPR   = IPP
C
      IF(ISELEP(1, 79) .EQ. 1  .OR.  ISELEP(1, 79) .EQ. 2) THEN
C
C                                 HIER: DER ENSEMBLE-MITTELWERT DES
C                                 DRUCKES EXISTIERT. ALS BEZUGSDRUCK
C                                 WIRD DAHER <P(...)> VERWENDET
C
         IF(KMXA .EQ. 1) THEN
            KPPR   = 1
            KPPRP  = 1
         END IF
         IF(JMXA .EQ. 1) THEN
            JPPR   = 1
         END IF
         IF(IMXA .EQ. 1) THEN
            IPPR   = 1
         END IF
C
         APREF  = 0.5 * (AP (KPPR , JPPR, IPPR)
     $          +        AP (KPPRP, JPPR, IPPR))
C
C                                 REDUZIERUNG DES DRUCKNIVEAUS
C
         DO 100 I = 1,IMXA
            DO 100 J = 1,JMXA
               DO 100 K = 1,KMXA
  100             AP (K,J,I) = AP (K,J,I) - APREF
C
      ELSE
C
C                                 DA <P> NICHT EXISTIERT WIRD ALS
C                                 REFERENZDRUCK  P(...) GEWAEHLT
C
         APREF  = 0.5 * ( P (KPPR , JPPR, IPPR)
     $          +         P (KPPRP, JPPR, IPPR))
      END IF
C
C                                 REDUZIERUNG DES DRUCKNIVEAUS
C
      DO 200 I = 1,IMX
         DO 200 J = 1,JMX
            DO 200 K = 1,KMX
  200          P (K,J,I) = P (K,J,I) - APREF
C
      RETURN
      END
