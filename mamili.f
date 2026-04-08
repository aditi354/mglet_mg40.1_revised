










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
      SUBROUTINE MAMILI  (PHI,JEB,KKNL,JJNL,ILIMXP,
     $                    KBEG,KEND,JBEG,JEND,IBEG,IEND)
C*STARLET***************************************************************
C        M A M I L I      FOLGENDE MAXIMA UND MINIMA FUER DIE "LINIEN"-
C                         FELDER WERDEN BESTIMMT:
C                         - MAXIMUM DES ORDINATENWERTES (ORDMAX) UND
C                           DER ZUGEHOERIGE ABSZISSENWERT (PEAKMA)
C                         - MINIMUM DES ORDINATENWERTES (ORDMIN) UND
C                           DER ZUGEHOERIGE ABSZISSENWERT (PEAKMI)
C                         - MAXIMUM DES ABSZISSENWERTES (ABSMAX)
C                         - MINIMUM DES ABSZISSENWERTES (ABSMIN)
C*STARLET***************************************************************
C
C PARAM: PHI (KKNL,     + BELIEBIGES "LINIEN"-FELD
C        JJNL,ILIMXP)
C        JEB            - GIBT AN, IN WELCHER J-EBENE DIE ABSZISSENWERTE
C                         STEHEN
C        KKNL,JJNL,     - ARRAYDIMENSIONEN
C        ILIMXP
C        KBEG           - BEGINN D. DO-SCHLEIFEN F. K-RICHTUNG
C        KEND           + ENDE D. DO-SCHLEIFEN F. K-RICHTUNG
C        JBEG,JEND      - GRENZEN D. DO-SCHLEIFEN F. J-RICHTUNG
C        IBEG,IEND      - GRENZEN D. DO-SCHLEIFEN F. I-RICHTUNG
C
C DEFINE DIREKTIVEN     : KEINE
C
C UPROG                 : ERRR
C
C        19.12.88 (HW)  : ORIGINAL
C
C*STARLET***************************************************************
C

      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
      REAL            PHI (KKNL,JJNL,ILIMXP)
C
      IF(JEB .NE. 1  .AND.  JEB .NE. 2) CALL ERRR (501,' MAMILI   ')
C
      DO 100 IL = IBEG,IEND
         KEND   = IFIX (PHI (KKNL-12, 1, IL))
         DO 100 J = JBEG,JEND
            ORDMAX = - (ABS (GREAT))
            ABSMAX = - (ABS (GREAT))
            PEAKMA = - (ABS (GREAT))
            ORDMIN =    ABS (GREAT)
            ABSMIN =    ABS (GREAT)
            PEAKMI =    ABS (GREAT)
            DO 200 K = KBEG,KEND
C
C                                 ZERO.. = 1.0 : NEUES MAXIMUM BZW.
C                                                MINIMUM
C                                 ZERO.. = 0.0 : BETRACHTETER WERT IST
C                                                UNINTERESSANT FUER MAX.
C                                                BZW. MIN.
C
               ZEROMA = 0.5 + SIGN (0.5,(PHI(K,J,IL) - ORDMAX))
               ZEROMI = 0.5 + SIGN (0.5,(ORDMIN - PHI(K,J,IL)))
               ORDINA = PHI (K,J  ,      IL        )
               ABSZIS = PHI (K,JEB,2-JEB+IL*(JEB-1))
               ORDMAX = ORDMAX * (1.0 - ZEROMA)
     $                + ORDINA *        ZEROMA
               PEAKMA = PEAKMA * (1.0 - ZEROMA)
     $                + ABSZIS *        ZEROMA
C
               ORDMIN = ORDMIN * (1.0 - ZEROMI)
     $                + ORDINA *        ZEROMI
               PEAKMI = PEAKMI * (1.0 - ZEROMI)
     $                + ABSZIS *        ZEROMI
               ABSMAX = AMAX1 (ABSMAX, ABSZIS)
  200          ABSMIN = AMIN1 (ABSMIN, ABSZIS)
C
C                                 DIE GEFUNDENEN MAXIMA U. MINIMA WERDEN
C                                 IM PHI-FELD GESPEICHERT (Z.B. FUER
C                                 DIE GRAPHIK)
C
            PHI (KKNL-16, J,IL) = ORDMAX
            PHI (KKNL-17, J,IL) = ORDMIN
            PHI (KKNL-18, J,IL) = ABSMAX
            PHI (KKNL-19, J,IL) = ABSMIN
            PHI (KKNL-20, J,IL) = PEAKMA
            PHI (KKNL-21, J,IL) = PEAKMI
  100 CONTINUE
C
      RETURN
      END
