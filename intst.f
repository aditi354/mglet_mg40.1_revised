










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
      SUBROUTINE INTST   (PHIA,KSTAGA,JSTAGA,ISTAGA,
     $                    PHIB,KSTAGB,JSTAGB,ISTAGB,
     $                    PHIAIN,PHIBIN,K, J, I, KK, JJ, II)
C*STARLET***************************************************************
C        I N T S T        ZWEI IM DREIDIMENSIONALEN FELD BELIEBIG VER-
C                         SCHOBENE GROESSEN 'PHIA' UND 'PHIB' WERDEN
C                         AUF EINEN GEMEINSAMEN PUNKT INTERPOLIERT,
C                         DER AUS DER VERSCHIEBUNG VON 'PHIA' UND
C                         'PHIB' BESTIMMT WIRD. NSTAGA UND NSTAGB
C                         DUERFEN BELIEBIG VERSCHIEDEN SEIN.
C*STARLET***************************************************************
C
C PARAM: PHIA(KK,JJ,II) - ENTHAELT DIE FLUKTUATIONEN DER GROESSE PHIA
C        K-, J-, ISTAGA - KENNZEICHNET DIE VERSCHIEBUNG DER GROESSE 'A'
C                         NSTAGA = 0: KEINE VERSCHIEBUNG IN N-RICHTUNG
C                         NSTAGA = 1: DIE GROESSE 'A' IST IN POSITIVER
C                                     N-RICHTUNG IM MASCHENGITTER VER-
C                                     SCHOBEN
C        PHIB(KK,JJ,II) - ENTHAELT DIE FLUKTUATIONEN DER GROESSE PHIB
C        K-, J-, ISTAGB - KENNZEICHNET DIE VERSCHIEBUNG DER GROESSE 'B'
C                         NSTAGB = 0: KEINE VERSCHIEBUNG IN N-RICHTUNG
C                         NSTAGB = 1: DIE GROESSE 'B' IST IN POSITIVER
C                                     N-RICHTUNG IM MASCHENGITTER VER-
C                                     SCHOBEN
C        PHIAIN         + INTERPOLIERTER WERT DER GROESSE 'A' AM
C                         GITTERPUNKT K, J, I
C        PHIBIN         + INTERPOLIERTER WERT DER GROESSE 'B' AM
C                         GITTERPUNKT K, J, I
C        K, J, I        - MOMENTAN BETRACHTETER MASCHENGITTERPUNKT
C        KK, JJ, II     - ARRAYDIMENSIONEN
C
C UPROG                 : ERRR
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        14.09.88 (HW)  : ORIGINAL
C
C*STARLET***************************************************************
C
      REAL     PHIA (KK, JJ, II ),   PHIB (KK, JJ, II )
C
C                                 UEBERPRUEFUNG, OB DIE GROESSEN
C                                 IN ZULAESSIGER WEISE VERSCHOBEN SIND
C
      LTNULL = MIN0 (KSTAGA,JSTAGA,ISTAGA,KSTAGB,JSTAGB,ISTAGB)
      MTEINS = MAX0 (KSTAGA,JSTAGA,ISTAGA,KSTAGB,JSTAGB,ISTAGB)
C
      IF(LTNULL .LT. 0 .OR. MTEINS .GT. 1) CALL ERRR (501,' INTST    ')
C
C                                 BETRAG DER DIFFERENZ ZWISCHEN
C                                 NSTAGA UND NSTAGB. (IST ENTWEDER
C                                 0.0 ODER 1.0)
C
      IDKS   = IABS (KSTAGA - KSTAGB)
      IDJS   = IABS (JSTAGA - JSTAGB)
      IDIS   = IABS (ISTAGA - ISTAGB)
      DKS    = FLOAT (IDKS)
      DJS    = FLOAT (IDJS)
      DIS    = FLOAT (IDIS)
C
      FRACT  = 1.0 / FLOAT(2**(IDKS + IDJS + IDIS))
C
C                                 INTERPOLATION DER VARIABLEN "A"
C
      PHIAIN = FRACT * (
     $         PHIA (K        ,J        ,I        )
     $       + PHIA (K +KSTAGB,J        ,I        ) * DKS
     $       + PHIA (K        ,J +JSTAGB,I        )       * DJS
     $       + PHIA (K        ,J        ,I +ISTAGB)             * DIS
     $       + PHIA (K        ,J +JSTAGB,I +ISTAGB)       * DJS * DIS
     $       + PHIA (K +KSTAGB,J        ,I +ISTAGB) * DKS       * DIS
     $       + PHIA (K +KSTAGB,J +JSTAGB,I        ) * DKS * DJS
     $       + PHIA (K +KSTAGB,J +JSTAGB,I +ISTAGB) * DKS * DJS * DIS
     $       )
C
C                                 INTERPOLATION DER VARIABLEN "B"
C
      PHIBIN = FRACT * (
     $         PHIB (K        ,J        ,I        )
     $       + PHIB (K +KSTAGA,J        ,I        ) * DKS
     $       + PHIB (K        ,J +JSTAGA,I        )       * DJS
     $       + PHIB (K        ,J        ,I +ISTAGA)             * DIS
     $       + PHIB (K        ,J +JSTAGA,I +ISTAGA)       * DJS * DIS
     $       + PHIB (K +KSTAGA,J        ,I +ISTAGA) * DKS       * DIS
     $       + PHIB (K +KSTAGA,J +JSTAGA,I        ) * DKS * DJS
     $       + PHIB (K +KSTAGA,J +JSTAGA,I +ISTAGA) * DKS * DJS * DIS
     $       )
C
      RETURN
      END
