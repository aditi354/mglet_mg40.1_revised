










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
      SUBROUTINE KJIDCO  (ICODE,CDIRC,IDIRC,KAUFP,JAUFP,IAUFP)
C*STARLET***************************************************************
C        K J I D C O      IN KJIDCO WIRD DIE IM ISLINP-FELD KODIERT
C                         GESPEICHERTE INFORMATION UEBER DEN AUFPUNKT
C                         WIEDER DEKODIERT : ICODE --> CD, KAUFP, JAUFP,
C                         IAUFP
C*STARLET***************************************************************
C
C PARAM: ICODE          - NORMALERWEISE 10-STELLIGE INTEGERVARIABLE,
C                         DIE DEN AUFPUNKT IN KODIERTER FORM ENTHAELT
C        CDIRC          + CHARACTER (LEN=1) VARIABLE FUER DIE RICHTUNGS-
C                         INFORMATION
C        IDIRC          + INTEGER-VARIABLE FUER DIE RICHTUNGS-
C                         INFORMATION
C        K-, J-, IAUFP  + INDIZES DES AUFPUNKTES
C                         FALLS DIE N-RICHTUNG (N = K, J, I) EINE HOMOGE
C                         RICHTUNG IST, WIRD DAS VORZEICHEN NEGATIV !!
C                         (NOTWENDIG FUER DIE GRAPHIK ALS INFO)
C
C UPROG                 : KEINE
C
C DEFINE-DIREKTIVEN     : XHOMOG, YHOMOG, ZHOMOG
C
C        11.10.88 (HW)  : ORIGINAL
C
C*STARLET***************************************************************
C
      CHARACTER (LEN=1)  CDIRC
C
C                                 NSTELL MUSS UNBEDINGT MIT NSTELL VON
C                                 SUBR. KJIECO UEBEREINSTIMMEN !!!!
      NSTELL = 3
C
      IREST  = ICODE
      IDIRC  = IREST         / (10**(3*NSTELL))
      IREST  = IREST - IDIRC *  10**(3*NSTELL)
      KAUFP  = IREST         / (10**(2*NSTELL))
      IREST  = IREST - KAUFP *  10**(2*NSTELL)
      JAUFP  = IREST         / (10**(1*NSTELL))
      IREST  = IREST - JAUFP *  10**(1*NSTELL)
      IAUFP  = IREST
C
      IF(IDIRC .EQ. 0) CDIRC = '0'
      IF(IDIRC .EQ. 1) CDIRC = 'X'
      IF(IDIRC .EQ. 2) CDIRC = 'Y'
      IF(IDIRC .EQ. 3) CDIRC = 'Z'
C
      RETURN
      END
