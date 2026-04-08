










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
      SUBROUTINE BODYST(KK,JJ,II,KMX,JMX,IMX,NGRDMX,
     1                  MAXCBP,B,NBND,IGRDMX,
     1                  LALFAX,LGAMAX,LALFAY,LGAMAY,LALFAZ,LGAMAZ,
     2                  ICUBPT,JCUBPT,KCUBPT,IDXGRD,LCUB)

      INTEGER II,JJ,KK,NBND,JUMP,NGRDMX,MAXCBP,GRID,IDXAKT,IGRDMX
      INTEGER ICUBPT(MAXCBP),JCUBPT(MAXCBP),KCUBPT(MAXCBP)
      INTEGER IDXGRD(9,NGRDMX),IMX,JMX,KMX
      REAL    LCUB(8,MAXCBP)
      REAL LALFAX(II,NGRDMX),LGAMAX(II,NGRDMX),LALFAY(JJ,NGRDMX),
     1     LGAMAY(JJ,NGRDMX),LALFAZ(KK,NGRDMX),LGAMAZ(KK,NGRDMX)
      REAL B(KK,JJ,II)

C     BESTIMMUNG DER START- UND ENDWERTE DER EINZELNEN
C     SCHLEIFEN

      INTEGER ISTA0,IEND0,ISTA1,IEND1,JSTA0,JEND0,JSTA1,JEND1,KSTA0,
     1        KEND0,KSTA1,KEND1,IMGSTA,IMGEND,JMGSTA,JMGEND,KMGSTA,
     2        KMGEND

      JUMP = 1

      ISTA0 = NBND + 1
      ISTA1 = ISTA0 + JUMP
      IMGSTA = NBND + 1

      IEND0 = IMX - NBND
      IEND1 = IEND0
      IMGEND = IMX - NBND


      JSTA0 = NBND + 1
      JSTA1 = JSTA0 + JUMP
      JMGSTA = NBND + 1

      JEND0 = JMX - NBND
      JEND1 = JEND0
      JMGEND = JMX - NBND


      KSTA0 = NBND + 1
      KSTA1 = KSTA0 + JUMP
      KMGSTA = NBND + 1

      KEND0 = KMX - NBND
      KEND1 = KEND0
      KMGEND = KMX - NBND

      IDXAKT = 0
      JUMP = 1

      DO GRID = IGRDMX,1,-1

        IDXGRD(1,GRID) = IDXAKT

C       GRID(0,0,0)
        CALL SETIDX(KK,JJ,II,KMX,JMX,IMX,B,JUMP,NBND,MAXCBP,
     1              ISTA0,IEND0,JSTA0,JEND0,KSTA0,KEND0,
     2              LALFAX(1,GRID),LGAMAX(1,GRID),
     3              LALFAY(1,GRID),LGAMAY(1,GRID),
     4              LALFAZ(1,GRID),LGAMAZ(1,GRID),
     5              IDXAKT,ICUBPT,JCUBPT,KCUBPT,IDXGRD(2,GRID),LCUB)


C     GRID(0,1,1)

        CALL SETIDX(KK,JJ,II,KMX,JMX,IMX,B,JUMP,NBND,MAXCBP,
     1              ISTA1,IEND1,JSTA1,JEND1,KSTA0,KEND0,
     2              LALFAX(1,GRID),LGAMAX(1,GRID),
     3              LALFAY(1,GRID),LGAMAY(1,GRID),
     4              LALFAZ(1,GRID),LGAMAZ(1,GRID),
     5              IDXAKT,ICUBPT,JCUBPT,KCUBPT,IDXGRD(3,GRID),LCUB)


C     GRID(1,0,1)

        CALL SETIDX(KK,JJ,II,KMX,JMX,IMX,B,JUMP,NBND,MAXCBP,
     1              ISTA1,IEND1,JSTA0,JEND0,KSTA1,KEND1,
     2              LALFAX(1,GRID),LGAMAX(1,GRID),
     3              LALFAY(1,GRID),LGAMAY(1,GRID),
     4              LALFAZ(1,GRID),LGAMAZ(1,GRID),
     5              IDXAKT,ICUBPT,JCUBPT,KCUBPT,IDXGRD(4,GRID),LCUB)

C     GRID(1,1,0)

        CALL SETIDX(KK,JJ,II,KMX,JMX,IMX,B,JUMP,NBND,MAXCBP,
     1              ISTA0,IEND0,JSTA1,JEND1,KSTA1,KEND1,
     2              LALFAX(1,GRID),LGAMAX(1,GRID),
     3              LALFAY(1,GRID),LGAMAY(1,GRID),
     4              LALFAZ(1,GRID),LGAMAZ(1,GRID),
     5              IDXAKT,ICUBPT,JCUBPT,KCUBPT,IDXGRD(5,GRID),LCUB)

C     RED-SCHRITT ABGESCHLOSSEN 

C     BLACK POINTS !
C     GRID(1,0,0)

        CALL SETIDX(KK,JJ,II,KMX,JMX,IMX,B,JUMP,NBND,MAXCBP,
     1              ISTA0,IEND0,JSTA0,JEND0,KSTA1,KEND1,
     2              LALFAX(1,GRID),LGAMAX(1,GRID),
     3              LALFAY(1,GRID),LGAMAY(1,GRID),
     4              LALFAZ(1,GRID),LGAMAZ(1,GRID),
     5              IDXAKT,ICUBPT,JCUBPT,KCUBPT,IDXGRD(6,GRID),LCUB)

C     GRID(1,1,1)

        CALL SETIDX(KK,JJ,II,KMX,JMX,IMX,B,JUMP,NBND,MAXCBP,
     1              ISTA1,IEND1,JSTA1,JEND1,KSTA1,KEND1,
     2              LALFAX(1,GRID),LGAMAX(1,GRID),
     3              LALFAY(1,GRID),LGAMAY(1,GRID),
     4              LALFAZ(1,GRID),LGAMAZ(1,GRID),
     5              IDXAKT,ICUBPT,JCUBPT,KCUBPT,IDXGRD(7,GRID),LCUB)

C     GRID(0,0,1)

        CALL SETIDX(KK,JJ,II,KMX,JMX,IMX,B,JUMP,NBND,MAXCBP,
     1              ISTA1,IEND1,JSTA0,JEND0,KSTA0,KEND0,
     2              LALFAX(1,GRID),LGAMAX(1,GRID),
     3              LALFAY(1,GRID),LGAMAY(1,GRID),
     4              LALFAZ(1,GRID),LGAMAZ(1,GRID),
     5              IDXAKT,ICUBPT,JCUBPT,KCUBPT,IDXGRD(8,GRID),LCUB)


C     GRID(0,1,0)

        CALL SETIDX(KK,JJ,II,KMX,JMX,IMX,B,JUMP,NBND,MAXCBP,
     1              ISTA0,IEND0,JSTA1,JEND1,KSTA0,KEND0,
     2              LALFAX(1,GRID),LGAMAX(1,GRID),
     3              LALFAY(1,GRID),LGAMAY(1,GRID),
     4              LALFAZ(1,GRID),LGAMAZ(1,GRID),
     5              IDXAKT,ICUBPT,JCUBPT,KCUBPT,IDXGRD(9,GRID),LCUB)

C     BLACK-SCHRITT ABGESCHLOSSEN 

C     NEUES GITTER
        JUMP = 2 * JUMP

        ISTA1 = ISTA0 + JUMP



        JSTA1 = JSTA0 + JUMP



        KSTA1 = KSTA0 + JUMP


      ENDDO

      RETURN
      END

C23456
