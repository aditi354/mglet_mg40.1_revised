










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
      SUBROUTINE SETGRD(LL,LMX,DDS,DS,S,STOT,SFRO,SCUB,LB1,LB2,MB1,MB2,
     $                  SFAK)
C*STAR*****************************************************************
C*STAR*  S E T G R D    ERZEUGT GITTERKOORDINATEN
C*STAR*****************************************************************
C
C FUNKT.:ERMITTLUNG DER GITTERKOORDINATEN BEI INAEQUIDIST. GITTER
C
C PARAM :  S(LL)        + KOORDINATENFELD
C         DS(LL)        + ABSTAND DER GITTERPUNKTE
C        DDS(LL)        + KANTENLAENGE DER KONTROLLVOLUMINA
C        STOT           - GESAMTLAENGE BERECH.- GEBIET
C        SFRO           - LAENGE VOR CUBUS
C        SCUB           - LAENGE CUBUS
C        LL             - FELDDIMENSIONEN
C        LMX            - GRENZE D. BERCHNUNGSGEBIETES(MIT BOUND)
C        LB1            - CUBUS -  ANFANG
C        LB2            - CUBUS -  ENDE
C        MB1            - ZAHL AEQUIDIST. ZELLEN VOR CUBUS
C        MB2            - ZAHL AEQUIDIST. ZELLEN NACH CUBUS
C        SFAK           - DIE ZELLEN AN DEN KANTEN DES KUBUS HABEN
C                         DIE ABMESSUNG  SFAK*DX(AEQUIDISTANTES GITTER)
C                         = SFAK * SCUB / (ANZAHL DER ZELLEN IM KUBUS)
C
C UPROG                 : ERRR,  ZBRENT
C
C DEFINE-DIREKTIVEN     : KEINE
C
C VERS:  10.04.84 (SF)  - ORIGINAL
C        19.12.84 (HW)  - "HINTERER BEREICH" KANN NULL ZELLEN BEINHALTEN
C        04.01.85 (HW)  - FAKTOR "B" AUF 3.0 GESETZT (FRUEHER : B=2.0)
C        26.02.85 (HW)  : BEREICHSWEISE AEQUIDISTANTES GITTER AM KUBUS
C                         MOEGLICH (FAKTOR "B" AUF 5.0 GESETZT)
C        11.03.85 (HW)  : IN Z-RICHTUNG (GENAUER FALLS LB1 .EQ. LSTART)
C                         KANN DS IM BEREICH VON IB1 VON DS IM BEREICH
C                         IB2 UNTERSCHIEDLICH SEIN. DS(IB1) = SF12*
C                         DS(IB2)
C        14.11.88 (HW)  : SFAK WIRD JETZT UEBERGEBEN
C
C*STAR*****************************************************************
C
      COMMON /DBLOCK/ DSC2,SDIST,NDIST,MDIST
      EXTERNAL FGRID
C
C
      REAL S(LL),  DS(LL),  DDS(LL)
C
      DATA LSTART /3/
      DATA EPS,NSIG /0.0,13/
C
      WRITE (6,100)
  100 FORMAT(1H1,9(1H-),13H S E T G R D ,110(1H-))
C
      LSTOP = LMX-2
      SBAC = STOT-SFRO-SCUB
C
C                                  BEREICHSWEISE AEQUIDISTANTES GITTER
C                                  AM KUBUS
C
C                                  DS(IB1) = SF12 * DS(IB2)
      SF12   = 1.0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C                                  FUER DIE ZURUECKSPRINGENDE STUFE
C
*     IF(LB1 .EQ. LSTART  .AND.  SFAK .NE. 1.0) SF12 = 2.0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DSC2   = SCUB/FLOAT(LB2-LB1+1) * SFAK
      IF(DSC2 .LE. 0.0) DSC2 = STOT/FLOAT(LMX-4)
      DSC1   = DSC2 * SF12
      IF(LB2-LB1-MB1-MB2 .GT. 0) GOTO 2001
         DSCL = DSC1
         GOTO 2002
 2001 DSCL   = (SCUB - DSC1*(FLOAT(MB1)+0.5) - DSC2*(FLOAT(MB2)+0.5))
     $       / FLOAT(LB2-LB1-MB1-MB2)
C
 2002 LAEQ1  = LB1-MB1-1
      IF((LB1 .EQ. LSTART) .AND. (LAEQ1 .LE. 0)) LAEQ1 = 1
      LAEQ2  = LB1+MB1-1
      DO 10 L= LAEQ1,LAEQ2
   10    DS(L) = DSC1
C
      LAEQ1  = LAEQ2+1
      LAEQ2  = LB2-MB2-1
      IF(LAEQ2 .LT. LAEQ1) GOTO 2100
      DO 11 L= LAEQ1,LAEQ2
   11    DS(L) = DSCL
C
 2100 IF(LAEQ2 .EQ. LSTOP) GOTO 3000
C
C
      LAEQ1  = LAEQ2+1
      LAEQ2  = LB2+MB2
      IF(LAEQ2 .LT. LAEQ1) GOTO 2200
      DO 12 L= LAEQ1,LAEQ2
   12    DS(L) = DSC2
C
C                                  INAEQUIDIST. GITTER FUER HINTEREN
C                                  BEREICH
 2200 SDIST = SBAC
      NDIST = LSTOP-LB2
C
      IF(NDIST .LE. 0) GOTO 3000
C
      MDIST = MB2
      IF(NDIST-MDIST .LT. 3) CALL ERRR(600,' SETGRD3  ')
C
      MAXFN = 1000
      A = .49
      B = 5.0
      CALL ZBRENT(FGRID,EPS,NSIG,A,B,MAXFN,IER)
      IF(IER .NE. 0) CALL ERRR(IER+600,' SETGRD3  ')
      WRITE(6,110)A,B,MAXFN,NSIG
  110 FORMAT(1H0,10X,17HHINTERER BEREICH:,//11X,3HA: ,F9.6,3X,3HB: ,
     $       F9.6,3X,7HMAXFN: ,I4,4X,6HNSIG: ,I2)
C
C                                  DS IN POS. RICHTG. ERMITTELN
C
      LVAR1 = LAEQ2+1
      LVAR2 = LSTOP-1
      DO 20 L=LVAR1,LVAR2
   20    DS(L) = DS(L-1)*B
C
 3000 DS(LSTOP) = DS(LSTOP-1)
      DS(LSTOP+1) = DS(LSTOP-1)
      DS(LSTOP+2) = DS(LSTOP-1)
C
C                                  INAEQUIDIST. GITTER FUER VORDEREN
C                                  BEREICH
      IF(LB1 .LE. LSTART) GOTO 4000
C
      SDIST = SFRO
      NDIST = LB1-LSTART
      MDIST = MB1
      IF(NDIST-MDIST .LT. 3) CALL ERRR(800,' SETGRD3  ')
C
      MAXFN = 1000
      A = .49
      B = 5.0
      CALL ZBRENT(FGRID,EPS,NSIG,A,B,MAXFN,IER)
      IF(IER .NE. 0) CALL ERRR(IER+800,' SETGRD3  ')
      WRITE(6,120)A,B,MAXFN,NSIG
  120 FORMAT(1H0,10X,17HVORDERER BEREICH:,//11X,3HA: ,F9.6,3X,3HB: ,
     $       F9.6,3X,7HMAXFN: ,I4,4X,6HNSIG: ,I2)
C
C                                  DS IN NEG. RICHTG. ERMITTELN
C
      LVAR1 = LSTART
      LVAR2 = LB1-MB1-2
      DO 30 LB=LVAR1,LVAR2
         L = LVAR2-LB+LVAR1
   30    DS(L) = DS(L+1)*B
C
 4000 CONTINUE
      DS(LSTART-1) = DS(LSTART)
      DS(LSTART-2) = DS(LSTART)
C
C                                  BESTIMMUNG V. DDS
C
      DDS(1) = DS(1)
      DO 50 L=1,LSTOP
   50    DDS(L+1) = (DS(L)+DS(L+1))*0.5
      DDS(LMX) = DS(LSTOP)
C
C                                  BESTIMMUNG V. S
C
      S(1) = 0.0
      DO 60 L=2,LMX
   60    S(L) = S(L-1)+DS(L-1)
      SST = -((S(LB1)+S(LB1-1))/2.0)
      IF(SCUB .EQ. 0.0) SST = -(S(3)+S(2))*0.5
      IF(LB1 .LE. LSTART) GOTO 5000
      SST = -(S(LB1)+(S(LB2)-S(LB1))/2.0)
 5000 CONTINUE
      DO 70 L=1,LMX
   70    S(L) = SST+S(L)
C
C                                  AUSGABE ALLER WERTE
C
      WRITE(6,130)(K,K=1,9)
  130 FORMAT(1H0,2X,9(11X,I1),11X,1H0)
      WRITE(6,140)
  140 FORMAT(1H0,4X,2HS:)
      WRITE(6,190)S
      WRITE(6,160)
  160 FORMAT(1H0,3X,3HDS:)
      WRITE(6,190)DS
      WRITE(6,180)
  180 FORMAT(1H0,2X,4HDDS:)
      WRITE(6,190)DDS
  190 FORMAT(1H+,5(10X,10(F9.6,3X),/,1X))
      RETURN
      END
      REAL FUNCTION FGRID(A)
C*STAR*****************************************************************
C*STAR*  F G R I D     LIEFERT 0.0 FUER RICHTIGES A
C*STAR*****************************************************************
C
C  PARAM:   A           - ARGUMENT VON FGRID
C        // DSC2        - ABSTAND DER GITTERPUNKTE AM CUBUS
C                         IM BEREICH DER VORDER- BZW. HINTERKANTE
C        // SDIST       - LAENGE VOR OD. HINTER CUBUS
C        // NDIST       - N PUNKTE AUF SDIST
C        // MDIST       - ZAHL AEQUIDIST. ZELLEN VOR/NACH CUBUS
C
C*STAR*****************************************************************
C
      COMMON /DBLOCK/ DSC2,SDIST,NDIST,MDIST
C
C
      NPOTMX = NDIST-MDIST-1
      SUM = (0.5+FLOAT(MDIST))*DSC2-SDIST
C
      DO 10 M=1,NPOTMX
   10    SUM = SUM+DSC2*A**M
C
      FGRID = SUM+0.5*DSC2*A**NPOTMX
C
      RETURN
      END
