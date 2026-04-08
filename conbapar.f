










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
       SUBROUTINE CONBAPAR (KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                   U,V,W,P,G,H2D1,H2D2,UBA,VBA,WBA,PBA,GBA,
     $                   KMXN,JMXN,IMXN,
     $                   XN,YN,ZN,DXN,DYN,DZN,DDXN,DDYN,DDZN,
     $                   UN,VN,WN,PN,GN,HN,
     $                   JSTA,JSTO,KSTA,KSTO,
     $                   IPOS,JPOS,KPOS,ITYP,IRB,JRB,KRB,
     $                   IPROC,IPROCNBR
     $                  )      
C*MGLET*****************************************************************
C     C O N B A P A R    SETZEN DER RANDBEDINGUNGEN FUER DIE
C                         GESCHWINDIGKEITSFELDER. AUS PARENT-Gitter
C*MGLET*****************************************************************
C
C PARAM: 
C        KMX, JMX, IMX  - DIMENSIONEN DES GITTERS 
C        U(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        V(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        W(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        T(KK,JJ,II)    + SCALAR T FIELD
C        P(KMX,JMX,IMX) - DRUCKFELD
C        G(KMX,JMX,IMX) - EFFEKTIVE DYNAMISCHE VISKOSITAET (= MUE)
C        ...N           - ENTSPRECHENDE GROESSE AUS DEM NACHBARGITTER
C        JSTART,JSTOP
C        KSTART,KSTOP   - BEREICH IN DEM RANDBEDINGUNG GESETZT WIRD
C        IPOS,JPOS,KPOS - POSITION IM NACHBARGITTER, AUF DER DER PUNKT
C                         MIT DEN INDIZES (3,3,3) ZU LIEGEN KOMMT
C        ITYP           - HOLLERITH-KONSTANTE:
C                         'P'  VOR DER DRUCKKORREKTUR
C                         'T'  VOR DEM ZEITSCHRITT
C            IRB                  - INTEGER-KENNZAHL I-RICHTUNG RED-BLACK
C            JRB                  - INTEGER-KENNZAHL J-RICHTUNG RED-BLACK
C            KRB                  - INTEGER-KENNZAHL K-RICHTUNG RED-BLACK
C
C VERS:  27.10.97 (MM)  : ORIGINAL AUS BBOPAR ABGELEITET
C        28.02.03 (TB)  : SCALAR BOUNDARY TREATMENT IMPLEMENTED
C
C DEFINE-DIREKTIVEN     : QUD 
C
C*STARLET***************************************************************

C
      CHARACTER (LEN=1)  ITYP
C
      REAL
     $     U(KMX,JMX,IMX),  V(KMX,JMX,IMX),    W(KMX,JMX,IMX),
     $     P(KMX,JMX,IMX),  G(KMX,JMX,IMX),      
     $       H2D1(KMX,JMXN),  H2D2(JMX,IMX),
     $       UBA(KMX,JMX,2),    VBA(KMX,JMX,2),    WBA(KMX,JMX,2),
     $       PBA(KMX,JMX,2),    GBA(KMX,JMX,2),
     $       X(IMX),            Y(JMX),            Z(KMX),
     $      DX(IMX),           DY(JMX),           DZ(KMX),
     $     DDX(IMX),          DDY(JMX),          DDZ(KMX)

      REAL
     $    UN(KMXN,JMXN,IMXN),  VN(KMXN,JMXN,IMXN),  WN(KMXN,JMXN,IMXN),
     $    PN(KMXN,JMXN,IMXN),  GN(KMXN,JMXN,IMXN),  HN(KMXN,JMXN     ),
     $       XN(IMXN),            YN(JMXN),            ZN(KMXN),
     $      DXN(IMXN),           DYN(JMXN),           DZN(KMXN),
     $     DDXN(IMXN),          DDYN(JMXN),          DDZN(KMXN)

      REAL I3DO2,I2DO2
       CHARACTER (LEN=16) CIDREC
C
C

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                 TEST OF DEFINE-DIRECTIVES

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      NBND = 2
C
C
C                         FINE-GRID: I-INDEX=IMX -1  (BAC)
              IIF = IMX - 1
C                         I-INDEX OF COARSE-GRID
C                                   (=PARENT, CALLED NEIGHBOUR)
            IC = IPOS - 1 + (IF-1)/2
C
C
C                                     ORDNUNG DER INTERPOLATION
C
          IORDER = 1
C
C
          IORDER_P = 1
          IORDER = 1

C
C                                 FUER ZEITSCHRITT
      IF (ITYP.EQ.'T') THEN

      CALL XTRACTYZ (KMXN,JMXN,IMXN,IMX,X,DX,DDX,XN,DXN,DDXN,
     $               UN,HN,
     $               JSTA,JSTO,KSTA,KSTO,
     $               IPOS,JPOS,KPOS, 1,1)
      CALL PROLONG1 (KMX,JMXN,KMXN,NBND,Z,DZ,DDZ,ZN,DZN,DDZN,
     $               HN,H2D1,H2D2,
     $               JSTA,JSTO,KSTA,KSTO,
     $               JPOS,KPOS,0,IORDER)
      CALL PROLONG2 (KMX,JMX,KMXN,NBND,Y,DY,DDY,YN,DYN,DDYN,
     $               H2D1,UBA,H2D2,
     $               JSTA,JSTO,KSTA,KSTO,
     $               JPOS,KPOS,0,IORDER)
      CALL XTRACTYZ (KMXN,JMXN,IMXN,IMX,X,DX,DDX,XN,DXN,DDXN,
     $               VN,HN,
     $               JSTA,JSTO,KSTA,KSTO,
     $               IPOS,JPOS,KPOS,0,1)
      CALL PROLONG1 (KMX,JMXN,KMXN,NBND,Z,DZ,DDZ,ZN,DZN,DDZN,
     $               HN,H2D1,H2D2,
     $               JSTA,JSTO,KSTA,KSTO,
     $               JPOS,KPOS,0,IORDER)
      CALL PROLONG2 (KMX,JMX,KMXN,NBND,Y,DY,DDY,YN,DYN,DDYN,
     $               H2D1,VBA,H2D2,
     $               JSTA,JSTO,KSTA,KSTO,
     $               JPOS,KPOS,1,IORDER)


      CALL XTRACTYZ (KMXN,JMXN,IMXN,IMX,X,DX,DDX,XN,DXN,DDXN,
     $               WN,HN,
     $               JSTA,JSTO,KSTA,KSTO,
     $               IPOS,JPOS,KPOS,0,1)
      CALL PROLONG1 (KMX,JMXN,KMXN,NBND,Z,DZ,DDZ,ZN,DZN,DDZN,
     $               HN,H2D1,H2D2,
     $               JSTA,JSTO,KSTA,KSTO,
     $               JPOS,KPOS,1,IORDER)
      CALL PROLONG2 (KMX,JMX,KMXN,NBND,Y,DY,DDY,YN,DYN,DDYN,
     $               H2D1,WBA,H2D2,
     $               JSTA,JSTO,KSTA,KSTO,
     $               JPOS,KPOS,0,IORDER)

C-----------------------------------------------------------------------
C                                  SCALAR T
C-----------------------------------------------------------------------
C
C---- ------------------------------------------------------------------
C                                  Druck ausserhalb des lokalen Gitters


      CALL XTRACTYZ (KMXN,JMXN,IMXN,IMX,X,DX,DDX,XN,DXN,DDXN,
     $               PN,HN,
     $               JSTA,JSTO,KSTA,KSTO,
     $               IPOS,JPOS,KPOS,0,IORDER_P)
      CALL PROLONG1 (KMX,JMXN,KMXN,NBND,Z,DZ,DDZ,ZN,DZN,DDZN,
     $               HN,H2D1,H2D2,
     $               JSTA,JSTO,KSTA,KSTO,
     $               JPOS,KPOS,0,IORDER_P)
      CALL PROLONG2 (KMX,JMX,KMXN,NBND,Y,DY,DDY,YN,DYN,DDYN,
     $               H2D1,PBA,H2D2,
     $               JSTA,JSTO,KSTA,KSTO,
     $               JPOS,KPOS,0,IORDER_P)

      CALL XTRACTYZ (KMXN,JMXN,IMXN,IMX,X,DX,DDX,XN,DXN,DDXN,
     $               GN,HN,
     $               JSTA,JSTO,KSTA,KSTO,
     $               IPOS,JPOS,KPOS,0,1)
      CALL PROLONG1 (KMX,JMXN,KMXN,NBND,Z,DZ,DDZ,ZN,DZN,DDZN,
     $               HN,H2D1,H2D2,
     $               JSTA,JSTO,KSTA,KSTO,
     $               JPOS,KPOS,0,1)
      CALL PROLONG2 (KMX,JMX,KMXN,NBND,Y,DY,DDY,YN,DYN,DDYN,
     $               H2D1,GBA,H2D2,
     $               JSTA,JSTO,KSTA,KSTO,
     $               JPOS,KPOS,0,1)


C                                 FUER DRUCKKORREKTUR
C                                 NUR NORMALKOMPONENTE MUSS GESETZT
C                                 WERDEN

       ELSEIF (ITYP.EQ.'P') THEN

C                   NORMAL-KOMPONENT MUST NOT BE SET
C                   BECAUSE PRESSURE IN THE BOUNDARY-CELLS
C                   IS NOT MODIFIED, SO NORMAL-KOMPONENTS MUST
C                   BE ABLE TO DEVELOP FREELY

               RETURN


C
       ELSE
                CALL ERRR (501,' CONBAPAR')
       ENDIF
C
C	write (0,*)'CONBAPAR.................'
C	do k=1,kmx
C	write (0,*) 'k: ',k,'   WBA: ',WBA(k,3,1)
C	enddo

      RETURN
      END











