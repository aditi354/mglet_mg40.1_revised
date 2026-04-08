










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
       SUBROUTINE CONTOPAR (KMX,JMX,IMX,X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                   U,V,W,P,G,H2D1,H2D2,UTO,VTO,WTO,PTO,GTO,
     $                   KMXN,JMXN,IMXN,
     $                   XN,YN,ZN,DXN,DYN,DZN,DDXN,DDYN,DDZN,
     $                   UN,VN,WN,PN,GN,HN,
     $                   ISTA,ISTO,JSTA,JSTO,
     $                   IPOS,JPOS,KPOS,ITYP,IRB,JRB,KRB,
     $                   IPROC,IPROCNBR
     $                  )     
C*MGLET*****************************************************************
C     C O N T O P A R    SETZEN DER RANDBEDINGUNGEN FUER DIE
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
C VERS:  27.07.95 (MM)  : ORIGINAL AUS BBOPAR ABGELEITET
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
     $       H2D1(JMX,IMX),  H2D2(JMX,IMX),
     $       UTO(JMX,IMX,2),    VTO(JMX,IMX,2),    WTO(JMX,IMX,2),
     $       PTO(JMX,IMX,2),    GTO(JMX,IMX,2),
     $       X(IMX),            Y(JMX),            Z(KMX),
     $      DX(IMX),           DY(JMX),           DZ(KMX),
     $     DDX(IMX),          DDY(JMX),          DDZ(KMX)

      REAL
     $    UN(KMXN,JMXN,IMXN),  VN(KMXN,JMXN,IMXN),  WN(KMXN,JMXN,IMXN),
     $    PN(KMXN,JMXN,IMXN),  GN(KMXN,JMXN,IMXN),  HN(     JMXN,IMXN),
     $       XN(IMXN),            YN(JMXN),            ZN(KMXN),
     $      DXN(IMXN),           DYN(JMXN),           DZN(KMXN),
     $     DDXN(IMXN),          DDYN(JMXN),          DDZN(KMXN)


      REAL I3DO2,I2DO2
       CHARACTER (LEN=16) CIDREC
C

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                 TEST OF DEFINE-DIRECTIVES

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      NBND = 2
CRED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-
      IF (IRB .EQ. 0) THEN
         ISTART = 1
         ISTOP  = IMX
      ELSEIF (IRB .EQ. 1) THEN
         ISTART = NBND + 1
         ISTOP  = IMX/2
      ELSEIF (IRB .EQ. 2) THEN
         ISTART = IMX/2 + 1
         ISTOP  = IMX - NBND
      ELSE
         CALL ERRR (504,'IRB BTOPAR')
      ENDIF   
C
      IF (JRB .EQ. 0) THEN
         JSTART = 1
         JSTOP  = JMX
      ELSEIF (JRB .EQ. 1) THEN
         JSTART =   2  + 1
         JSTOP  = JMX/2
      ELSEIF (JRB .EQ. 2) THEN
         JSTART = JMX/2 + 1
         JSTOP  = JMX -   2 
      ELSE
         CALL ERRR (504,'JRB BTOPAR')
      ENDIF   
CRED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-RED-BLACK-

C
C                         FINE-GRID: K-INDEX=KMX -1  (TOP)
               KF = KMX - 1
C                         K-INDEX OF COARSE-GRID
C                                   (=PARENT, CALLED NEIGHBOUR)
            KC = KPOS - 1 + (KF-1)/2
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

      CALL XTRACTXY (KMXN,JMXN,IMXN,KMX,Z,DZ,DDZ,ZN,DZN,DDZN,
     $               UN,HN,
     $               ISTA,ISTO,JSTA,JSTO,
     $               IPOS,JPOS,KPOS,0,1)
      CALL PROLONG1 (JMX,IMXN,JMXN,NBND,Y,DY,DDY,YN,DYN,DDYN,
     $               HN,H2D1,H2D2,
     $               ISTA,ISTO,JSTA,JSTO,
     $               IPOS,JPOS,0,IORDER)
      CALL PROLONG2 (JMX,IMX,IMXN,NBND,X,DX,DDX,XN,DXN,DDXN,
     $               H2D1,UTO,H2D2,
     $               ISTA,ISTO,JSTA,JSTO,
     $               IPOS,JPOS,1,IORDER)
      CALL XTRACTXY (KMXN,JMXN,IMXN,KMX,Z,DZ,DDZ,ZN,DZN,DDZN,
     $               VN,HN,
     $               ISTA,ISTO,JSTA,JSTO,
     $               IPOS,JPOS,KPOS,0,1)
      CALL PROLONG1 (JMX,IMXN,JMXN,NBND,Y,DY,DDY,YN,DYN,DDYN,
     $               HN,H2D1,H2D2,
     $               ISTA,ISTO,JSTA,JSTO,
     $               IPOS,JPOS,1,IORDER)
      CALL PROLONG2 (JMX,IMX,IMXN,NBND,X,DX,DDX,XN,DXN,DDXN,
     $               H2D1,VTO,H2D2,
     $               ISTA,ISTO,JSTA,JSTO,
     $               IPOS,JPOS,0,IORDER)


      CALL XTRACTXY (KMXN,JMXN,IMXN,KMX,Z,DZ,DDZ,ZN,DZN,DDZN,
     $               WN,HN,
     $               ISTA,ISTO,JSTA,JSTO,
     $               IPOS,JPOS,KPOS, 1,1)
      CALL PROLONG1 (JMX,IMXN,JMXN,NBND,Y,DY,DDY,YN,DYN,DDYN,
     $               HN,H2D1,H2D2,
     $               ISTA,ISTO,JSTA,JSTO,
     $               IPOS,JPOS,0,1)
      CALL PROLONG2 (JMX,IMX,IMXN,NBND,X,DX,DDX,XN,DXN,DDXN,
     $               H2D1,WTO,H2D2,
     $               ISTA,ISTO,JSTA,JSTO,
     $               IPOS,JPOS,0,1)

C-----------------------------------------------------------------------
C                                  SCALAR T
C-----------------------------------------------------------------------
C
C---- ------------------------------------------------------------------
C                                  Druck ausserhalb des lokalen Gitters


      CALL XTRACTXY (KMXN,JMXN,IMXN,KMX,Z,DZ,DDZ,ZN,DZN,DDZN,
     $               PN,HN,
     $               ISTA,ISTO,JSTA,JSTO,
     $               IPOS,JPOS,KPOS,0,IORDER_P)
      CALL PROLONG1 (JMX,IMXN,JMXN,NBND,Y,DY,DDY,YN,DYN,DDYN,
     $               HN,H2D1,H2D2,
     $               ISTA,ISTO,JSTA,JSTO,
     $               IPOS,JPOS,0,IORDER_P)
      CALL PROLONG2 (JMX,IMX,IMXN,NBND,X,DX,DDX,XN,DXN,DDXN,
     $               H2D1,PTO,H2D2,
     $               ISTA,ISTO,JSTA,JSTO,
     $               IPOS,JPOS,0,IORDER_P)


      CALL XTRACTXY (KMXN,JMXN,IMXN,KMX,Z,DZ,DDZ,ZN,DZN,DDZN,
     $               GN,HN,
     $               ISTA,ISTO,JSTA,JSTO,
     $               IPOS,JPOS,KPOS,0,1)
      CALL PROLONG1 (JMX,IMXN,JMXN,NBND,Y,DY,DDY,YN,DYN,DDYN,
     $               HN,H2D1,H2D2,
     $               ISTA,ISTO,JSTA,JSTO,
     $               IPOS,JPOS,0,1)
      CALL PROLONG2 (JMX,IMX,IMXN,NBND,X,DX,DDX,XN,DXN,DDXN,
     $               H2D1,GTO,H2D2,
     $               ISTA,ISTO,JSTA,JSTO,
     $               IPOS,JPOS,0,1)


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
                CALL ERRR (501,' BTOPAR')
       ENDIF
C
C	write (0,*)'CONTOPAR.................'
C	do i=1,imx
C	write (0,*) 'i: ',i,'   UTO: ',UTO(3,i,1)
C	enddo
      RETURN
      END











