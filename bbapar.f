










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
      SUBROUTINE BBAPAR (KMX,JMX,IMX,
     $                   U,V,W,P,G,UBA,VBA,WBA,PBA,GBA,
     $                   JSTART,JSTOP,KSTART,KSTOP,
     $                   ITYP,IRB,JRB,KRB,DX,DY,DZ,DDX,DDY,DDZ
     $                   )      
C*MGLET*****************************************************************
C        B B A P A R    SETZEN DER RANDBEDINGUNGEN FUER DIE
C                GESCHWINDIGKEITSFELDER AND SCALAR T. AUS PARENT-GITTERN
C*MGLET*****************************************************************
C
C PARAM: 
C        KMX, JMX, IMX  - DIMENSIONEN DES GITTERS 
C        U(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        V(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        W(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        T(KMX,JMX,IMX) + SCALAR T FIELD
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
C        IRB            - INTEGER-KENNZAHL I-RICHTUNG RED-BLACK
C        JRB            - INTEGER-KENNZAHL J-RICHTUNG RED-BLACK
C        KRB            - INTEGER-KENNZAHL K-RICHTUNG RED-BLACK
C        IGRID          - GITTER, FUER DAS DIE RANDBEDINGUNG GESETZT WIRD
C        INBR           - GITTER, AUS DEM DIE RANDBEDINGUNG KOMMT
C
C VERS:  27.07.95 (MM)  : ORIGINAL AUS BBACON ABGELEITET
C        27.10.97 (MM)    SETZEN DER WERTE AUS BUFFERN UBA,VBA...
C        28.02.03 (TB)  : SCALAR BOUNDARY TREATMENT IMPLEMENTED
C
C
C*MGLET*****************************************************************
C

      COMMON /COSTRLES/
     $                  VERS,   NRRUN,
     $                  DREAD,  DWRITE, DCONT,  DCUB,   
     $                  LBODYINSERT,
     $                  LGRIDNEW,    LGRIDREMOVE,
     $                  NPRNEU, FPRNEU, MTSTEP,
     $                  DT,     ITPRIN, IPINF,  ITINT,
     $                  IPRINT_WSS, IPRINT_WNS, ITFLUC, ITMIT,  
     $                  MPCORR, EPCORR, EPFAK,  MPCVOR, MPCNACH,
     $                  IVPINF, 
     $                  OMG,    LDIMLO, ISETRE,
     $                  LINPRN, IWRB,   LREC

      INTEGER 
     $                  NRRUN,
     $                  MTURB,  NPRNEU, MTSTEP,
     $                          ITPRIN, IPINF,  ITINT,
     $                  ITFLUC, ITMIT,  MPCORR,        
     $                                          ISETRE,
     $                          MSLIN,  IWRB          

      LOGICAL
     $                  DREAD,  DWRITE, DCONT,  DCUB,   
     $                  LBODYINSERT,
     $                  LGRIDNEW,       LGRIDREMOVE,
     $                  LDIMLO, LINPRN, LREC

      REAL
     $                  DT,     EPCORR,  OMG, FPRNEU

      CHARACTER (LEN=8)       VERS


      PARAMETER (NPPHYS_MAX=1)
      COMMON /COPHYSPAR/
     $                  MTURB,  TU_LEVEL,   RHO,    GMOL,   UGRID,
     $                  VREF,   EXPON,  UTAUX,
     $                  CIDUFR, UFRCON, UFRFREQ, DELTA,  XREF,
     $                  IPP,    JPP,    KPP,
     $                  XPER,   YPER,   ZPER, 
     $                  NXPER,  NYPER,  NZPER,
     $                  NPPHYS, PPHYS, XPPHYS

      INTEGER           IPP,  JPP,  KPP,   NPPHYS
     
      REAL              TU_LEVEL,  RHO,   GMOL,   UGRID,  VREF,  
     $                  EXPON, UFRCON, 
     $                  PPHYS(NPPHYS_MAX), XPPHYS(NPPHYS_MAX)

      CHARACTER (LEN=16)      CIDUFR

C
      CHARACTER (LEN=1)  ITYP
C
      REAL
     $    U(KMX,JMX,IMX),   V(KMX,JMX,IMX),   W(KMX,JMX,IMX),
     $    P(KMX,JMX,IMX),   G(KMX,JMX,IMX)
      REAL
     $     UBA(KMX,JMX,2),  VBA(KMX,JMX,2),  WBA(KMX,JMX,2),
     $     PBA(KMX,JMX,2),  GBA(KMX,JMX,2)
      REAL  DX(IMX), DY(JMX), DZ(KMX)
      REAL DDX(IMX),DDY(JMX),DDZ(KMX)
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IM1 = IMX-1
      IM2 = IMX-2
      IM3 = IMX-3
C
C
C                         FINE-GRID: I-INDEX= IMX - 1 (BACK)
               IIF = IMX - 1
C                         I-INDEX OF COARSE-GRID
C                                   (=PARENT, CALLED NEIGHBOUR)
C INDEX OF PARENT:               IC = IPOS - 1 + (IF-1)/2
C INDEX IN BUFFER:
                         IC = 1
C
C
C                                 FUER ZEITSCHRITT
C                                 NORMALKOMPONENTE WIRD EINEN INDEX
C                                 AUSSERHALB  GESETZT
      IF (ITYP.EQ.'T') THEN

         DO J=JSTART,JSTOP
         DO K=KSTART,KSTOP

               U(K,J,IIF-1) = UBA(K,J,1)
               U(K,J,IIF)   = UBA(K,J,1)
               V(K,J,IIF)   = VBA(K,J,1)
               W(K,J,IIF)   = WBA(K,J,1)
               P(K,J,IIF)   = PBA(K,J,1)
               G(K,J,IIF)   = GBA(K,J,1)


         ENDDO
         ENDDO


C
C                                 FUER DRUCKKORREKTUR
C                                 NUR NORMALKOMPONENTE MUSS GESETZT
C                                 WERDEN
       ELSEIF (ITYP.EQ.'P') THEN

C                   NORMAL-KOMPONENT MUST NOT BE SET
C                   BECAUSE PRESSURE IN THE BOUNDARY-CELLS
C                   IS NOT MODIFIED, SO NORMAL-KOMPONENTS MUST
C                   BE ABLE TO DEVELOP FREELY
               RETURN

         DO J=JSTART,JSTOP
         DO K=KSTART,KSTOP

                U(K,J,IIF-1) = UBA(K,J,1)

         ENDDO
         ENDDO
C
       ELSE
                CALL ERRR (501,' BBAPAR')
       ENDIF


      RETURN
      END

