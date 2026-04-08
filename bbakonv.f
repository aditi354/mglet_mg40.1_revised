










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
      SUBROUTINE BBAKONV (KMX,JMX,IMX,
     $                    U,V,W,P,G,
     $                    JSTART,JSTOP,KSTART,KSTOP,ITYP,
     $                    IDIM3D,IDIM2D,IDIM1D,
     $                    X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,
     $                    UBA,VBA,WBA
     $                   )      
C*STARLET***************************************************************
C        B B A K O N V  SETZEN DER RANDBEDINGUNGEN FUER DIE
C                         GESCHWINDIGKEITSFELDER AND SCALAR T.
C                         (LARGE-EDDY-SIMULATIO/DNS)
C                    KONVECTIVE BOUNDARY CONDITION FOR THE OUTLET
C  VERS: 10.09.1995 ORIGINAL (AO)
C*STARLET***************************************************************
C
C PARAM: 
C        KMX, JMX, IMX  - DIMENSIONEN DES GITTERS 
C        U(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        V(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        W(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        T(KMX,JMX,IMX) + SCALAR T FIELD
C        P(KMX,JMX,IMX) - DRUCKFELD
C        G(KMX,JMX,IMX) - EFFEKTIVE DYNAMISCHE VISKOSITAET (= MUE)
C        ITYP           - HOLLERITH-KONSTANTE:
C                         'P'  VOR DER DRUCKKORREKTUR
C                         'T'  VOR DEM ZEITSCHRITT
C
C        28.02.03 (TB)  : SCALAR BOUNDARY TREATMENT IMPLEMENTED
C
C DEFINE-DIREKTIVEN     : QUD 
C
C*STARLET***************************************************************
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

C
      CHARACTER (LEN=1)  ITYP
C
      REAL
     $    U(KMX,JMX,IMX),   V(KMX,JMX,IMX),   W(KMX,JMX,IMX),
     $    P(KMX,JMX,IMX),   G(KMX,JMX,IMX),
     $    UBA(KMX,JMX,3), VBA(KMX,JMX,3), WBA(KMX,JMX,3),   
     $    DX(IMX), DY(JMX),DZ(KMX),
     $    DDX(IMX),DDY(JMX),DDZ(KMX) 
C
C
      IM1 = IMX-1
      IM2 = IMX-2
      IM3 = IMX-3
C
C
C                                 **************************************
C                                 KONVECTIVE BOUNDARY CONDITION FROM MOIN (AO)
C                                 **************************************
C
C
C
C                                 VOR DER DRUCKKORREKTUR
C                                 ----------------------
      IF(ITYP .EQ. 'P') THEN

         DO J=JSTART+1,JSTOP
         DO K=KSTART+1,KSTOP

C              P(K,J,IM1) = 0.0
C              P(K,J,IM1) = P(K,J,IM2)

         ENDDO
         ENDDO

         RETURN

      ENDIF
         UIN = 1.0

         DO J=JSTART+1,JSTOP
         DO K=KSTART+1,KSTOP
         
C           UIN = UBA(K,J,2)
C           UIN = MAX(0.0,UIN)

            U(K,J,IM1) = UBA(K,J,1) - (DT / DDX(IM1)) * UIN * 
     $                 ( UBA(K,J,1) - UBA(K,J,2) )
            V(K,J,IM1) = VBA(K,J,1) - (DT / DX(IM2))  * UIN * 
     $                 ( VBA(K,J,1) - VBA(K,J,2) )
            W(K,J,IM1) = WBA(K,J,1) - (DT / DX(IM2))  * UIN * 
     $                 ( WBA(K,J,1) - WBA(K,J,2) )
C           G(K,J,IM1) = GBA(K,J,1) - (DT / DX(IM2))  * UBA(K,J,3) *
C    $                 ( GBA(K,J,1) - GBA(K,J,2) )
            G(K,J,IM1) = G(K,J,IM2)
            P(K,J,IM1) = 0.0

C           IF THE MASS FLOW IS NEGATIVE THEN GIVE
C           INFORMATION ABOUT IN OUTFILE


         ENDDO
         ENDDO


      RETURN
      END

