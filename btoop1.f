










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
      SUBROUTINE BTOOP1 (KMX,JMX,IMX,X,
     $                    U,V,W,P,G,
     $                    ISTART,ISTOP,JSTART,JSTOP,ITYP
     $                   )     
C*STARLET***************************************************************
C        B T O O P 1   SETZEN DER RANDBEDINGUNGEN FUER DIE
C                         GESCHWINDIGKEITSFELDER AND SCALAR T.
C                         (LARGE-EDDY-SIMULATION)
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
C VERS:   8. 3.93 (MM)  : ORIGINAL
C                         AUS BTOPLE (BTOPL) ABGELEITET
C        28.02.03 (TB)  : SCALAR BOUNDARY TREATMENT IMPLEMENTED
C
C DEFINE-DIREKTIVEN     : QUD 
C
C*STARLET***************************************************************
C

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


      COMMON /KONSTA/  GREAT,SMALL,RINDEF,SMAONE,PRESET
      SAVE   /KONSTA/
C
      CHARACTER (LEN=1)  ITYP
C
      REAL
     $            U(KMX,JMX,IMX),   V(KMX,JMX,IMX),   W(KMX,JMX,IMX),
     $            P(KMX,JMX,IMX),   G(KMX,JMX,IMX)
      REAL        X( IMX )
C
C
      KM1 = KMX-1
      KM2 = KMX-2
      KM3 = KMX-3
      KM4 = KMX-4

C
C       do i=1,npphys
C         write(6,*)'PPHYS:',PPHYS(i),XPPHYS(i)
C       enddo
C
C                                 START- AND STOP-INDICES FOR PRESSURE
C                                 ON TOP
      IX_START = 1
      DO I=2,NPPHYS
         IF ( X(1)   .LE. XPPHYS(I) ) THEN
            IX_START = I-1
            GOTO 30
         ENDIF
      ENDDO
   30 CONTINUE

      IX_STOP = NPPHYS
      DO I=NPPHYS-1,1,-1
         IF ( X(IMX) .GE. XPPHYS(I) ) THEN
            IX_STOP  = I+1
            GOTO 31
         ENDIF
      ENDDO
   31 CONTINUE
C
C
C                                 **************************************
C                                 OPEN-WALL, GRADIENT = 0
C                                 **************************************
C
C
C
       IF(ITYP .EQ. 'P') RETURN

       DO I=ISTART,ISTOP
C                                  CALCULATION OF PINFTY

          IF ( X(I) .LE. XPPHYS(1     )) THEN
             PINFTY = PPHYS(1)
          ELSEIF ( X(I) .GE. XPPHYS(NPPHYS)) THEN
             PINFTY = PPHYS(NPPHYS)
          ELSE

             DO IX=IX_START,IX_STOP
                IF ((X(I)-XPPHYS(IX))*(X(I)-XPPHYS(IX+1)) 
     $               .LE. 0.0) THEN
                   PINFTY = PPHYS(IX) + 
     $                  (PPHYS(IX+1)-PPHYS(IX))*
     $                  (X(I)       -XPPHYS(IX))/
     $                  (XPPHYS(IX+1)-XPPHYS(IX))
                ENDIF
             ENDDO
          ENDIF
C
C---------------------------------------------------------------
C          write(6,*)'pinfty:',i,x(i),pinfty

          DO J=JSTART,JSTOP

               U(KM1,J,I) = U(KM2,J,I)
               V(KM1,J,I) = V(KM2,J,I)
C_030297               W(KM2,J,I) = W(KM3,J,I)
               W(KM1,J,I) = W(KM2,J,I)
               P(KM1,J,I) = PINFTY
               G(KM1,J,I) = G(KM2,J,I)

         ENDDO
         ENDDO


      RETURN
      END

