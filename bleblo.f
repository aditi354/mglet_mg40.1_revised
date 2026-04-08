










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
      SUBROUTINE BLEBLO (KMX,JMX,IMX,    
     $                      U,V,W,P,G,
     $                      ISTART,ISTOP,KSTART,KSTOP,
     $                      FREQ,AU,AV,AW,ANIV,TIMEPHYS,     
     $                      X,Y,Z,FLOWTYP,UFR,VFR,WFR
     $                   )     
C*STARLET***************************************************************
C        B F R B L O   SETZEN DER RANDBEDINGUNGEN FUER DIE STROEMUNGS-
C                      MANIPULATION DURCH EINEN JETSTRAHL
C
C*STARLET***************************************************************
C
C PARAM: 
C        U(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        V(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        W(KMX,JMX,IMX) + GESCHWINDIGKEITSFELD
C        T(KMX,JMX,IMX) + SCALAR T FIELD
C        P(KMX,JMX,IMX) - DRUCKFELD
C        G(KMX,JMX,IMX) - EFFEKTIVE DYNAMISCHE VISKOSITAET (= MUE)
C        AUB            - AMPLITUDE DES STRAHLS IN X-RICHTUNG
C        AVB            - AMPLITUDE DES STRAHLS IN Y-RICHTUNG
C        AWB            - AMPLITUDE DES STRAHLS IN Z-RICHTUNG
C        ANIV      - GIBT DAS NULLNIVEAU AN UM DAS DER COSINUS SCHWINGT
C        NUMJET         - NUMBER OF JETS
C        FREQB          - EINSTROEMFREQUENZ 
C        ITYP           - HOLLERITH-KONSTANTE:
C                         'P'  VOR DER DRUCKKORREKTUR
C                         'T'  VOR DEM ZEITSCHRITT
C
C VERS:  25.02.94 (RK)  : ORIGINAL
C VERS:  28.09.94 (AO)  : TOTALLY NEW VERSION AND ADAPTED TO MGLET
C VERS:  04.02.97 (AM)  : NEU INLET BOUNDARY CONDITION
C        28.02.03 (TB)  : SCALAR BOUNDARY TREATMENT IMPLEMENTED
C
C
C*STARLET***************************************************************
C
C
C
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

C
      REAL
     $            U(KMX,JMX,IMX),   V(KMX,JMX,IMX),   W(KMX,JMX,IMX),
     $            P(KMX,JMX,IMX),   G(KMX,JMX,IMX),   TIMEPHYS,
     $            X(IMX),           Y(JMX),           Z(KMX),
     $        UFR(KMX,JMX,2),     VFR(KMX,JMX,2),      WFR(KMX,JMX,2) 
      REAL RANF,RWORK
      EXTERNAL RANF

      CHARACTER (LEN=16) FLOWTYP
C
C

C                               QUADRATIC UPSTREAM IST NICHT DEFINIERT
        DATA NQUD /0/
        JSTART = JMX-2
        JSTOP  = JMX-2

C------------------RUNDER JET, PARABOLISCHES PROFIL----------------------------

      IF (FLOWTYP.EQ.'ROUNDJET') THEN


         VALUEK = 2.0 * 3.1415927 * FREQ * TIMEPHYS
         COSUS= COS(VALUEK)

         VALUEU = ANIV + AU * COSUS
         VALUEV = ANIV + AV * COSUS
         VALUEW = ANIV + AW * COSUS
         VALUEG = GMOL
  
         XA = X(ISTART)
         XB= X(ISTOP)
         ZA = Z(KSTART)
         ZB= Z(KSTOP)

         ZCENT = (ZB + ZA)/2.0
         ZR    = (ZB - ZA)/2.0
         XCENT = (XB + XA)/2.0
         XR    = (XB - XA)/2.0

         DO I = ISTART,ISTOP
           DO J = JSTART, JSTOP
              DO K = KSTART, KSTOP

                 R = SQRT((Z(K) - ZCENT)**2 + (X(I) - XCENT)**2)

                 R = R/ZR

                 V(K,J,I) = V(K,J,I) + MIN(VALUEV * (1.0 - R**2), 0.0)
              ENDDO
           ENDDO
        ENDDO
C        CALL WRITE3DDIAGY(KMX,JMX,IMX,V,20,JSTOP)
C        STOP 'BLEBLO'
      ENDIF

C------------------ECKIGER JET, DUCT-ROFIL----------------------------
      IF (FLOWTYP.EQ.'SQUAREJET') THEN

         XA = X(ISTART)
         XB= X(ISTOP)
         ZA = Z(KSTART)
         ZB= Z(KSTOP)

         DO I = ISTART,ISTOP
           DO J = JSTART,JSTOP
              DO K = KSTART, KSTOP
                 RWORK = 0.0 
                 DO L = 0,4
                    DO M = 0,4
                       RWORK = RWORK+(SIN((2*L+1)*3.14*X(I)/(XB-XA))*
     $                         SIN((2*M+1)*3.14*Z(K)/(ZB-ZA)))/
     $                         ((2*L+1)*(2*M+1)*((2*L+1)**2+(2*M+1)**2))
                    ENDDO
                 ENDDO

                 V(K,J,I) = AV/0.2138*RWORK
              ENDDO
           ENDDO
        ENDDO
C        OPEN(70)
C        DO I=ISTART,ISTOP
C           DO K=KSTART,KSTOP
C           WRITE(70,*)'T_PROF(Y):',Z(K),X(I),T(K,JSTOP,I),JSTOP,JSTART
C        ENDDO
C        WRITE(70,*)
C      ENDDO
C      CLOSE(70)
C        STOP 'HALT'      
      ENDIF

      RETURN
      END

