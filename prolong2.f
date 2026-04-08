










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
      SUBROUTINE PROLONG2 (JMX,IMX,IMXC,NBND,X,DX,DDX,XC,DXC,DDXC,
     $                   IN,OUT,HILF,
     $                   ISTA,ISTO,JSTA,JSTO,
     $                   IPOS,JPOS,ISTAG,IORDER)
C*MGLET*****************************************************************
C     P R O L O N G 2    PROLONGATES THE SECOND DIRECTION IN A 2D FIELD
C*MGLET*****************************************************************
C
C PARAM: 
C             JMX, IMX  - DIMENSIONS OF FINE GRID
C                 IMXC  - DIMENSION OF COARSE GRID
C       IN(KMX,JMX,IMX) - INPUT-2D-FIELD COARSE
C      OUT(KMX,JMX,IMX) + OUTPUT-2D-FIELD FINE
C        ISTA,ISTO,JSTA,JSTO - INDEX-REGION
C        IPOS,JPOS,KPOS - POSITIONS OF FINE GRID IN COARSE GRID, 
C                         AUF DER DER PUNKT
C                         MIT DEN INDIZES (3,3,3) ZU LIEGEN KOMMT
C               ISTAG   - INIDCATOR, IF VARIABLE IS STAGGERED IN FIRST DIM.
C                         0: NON-STAGGERED; 1: STAGGERED
C           IORDER      - INTERPOLATION ORDER
C
C VERS:  16.09.96 (MM)  : ORIGINAL AUS CONTOPAR ABGELEITET
C
C DEFINE-DIREKTIVEN     : QUD 
C
C*STARLET***************************************************************

C
      CHARACTER (LEN=1)  ITYP
C
      REAL
     $    IN(JMX,IMXC),  OUT(JMX,IMX),  HILF(JMX,IMX),
     $       X(JMX),           XC(JMX),
     $      DX(JMX),          DXC(JMX),
     $     DDX(JMX),         DDXC(JMX)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PI2 = ATAN(1.0)*8.0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  VARIABLE NON-STAGGERED IN SECOND DIR.
      IF (ISTAG .EQ. 0) THEN
C                                Interpolation of first order
         IF (IORDER .EQ. 1) THEN

            DO IF=ISTA,ISTO,2
               IC = IPOS - 1 + (IF-1)/2

               DO JF = JSTA,JSTO

                  OUT(JF,IF  ) = IN(JF,IC)
                  OUT(JF,IF+1) = IN(JF,IC)

               ENDDO
            ENDDO
            
         ENDIF
C                                Interpolation of second order
         IF (IORDER .EQ. 2) THEN
            DO IF=ISTA,ISTO,2
               IC = IPOS - 1 + (IF-1)/2

               DO JF = JSTA,JSTO
                  
                  OUT(JF,IF  ) = IN(JF,IC  )*0.75  +
     +                           IN(JF,IC-1)*0.25

                  OUT(JF,IF+1) = IN(JF,IC  )*0.75  +
     +                           IN(JF,IC+1)*0.25
                   
                ENDDO
            ENDDO
         ENDIF
C                                Interpolation of third order, 
C                                Lagrange-Polynomial (Carnahan et al. 1969)
         IF (IORDER .EQ. 3) THEN
            DO IF=ISTA,ISTO,2
               IC = IPOS - 1 + (IF-1)/2

               DO JF = JSTA,JSTO

                  D1X0 =    0.75*DXC(IC-1)
                  D1X1 = -  0.25*DXC(IC-1)
                  D1X2 = - (0.25*DXC(IC-1) + DXC(IC))

                  D2X0 =   (0.25*DXC(IC) + DXC(IC-1))
                  D2X1 =    0.25*DXC(IC)
                  D2X2 =  - 0.75*DXC(IC)

                  R1_0 = D1X1*D1X2 / (DXC(IC-1)*(DXC(IC-1)+DXC(IC)))
                  R1_1 = D1X0*D1X2 / (DXC(IC-1)*( - DXC(IC)))
                  R1_2 = D1X0*D1X1 / (DXC(IC  )*(DXC(IC-1)+DXC(IC)))

                  R2_0 = D2X1*D2X2 / (DXC(IC-1)*(DXC(IC-1)+DXC(IC)))
                  R2_1 = D2X0*D2X2 / (DXC(IC-1)*( - DXC(IC)))
                  R2_2 = D2X0*D2X1 / (DXC(IC  )*(DXC(IC-1)+DXC(IC)))

                  OUT(JF,IF  ) = IN(JF,IC-1)*R1_0  +
     +                           IN(JF,IC  )*R1_1  +
     +                           IN(JF,IC+1)*R1_2

                  OUT(JF,IF+1) = IN(JF,IC-1)*R2_0  +
     +                           IN(JF,IC  )*R2_1  +
     +                           IN(JF,IC+1)*R2_2

                ENDDO
            ENDDO
         ENDIF
C                                Interpolation of third order, 
C                                Lagrange-Polynomial (Carnahan et al. 1969)
C                                PLUS CONTINUITY!!!!!!!!!!!!!!
         IF (IORDER .EQ.13) THEN
            DO IF=ISTA,ISTO,2
               IC = IPOS - 1 + (IF-1)/2

               DO JF = JSTA,JSTO

                  D1X0 =    0.75*DXC(IC-1)
                  D1X1 = -  0.25*DXC(IC-1)
                  D1X2 = - (0.25*DXC(IC-1) + DXC(IC))

                  D2X0 =   (0.25*DXC(IC) + DXC(IC-1))
                  D2X1 =    0.25*DXC(IC)
                  D2X2 =  - 0.75*DXC(IC)

                  R1_0 = D1X1*D1X2 / (DXC(IC-1)*(DXC(IC-1)+DXC(IC)))
                  R1_1 = D1X0*D1X2 / (DXC(IC-1)*( - DXC(IC)))
                  R1_2 = D1X0*D1X1 / (DXC(IC  )*(DXC(IC-1)+DXC(IC)))

                  R2_0 = D2X1*D2X2 / (DXC(IC-1)*(DXC(IC-1)+DXC(IC)))
                  R2_1 = D2X0*D2X2 / (DXC(IC-1)*( - DXC(IC)))
                  R2_2 = D2X0*D2X1 / (DXC(IC  )*(DXC(IC-1)+DXC(IC)))

                  OUT(JF,IF  ) = IN(JF,IC-1)*R1_0  +
     +                           IN(JF,IC  )*R1_1  +
     +                           IN(JF,IC+1)*R1_2

                  OUT(JF,IF+1) = IN(JF,IC-1)*R2_0  +
     +                           IN(JF,IC  )*R2_1  +
     +                           IN(JF,IC+1)*R2_2

                  DIFF = IN(JF,IC) - 0.5*(OUT(JF,IF  )+OUT(JF,IF+1))
                  OUT(JF,IF  ) = OUT(JF,IF  ) + DIFF
                  OUT(JF,IF+1) = OUT(JF,IF+1) + DIFF
                ENDDO
            ENDDO
         ENDIF
C                                Interpolation of fourth order,
C                                Lagrange-Polynomial (Carnahan et al. 1969)
         IF (IORDER .EQ. 4) THEN
            DO IF=MAX(3,ISTA),MIN(IMX-NBND,ISTO),2
               IC = IPOS - 1 + (IF-1)/2

               DO JF = JSTA,JSTO

                  D1X0 =    0.75*DXC(IC-1) + DXC(IC-2)
                  D1X1 =    0.75*DXC(IC-1)
                  D1X2 = -  0.25*DXC(IC-1)
                  D1X3 = - (0.25*DXC(IC-1) + DXC(IC))

                  D2X0 =   (0.25*DXC(IC) + DXC(IC-1))
                  D2X1 =    0.25*DXC(IC)
                  D2X2 =  - 0.75*DXC(IC)
                  D2X3 =  - 0.75*DXC(IC) - DXC(IC+1)

                  R1_0 = D1X1*D1X2*D1X3 / 
     $                  (-1.0*(DXC(IC-2))*
     $                        (DXC(IC-2)+DXC(IC-1))*
     $                        (DXC(IC-2)+DXC(IC-1)+DXC(IC)))

                  R1_1 = D1X0*D1X2*D1X3 / 
     $                  (     (DXC(IC-2))*
     $                        (DXC(IC-1))*
     $                        (DXC(IC-1)+DXC(IC)))

                  R1_2 = D1X0*D1X1*D1X3 / 
     $                  (-1.0*(DXC(IC-2)+DXC(IC-1))*
     $                        (DXC(IC-1))*
     $                        (DXC(IC)))

                  R1_3 = D1X0*D1X1*D1X2 / 
     $                  (     (DXC(IC-2)+DXC(IC-1)+DXC(IC))*
     $                        (DXC(IC-1)+DXC(IC  ))*
     $                        (DXC(IC  )))

                  R2_0 = D2X1*D2X2*D2X3 / 
     $                  (-1.0*(DXC(IC-1))*
     $                        (DXC(IC-1)+DXC(IC  ))*
     $                        (DXC(IC-1)+DXC(IC  )+DXC(IC+1)))

                  R2_1 = D2X0*D2X2*D2X3 / 
     $                  (     (DXC(IC-1))*
     $                        (DXC(IC  ))*
     $                        (DXC(IC  )+DXC(IC+1)))

                  R2_2 = D2X0*D2X1*D2X3 / 
     $                  (-1.0*(DXC(IC-1)+DXC(IC  ))*
     $                        (DXC(IC  ))*
     $                        (DXC(IC+1)))

                  R2_3 = D2X0*D2X1*D2X2 / 
     $                  (     (DXC(IC-1)+DXC(IC  )+DXC(IC+1))*
     $                        (DXC(IC  )+DXC(IC+1))*
     $                        (DXC(IC+1)))

                  OUT(JF,IF  ) = IN(JF,IC-2)*R1_0  +
     +                           IN(JF,IC-1)*R1_1  +
     +                           IN(JF,IC  )*R1_2  +
     +                           IN(JF,IC+1)*R1_3

                  OUT(JF,IF+1) = IN(JF,IC-1)*R2_0  +
     +                           IN(JF,IC  )*R2_1  +
     +                           IN(JF,IC+1)*R2_2  +
     +                           IN(JF,IC+2)*R2_3

                ENDDO
            ENDDO
         ENDIF
C                                  FOURIER-INTERPOLATION
         IF (IORDER .EQ.99) THEN
C                                  INITIALISIEREN DES OUT-FELDES
            DO I=1,IMX
               DO J=1,JMX
                  HILF(J,I)=0.0
               ENDDO
               DO J=1,JMX
                  OUT(J,I)=0.0
               ENDDO
            ENDDO
C                                  BELEGUNG DES 2D-FELDES
            DO JC = NBND+1,JMX-NBND
               DO IC = NBND+1,IMXC-NBND
                  
                  HILF(JC,IC) = IN(JC,IC)
                  
               ENDDO
            ENDDO
C                                  FFT IN ZWEITEM INDEX


C                                  SHIFTEN DER WELLENZAHLEN
            DO IK = 1,(IMXC-(2*NBND))/2+NBND+1
               IC = 2*IK+NBND

               RCOS=COS( PI2 * FLOAT(IK-1) * 0.5 * DX(NBND+1) 
     $              / (X(IMX-NBND)-X(NBND)))
               RSIN=SIN( PI2 * FLOAT(IK-1) * 0.5 * DX(NBND+1) 
     $              / (X(IMX-NBND)-X(NBND)))

               DO JC = NBND+1,JMX-NBND
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC HIER: REALTEIL DER WELLENZAHL JK
                  OUT(JC,IC  )=HILF(JC,IC  ) * RCOS * 2.0 - 
     $                         HILF(JC,IC-1) * RSIN * 2.0


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC HIER: IMAGINAERTEIL DER WELLENZAHL JK
                  OUT(JC,IC-1)=HILF(JC,IC-1) * RCOS * 2.0 + 
     $                         HILF(JC,IC  ) * RSIN * 2.0

               ENDDO
            ENDDO
C                                         WELLENZAHLEN > (IMXC-(2*NBND))
C                                         WERDEN MIT NULL BELEGT
            DO JC = NBND+1,JMX-NBND
               DO IC = (IMXC-(2*NBND))+NBND+1,IMX
               
                  OUT(JC,IC)=0.0

               ENDDO
            ENDDO
               

C                                  RUECKTRANSFORMATION

C                                  PERIODIC BOUNDARY CONDITIONS
            DO J=1,JMX
               OUT( J ,    NBND  ) = OUT( J ,IMX-NBND  )
               OUT( J ,IMX-NBND+1) = OUT( J ,    NBND+1)
            ENDDO

         ENDIF

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  VARIABLE STAGGERED IN SECOND DIR.
      ELSE

C                                Interpolation of first order
         IF (IORDER .EQ. 1) THEN
           DO IF=ISTA,ISTO,2
              IC = IPOS - 1 + (IF-1)/2

              DO JF = JSTA,JSTO

                 OUT(JF,IF  ) = 0.5*( IN(JF,IC) + IN(JF,IC-1))
                 OUT(JF,IF+1) = IN(JF,IC)

              ENDDO
           ENDDO


         ENDIF
C                                Interpolation of second order
         IF (IORDER .EQ. 2) THEN
           DO IF=ISTA,ISTO,2
              IC = IPOS - 1 + (IF-1)/2

              DO JF = JSTA,JSTO
                  
                  OUT(JF,IF  ) = 1./DDXC(IC) *
     +                          (IN(JF,IC  )*DDX(IF  ) +
     +                           IN(JF,IC-1)*DDX(IF+1))

                  OUT(JF,IF+1) = IN(JF,IC)
                   
                ENDDO
            ENDDO
         ENDIF
C                                Interpolation of third order, 
C                                Lagrange-Polynomial (Carnahan et al. 1969)
         IF (IORDER .EQ. 3) THEN
           DO IF=ISTA,ISTO,2
              IC = IPOS - 1 + (IF-1)/2

              DO JF = JSTA,JSTO

                  D1X0 =    0.25*(DXC(IC-1) + DDXC(IC))
                  D1X1 = -  0.25*(DXC(IC  ) + DDXC(IC))
                  D1X2 = -  0.25*(DXC(IC  ) + DDXC(IC)) - DDXC(IC+1)

                  R1_0 = D1X1*D1X2 / (DDXC(IC)*(DDXC(IC)+DDXC(IC+1)))
                  R1_1 = D1X0*D1X2 / (DDXC(IC)*( - DDXC(IC+1)))
                  R1_2 = D1X0*D1X1 / (DDXC(IC+1)*(DDXC(IC+1)+DDXC(IC)))

                  OUT(JF,IF  ) = IN(JF,IC-1)*R1_0  +
     +                           IN(JF,IC  )*R1_1  +
     +                           IN(JF,IC+1)*R1_2

                  OUT(JF,IF+1) = IN(JF,IC)

                ENDDO
            ENDDO
         ENDIF

      ENDIF


      RETURN
      END

