










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
      SUBROUTINE  SWCLE3D(KK,JJ,II,KMX,JMX,IMX,DX,DY,DZ,DDX,DDY,DDZ,
     $                    U,V,W,P,G,B,
     $                    WCU,WCV,WCW,GMOL,RHO,UGRID,
     $                    MTURB,ZTOT,GRADPX,
     $                    IC1,IC2,JC1,JC2,KC1,KC2,
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,
     $                    FDUK,FDVK,FDUJ,FDWJ,FDVI,FDWI,
     $                    RDDX,RDDY,RDDZ)
C*STARLET***************************************************************
C        S W C L E 3 D     WALL-CORRECTION (LAM-TUR)
C*STARLET***************************************************************
C
C PARAM: GMOL             - MOLEKULARE DYNAM. VISKOSITAET ( = MUE)
C        CAPPA            - VON KARMAN KONSTANTE
C        CMUE*RHO             - KONSTANTE (G = CMUE*RHO*TK**2/TE)
C        ECONST           - RAUHIGKEITSPARAMETER FUER WANDGESETZ
C        MTURB            - 0 = LAMINAR, 1 = TURBULENT
C        IC1,JC1,KC1    - LINKE  RAENDER DER BOUNDING BOX
C        IC2,JC2,KC2    - RECHTE RAENDER DER BOUNDING BOX
C
C WARNG: IN WANDNORMALRICHTUNG WIRD D = DD VORAUSGESETZT
C        ALLE WC''S  HABEN NEGATIVES VORZEICHEN
C
C DEFINE-DIREKTIVEN     : FRPER, RIPER, TONOS
C
C VERS:  12.05.86 (HW)  : SWCLEC AUS SWCLE1 ABGELEITET. SWCLEC IST FUER
C                         EINBAUTEN, DIE AUS SAEULEN AUFGEBAUT SIND,
C                         GEEIGNET. DIE WANDSCHUBSPANNUNG WIRD AUS DEM
C                         1/7-POTENZGESETZ BERECHNET.
C        02.04.92 (MM)  : WANDERKENNUNG JETZT MIT G-FELD               
C         4. 4.92 (MM)  : SCHLEIFEN ZUR WANDKORREKTUR AM KOERPER
C                         LAUFEN NUR INNERHALB DER
C                         BOUNDING BOX (IC1,IC2,JC1,JC2,KC1,KC2)
C        14.08.1996 (AO): RDDX(II),RDDY(JJ) UND RDDZ(KK) INTRODUCED;
C                        IF DNS THE WALL-SHEAR-STRESS WILL BE CALCULATED
C                         DIRECTLY WITHOUT CALLING THE FUNCTION TAUWIN
C        02.07.03 (FS)  : WANDSCHUPSPANNUNG MIT 4TER ORDNUNG ENIGEFÜHRT
C*STARLET***************************************************************

      IMPLICIT NONE
      INTEGER KMX,JMX,IMX,KK, JJ, II,
     $        KC1,JC1,IC1,KC2,JC2,IC2,
     $        KM1,JM1,IM1,KM2,JM2,IM2,
     $        KM3,JM3,IM3,K,  J,  I,
     $        NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,MTURB

      COMMON /CONLES/  CONV2S,CMUE,CAPPA,ECONST
      SAVE   /CONLES/

      COMMON /CONPOT/  CWA,CWB,CPO1,CPO2,CPO3,CPO4,CPO5,CPO6,CPO7,CPO8,
     $                 CPO9,CPO10,CPO11,CPO12
      SAVE   /CONPOT/

      REAL    U  (KK,JJ,II),V  (KK,JJ,II),W  (KK,JJ,II),P(KK,JJ,II),
     $        G  (KK,JJ,II),B  (KK,JJ,II), 
     $        WCU(KK,JJ,II),WCV(KK,JJ,II),WCW(KK,JJ,II)


      REAL    DX(II),      DY(JJ),      DZ(KK),
     $        DDX(II),     DDY(JJ),     DDZ(KK),
     $        RDDX(II),    RDDY(JJ),    RDDZ(KK)

      REAL RHO,GMOL,UGRID,ZTOT,GRADPX,CONV2S,CMUE,CAPPA,ECONST,CWA,CWB,
     $     CPO1,CPO3,CPO2,CPO4,CPO5,CPO6,CPO7,CPO8,CPO9,
     $     CPO10,CPO11,CPO12,TAUWIN

      REAL UXY,UXZ,VXY,VYZ,WXY,WXZ,WYZ,
     $     WCUY,WCUZ,WCVX,WCVZ,WCWX,WCWY,
     $     FXNE,FXPO,FYNE,FYPO,FZNE,FZPO

      REAL FDUK(1),FDVK(1)
      REAL FDUJ(1),FDWJ(1)
      REAL FDVI(1),FDWI(1)
C         WRITE(6,*)'IN SWCLE3D1',IMX,JMX,KMX,II,JJ,KK
         KM1  = KMX-1
         KM2  = KMX-2
         JM1  = JMX-1
         JM2  = JMX-2
         JM3  = JMX-3
         IM1  = IMX-1
         IM2  = IMX-2
         IM3  = IMX-3

C                                 WALL-CORRECTION AUF NULL
         DO I = 1,IMX
            DO J=1,JMX
               DO K=1,KMX
                  WCU(K,J,I)  = 0.0
                  WCV(K,J,I)  = 0.0
                  WCW(K,J,I)  = 0.0
               ENDDO
            ENDDO
         ENDDO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  BOTTOM NOSLIP
C
      IF ( NBOT .EQ. 5 ) THEN


         K = 3

         DO I=2,IM2
            DO J=2,JM1

C                                 **************************************
C                                 UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C                                 **************************************
         UXZ     = DX(I)*DDZ(K)
         UXY     = DX(I)*DDY(J)

C                                 WANDFUNKTIONEN

         WCUZ  = UXY * 
     &   2.0*GMOL*ABS(U(K,J,I)+UGRID)*RDDZ(K)
     &   *SIGN(1.0,(U(K,J,I)+UGRID))
C
         WCU(K,J,I) = - WCUZ
C
C                                 **************************************
C                                 VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
C                                 **************************************
C
         VYZ   = DY(J)*DDZ(K)
         VXY   = DY(J)*DDX(I)
C
C                                 WANDFUNKTIONEN
C
         WCVZ  = VXY * 
     &   2.0*GMOL*ABS(V(K,J,I))*RDDZ(K)*SIGN(1.0,V(K,J,I))

         WCV(K,J,I) = - WCVZ

         ENDDO
      ENDDO
C      STOP 'DEBUG'
      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  TOP NOSLIP

      IF ( NTOP .EQ. 5 ) THEN


         K = KMX-2

         DO I=2,IM1
            DO J=2,JM1

C                                 **************************************
C                                 UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C                                 **************************************
         UXZ     = DX(I)*DDZ(K)
         UXY     = DX(I)*DDY(J)
C
C                                 WANDFUNKTIONEN
C
         WCUZ  = UXY * 
     &   2.0*GMOL*ABS(U(K,J,I)+UGRID)*RDDZ(K)*SIGN(1.0,U(K,J,I)+UGRID)
C
         WCU(K,J,I) = - WCUZ
C
C                                 **************************************
C                                 VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
C                                 **************************************
C
         VYZ   = DY(J)*DDZ(K)
         VXY   = DY(J)*DDX(I)
C
C                                 WANDFUNKTIONEN
C
         WCVZ  = VXY * 
     &   2.0*GMOL*ABS(V(K,J,I))*RDDZ(K)*SIGN(1.0,V(K,J,I))

         WCV(K,J,I) = - WCVZ

         ENDDO
      ENDDO

      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  LEFT NOSLIP

      IF ( NLFT .EQ. 5 ) THEN


         J = JMX-2

         DO I=2,IM1
            DO K=2,KM1

C                                 **************************************
C                                 UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C                                 **************************************
         UXZ     = DX(I)*DDZ(K)
         UXY     = DX(I)*DDY(J)
C
C                                 WANDFUNKTIONEN
C
         WCUY  = UXZ *
     &   2.0*GMOL*ABS(U(K,J,I)+UGRID)*RDDY(J)*SIGN(1.0,U(K,J,I)+UGRID)

         WCU(K,J,I) = - WCUY

C                                 **************************************
C                                 WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
C                                 **************************************

         WXZ   = DZ(K)*DDX(I)
         WYZ   = DZ(K)*DDY(J)

C                                 WANDFUNKTIONEN

         WCWY  = WXZ *
     &   2.0*GMOL*ABS(W(K,J,I))*RDDY(J)*SIGN(1.0,W(K,J,I))

         WCW(K,J,I) = - WCWY

         ENDDO
      ENDDO

      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  RIGHT NOSLIP

      IF ( NRGT .EQ. 5 ) THEN


         J = 3

         DO I=2,IM1
            DO K=2,KM1

C                                 **************************************
C                                 UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C                                 **************************************
         UXZ     = DX(I)*DDZ(K)
         UXY     = DX(I)*DDY(J)
C
C                                 WANDFUNKTIONEN
C
         WCUY  = UXZ * 
     &   2.0*GMOL*ABS(U(K,J,I)+UGRID)*RDDY(J)*SIGN(1.0,U(K,J,I)+UGRID)
C
         WCU(K,J,I) = - WCUY
C
C                                 **************************************
C                                 WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
C                                 **************************************
C
         WXZ   = DZ(K)*DDX(I)
         WYZ   = DZ(K)*DDY(J)
C
C                                 WANDFUNKTIONEN
C
         WCWY  = WXZ * 
     &   2.0*GMOL*ABS(W(K,J,I))*RDDY(J)*SIGN(1.0,W(K,J,I))

         WCW(K,J,I) = - WCWY

         ENDDO
      ENDDO

      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  FRONT NOSLIP

      IF ( NFRO .EQ. 5 ) THEN

         I = 3

         DO J=2,JM1
            DO  K=2,KM1
C                                 **************************************
C                                 VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
C                                 **************************************
               VYZ     = DY(J)*DDZ(K)
C
C                                 WANDFUNKTIONEN
C
         WCVX  = VYZ * 
     &   2.0*GMOL*V(K,J,I)*RDDX(I)
C
         WCV(K,J,I) = - WCVX + WCV(K,J,I)
C
C                                 **************************************
C                                 WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
C                                 **************************************
C
         WYZ   = DZ(K)*DDY(J)
C
C                                 WANDFUNKTIONEN
C
         WCWX  = WYZ * 
     &   2.0*GMOL*W(K,J,I)*RDDX(I)

         WCW(K,J,I) = - WCWX + WCW(K,J,I)

         ENDDO
      ENDDO

      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  BACK NOSLIP

      IF ( NBAC .EQ. 5 ) THEN

            I = IM2

         DO J=2,JM1
            DO  K=2,KM1
C                                 **************************************
C                                 VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
C                                 **************************************
               VYZ     = DY(J)*DDZ(K)
C
C                                 WANDFUNKTIONEN
C
         WCVX  = VYZ *
     &   2.0*GMOL*V(K,J,I)*RDDX(I)
C
         WCV(K,J,I) = - WCVX + WCV(K,J,I)
C
C                                 **************************************
C                                 WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
C                                 **************************************
C
         WYZ   = DZ(K)*DDY(J)
C
C                                 WANDFUNKTIONEN
C
         WCWX  = WYZ *
     &   2.0*GMOL*W(K,J,I)*RDDX(I)

         WCW(K,J,I) = - WCWX + WCW(K,J,I)

         ENDDO
      ENDDO

      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  KOERPER NOSLIP

      IF ( NCUB .EQ. 5 ) THEN

         DO I = 2,IM1

          IF ((I.GE.IC1).AND.(I.LE.IC2)) THEN

          DO J=JC1,JC2
          DO K=KC1,KC2
C
C                                 **************************************
C                                 UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C                                 **************************************
C
          UXZ     = DX(I)*DDZ(K)
          UXY     = DX(I)*DDY(J)
C
C                                 SCHALTFAKTOREN BERECHNUNG
C
          FYPO  = (SIGN(0.25,G(K,J,I+1)*G(K,J+1,I+1))-0.25)
     $          + (SIGN(0.25,G(K,J,I  )*G(K,J+1,I  ))-0.25)
          FYNE  = (SIGN(0.25,G(K,J,I+1)*G(K,J-1,I+1))-0.25)
     $          + (SIGN(0.25,G(K,J,I  )*G(K,J-1,I  ))-0.25)
          FZPO  = (SIGN(0.25,G(K,J,I+1)*G(K+1,J,I+1))-0.25)
     $          + (SIGN(0.25,G(K,J,I  )*G(K+1,J,I  ))-0.25)
          FZNE  = (SIGN(0.25,G(K,J,I+1)*G(K-1,J,I+1))-0.25)
     $          + (SIGN(0.25,G(K,J,I  )*G(K-1,J,I  ))-0.25)
C
C                                 WANDFUNKTIONEN
C
          WCUY  = UXZ * 
     &    2.0*GMOL*ABS(U(K,J,I)+UGRID)*RDDY(J)*SIGN(1.0,U(K,J,I)+UGRID)
          WCUZ  = UXY *
     &    2.0*GMOL*ABS(U(K,J,I)+UGRID)*RDDZ(K)*SIGN(1.0,U(K,J,I)+UGRID)
C
          WCU(K,J,I) = FYPO*WCUY + FYNE*WCUY + FZPO*WCUZ + FZNE*WCUZ

          ENDDO
          ENDDO

          DO J=JC1,JC2
          DO K=KC1,KC2
C                                 **************************************
C                                 VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
C                                 **************************************
C
          VYZ   = DY(J)*DDZ(K)
          VXY   = DY(J)*DDX(I)
C
C                                 SCHALTFAKTOREN BERECHNUNG
C
          FXPO  = (SIGN(0.25,G(K,J+1,I)*G(K,J+1,I+1))-0.25)
     $          + (SIGN(0.25,G(K,J,I  )*G(K,J,I+1  ))-0.25)
          FXNE  = (SIGN(0.25,G(K,J+1,I)*G(K,J+1,I-1))-0.25)
     $          + (SIGN(0.25,G(K,J,I  )*G(K,J,I-1  ))-0.25)
          FZPO  = (SIGN(0.25,G(K,J+1,I)*G(K+1,J+1,I))-0.25)
     $          + (SIGN(0.25,G(K,J,I  )*G(K+1,J,I  ))-0.25)
          FZNE  = (SIGN(0.25,G(K,J+1,I)*G(K-1,J+1,I))-0.25)
     $          + (SIGN(0.25,G(K,J,I  )*G(K-1,J,I  ))-0.25)
C
C                                 WANDFUNKTIONEN
C
          WCVX  = VYZ *
     &    2.0*GMOL*ABS(V(K,J,I))*RDDX(I)*SIGN(1.0,V(K,J,I))
          WCVZ  = VXY *
     &    2.0*GMOL*ABS(V(K,J,I))*RDDZ(K)*SIGN(1.0,V(K,J,I))

          WCV(K,J,I) = FXPO*WCVX + FXNE*WCVX + FZPO*WCVZ + FZNE*WCVZ

          ENDDO
          ENDDO

          DO J=JC1,JC2
          DO K=KC1,KC2
C                                 **************************************
C                                 WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
C                                 **************************************
C
          WXZ   = DZ(K)*DDX(I)
          WYZ   = DZ(K)*DDY(J)
C
C                                 SCHALTFAKTOREN BERECHNUNG
C
          FXPO  = (SIGN(0.25,G(K+1,J,I)*G(K+1,J,I+1))-0.25)
     $          + (SIGN(0.25,G(K,J,I  )*G(K,J,I+1  ))-0.25)
          FXNE  = (SIGN(0.25,G(K+1,J,I)*G(K+1,J,I-1))-0.25)
     $          + (SIGN(0.25,G(K,J,I  )*G(K,J,I-1  ))-0.25)
          FYPO  = (SIGN(0.25,G(K+1,J,I)*G(K+1,J+1,I))-0.25)
     $          + (SIGN(0.25,G(K,J,I  )*G(K,J+1,I  ))-0.25)
          FYNE  = (SIGN(0.25,G(K+1,J,I)*G(K+1,J-1,I))-0.25)
     $          + (SIGN(0.25,G(K,J,I  )*G(K,J-1,I  ))-0.25)
C
C                                 WANDFUKTIONEN
C
          WCWX  = WYZ *
     &    2.0*GMOL*ABS(W(K,J,I))*RDDX(I)*SIGN(1.0,W(K,J,I))
          WCWY  = WXZ *
     &    2.0*GMOL*ABS(W(K,J,I))*RDDY(J)*SIGN(1.0,W(K,J,I))

          WCW(K,J,I) = FXPO*WCWX + FXNE*WCWX + FYPO*WCWY + FYNE*WCWY

          ENDDO
          ENDDO
          ENDIF
      ENDDO
      ENDIF


      RETURN
      END
      FUNCTION TAUWIN   (UQUER,DDS)
C*STARLET***************************************************************
C        T A U W P      BERECHNUNG DER WANDSCHUBSPANNUNG MITTELS DES
C                       1/7 POTENZGESETZES. DAS VORZEICHEN DER WAND-
C                       SCHUBSPG. ENTSPRICHT DEM VON UQUER !
C*STARLET***************************************************************
C
C PARAM: UQUER          - GESCHWINDIGKEITSKOMPONENTE AM WANDNAECHSTEN
C                         GITTERPUNKT (UQUER IST EIN MITTELWERT UEBER
C                         DIE MASCHENFLAECHE SENKRECHT ZU UQUER)
C        DDS            - KANTENLAENGE DER WANDNAECHSTEN ZELLE
C                         (SENKRECHT ZU UQUER UND ZUR WAND )
C
C VERS:  10.09.85 (HW)  : ORIGINAL
C                         IDENTISCH MIT TAUWP, NUR FUER INLINING
C                         AUF VPP IN SWCLE1 EINGEFUEGT
C
C DEFINE-DIREKTIVEN     : NATWBC
C
C*STARLET***************************************************************
C

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


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


      COMMON /CONPOT/  CWA,CWB,CPO1,CPO2,CPO3,CPO4,CPO5,CPO6,CPO7,CPO8,
     $                 CPO9,CPO10,CPO11,CPO12
      SAVE   /CONPOT/
C
C
      VZ     = SIGN(1.0,UQUER)
      UQUERN = ABS(UQUER)
C
C
C                                 HIER: UQUER LIEGT IN DER VISKOSEN
C                                 UNTERSCHICHT
C
 2010 TAUWIN  = 2.0*GMOL*UQUERN/DDS * VZ
      RETURN
C
      END
