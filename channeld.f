










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
      SUBROUTINE CHANNELD(KK,JJ,II,NBND,
     $     DX,DY,DZ,DDX,DDY,DDZ,U,V,W,
     $     WALLSSX,WALLSSY,WALLSSZ,GMOL,RHO,UGRID,
     $     NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,IGRID
     $     )

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     BERECHNUNG DER GEMITTELTEN WANDSCHUBSPANNUNG
C
C     VERSION VOM 25.1.1998 (MM) 
C                           (NICHT ALLE SCHUBSPANNUNGEN IMPLEMENTIERT!)
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )

C
C     VERS:   28.09.95 (AO) JET VARIABLES INTRODUCED
C             12.12.02 (TB) SCALAR PARAMETERS (FIX VALUE, GRADIENT) ADDED
C
C     XMPOS1,XMPOS2,YMPOS1,YMPOS2,ZMPOS1,ZMPOS2: BEREICH DER EFFEKTMESSUNG 
C                                                AUFGRUND DER MANIPULATION

      COMMON /COBOUND/
     &                 NBOCD,
     &                NBOCONDS,     LARBOCONDS,     ITYPBOCONDS,
     &                LBOGRIDS,    LPOSBOGRIDS,
     &                 FRONT,    BACK,     RIGHT,     LEFT,
     &                 BOTTOM,   TOP,      CUBE,
     &                 IBPOS,   JBPOS,    KBPOS,
     &                 IBANF,   JBANF,    KBANF,
     &                 IBEND,   JBEND,    KBEND,
     &                 XBANF,   YBANF,    ZBANF,
     &                 XBEND,   YBEND,    ZBEND,
     &                 ANIVEAU,
     &                 AUB, AVB, AWB,
     &                 FREQB,  FLOWTYP, WAVENUMBER,
     &                 LOPT,NTOPT1,NTOPT2,
     &                 XMPOS1,YMPOS1,ZMPOS1,
     &                 XMPOS2,YMPOS2,ZMPOS2,
     &                 TIMEALT,XRTALT1,XRTALT2,FREQOPT,PERIODE,ITALT,
     &                 FREQALT1,FREQALT2,XRMIN,FXRMIN,NXRMIN,
     &                 RANNUM,PHASE


      INTEGER
     &       NBOCD   (                9, MAXGRIDS),
     &       NBOCONDS(                9, MAXGRIDS),
     &     LARBOCONDS( 6, MAXBOCONDS, 9, MAXGRIDS),
     &    ITYPBOCONDS(    MAXBOCONDS, 9, MAXGRIDS),
     &       LBOGRIDS(    MAXBOCONDS, 9, MAXGRIDS),
     &    LPOSBOGRIDS( 3, MAXBOCONDS, 9, MAXGRIDS),
     &   IBPOS(MAXBOCONDS,9,MAXGRIDS),  JBPOS(MAXBOCONDS,9,MAXGRIDS),
     &   KBPOS(MAXBOCONDS,9,MAXGRIDS),
     &   IBANF(MAXBOCONDS,9,MAXGRIDS),  IBEND(MAXBOCONDS,9,MAXGRIDS),
     &   JBANF(MAXBOCONDS,9,MAXGRIDS),  JBEND(MAXBOCONDS,9,MAXGRIDS),
     &   KBANF(MAXBOCONDS,9,MAXGRIDS),  KBEND(MAXBOCONDS,9,MAXGRIDS),
     &   LOPT(MAXBOCONDS,MAXGRIDS),
     &   NTOPT1(MAXBOCONDS,MAXGRIDS),NTOPT2(MAXBOCONDS,MAXGRIDS)



      CHARACTER (LEN=16)
     &      FRONT(MAXBOCONDS,MAXGRIDS),   BACK(MAXBOCONDS,MAXGRIDS),
     &      RIGHT(MAXBOCONDS,MAXGRIDS),   LEFT(MAXBOCONDS,MAXGRIDS),
     &     BOTTOM(MAXBOCONDS,MAXGRIDS),    TOP(MAXBOCONDS,MAXGRIDS),
     &       CUBE(MAXBOCONDS,MAXGRIDS),
     &      FLOWTYP(MAXBOCONDS,9,MAXGRIDS)

      REAL  ANIVEAU(MAXBOCONDS,9,MAXGRIDS), AUB(MAXBOCONDS,9,MAXGRIDS),
     &      AVB(MAXBOCONDS,9,MAXGRIDS), AWB(MAXBOCONDS,9,MAXGRIDS),
     &      FREQB(MAXBOCONDS,9,MAXGRIDS),
     &      XBANF(MAXBOCONDS,9,MAXGRIDS),XBEND(MAXBOCONDS,9,MAXGRIDS),
     &      YBANF(MAXBOCONDS,9,MAXGRIDS),YBEND(MAXBOCONDS,9,MAXGRIDS),
     &      ZBANF(MAXBOCONDS,9,MAXGRIDS),ZBEND(MAXBOCONDS,9,MAXGRIDS),
     &      XM1(MAXBOCONDS,9,MAXGRIDS),XM2(MAXBOCONDS,9,MAXGRIDS),
     &      XM3(MAXBOCONDS,9,MAXGRIDS),
     &      YM1(MAXBOCONDS,9,MAXGRIDS),YM2(MAXBOCONDS,9,MAXGRIDS),
     &      YM3(MAXBOCONDS,9,MAXGRIDS),
     &      ZM1(MAXBOCONDS,9,MAXGRIDS),ZM2(MAXBOCONDS,9,MAXGRIDS),
     &      ZM3(MAXBOCONDS,9,MAXGRIDS),
     &      WAVENUMBER(MAXBOCONDS,9,MAXGRIDS),
     &      XMPOS1(MAXBOCONDS,MAXGRIDS),
     &      XMPOS2(MAXBOCONDS,MAXGRIDS),
     &      YMPOS1(MAXBOCONDS,MAXGRIDS),
     &      YMPOS2(MAXBOCONDS,MAXGRIDS),
     &      ZMPOS1(MAXBOCONDS,MAXGRIDS),
     &      ZMPOS2(MAXBOCONDS,MAXGRIDS),
     &      RANNUM,
     &      PHASE(MAXBOCONDS,9,MAXGRIDS)

      REAL U(KK,JJ,II), V(KK,JJ,II), W(KK,JJ,II)

      REAL DDX(II),     DDY(JJ),     DDZ(KK)
      REAL  DX(II),      DY(JJ),      DZ(KK)

      REAL WALLSSX(6), WALLSSY(6), WALLSSZ(6)

      REAL LX,LY,LZ
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
      IMX = II
      JMX = JJ
      KMX = KK
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

      DO I=1,6
         WALLSSX(I) = 0.0
         WALLSSY(I) = 0.0
         WALLSSZ(I) = 0.0
      ENDDO

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C     BOTTOM- FLAECHE

      IF (NBOT .EQ. 5) THEN

         WALLSSX(5) = 0.0
         WALLSSY(5) = 0.0
         LX   = 0.0
         LY   = 0.0
      
         DO I=NBND+1,IMX-NBND
            LX = LX+DDX(I)
         ENDDO
         DO J=NBND+1,JMX-NBND
            LY = LY+DDY(J)
         ENDDO


         DO I=NBND+1,IMX-NBND
            DO J=NBND+1,JMX-NBND
               WALLSSX(5) = WALLSSX(5) + 
     $              (U(NBND+1,J,I)+UGRID)*DDY(J)*DX(I)
               WALLSSY(5) = WALLSSY(5) + 
     $              (V(NBND+1,J,I))*DY(J)*DDX(I)
            ENDDO
         ENDDO

         WALLSSX(5) = WALLSSX(5)*GMOL / (LX*LY* 0.5*DDZ(1+NBND) )
         WALLSSY(5) = WALLSSY(5)*GMOL / (LX*LY* 0.5*DDZ(1+NBND) )

      ENDIF
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C     TOP- FLAECHE

      IF (NTOP .EQ. 5) THEN

         WALLSSX(6) = 0.0
         WALLSSY(6) = 0.0
         LX   = 0.0
         LY   = 0.0
      
         DO I=NBND+1,IMX-NBND
            LX = LX+DDX(I)
         ENDDO
         DO J=NBND+1,JMX-NBND
            LY = LY+DDY(J)
         ENDDO


         DO I=NBND+1,IMX-NBND
            DO J=NBND+1,JMX-NBND
               WALLSSX(6) = WALLSSX(6) + 
     $              (U(KMX-NBND,J,I)+UGRID)*DDY(J)*DX(I)
               WALLSSY(6) = WALLSSY(6) + 
     $              (V(KMX-NBND,J,I))*DY(J)*DDX(I)

            ENDDO
         ENDDO

         WALLSSX(6) = WALLSSX(6)*GMOL / (LX*LY* 0.5*DDZ(KMX-NBND) )
         WALLSSY(6) = WALLSSY(6)*GMOL / (LX*LY* 0.5*DDZ(KMX-NBND) )

      ENDIF

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C     RIGHT- FLAECHE

      IF (NRGT .EQ. 5) THEN

         WALLSSX(3) = 0.0
         WALLSSZ(3) = 0.0

         LX   = 0.0
         LZ   = 0.0
      
         DO I=NBND+1,IMX-NBND
            LX = LX+DDX(I)
         ENDDO
         DO K=NBND+1,KMX-NBND
            LZ = LZ+DDZ(K)
         ENDDO


         DO I=NBND+1,IMX-NBND
            DO K=NBND+1,KMX-NBND
               WALLSSX(3) = WALLSSX(3) + 
     $              (U(K,1+NBND,I)+UGRID)*DDZ(K)*DDX(I)
               WALLSSZ(3) = WALLSSZ(3) +
     $              W(K,1+NBND,I) *DDZ(K)*DDX(I)
            ENDDO
         ENDDO

         WALLSSX(3) = WALLSSX(3)*GMOL / (LX*LZ* 0.5*DDY(1+NBND) )
         WALLSSZ(3) = WALLSSZ(3)*GMOL / (LX*LZ* 0.5*DDY(1+NBND) )
      ENDIF
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C     LEFT- FLAECHE

      IF (NLFT .EQ. 5) THEN

         WALLSSX(4) = 0.0
         WALLSSZ(4) = 0.0
         LX   = 0.0
         LZ   = 0.0
      
         DO I=NBND+1,IMX-NBND
            LX = LX+DDX(I)
         ENDDO
         DO K=NBND+1,KMX-NBND
            LZ = LZ+DDZ(K)
         ENDDO


         DO I=NBND+1,IMX-NBND
            DO K=NBND+1,KMX-NBND
               WALLSSX(4) = WALLSSX(4) + 
     $              (U(K,JMX-NBND,I)+UGRID)*DDZ(K)*DDX(I)
               WALLSSZ(4) = WALLSSZ(4) +
     $              W(K,JMX-NBND,I) *DDZ(K)*DDX(I)
            ENDDO
         ENDDO

         WALLSSX(4) = WALLSSX(4)*GMOL / (LX*LZ* 0.5*DDY(JMX-NBND) ) 
         WALLSSZ(4) = WALLSSZ(4)*GMOL / (LX*LZ* 0.5*DDY(JMX-NBND) ) 

      ENDIF


CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C     FRONT- FLAECHE

      IF (NFRO .EQ. 5) THEN

         WALLSSZ(1) = 0.0
         LY   = 0.0
         LZ   = 0.0
      
         DO J=NBND+1,JMX-NBND
            LY = LY+DDY(J)
         ENDDO
         DO K=NBND+1,KMX-NBND
            LZ = LZ+DDZ(K)
         ENDDO


         DO J=NBND+1,JMX-NBND
            DO K=NBND+1,KMX-NBND
               WALLSSZ(1) = WALLSSZ(1) + 
     $              W(K,J,1+NBND)*DDZ(K)*DDY(J)
            ENDDO
         ENDDO

         WALLSSZ(1) = WALLSSZ(1)*GMOL / (LY*LZ* 0.5*DDX(1+NBND) ) 
      ENDIF


CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C     BACK - FLAECHE

      IF (NBAC .EQ. 5) THEN

         WALLSSZ(2) = 0.0
         LY   = 0.0
         LZ   = 0.0
      
         DO J=NBND+1,JMX-NBND
            LY = LY+DDY(J)
         ENDDO
         DO K=NBND+1,KMX-NBND
            LZ = LZ+DDZ(K)
         ENDDO


         DO J=NBND+1,JMX-NBND
            DO K=NBND+1,KMX-NBND
               WALLSSZ(2) = WALLSSZ(2) + 
     $              W(K,J,IMX-NBND)*DDZ(K)*DDY(J)
            ENDDO
         ENDDO

         WALLSSZ(2) = WALLSSZ(2)*GMOL / (LY*LZ* 0.5*DDX(IMX-NBND) ) 
      ENDIF


      RETURN
      END

