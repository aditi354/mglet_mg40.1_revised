










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
      SUBROUTINE BBOBLO (KMX,JMX,IMX,
     $                      U,V,W,P,G,
     $                      ISTART,ISTOP,JSTART,JSTOP,KSTART,KSTOP,
     $                      FREQ,AU,AV,AW,ANIV,TIMEPHYS,
     $                      X,Y,Z,FLOWTYP,UBO,VBO,WBO,
     $                      WAVENUMBER,RANNUM,PHASE
     $                   )     
C*STARLET***************************************************************
C        B B O B L O   SETZEN DER RANDBEDINGUNGEN FUER DIE STROEMUNGS-
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
C        ANIV        - GIBT DAS NULLNIVEAU AN UM DAS DER COCOSUS SCHWINGT
C        NUMJET         - NUMBER OF JETS
C        FREQB          - EINSTROEMFREQUENZ 
C        ITYP           - HOLLERITH-KONSTANTE:
C                         'P'  VOR DER DRUCKKORREKTUR
C                         'T'  VOR DEM ZEITSCHRITT
C        WAVENUMBER       ANZAHL DER PERIODEN IN Y-RICHTUNG (LAMBDAWIRBEL) 
C 
C VERS:  25.02.94 (RK)  : ORIGINAL
C VERS:  28.09.94 (AO)  : TOTALLY NEW VERSION AND ADAPTED TO MGLET
C VERS:  21.04.98 (AO)  : KSTART, KSTOP ERMOEGLICHT MANIPULATION
C                         AN STUFENKANTEN ODER IM RAUM SELBER
C VERS:  29.10.98 (AO)  : IN WELCHE RICHTUNG WIRD MANIPULIERT ?
C                         IBLOW=1 FALLS BLOWING IN X-RICHTUNG
C VERS:  04.06.02 (TB)  : IBLOW-ELSE PART REMOVED  
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
     $        UBO(2,JMX,IMX),     VBO(2,JMX,IMX),      WBO(2,JMX,IMX) 
      REAL RANF, WAVENUMBER
      EXTERNAL RANF


      CHARACTER (LEN=16) FLOWTYP
C


C                               QUADRATIC UPSTREAM IST NICHT DEFINIERT
        DATA NQUD /0/



C                                 ***********************************
C                                 RECHTECKIGES EINSTROEMGEBIET BOTTOM
C                                 ***********************************
       IBLOW=0.0
C          IST DIE BLOWING-RICHTUNG in X ?
C     IF(ISTOP-ISTART.EQ.1) THEN
C        IBLOW=1
C        WRITE(*,*) 'BLOWING IN X-DIRECTION'
C     ENDIF

C                               LOOP OVER EVERY FREQUENCY OF THE JET
C
	 IF (FLOWTYP.EQ.'UNIFORM') THEN

         VALUEK = 2.0 * 3.1415927 * FREQ * TIMEPHYS
C
C           AUF ORIGIN COMPILERFEHLER !!!
C           COS(VALUEK) WIRD IMMER NULL
C           DESWEGEN NUN MIT SINUS BELEGT  (AO 02.04.1998)
C           SINUS + pi/2 ergibt wieder cosinus
C A1_TB240402 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         COSUS= -SIN(VALUEK)
C         COSUS= SIN(PHASE+VALUEK+1.5707963)
         VALUEU = ANIV + AU * COSUS
         VALUEV = ANIV + AV * COSUS
         VALUEW = ANIV + AW * COSUS
C         WRITE(6,*)'BBOBLO(ANIV,AW,VALUEW)  :',ANIV,AW,VALUEW
         VALUEG = GMOL
	 YA = Y(JSTART)
	 YB= Y(JSTOP)


         IF (IBLOW.EQ.1) THEN
C             MANIPULATION IN X-RICHTUNG !!!!
C----------- U - KOMPONENTE ------------------------------------------
         DO I = ISTART, ISTOP
         DO J = JSTART, JSTOP
         DO K = KSTART, KSTOP
         VWY = 1.0
         U(K,J,I) = VALUEU*VWY
         ENDDO
         ENDDO
         ENDDO

C----------- V - KOMPONENTE ------------------------------------------
         DO I=ISTART,ISTOP
         DO J=JSTART,JSTOP-1
         DO K = KSTART, KSTOP
         V(K,J,I) = -V(K,J,I+1) + VALUEV*VWY
         ENDDO
         ENDDO
         ENDDO

C----------- W - UND G - KOMPONENTE ----------------------------------
      ELSE
C C9_TB260502 Achung: Doppelbelegung CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         VALUEK = 2.0 * 3.1415927 * FREQ * TIMEPHYS
C         DO I=ISTART,ISTOP
C         DO J=JSTART,JSTOP
C         DO K = KSTART, KSTOP
C        W(K,J,I)     = -W(K,J,I+1)+ VALUEW*VWY
C        G(K,J,I)     =                   VALUEG
C        ENDDO
C        ENDDO
C        ENDDO
C---------------------------------------------------------------------
      ENDIF
C----------- U - KOMPONENTE ------------------------------------------

         DO I = ISTART, ISTOP
         DO J = JSTART, JSTOP
         DO K = KSTART, KSTOP
         VWY = 1.0

C        U(K,J  ,I  ) = -U(K+1,J  ,I  ) +  VALUEU*VWY
         U(K,J  ,I  ) =  VALUEU*VWY + U(K,J,I)

         ENDDO
         ENDDO
         ENDDO

C----------- V - KOMPONENTE ------------------------------------------

         DO I=ISTART,ISTOP
         DO J=JSTART,JSTOP-1
         DO K = KSTART, KSTOP

C        V(K,J  ,I  ) = -V(K+1,J  ,I  ) + VALUEV*VWY
         V(K,J  ,I  ) = VALUEV*VWY + V(K,J,I)

         ENDDO
         ENDDO
         ENDDO

C----------- W - UND G - T - KOMPONENTE --------------------------------


	 VALUEK = 2.0 * 3.1415927 * FREQ * TIMEPHYS
         DO I=ISTART,ISTOP
         DO J=JSTART,JSTOP
         DO K = KSTART, KSTOP

C        WRITE(6,*)'BBOBLO(K,J,I,W)  :',K, J, I, W(K,J  ,I  )
        W(K,J  ,I  )     =   VALUEW*VWY + W(K,J,I)
        G(K,J  ,I  )     =                   VALUEG
C        WRITE(6,*)'BBOBLO(K,J,I,W)  :',K, J, I, W(K,J  ,I  )
        ENDDO
        ENDDO
        ENDDO

      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	 IF (FLOWTYP.EQ.'UNIFORM-1') THEN

         VALUEK = 2.0 * 3.1415927 * FREQ * TIMEPHYS
C
         COSUS= -SIN(VALUEK+1.5707963)
         VALUEU = ANIV + AU * COSUS
         VALUEV = ANIV + AV * COSUS
         VALUEW = ANIV + AW * COSUS
         VALUEG = GMOL
	 YA = Y(JSTART)
	 YB= Y(JSTOP)


         IF (IBLOW.EQ.1) THEN
C             MANIPULATION IN X-RICHTUNG !!!!
C----------- U - KOMPONENTE ------------------------------------------
         DO I = ISTART, ISTOP
         DO J = JSTART, JSTOP
         DO K = KSTART, KSTOP
         VWY = 1.0
         U(K,J,I) = VALUEU*VWY
         ENDDO
         ENDDO
         ENDDO

C----------- V - KOMPONENTE ------------------------------------------
         DO I=ISTART,ISTOP
         DO J=JSTART,JSTOP-1
         DO K = KSTART, KSTOP
         V(K,J,I) = -V(K,J,I+1) + VALUEV*VWY
         ENDDO
         ENDDO
         ENDDO

C----------- W - UND G - KOMPONENTE ----------------------------------
      ELSE

	 VALUEK = 2.0 * 3.1415927 * FREQ * TIMEPHYS
         DO I=ISTART,ISTOP
         DO J=JSTART,JSTOP
         DO K = KSTART, KSTOP
        W(K,J,I)     = -W(K,J,I+1)+ VALUEW*VWY
        G(K,J,I)     =                   VALUEG
        ENDDO
        ENDDO
        ENDDO
C---------------------------------------------------------------------
      ENDIF
C----------- U - KOMPONENTE ------------------------------------------

         DO I = ISTART, ISTOP
         DO J = JSTART, JSTOP
         DO K = KSTART, KSTOP
         VWY = 1.0

C        U(K,J  ,I  ) = -U(K+1,J  ,I  ) +  VALUEU*VWY
         U(K,J  ,I  ) =  VALUEU*VWY

         ENDDO
         ENDDO
         ENDDO

C----------- V - KOMPONENTE ------------------------------------------

         DO I=ISTART,ISTOP
         DO J=JSTART,JSTOP-1
         DO K = KSTART, KSTOP

C        V(K,J  ,I  ) = -V(K+1,J  ,I  ) + VALUEV*VWY
         V(K,J  ,I  ) = VALUEV*VWY

         ENDDO
         ENDDO
         ENDDO

C----------- W - UND G -, T- KOMPONENTE --------------------------------


	 VALUEK = 2.0 * 3.1415927 * FREQ * TIMEPHYS
         DO I=ISTART,ISTOP
         DO J=JSTART,JSTOP
         DO K = KSTART, KSTOP

        W(K,J  ,I  )     =   VALUEW*VWY
        G(K,J  ,I  )     =                   VALUEG
        ENDDO
        ENDDO
        ENDDO

      ENDIF


      IF (FLOWTYP.EQ.'WAVE') THEN
C     DO IFREQ = 1 , NUMFREQ

	 XA = X(ISTART)
	 XB= X(ISTOP)
	 YA = Y(JSTART)
	 YB= Y(JSTOP)


C----------- U - KOMPONENTE ------------------------------------------

	 VALUEK = 2.0 * 3.1415927 * FREQ * TIMEPHYS
         COSUS= COS(VALUEK)
         DO J = JSTART, JSTOP
         DO I = ISTART, ISTOP
	 ZETAX = 2.0 * 3.1415927* (X(I)-XA)/(XB-XA)
	 ZETAY =2.0* WAVENUMBER * 3.1415927*(Y(J)-YA)/(YB-YA) 
C	 VWX = SIN(ZETAX) * (1-COS(ZETAX))
	 VWX = SIN(ZETAX)
	 VWY = SIN(ZETAY)
         VALUEU = ANIV + AU * COSUS

            U(2,J  ,I  ) = U(2,J  ,I ) + VALUEU * VWX*VWY

         ENDDO
         ENDDO

C----------- V - KOMPONENTE ------------------------------------------

	 VALUEK = 2.0 * 3.1415927 * FREQ * TIMEPHYS
         COSUS= COS(VALUEK)
         DO J=JSTART,JSTOP
         DO I=ISTART,ISTOP
	 ZETAX = 2.0 * 3.1415927* (X(I)-XA)/(XB-XA)
	 ZETAY =2.0* WAVENUMBER * 3.1415927*(Y(J)-YA)/(YB-YA) 
C	 VWX = SIN(ZETAX) * (1-COS(ZETAX))
	 VWX = SIN(ZETAX)
	 VWY = SIN(ZETAY)
C	 VWX = SIN(ZETAY)*(1-COS(ZETAY)) 
         VALUEV = ANIV + AV * COSUS
C         DO J=JSTART,JSTOP

            V(2,J  ,I  ) = V(2,J  ,I ) + VALUEV * VWX*VWY

         ENDDO
         ENDDO
 
C----------- W - UND G -, T - KOMPONENTE -------------------------------

	 VALUEK = 2.0 * 3.1415927 * FREQ * TIMEPHYS
         COSUS= COS(VALUEK)
	 DO J=JSTART,JSTOP
         DO I=ISTART,ISTOP
	 ZETAX = 2.0 * 3.1415927* (X(I)-XA)/(XB-XA)
	 ZETAY = 2.0*WAVENUMBER * 3.1415927*(Y(J)-YA)/(YB-YA) 
C	 VWX = SIN(ZETAX) * (1-COS(ZETAX))
	 VWX = SIN(ZETAX)
	 VWY = SIN(ZETAY)
         VALUEW = ANIV + AW * COSUS
C         VALUEW = ANIV + AW *  VWX 
         VALUEG = GMOL

            W(2,J  ,I  )     =W(2,J  ,I  ) +  VALUEW*VWX*VWY
            G(2,J  ,I  )     =                   VALUEG

         ENDDO
         ENDDO
 
C     ENDDO
      ENDIF
C------------------PARABELFOERMIGER JET----------------------------

      IF (FLOWTYP.EQ.'PARABEL') THEN
C     DO IFREQ = 1 , NUMFREQ

	 XA = X(ISTART)
	 XB= X(ISTOP)
	 YA = Y(JSTART)
	 YB= Y(JSTOP)


C----------- U - KOMPONENTE ------------------------------------------

	 VALUEK = 2.0 * 3.1415927 * FREQ * TIMEPHYS
         COSUS= COS(VALUEK)
         DO I = ISTART, ISTOP
           DO J = JSTART, JSTOP
	 ZETAX = 2.0 * 3.1415927* (X(I)-XA)/(XB-XA)
	 ZETAY = 2.0*WAVENUMBER * 3.1415927*(Y(J)-YA)/(YB-YA) 
C	 VWX = SIN(ZETA) * (1-COS(ZETA))
	 VWX = SIN(ZETAX)
	 VWY = SIN(ZETAY)
         VALUEU = ANIV + AU * COSUS
           U(2,J  ,I  ) = U(2,J  ,I  ) + VALUEU*VWX * VWY
           ENDDO
         ENDDO


C----------- V - KOMPONENTE ------------------------------------------

	 VALUEK = 2.0 * 3.1415927 * FREQ * TIMEPHYS
         COSUS= COS(VALUEK)
         VALUEV = ANIV + AV * COSUS
         DO I=ISTART,ISTOP
         DO J=JSTART,JSTOP
	 ZETAX = 2.0 * 3.1415927* (X(I)-XA)/(XB-XA)
	 ZETAY = 2.0*WAVENUMBER * 3.1415927*(Y(J)-YA)/(YB-YA) 
C	 VWX = SIN(ZETA) * (1-COS(ZETA))
	 VWX = SIN(ZETAX)
	 VWY = SIN(ZETAY)
         VALUEV = ANIV + AV * COSUS
            V(2,J  ,I  ) = V(2,J  ,I  ) + VALUEV*VWX*VWY
         ENDDO
         ENDDO
 
C----------- W - UND G -, T - KOMPONENTE -------------------------------

	 VALUEK = 2.0 * 3.1415927 * FREQ * TIMEPHYS
         COSUS= COS(VALUEK)
         DO I=ISTART,ISTOP
         DO J=JSTART,JSTOP
	 ZETAY = 2.0*WAVENUMBER * 3.1415927*(Y(J)-YA)/(YB-YA) 
	 VWY = SIN(ZETAY)
	 VWX =(((XA-XB)/2)**2-(X(I)-(XB+XA)/2)**2)/((XA-XB)/2)**2
         VALUEW = ANIV + AW * COSUS
         VALUEG = GMOL
C         DO J=JSTART,JSTOP

            W(2,J  ,I  )     = W(2,J  ,I  ) + VALUEW*VWX*VWY
            G(2,J  ,I  )     =                   VALUEG
         ENDDO
         ENDDO
 
      ENDIF

C--------RANDOM GENERATOR UNIFORM OVER A SPECIFIED SCALE ----------

      IF (FLOWTYP.EQ.'RANDOM') THEN

C        VALUEU = ANIV + AU * RANF()
         VALUEU = ANIV + AU * RANNUM
         DO I = ISTART, ISTOP
         DO J = JSTART, JSTOP
         DO K = KSTART, KSTOP

         U(K,J  ,I  ) =  VALUEU

         ENDDO
         ENDDO
         ENDDO

C----------- V - KOMPONENTE ------------------------------------------

C        VALUEV = ANIV + AV * RANF()
         VALUEV = ANIV + AV * RANNUM
         DO I=ISTART,ISTOP
         DO J=JSTART,JSTOP-1
         DO K = KSTART, KSTOP

         V(K,J  ,I  ) = VALUEV

         ENDDO
         ENDDO
         ENDDO

C----------- W - UND G -, T - KOMPONENTE -------------------------------


C        VALUEW = ANIV + AW * RANF()
         VALUEW = ANIV + AW * RANNUM
         VALUEG = GMOL
         DO I=ISTART,ISTOP
         DO J=JSTART,JSTOP
         DO K = KSTART, KSTOP

        W(K,J  ,I  )     =   VALUEW
        G(K,J  ,I  )     =                   VALUEG
        ENDDO
        ENDDO
        ENDDO



      ENDIF
      

	 IF (FLOWTYP.EQ.'KILLV') THEN


         DO I = ISTART, ISTOP
         DO J = JSTART, JSTOP
         DO K = KSTART, KSTOP
               V(K,J,I) = 0.0
         ENDDO
         ENDDO
         ENDDO

         ENDIF

	 IF (FLOWTYP.EQ.'KILLVWRAND') THEN


         DO I = ISTART, ISTOP
         DO J = JSTART, JSTOP
         DO K = KSTART, KSTOP
               V(K,J,I) = -V(K+1,J,I)
               W(K,J,I) = -W(K+1,J,I)
           VMAX = MAX(V(K,J,I),VMAX)
           WMAX = MAX(W(K,J,I),WMAX)
         ENDDO
         ENDDO
         ENDDO
 
      WRITE(*,*) 'MANIPULATION','VMAX=',VMAX,'WMAX=',WMAX

         ENDIF

      RETURN
      END

