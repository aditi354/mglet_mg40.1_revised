










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
         SUBROUTINE MGNBRCHECK (IG1,IBOC,IDIR,IG2,
     $                          IDIM1D,IDIM2D,NBND,X,Y,Z,
     $                          XSHIFT,YSHIFT,ZSHIFT,
     $                          BUF,IFOUND,ICOMPLETE)
C
C--MGLET----------------------------------------------------------------
C
C                    PRUEFT NACH, OB DIE ZWEI GITTER IG1 UND IG2
C                    IN DER RICHTUNG IDIR NACHBARN SIND
C                    UND OB SIE SICH UEBERLAPPEN
C                    UND OB DIE RANDFLAECHE VON GITTER 1 VOLLSTAENDIG
C                    VON NACHBARN UEBERLAPPT IST
C                  
C
C
C        17. 2.94 (MM)  : ORIGINAL
C
C--MGLET----------------------------------------------------------------
C

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


      COMMON /COMGRID/
     &               NGRID,    NGRDOLD,   NGRDSET, NGRDDFD,
     &               NAUFP,    NAUFOLD,
     &                 KMX,  JMX,  IMX,
     &                 KMXA, JMXA, IMXA,
     &                IP3D, IP2D, IP1D, IPBB,IPB3, IPBU,
     &                NOF3D,NOF2D,NOF1D,NOFBB,NOFB3,NOFBU,
     &                IPA,   IP1L,  IP2L,
     &                NOFA,  NOF1L, NOF2L,  
     &                 IC1,  IC2,  JC1,  JC2,  KC1,  KC2

      INTEGER
     &         KMX(MAXGRIDS),     JMX(MAXGRIDS),       IMX(MAXGRIDS),
     &        KMXA(MAXGRIDS),    JMXA(MAXGRIDS),      IMXA(MAXGRIDS),
     &        IP3D(MAXGRIDS),    IP2D(MAXGRIDS),      IP1D(MAXGRIDS),
     &        IPBB(MAXGRIDS),    IPB3(MAXGRIDS),      IPBU(MAXGRIDS),
     &        NOF3D,   NOF2D,   NOF1D,    NOFBB,   NOFBU,
     &         IPA(MAXGRIDS),    IP1L(MAXGRIDS),      IP2L(MAXGRIDS),
     &         NOFA,   NOF1L,   NOF2L,    NAUFP,
     &         IC1(MAXGRIDS),     IC2(MAXGRIDS),
     &         JC1(MAXGRIDS),     JC2(MAXGRIDS),
     &         KC1(MAXGRIDS),     KC2(MAXGRIDS)
 

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
C
      REAL    
     $            X(IDIM1D),     Y(IDIM1D),     Z(IDIM1D)

      REAL BUF(IDIM2D)
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      IFOUND = 0
      ICOMPLETE = 0
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                                      DIMENSIONIERUNGEN UND POINTER 
C
      CALL MGDPB (KK1,JJ1,II1,IP31,IP21,IP11,IBB1,IBU1,
     $                            NFR,NBA,NRI,NLE,NBO,NTO,NCU,IG1)
C
      CALL MGDPB (KK2,JJ2,II2,IP32,IP22,IP12,IBB2,IBU2,
     $                            NFR,NBA,NRI,NLE,NBO,NTO,NCU,IG2)
C
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C                                UEBERLAPPEN SICH DIE GITTER?
C
        CALL MGOVERLAP (II1,X(IP11),II2,X(IP12),NBND,XSHIFT,
     $                  IA1,IE1,IA2,IE2)
        CALL MGOVERLAP (JJ1,Y(IP11),JJ2,Y(IP12),NBND,YSHIFT,
     $                  JA1,JE1,JA2,JE2)
        CALL MGOVERLAP (KK1,Z(IP11),KK2,Z(IP12),NBND,ZSHIFT,
     $                  KA1,KE1,KA2,KE2)
C
CC      WRITE(6,*) 'MELDUNG AUS mgnbrcheck.sr35 '
CC      WRITE(6,'(A,20I4)') 'UEBERLAPPEN:',IG1,IG2,IA1,IE1,IA2,IE2
CC     $                    ,JA1,JE1,JA2,JE2,KA1,KE1,KA2,KE2
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C                                         FRONT-FLAECHE VON GITTER 1
      IF (IDIR.EQ.1) THEN

C170795        IF ((IA1 .EQ. 1) .AND. (IE1 .EQ. 4)) THEN
        IF ( IE1-NBND .EQ. IE2-(II2-NBND) ) THEN

          IF  (( JE1-NBND .GT. NBND ) .AND. ( JA1+NBND .LE. JJ1-NBND )
     $   .AND. ( KE1-NBND .GT. NBND ) .AND. ( KA1+NBND .LE. KK1-NBND ))
     $    THEN

              IFOUND = 1

              CALL MGNBRBUF (JJ1,KK1,BUF,JA1,JE1,KA1,KE1,ICOMPLETE)
              CALL SETCOBONE (IG1, 1 ,IBOC,FRONT(IBOC,IG1),
     $                        0,0,JA1,JE1,KA1,KE1,IG2,X,Y,Z,IDIM1D)

           ENDIF

        ELSE

           IFOUND = 0

        ENDIF

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C                                         BACK-FLAECHE VON GITTER 1
      ELSEIF (IDIR.EQ.2) THEN

C170795        IF ((IA1 .EQ. II1-NBND-1) .AND. (IE1 .EQ. II1)) THEN
        IF ( IE2-NBND .EQ. IE1-(II1-NBND) ) THEN

          IF  (( JE1-NBND .GT. NBND ) .AND. ( JA1+NBND .LE. JJ1-NBND )
     $   .AND. ( KE1-NBND .GT. NBND ) .AND. ( KA1+NBND .LE. KK1-NBND ))
     $    THEN

              IFOUND = 1

              CALL MGNBRBUF (JJ1,KK1,BUF,JA1,JE1,KA1,KE1,ICOMPLETE)
              CALL SETCOBONE (IG1, 2 ,IBOC,BACK(IBOC,IG1),
     $                        0,0,JA1,JE1,KA1,KE1,IG2,X,Y,Z,IDIM1D)

           ENDIF

        ELSE

           IFOUND = 0

        ENDIF

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C                                         RIGHT-FLAECHE VON GITTER 1
      ELSEIF (IDIR.EQ.3) THEN

C170795        IF ((JA1 .EQ. 1) .AND. (JE1 .EQ. 4)) THEN
        IF ( JE1-NBND .EQ. JE2-(JJ2-NBND) ) THEN

          IF  (( IE1-NBND .GT. NBND ) .AND. ( IA1+NBND .LT. II1-NBND )
     $   .AND. ( KE1-NBND .GT. NBND ) .AND. ( KA1+NBND .LT. KK1-NBND ))
     $    THEN

              IFOUND = 1

              CALL MGNBRBUF (II1,KK1,BUF,IA1,IE1,KA1,KE1,ICOMPLETE)
              CALL SETCOBONE (IG1, 3 ,IBOC,RIGHT(IBOC,IG1),
     $                        IA1,IE1,0,0,KA1,KE1,IG2,X,Y,Z,IDIM1D)

           ENDIF

        ELSE

           IFOUND = 0

        ENDIF

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C                                         LEFT-FLAECHE VON GITTER 1
      ELSEIF (IDIR.EQ.4) THEN

C170795        IF ((JA1 .EQ. JJ1-NBND-1) .AND. (JE1 .EQ. JJ1)) THEN
        IF ( JE2-NBND .EQ. JE1-(JJ1-NBND) ) THEN

          IF  (( IE1-NBND .GT. NBND ) .AND. ( IA1+NBND .LT. II1-NBND )
     $   .AND. ( KE1-NBND .GT. NBND ) .AND. ( KA1+NBND .LT. KK1-NBND ))
     $    THEN

              IFOUND = 1

              CALL MGNBRBUF (II1,KK1,BUF,IA1,IE1,KA1,KE1,ICOMPLETE)
              CALL SETCOBONE (IG1, 4 ,IBOC,LEFT(IBOC,IG1),
     $                        IA1,IE1,0,0,KA1,KE1,IG2,X,Y,Z,IDIM1D)

           ENDIF

        ELSE

           IFOUND = 0

        ENDIF


C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C                                         BOTTOM-FLAECHE VON GITTER 1
      ELSEIF (IDIR.EQ.5) THEN

C170795        IF ((KA1 .EQ. 1) .AND. (KE1 .EQ. 4)) THEN
        IF ( KE1-NBND .EQ. KE2-(KK2-NBND) ) THEN

          IF  (( IE1-NBND .GT. NBND ) .AND. ( IA1+NBND .LT. II1-NBND )
     $   .AND. ( JE1-NBND .GT. NBND ) .AND. ( JA1+NBND .LT. JJ1-NBND ))
     $    THEN

              IFOUND = 1

              CALL MGNBRBUF (II1,JJ1,BUF,IA1,IE1,JA1,JE1,ICOMPLETE)
              CALL SETCOBONE (IG1, 5 ,IBOC,BOTTOM(IBOC,IG1),
     $                        IA1,IE1,JA1,JE1,0,0,IG2,X,Y,Z,IDIM1D)

           ENDIF

        ELSE

           IFOUND = 0

        ENDIF

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
C                                         TOP-FLAECHE VON GITTER 1
      ELSEIF (IDIR.EQ.6) THEN

C170795        IF ((KA1 .EQ. KK1-NBND-1) .AND. (KE1 .EQ. KK1)) THEN
        IF ( KE2-NBND .EQ. KE1-(KK1-NBND) ) THEN

          IF  (( IE1-NBND .GT. NBND ) .AND. ( IA1+NBND .LT. II1-NBND )
     $   .AND. ( JE1-NBND .GT. NBND ) .AND. ( JA1+NBND .LT. JJ1-NBND ))
     $    THEN

              IFOUND = 1

              CALL MGNBRBUF (II1,JJ1,BUF,IA1,IE1,JA1,JE1,ICOMPLETE)
              CALL SETCOBONE (IG1, 6 ,IBOC,TOP(IBOC,IG1),
     $                        IA1,IE1,JA1,JE1,0,0,IG2,X,Y,Z,IDIM1D)

           ENDIF

        ELSE

           IFOUND = 0

        ENDIF

C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C
      ELSE

         CALL ERRR (501,'MGNBRCHECK')

      ENDIF
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C

      RETURN
      END


