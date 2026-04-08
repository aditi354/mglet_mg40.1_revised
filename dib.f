










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
      SUBROUTINE DIB (CIDENT,IIDENT,RIDENT,KK,JJ,II,KMX,JMX,IMX,NBND,
     $                DDX,DDY,DDZ,DX,DY,DZ,X,Y,Z,NVAR,
     $                U,V,W,P,G,B,TK,TE,HILF,IGRID,IDIM3D)
C*STAR******************************************************************
C  D I B      EINLESEN VON DATEN (BINAER,KANAL1)
C*STAR********************************************** W.STEITZ 20.08.81 *
C                                         GEAENDERT: ST.FRANK 20.10.83 *
C                                                    ST.FRANK 23.05.84 *
C                                                    ST.FRANK 02.07.84 *
C                                                    H.WERNER 19.09.85
C                                                    M.MANHART 6. 4.92
C                                                    M.MANHART 7. 6.93
C                                                    T.BRUNNER 06.03.03
C     26. 5.95 (MM) : MPI EINGEFUEHRT, KOPF UM HILF ERWEITERT
C     06.02.03 (TB) : SCALAR TRANSPORT IMPLEMENTED
C
C  PARAMETER
C             KK, JJ, II     - ARRAYGRENZEN
C             KMX, JMX, IMX  - ANZAHL DER GITTERPKTE. IN Z-,Y-,X-DIR.
C             X(II)          - KOORDINATEN IN X-RICHTUNG
C             Y(JJ)          - KOORDINATEN IN Y-RICHTUNG
C             Z(KK)          - KOORDINATEN IN Z-RICHTUNG
C             DDX(II)        - X-KANTENLAENGE DER KONTROLLVOLUMINA
C             DDY(JJ)        - Y-KANTENLAENGE DER KONTROLLVOLUMINA
C             DDZ(KK)        - Z-KANTENLAENGE DER KONTROLLVOLUMINA
C             DX(II)         - X-ABSTAND DER GITTERPUNKTE
C             DY(JJ)         - Y-ABSTAND DER GITTERPUNKTE
C             DZ(KK)         - Z-ABSTAND DER GITTERPUNKTE
C             NVAR           - ANZAHL DER AUSZUGEBENDEN DATENFELDER
C                              FOR VERSION < 8 (older than 6.04.92)
C             U(KK,JJ,II)    - GESCHWINDIGKEITSFELD
C             V(KK,JJ,II)    - GESCHWINDIGKEITSFELD
C             W(KK,JJ,II)    - GESCHWINDIGKEITSFELD
C             T(KK,JJ,II)    - SCALAR FIELD
C             P(KK,JJ,II)    - DRUCKFELD
C             G(KK,JJ,II)    - GAMMAFELD
C             B(KK,JJ,II)    - WANDERKENNUNGSFELD
C             TK(KK,JJ,II)   - DISSIPATIONSFELD
C             TE(KK,JJ,II)   - TURBULENZENERGIEFELD
C
C  ANMERKUNG  CIDENT(10) ENTHAELT KENNUNG
C
C             CIDENT( 1) = "AERO WENGLE"
C             CIDENT( 2) = PROGRAMMNAME
C             CIDENT( 3) = LAUFENDE NUMMER
C             CIDENT( 4) = DATE()
C             CIDENT( 5) = CLOCK()
C
C             IIDENT(100) ENTHAELT INTEGER-KONSTANTEN
C
C             IIDENT( 1) = MTSTEP
C             IIDENT( 2) = MPCORR
C             IIDENT( 5) : VERSIONSNUMMER
C             IIDENT(10) = KB
C             IIDENT(11) = JB1
C             IIDENT(12) = JB2
C             IIDENT(13) = IB1
C             IIDENT(14) = IB2
C             IIDENT(50) = SCALAR FIELD IDENTIFIER
C                           0:= SCALAR FIELD UNSET, FIELD NOT READ
C                           1:= SCALAR FIELD SET, FIELD CAN BE READ & 
C                                                WILL BE WRITTEN)
C
C             RIDENT(100) ENTHAELT REAL-KONSTANTEN
C
C             RIDENT( 1) = RE
C             RIDENT( 2) = DT
C             RIDENT( 3) = EPCORR
C             RIDENT(10) = ZWIDTH
C             RIDENT(11) = YWIDTH
C             RIDENT(12) = XWIDTH
C             RIDENT(13) = BETA
C             ifdef _TSCAL_:
C             RIDENT(50) = CP              
C             RIDENT(51) = PRMOL           
C             RIDENT(52) = PRTURB          
C             RIDENT(53) = EXPONT          
C
C*STAR******************************************************************
      CHARACTER (LEN=8)    CIDENT(10)
      CHARACTER (LEN=80)  CT,  CH1680, TEXT
      INTEGER KK, JJ, II, KMX, JMX, IMX, NVAR, IIDENT(100)
      REAL    RIDENT(100), X(IMX), Y(JMX), Z(KMX),
     $        DX(IMX), DY(JMX), DZ(KMX), DDX(IMX), DDY(JMX), DDZ(KMX),
     $        U( IDIM3D ), V( IDIM3D ), W( IDIM3D ),
     $        P( IDIM3D ), G( IDIM3D ), 
     $        B( IDIM3D ), HILF( IDIM3D )
C
C
C
C                          
C
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC FALLS ERGEBNIS AUS MLET
C                                          WIRD FOLGENDES UEBERSPRUNGEN
          WRITE(6,*)'JETZT IN DIB:',IIDENT(5)
      IF (IIDENT(5) .GT. 8) THEN

        READ(1,END=10001) TEXT
10001   READ(1,END=10002) IGRIDOLD
        WRITE(6,*)'TEXT: ',TEXT
        WRITE(6,*)'IGRIDOLD: ',IGRIDOLD
10002   READ(1,END=10003) KKOLD,JJOLD,IIOLD,KMXOLD,JMXOLD,IMXOLD
C                        KKOLD,JJOLD,IIOLD SIND DUMMYVARIABLE
C                        UND WERDEN NICHT BENOETIGT
C
C                        PLAUSIBILITAETSTEST
C
10003   IF(IGRID.NE.IGRIDOLD) CALL ERRR(700,' DIB      ')
        IF(KMX.NE.KMXOLD) CALL ERRR(701,' DIB      ')
        IF(JMX.NE.JMXOLD) CALL ERRR(702,' DIB      ')
        IF(IMX.NE.IMXOLD) CALL ERRR(703,' DIB      ')        

      ENDIF
C
      READ(1,END=10004) 
     $     (X(I),I=1,IMXOLD),(DX(I),I=1,IMXOLD),(DDX(I),I=1,IMXOLD),
     $     (Y(J),J=1,JMXOLD),(DY(J),J=1,JMXOLD),(DDY(J),J=1,JMXOLD),
     $     (Z(K),K=1,KMXOLD),(DZ(K),K=1,KMXOLD),(DDZ(K),K=1,KMXOLD)
C                                VERLAENGERUNG DER K,J,I-RICHTUNG
      CALL MGBASB  (NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,IGRID)
C     ONLY EQUIDISTANT ELONGATION PERMITTET FOR K DIRECTION
      IF (KMX .GT. KMXOLD) THEN
         DO K=KMXOLD+1, KMX
              Z(K) =   Z(K-1) + DZ(KMXOLD)
             DZ(K) =  DZ(K-1)
            DDZ(K) = DDZ(K-1)
         ENDDO
      ENDIF
C     PERIODIC AND EQUIDISTANT ELONGATION PERMITTET FOR J DIRECTION
      IF (JMX .GT. JMXOLD) THEN
      	IF (NRGT .EQ. 1 .AND. NLFT .EQ. 1) THEN
C        PERIODIC ELONGATION
         DO J=JMXOLD+1, JMX
              Y(J) =   Y(J-1) + DY(J-JMXOLD+2)
             DY(J) =  DY(J-JMXOLD+2)
            DDY(J) = DDY(J-JMXOLD+2)
         ENDDO      		
        ELSE
C        EQUIDISTANT ELONGATION        	
         DO J=JMXOLD+1, JMX
              Y(J) =   Y(J-1) + DY(JMXOLD)
             DY(J) =  DY(J-1)
            DDY(J) = DDY(J-1)
         ENDDO
        ENDIF
      ENDIF
C     PERIODIC AND EQUIDISTANT ELONGATION PERMITTET FOR I DIRECTION
      IF (IMX .GT. IMXOLD) THEN
C       EQUIDISTANT ELONGATION        	
        DO I=IMXOLD+1, IMX
             X(I) =   X(I-1) + DX(IMXOLD)
            DX(I) =  DX(I-1)
           DDX(I) = DDX(I-1)
        ENDDO
      ENDIF            
C
C
10004 CONTINUE
C                                    **********************************
C
C                                    AB HIER VERSION VOR 6.4.92
C
C                                    **********************************
      IF (IIDENT(5).LT.8) THEN
C
      READ(1,END=10005) NVAR
      IF (NVAR .LT. 4 .OR. NVAR .GT. 7) THEN
          CALL ERRR(404,' DIB      ')
      ENDIF
10005 WRITE(6,*) ' NVAR= ',NVAR
C
            CALL READ3D(KK,JJ,II,KMX,JMX,IMX,U,'BINAER  ',1)
            CALL READ3D(KK,JJ,II,KMX,JMX,IMX,V,'BINAER  ',1)
            CALL READ3D(KK,JJ,II,KMX,JMX,IMX,W,'BINAER  ',1)
            CALL READ3D(KK,JJ,II,KMX,JMX,IMX,P,'BINAER  ',1)
C
      IF (NVAR .EQ. 4) RETURN
            CALL READ3D(KK,JJ,II,KMX,JMX,IMX,G,'BINAER  ',1)
C
CC      IF (NVAR .EQ. 5) RETURN
CC            CALL READ3D(KK,JJ,II,KMX,JMX,IMX,TK,'BINAER  ',1)
C
CC      IF (NVAR .EQ. 6) RETURN
CC            CALL READ3D(KK,JJ,II,KMX,JMX,IMX,TE,'BINAER  ',1)
C
C
C
C                                    **********************************
C
C                                    AB HIER VERSION NACH 6.4.92
C
C                                    **********************************

      ELSEIF (IIDENT(5).GE.8) THEN

            KANAL = 1

            CT = CH1680 (' U              ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,'BINAER  ',
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   U      ,CT,NBND,HILF,IGRID,IDIM3D)
            CT = CH1680 (' V              ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,'BINAER  ',
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   V      ,CT,NBND,HILF,IGRID,IDIM3D)
            CT = CH1680 (' W              ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,'BINAER  ',
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   W      ,CT,NBND,HILF,IGRID,IDIM3D)
            CT = CH1680 (' P              ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,'BINAER  ',
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   P      ,CT,NBND,HILF,IGRID,IDIM3D)
            CT = CH1680 (' G              ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,'BINAER  ',
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   G      ,CT,NBND,HILF,IGRID,IDIM3D)
            CT = CH1680 (' B              ')
            CALL DIBC   (KK ,JJ ,II ,KSTAG,JSTAG,ISTAG,KANAL,'BINAER  ',
     $                   ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                   B      ,CT,NBND,HILF,IGRID,IDIM3D)

      ELSE
      STOP 'UNBEKANNTE VERSIONSNUMMER IN DIB'
      ENDIF
C
C
      RETURN
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
C                END OF FILE:
10000 STOP 'EOF BEIM EINLESEN VON CIDENT, IIDENT ODER RIDENT'
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72
      END
      SUBROUTINE DIBC (KKA,JJA,IIA,KSTAG,JSTAG,ISTAG,KANAL,MODUS,
     $                 ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,TEEM,DTEM,
     $                 PHI,TEXT,NBND,HILF,IGRID,IDIM3D)
C*STAR******************************************************************
C  D I B C         EINLESEN VON DATEN (BINAER,KANAL NR. WIRD UEBERGEBEN
C                  DIBC ERLAUBT DAS EINLESEN VON FELDERN MIT UNTER-
C                  SCHIEDLICHSTEN FELDDIMENSIONEN. DIESE FELDER
C                  WERDEN IM ANSCHLUSS AN DIE "STANDARDFELDER"
C                  CIDENT ... U,V,W,T,P,G,TK,TE  UND ISELEP GESCHRIEBEN.
C                  BEIM EINLESEN WIRD VORAUSGESETZT, DASS DIESE FELDER
C                  BEREITS GELESEN WURDEN, DAS DATENFILE ALSO RICHTIG
C                  POSITIONIERT IST !
C***********************************************************************
C
C PARAM: KKA,JJA,IIA    - ARRAYDIMENSIONEN
C        KSTAG,J..,I..  + INTEGER KENNZAHLEN, DIE ANGEBEN OB DER
C                         DEFINITIONSPUNKT DER UEBERGEBENEN
C                         VARIABLEN 'PHI' IN Z-, Y-, ODER X-RICHTUNG
C                         VERSCHOBEN IST.
C                         0 : DEFINITIONPKT. IST ZELLMITTELPUNKT
C                         1 :       ""        IST DIE IN POS. KO-
C                                            ORDINATENRICHTUNG VER-
C                                            SCHOBENE GRENZFLAECHE
C        KANAL          - DIE DATEN WERDEN UEBER DIESEN KANAL
C                         EINGELESEN
C        MODUS          - CHARACTER-VARIABLE:
C                         'BINAER  ' : DATEN WERDEN UNFORMA-
C                                      TIERT GESCHRIEBEN
C                         'CODIERT ' : DATEN WERDEN FORMA-
C                                      TIERT GESCHRIEBEN
C        ZSLM, ZELM     + START- UND STOPKOORDINATE DER LINIEN-
C                         MITTELUNG IN Z-RICHTUNG (Z. ZT. KEINE
C                         LINIENMITTELUNG IN Z-RI. REALISIERT
C                         12.08.86)
C        YSLM, YELM     + WIE OBEN, JEDOCH IN Y-RICHTUNG
C        XSLM, XELM     + WIE OBEN, JEDOCH IN X-RICHTUNG
C        TSEM, TEEM     + PHYSIKAL. ZEIT, ZU DER DIE ERSTE BZW.
C                         LETZTE STICHPROBE FUER DIE BILDUNG VON
C                         ENSEMBLE-MITTELWERTEN GENOMMEN WURDE
C        DTEM           + IN DIESEN ZEITABSTAENDEN WURDEN JEWEILS
C                         STICHPROBEN GENOMMEN
C        PHI(KKA,JJA,IIA)+ ALLGEMEINE VARIABLE
C        TEXT           - CHARACTER VARIABLE ZUR IDENTIFIKATION
C
C  UPROG                 : ERRR
C
C  VERS:  03.07.86 (HW)  : KSTAG... EINGEFUEHRT
C         12.08.86 (HW)  : DIBC VOELLIG NEU UEBERARBEITET
C         09.12.86 (HW)  : ES KANN WAHLWEISE FORMATIERT ODER UNFORMA-
C                          TIERT GELESEN WERDEN. (VARIABLE 'MODUS'
C                          EINGEFUEHRT)
C         18.12.86 (HW)  : TIDENT IN CIDENT (CHARACTER (LEN=8) !!) UMGE-
C                          WANDELT WEGEN NOS/VE.
C     26. 5.95 (MM) : MPI EINGEFUEHRT, KOPF UM HILF UND IGRID ERWEITERT
C     06.02.03 (TB) : SPECIAL TREATMENT FOR SCALAR TRANSPORT VARIABLE T:
C                     WHEN RESTARTING FROM AN OLD RESULT WITHOUT SCALAR
C                     TRANSPORT BUT _TSCAL_ SET, A NEW T FIELD IS 
C                     INITIALISED AND WRITTEN TO VARIABLE T
C
C*STAR******************************************************************
C
      CHARACTER (LEN=8)   MODUS
      CHARACTER (LEN=80)  TEXT,TEXTIN
      REAL           PHI(IDIM3D),HILF(IDIM3D)
C

      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )


      INTEGER MAXLVLMAX,MINLVLMIN,MAXGRDSOFLVL

      PARAMETER ( MAXLVLMAX = 4,  MINLVLMIN = -4 )
      PARAMETER ( MAXGRDSOFLVL =260)

      COMMON /COLEVEL/ MAXLEVEL,MINLEVEL
      COMMON /COLEVEL/ NOFLEVEL,IGRDOFLEVEL,NOFPSDIR,IGRDOFPSDIR
      COMMON /COLEVEL/ NOFVPIT ,IGRDOFVPIT ,NOFTST,  IGRDOFTST
      COMMON /COLEVEL/ NOFSLCED,IGRDOFSLCED,NOFSLICE,IGRDOFSLICE
      COMMON /COLEVEL/ NOFSCAI, IGRDOFSCAI,NGRIDPERLEVEL,NOFTHISGRID

      INTEGER MAXLEVEL,MINLEVEL

      INTEGER NOFLEVEL   (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFLEVEL(MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFVPIT    (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFVPIT (MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFTST     (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFTST  (MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFSCAI    (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFSCAI (MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFSLCED   (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFSLCED(MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFSLICE   (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFSLICE(MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFPSDIR   (             MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFPSDIR(MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NGRIDPERLEVEL    
      INTEGER NOFTHISGRID(MAXGRDSOFLVL)
      
      COMMON /COGRDPRO/ LEVEL,LCHILD,XMIN,YMIN,ZMIN,XTOT,YTOT,ZTOT
      COMMON /COGRDPRO/ XHOMOG,YHOMOG,ZHOMOG,NXGRAE,NYGRAE,NZGRAE
      COMMON /COGRDPRO/ GRADPX,UBULKX,LTST,LVP,LSCAI,LPLEVEL,LPOISSONDIR
      COMMON /COGRDPRO/ LSLICE,NXSLICE,NYSLICE,NZSLICE,NVPGRIDS
      COMMON /COGRDPRO/ CONV1SANF,CONV1SEND,TRANSLES1,TRANSLES2
      COMMON /COGRDPRO/ GRADPXOLD

      INTEGER LEVEL(MAXGRIDS),NVPGRIDS(MAXGRIDS)
      INTEGER NXGRAE(MAXGRIDS),NYGRAE(MAXGRIDS),NZGRAE(MAXGRIDS)
      INTEGER NXSLICE(MAXGRIDS),NYSLICE(MAXGRIDS),NZSLICE(MAXGRIDS)

      REAL XTOT(MAXGRIDS),YTOT(MAXGRIDS),ZTOT(MAXGRIDS)
      REAL XMIN(MAXGRIDS),YMIN(MAXGRIDS),ZMIN(MAXGRIDS)
      REAL GRADPX(MAXGRIDS),UBULKX(MAXGRIDS)
      REAL CONV1SANF(MAXGRIDS),CONV1SEND(MAXGRIDS)
      REAL TRANSLES1(MAXGRIDS),TRANSLES2(MAXGRIDS),GRADPXOLD(MAXGRIDS)

      LOGICAL LCHILD(MAXGRIDS),LSLICE(MAXGRIDS),LPOISSONDIR(MAXGRIDS)
      LOGICAL XHOMOG(MAXGRIDS),YHOMOG(MAXGRIDS),ZHOMOG(MAXGRIDS)
      LOGICAL LTST(MAXGRIDS),LVP(MAXGRIDS),LSCAI(MAXGRIDS)
      LOGICAL LPLEVEL(MAXGRIDS)
      
      COMMON /COGRDCON/
     &                 IVPCHILD,
     &                 IPARENT,ISLPAR,
     &                 IPOSITION, JPOSITION, KPOSITION,
     &                 NOFSLCHILDS,IGRDOFSLCHILD,
     &                 ISLPOS, JSLPOS, KSLPOS,
     &                 IFRNBR, IBANBR, IRINBR, ILENBR,
     &                 IBONBR, ITONBR


      INTEGER
     &       IVPCHILD(MAXGRIDS),
     &       IPARENT(MAXGRIDS),ISLPAR(MAXGRIDS),
     & IPOSITION(MAXGRIDS), JPOSITION(MAXGRIDS), KPOSITION(MAXGRIDS),
     & NOFSLCHILDS(MAXGRIDS),IGRDOFSLCHILD(MAXGRIDS,MAXGRIDS),
     &    ISLPOS(MAXGRIDS), JSLPOS(MAXGRIDS), KSLPOS(MAXGRIDS),
     &       IFRNBR(MAXBOCONDS,MAXGRIDS), IBANBR(MAXBOCONDS,MAXGRIDS),
     &       IRINBR(MAXBOCONDS,MAXGRIDS), ILENBR(MAXBOCONDS,MAXGRIDS),
     &       IBONBR(MAXBOCONDS,MAXGRIDS), ITONBR(MAXBOCONDS,MAXGRIDS)
C


         TEXTIN(1:80) = ' '
         IFAIL = 0
C
C
      IF(MODUS .EQ. 'BINAER  ') THEN
C
C                                 **************************************
C                                 UNFORMATIERTES LESEN
C                                 **************************************
C
            WRITE (6,*) 'Trying DIBC',TEXT,kka,jja,iia
         READ (KANAL,END=10200) TEXTIN
            WRITE (6,*) 'DIBC',TEXT,TEXTIN
10200    READ (KANAL,END=10205) KKAI,JJAI,IIAI,KMXAI,JMXAI,IMXAI
            WRITE (6,*) 'DIBC',KKAI,JJAI,IIAI,KMXAI,JMXAI,IMXAI
10205    READ (KANAL,END=10206) KSTAG,JSTAG,ISTAG
            WRITE (6,*) 'DIBC',KSTAG,JSTAG,ISTAG
10206    READ (KANAL,END=10210) ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM,
     $                          TEEM,DTEM
            WRITE (6,*) 'DIBC',ZSLM,ZELM,YSLM,YELM,XSLM,XELM,TSEM
C
10210    IF(TEXTIN(1:16) .NE. TEXT(1:16)) THEN
C
C                                  DER DATENSATZ STIMMT NICHT MIT DEM
C                                  GEFORDERTEN UEBEREIN, DAHER MELDUNG
C
            WRITE (6,6005) TEXT,TEXTIN
            IFAIL = 999
         ELSE
         ENDIF
C
C                                 PLAUSIBILITAETSTEST
C
            IF((KMXAI .NE. KKA)
     $     .OR.(JMXAI .NE. JJA)     
     $     .OR.(IMXAI .NE. IIA)) THEN     	

               WRITE (6,*) 'MELDUNG AUS DIBC:'
               WRITE (6,*) ' IIA,  JJA,  KKA:',IIA,JJA,KKA
               WRITE (6,*) 'IMXAI,JMXAI,KMXAI',IMXAI,JMXAI,KMXAI
               WRITE (6,*) 'TEXT,TEXTIN:',TEXT,TEXTIN
               IFAIL = 714
            ENDIF
C


            CALL MGPOINT (IP3,IP2,IP1,IBB,IB3,IBU,IGRID)
C        write (6,*) 'dibc, mgpoint',myid,igrid,ip3,ip2,ip1
C------------------------------------------- IS IT A STATIST. ARRAY

      IF (TEXT(1:1) .EQ. "<" ) THEN

          CALL MGPOINA (IP3,I1L,I2L,IGRID)

      ENDIF

            CALL READ3D
     $           (KKA,JJA,IIA,KMXAI,JMXAI,IMXAI,PHI(IP3),MODUS,KANAL)

C
C
      ELSEIF(MODUS .EQ. 'CODIERT ') THEN
C
C                                 **************************************
C                                 FORMATIERTES LESEN
C                                 **************************************
C
C
         READ (KANAL,6010) TEXTIN
         READ (KANAL,6020) KKAI,JJAI,IIAI,KMXAI,JMXAI,IMXAI
         READ (KANAL,6030) KSTAG,JSTAG,ISTAG
         READ (KANAL,6040) ZSLM,ZELM,YSLM,YELM
         READ (KANAL,6040) XSLM,XELM,TSEM,TEEM
         READ (KANAL,6050) DTEM
C
         IF(TEXTIN(1:16) .NE. TEXT(1:16)) THEN
C
C                                  DER DATENSATZ STIMMT NICHT MIT DEM
C                                  GEFORDERTEN UEBEREIN, DAHER MELDUNG
C
            WRITE (6,6005) TEXT,TEXTIN
            IFAIL = 999
         ELSE
         ENDIF
C
C                                 PLAUSIBILITAETSTEST
C
            IF((KMXAI .NE. KKA)
     $     .OR.(JMXAI .NE. JJA)     
     $     .OR.(IMXAI .NE. IIA)) THEN     	

               WRITE (6,*) 'MELDUNG AUS DIBC:'
               WRITE (6,*) ' IIA,  JJA,  KKA:',IIA,JJA,KKA
               WRITE (6,*) 'IMXAI,JMXAI,KMXAI',IMXAI,JMXAI,KMXAI
               WRITE (6,*) 'TEXT,TEXTIN:',TEXT,TEXTIN
               IFAIL = 714
            ENDIF
C
            CALL MGPOINT (IP3,IP2,IP1,IBB,IB3,IBU,IGRID)
C------------------------------------------- IS IT A STATIST. ARRAY

      IF (TEXT(1:1) .EQ. "<" ) THEN

          CALL MGPOINA (IP3,I1L,I2L,IGRID)

      ENDIF
            CALL READ3D
     $           (KKA,JJA,IIA,KMXAI,JMXAI,IMXAI,PHI(IP3),MODUS,KANAL)

C
C
      ELSE
C
C                                 BINAER, ODER CODIERT?
C
         IFAIL = 501

      ENDIF
c     
C                                               ALLES KLAR, ODER NICHT?
      IF (IFAIL .GT. 0) CALL ERRR (IFAIL,' DIBC ')
C
      IF (LSLICE(IGRID)) THEN
         CALL SLICEFIELD (PHI,HILF,IDIM3D,IGRID)
      ENDIF
      
C
C
C
      RETURN
C
 6005 FORMAT (//,1X,'DAS EINZULESENDE FELD ',A80,/,1X,
     $        'STIMMT NICHT MIT DEM FELD AUF DEM DATENFILE',A80,
     $        /,1X,'UEBEREIN !!')
 6010 FORMAT (A80)
 6020 FORMAT (6(I9,1X))
 6030 FORMAT (3(I9,1X))
 6040 FORMAT (6(E12.5E3,1X))
 6050 FORMAT (E19.12E3)
      END
C-MGLET--------------------------------------------------------------

      SUBROUTINE READ3D(KK,JJ,II,KMXAI,JMXAI,IMXAI,PHI,MODUS,KANAL)

      CHARACTER (LEN=8)   MODUS
      REAL PHI(KK,JJ,II)

      IF(MODUS .EQ. 'BINAER  ') THEN

         DO 10220 I = 1,IMXAI
            READ (KANAL,END=10220) ((PHI(K,J,I),K=1,KMXAI),J=1,JMXAI)
10220    CONTINUE

      ELSE

         DO 700 I = 1,IMXAI
            READ (KANAL,6040) ((PHI(K,J,I),K=1,KMXAI),J=1,JMXAI)
  700    CONTINUE
         
      ENDIF

 6040 FORMAT (6(E12.5E3,1X))

      RETURN
      END
