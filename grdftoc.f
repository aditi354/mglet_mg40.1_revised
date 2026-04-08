










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
      SUBROUTINE GRDFTOC (IIF,XF,DXF,DDXF,IIC,XC,DXC,DDXC,CS,CSI)
C
C*MGLET*****************************************************************
C  G R D F T O C       ERZEUGEN EINES GITTERS FUER MULTIGRID
C                      ERZEUGT GROBGITTER AUS FEINGITTER
C                      DIE GITTERPUNKTE WERDEN IN MITTE DER ZELLE GELEGT
C
C   ORIGINAL:   17. 6. 93 (MM)
C
C*MGLET*****************************************************************
C
C  PARAMETER
C             IIF            - ANZAHL DER FEINGITTERPUNKTE
C              XF            - KOORDINATEN DES FEINGITTERS
C             DXF            - ABSTAND DER FEINGITTERPUNKTE
C            DDXF            - KANTENLAENGE DES FEINGITTERS
C
C             IIC            - ANZAHL DER GROBGITTERPUNKTE
C              XC            - KOORDINATEN DES GROBGITTERS
C             DXC            - ABSTAND DER GROBGITTERPUNKTE
C            DDXC            - KANTENLAENGE DES GROBGITTERS
C
C************************************************************************
C
C
      REAL XF(IIF),DXF(IIF),DDXF(IIF)
      REAL XC(IIC),DXC(IIC),DDXC(IIC)

      CHARACTER (LEN=1) CS,CSI
C
C-------------------------------------------------- CHECK AUF KONSISTENZ
C
      IF ( IIC-4 .NE. (IIF-4)/2 ) CALL ERRR (501,'GRDFTOC')
C
C----------------------------------------------- SETZEN DER GITTERPUNKTE
C
      DO IC = 2,IIC-1

         IF = (IC-2)*2 + 1

         XC( IC ) = XF( IF ) + DXF( IF )*0.5
       DDXC( IC ) = DDXF(IF) + DDXF(IF+1)

      ENDDO

C
C----------------------------------------------- UEBERKRAGENDE ZELLEN
C

      XC( 1 ) = XC( 2 )   - ( XC( 3   ) - XC( 2   ) )
      XC(IIC) = XC(IIC-1) + ( XC(IIC-1) - XC(IIC-2) )

      DDXC( 1 ) = DDXC( 2 )
      DDXC(IIC) = DDXC(IIC-1)
C
C----------------------------------------------- BERECHNUG DER ABSTAENDE
C

      DO I = 1,IIC-1

         DXC( I ) = XC(I+1) - XC( I )

      ENDDO

         DXC(IIC) = DXC(IIC-1)

C
C------------------------------------------------ LAENGE DER ZELLEN
C

C     DO I = 2,IIC-1

C        DDXC( I ) = 0.5*( DXC(I-1) + DXC( I ) )

C     ENDDO

C        DDXC( 1 ) = DDXC( 2 )
C        DDXC(IIC) = DDXC(IIC-1)

C
C------------------------------------------------------------
C
C
C                                 AUSGABE DER ERGEBNISSE
C
      WRITE (6,6110) CSI, CS, CSI, CS, CSI, CS, CSI, CS, CSI, CSI
      DO 350 I  = 1,IIC
         IDEL1I = (1/I) * (I/1)
         FAKTI  = DXC (I) / DXC (I-1+IDEL1I)
         SSTAG  = XC(I) + 0.5*DXC(I)
  350    WRITE (6,6120) I, XC(I), SSTAG, DXC(I), DDXC(I), FAKTI
C

 6110 FORMAT (/,'     ',A1,'         ',A1,'(',A1,')       ',A1,
     $        '-STAG (',A1,')       D',A1,'(',A1,')         DD',
     $        A1,'(',A1,')        F(',A1,')',/)

 6120 FORMAT (3X,I4,5(4X,F10.5))

C
C----------------------------------------------------------------
C
      RETURN
      END


