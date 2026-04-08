










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
      SUBROUTINE GRDCTOF (IIC,XC,DXC,DDXC,IIF,XF,DXF,DDXF,IPOS,CS,CSI)
C
C*MGLET*****************************************************************
C  G R D C T O F       ERZEUGEN EINES GITTERS FUER MULTIGRID
C                      ERZEUGT FEINGITTER AUS GROBGITTER
C                      DIE GITTERPUNKTE WERDEN IN MITTE DER ZELLE GELEGT
C
C   ORIGINAL:   12.10.93 (MM)    AUS GRDFTOC ABGELEITET
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
      REAL XF(IIF),DXF(IIF),DDXF(IIF)
      REAL XC(IIC),DXC(IIC),DDXC(IIC)

      CHARACTER (LEN=1) CS,CSI
C
C-------------------------------------------------- CHECK AUF KONSISTENZ
C
      IF ( IIC-4 .LT. (IIF-4)/2 ) CALL ERRR (501,'GRDCTOF')
C
C----------------------------------------------- GITTERABSTAENDE
C
      DO IF = 1,IIF-1,2

         IC = IPOS + (IF-3)/2

         DXF (IF  ) =  DDXC(IC) * 0.5 
         DXF (IF+1) =  DXC (IC) * 0.5

         DDXF(IF  ) =  0.25 * ( DXC(IC-1) + DDXC(IC  ))
         DDXF(IF+1) =  0.25 * ( DXC(IC  ) + DDXC(IC  ))

         XCL = XC(IC) - 0.5*DXC(IC-1)

         XF(IF  ) = XCL + 0.5*(DXF(MAX(1,IF-1)))
	     XF(IF+1) = XF(IF) + DXF(IF)

      ENDDO

C
C----------------------------------------------- GITTERPUNKTE
C

CTEST         XF(1) = XC(IPOS-1) - 0.25 *  DXC(IPOS-2)

CTEST      DO IF = 2,IIF

CTEST         XF(IF) = XF(IF-1) + DXF(IF-1)

CTEST      ENDDO
C
C------------------------------------------------------------
C
C
C                                 AUSGABE DER ERGEBNISSE
C
C
      WRITE (6,6110) CSI, CS, CSI, CS, CSI, CS, CSI, CS, CSI, CSI
      DO 350 I  = 1,IIF
         IDEL1I = (1/I) * (I/1)
         FAKTI  = DXF (I) / DXF (I-1+IDEL1I)
         SSTAG  = XF(I) + 0.5*DXF(I)
  350    WRITE (6,6120) I, XF(I), SSTAG, DXF(I), DDXF(I), FAKTI
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


