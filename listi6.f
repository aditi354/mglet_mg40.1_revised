










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
      SUBROUTINE LISTI6 (KK, JJ, II, CID, IID, RID, NOUT)
C**********************************************************************
C         L I S T I 6  LISTET IDENTS DER STAR-DATEN
C**********************************************************************
C
C  PARAMETER: KK,JJ,II  - ARRAYDIMENSIONEN
C             CID(10)   - CHARACTER (LEN=8) VARIABLE
C             IID(100)  - INTEGERFELD
C             RID(100)  - INTEGERFELD
C             NOUT      - KANALNUMMER FUER OUTPUT
C
C UPROG                 : KEINE
C
C DEFINE-DIREKTIVEN     : KEINE
C
C VERS:  20.06.89 (HW)  : LISTI6 AUS LISTID (MSTARP5) ABGELEITET
C
C*STARLET**************************************************************
C
      CHARACTER (LEN=8) CID(10)
      INTEGER      IID(100)
      REAL         RID(100)
C
      WRITE(NOUT,6000)
      WRITE(NOUT,6010)IID(3),CID(4),CID(5)
      WRITE(NOUT,6020)CID(2)
      WRITE(NOUT,6080)
      WRITE(NOUT,6030)'ZTOT  = ',RID(10),'YTOT  = ',RID(11),
     $                'XTOT   = ', RID(12)
      WRITE(NOUT,6040)'KK_PHY= ',KK-4,'JJ_PHY= ',JJ-4,
     $     'II_PHY = ',II-4
      WRITE(NOUT,6100)
      WRITE(NOUT,6050)IID(10),IID(11),IID(12),IID(13),IID(14)
      WRITE(NOUT,6110)
      WRITE(NOUT,6040)'MTSTEP= ',IID(1),'ITTOT = ',IID(30),
     *                'NPRTOT = ',IID(31)
      WRITE(NOUT,6060)RID(2),RID(30)
      WRITE(NOUT,6070)IID(2),RID(3)
      WRITE(NOUT,6030)'RHO   = ',RID(17),'GMOL  = ',RID(18),
     $                'GRADPX = ',RID(16)
      WRITE(NOUT,6090)
      WRITE(NOUT,6030)'U-REF = ',RID(75),'L-REF = ',RID(76),
     $                'T-REF  = ',RID(77)
      WRITE(NOUT,6120)
      RETURN
C
 6000 FORMAT(/,4X,10(1H*),'  EINIGE INFORMATIONEN ZU BEGINN DES',
     $       ' LAUFES (SUBR. LISTI6)',/)
 6010 FORMAT(4X,'NRRUN = ',I6,9X,'DATUM = ',A8,6X,'UHRZEIT= ',A8)
 6020 FORMAT(4X,'NAME  = ',A8)
 6030 FORMAT(4X,A8,1PE12.5,3X,A8,1PE12.5,3X,A9,1PE12.5)
 6040 FORMAT(4X,A8,I6,9X,A8,I6,9X,A9,I6)
 6050 FORMAT(4X,'KB    = ',I4,4X,'JB1 = ',I4,3X,'JB2 = ',I4,3X,'IB1 = '
     $       ,I4,3X,'IB2 = ',I4)
 6060 FORMAT(4X,'DT    = ',1PE12.5,3X,'T_GES = ',1PE12.5)
 6070 FORMAT(4X,'MPCORR= ',I6,9X,'EPCORR= ',1PE12.5)
 6080 FORMAT(/,4X,'ABMESSUNGEN DES BERECHNUNGSGEBIETES :',/)
 6090 FORMAT(/,4X,'BEZUGSGROESSEN :',/)
 6100 FORMAT(/,4X,'LAGE DES KUBUSSES :',/)
 6110 FORMAT(/,4X,'WESENTLICHE PARAMETER DES LAUFES :',/)
 6120 FORMAT(/)
      END
