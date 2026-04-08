










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
      SUBROUTINE DISSIPG (KK,JJ,II,KMX,JMX,IMX,DDX,DDY,DDZ,
     $                    DX,DY,DZ,X,Y,Z,UFG,VFG,WFG,B,DFG,DSUM,
     $                    NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,NCUB,ISUM)
C*MGLET*****************************************************************
C        D I S S I P G    IN DISSIPG WIRD DIE DISSIPATION DER GROBSTRUKTUR
C                         AUS DEN MOMENTAN VORHANDENEN FLUKTUATIONEN
C                         DES GESCHWINDIGKEITSFELDES GEBILDET.
C*MGLET*****************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        DDX,DDY,DDZ    - ABMESSUNGEN DER BASISZELLEN
C        DX,DY,DZ       - ABSTAND DER BASISZELLMITTELPUNKTE
C        X,Y,Z          - KOORDINATEN DER ZELLMITTELPUNKTE
C        UFG(KK,JJ,II)  - FLUKTUATIONEN DER U-KOMPONENTE
C                         UFG = U - <U>
C        VFG(KK,JJ,II)  - FLUKTUATIONEN DER V-KOMPONENTE
C                         VFG = V - <V>
C        WFG(KK,JJ,II)  - FLUKTUATIONEN DER W-KOMPONENTE
C                         WFG = W - <W>
C        B  (KK,JJ,II)  - WANDERKENNUNGSFELD
C        DFG(KK,JJ,II)  + ENTHAELT DIE DISSIPATION, DEFINITIONSPUNKT IST
C                         DER DRUCKPUNKT
C
C        ISUM           - SCHALTER, OB GESAMTENERGIE BERECHNET WERDEN
C                         SOLL:    0 --> NICHT BERECHNET
C                                  1 --> SUMME WIRD GEBILDET
C
C
C UPROG                 : KEINE
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        06.05.94 (MM)  : ORIGINAL
C        10.05.94 (MM)  : WANDBEHANDLUNG EINGEFUEHRT, NUR FUER DNS,
C                         NICHT FUER LES GUELTIG!!!!!
C
C*MGLET*****************************************************************
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
      REAL     UFG(KK,JJ,II),  VFG(KK,JJ,II),  WFG(KK,JJ,II),
     $           B(KK,JJ,II),  DFG(KK,JJ,II)
C
      REAL     DDX(II),   DDY(JJ),  DDZ(KK),
     $          DX(II),    DY(JJ),   DZ(KK),
     $           X(II),     Y(JJ),    Z(KK)
C
C                                 KONSTANTEN
C
      IM2    = IMX - 2
      JM2    = JMX - 2
      KM2    = KMX - 2
      IFRFIX = 0
      JRIFIX = 0
      KBOFIX = 0

      IF ( NFRO .EQ. 2) IFRFIX = 1
      IF ( NRGT .EQ. 2) JRIFIX = 1
      IF ( NBOT .EQ. 2) KBOFIX = 1
C
      DO 100 I = 3,IM2
         DDXI   = DDX(I)
         RDDX   = 1.0/DDXI
C                                             AN DER WAND AENDERUNG DES
C                                             ABSTANDES

         DXF    = 0.5*DX(I-1)*RDDX
         DXF2   = 2.0 *DXF

         DO 110 J = 3,JM2
            DDYJ   = DDY(J)
            RDDY   = 1.0/DDYJ
C           RDDYPL = 0.5/DDYJ
            DYF    = 0.5*DY(J-1)*RDDY
            DYF2   = 2.0 *DYF

            DO 120 K = 3,KM2
               DDZK   = DDZ(K)
               RDDZ   = 1.0/DDZK
C              RDDZPL = 0.5/DDZK
               DZF    = 0.5*DZ(K-1)*RDDZ
               DZF2   = 2.0 *DZF

               RDDXPL = 1.0/(DDXI
     $                      + DX(I+1)*(SIGN(0.25,B(K,J,I+1))+0.25)
     $                      + DX(I-1)*(SIGN(0.25,B(K,J,I-1))+0.25))

               RDDYPL = 1.0/(DDYJ
     $                      + DY(J+1)*(SIGN(0.25,B(K,J+1,I))+0.25)
     $                      + DY(J-1)*(SIGN(0.25,B(K,J-1,I))+0.25))

               RDDZPL = 1.0/(DDZK
     $                      + DZ(K  )*(SIGN(0.25,B(K+1,J,I))+0.25)
     $                      + DZ(K-1)*(SIGN(0.25,B(K-1,J,I))+0.25))

               DUDX = RDDX * (UFG(K  ,J  ,I  )-UFG(K  ,J  ,I-1))

               DUDY = RDDYPL *
     $               (((UFG(K  ,J+1,I  )-UFG(K  ,J-1,I  )) *      DXF
     $              +  (UFG(K  ,J+1,I-1)-UFG(K  ,J-1,I-1)) * (1.0-DXF)))

               DUDZ = RDDZPL *
     $               (((UFG(K+1,J  ,I  )-UFG(K-1,J  ,I  )) *      DXF
     $              +  (UFG(K+1,J  ,I-1)-UFG(K-1,J  ,I-1)) * (1.0-DXF)))

               DVDX = RDDXPL *
     $               (((VFG(K  ,J  ,I+1)-VFG(K  ,J  ,I-1)) *      DYF
     $              +  (VFG(K  ,J-1,I+1)-VFG(K  ,J-1,I-1)) * (1.0-DYF)))

               DVDY = RDDY   * (VFG(K  ,J  ,I  )-VFG(K  ,J-1,I  ))

               DVDZ = RDDZPL *
     $               (((VFG(K+1,J  ,I  )-VFG(K-1,J  ,I  )) *      DYF
     $              +  (VFG(K+1,J-1,I  )-VFG(K-1,J-1,I  )) * (1.0-DYF)))

               DWDX = RDDXPL *
     $               (((WFG(K  ,J  ,I+1)-WFG(K  ,J  ,I-1)) *      DZF
     $              +  (WFG(K-1,J  ,I+1)-WFG(K-1,J  ,I-1)) * (1.0-DZF)))

               DWDY = RDDYPL *
     $               (((WFG(K  ,J+1,I  )-WFG(K  ,J-1,I  )) *      DZF
     $              +  (WFG(K-1,J+1,I  )-WFG(K-1,J-1,I  )) * (1.0-DZF)))

               DWDZ = RDDZ   * (WFG(K,J,I) - WFG(K-1,J,I))

               DFG(K,J,I) = GMOL*
     $              ((2.0*(DUDX*DUDX    + DVDY*DVDY    + DWDZ*DWDZ  )
     $              +     (DUDY+DVDX)*(DUDY+DVDX)
     $              +     (DUDZ+DWDX)*(DUDZ+DWDX)
     $              +     (DVDZ+DWDY)*(DVDZ+DWDY)))


  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      IF (ISUM.EQ.1) THEN

        DO I = 3,IM2
        DO J = 3,JM2
        DO K = 3,KM2

           DSUM = DSUM + DFG(K,J,I)*DDX(I)*DDY(J)*DDZ(K)
  
        ENDDO
        ENDDO
        ENDDO

           XL = X(  IM2   ) + DX(  IM2     )*0.5
     $        -(X(3-IFRFIX) - DX(3-IFRFIX-1)*0.5)

           YL = Y(  JM2   ) + DY(  JM2     )*0.5
     $        -(Y(3-JRIFIX) - DY(3-JRIFIX-1)*0.5)

           ZL = Z(  KM2   ) + DZ(  KM2     )*0.5
     $        -(Z(3-KBOFIX) - DZ(3-KBOFIX-1)*0.5)

           DSUM = DSUM / (XL*YL*ZL)

      ENDIF

      RETURN
      END
