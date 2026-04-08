










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
      SUBROUTINE DISSIPGA (KKA,JJA,IIA,KMX,JMX,IMX,DDX,DDY,DDZ,
     $                    DX,DY,DZ,X,Y,Z,UAV,VAV,WAV,B,DFG,
     $                    XHOMOG,YHOMOG,ZHOMOG,KK,JJ,II)
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
C        01.10.03 (FS)  : ORIGINAL AUS DISSIPG UEBERNOMMEN
C                       : BERECHNUNG VON <SIJ><SIJ> G aus den mittleren
C                       : Geschwindigkeiten fuer belibige Kombination
C                       : von homogenen Richtngen
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

      INTEGER  KK,JJ,II, HOI,HOJ,HOK,HOJJ,HOII,HOKK

      REAL     UAV(KKA,JJA,IIA),VAV(KKA,JJA,IIA),WAV(KKA,JJA,IIA),
     $         B(KK,JJ,II),DFG(KKA,JJA,IIA)

      LOGICAL XHOMOG,YHOMOG,ZHOMOG
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

C Geometriefaktoren haben immer die dimensionierung KK,JJ,II
C Die Ergebnisfelder hingegen KKA,JJA,IIA. HO* steuern indiziierung
C wenn homogene richtungen Vorliegen (-> aequidistantes gitter) 
C Geometriefaktoren immer die gleichen. 
      HOI = 0
      HOJ = 0
      HOK = 0
      HOJJ = 1
      HOII = 1
      HOKK = 1
      IF (XHOMOG) THEN
         IA = 1
         IE = 1
         HOI = 2
         HOII = 0
      ELSE 
         IA = 3
         IE = IIA-2
      ENDIF
      IF (YHOMOG) THEN
         JA = 1
         JE = 1
         HOJ = 2
         HOJJ = 0
      ELSE 
         JA = 3
         JE = JJA-2
      ENDIF
      IF (ZHOMOG) THEN
         KA = 1
         KE = 1
         HOK = 2
         HOKK = 0
      ELSE 
         KA = 3
         KE = KKA-2
      ENDIF
      DO I=IA,IE
         IF(XHOMOG) THEN
            DDXI   = DDX(3)
            RDDX   = 1.0/DDXI
            DXF    = 0.5*DX(3-1)*RDDX
            DXF2   = 2.0 *DXF
         ELSE
            DDXI   = DDX(I)
            RDDX   = 1.0/DDXI
            DXF    = 0.5*DX(I-1)*RDDX
            DXF2   = 2.0 *DXF
         ENDIF
         DO J=JA,JE
            IF(YHOMOG) THEN
               DDYJ   = DDY(3)
               RDDY   = 1.0/DDYJ
               DYF    = 0.5*DY(3-1)*RDDY
               DYF2   = 2.0 *DYF
            ELSE
               DDYJ   = DDY(J)
               RDDY   = 1.0/DDYJ
               DYF    = 0.5*DY(J-1)*RDDY
               DYF2   = 2.0 *DYF
            ENDIF
            DO K=KA,KE
               IF(ZHOMOG) THEN
               DDZK   = DDZ(3)
               RDDZ   = 1.0/DDZK
               DZF    = 0.5*DZ(3-1)*RDDZ
               DZF2   = 2.0 *DZF
            ELSE
               DDZK   = DDZ(K)
               RDDZ   = 1.0/DDZK
               DZF    = 0.5*DZ(K-1)*RDDZ
               DZF2   = 2.0 *DZF
            ENDIF
            
               RDDXPL = 1.0/(DDXI
     $           + DX(I+HOI+1)*(SIGN(0.25,B(K+HOK,J+HOJ,I+HOI+1))+0.25)
     $           + DX(I+HOI-1)*(SIGN(0.25,B(K+HOK,J+HOJ,I+HOI-1))+0.25))

               RDDYPL = 1.0/(DDYJ
     $           + DY(J+HOJ+1)*(SIGN(0.25,B(K+HOK,J+HOJ+1,I+HOI))+0.25)
     $           + DY(J+HOJ-1)*(SIGN(0.25,B(K+HOK,J+HOJ-1,I+HOI))+0.25))

               RDDZPL = 1.0/(DDZK
     $           + DZ(K+HOK  )*(SIGN(0.25,B(K+HOK+1,J+HOJ,I+HOI))+0.25)
     $           + DZ(K+HOK-1)*(SIGN(0.25,B(K+HOK-1,J+HOJ,I+HOI))+0.25))

               IF(XHOMOG) THEN
                  DUDX = 0.0
                  DWDX = 0.0
                  DVDX = 0.0
               ELSE
                  DUDX = RDDX * (UAV(K  ,J  ,I  )-UAV(K  ,J  ,I-1))
                  DVDX = RDDXPL *
     $      (((VAV(K  ,J  ,I+1)-VAV(K,J  ,I-1))*DYF
     $     +  (VAV(K  ,J-1*HOJJ,I+1)-VAV(K,J-1*HOJJ,I-1))*(1.0-DYF)))
                  DWDX = RDDXPL *
     $      (((WAV(K  ,J  ,I+1)-WAV(K  ,J  ,I-1)) *      DZF
     $   +  (WAV(K-1*HOKK,J  ,I+1)-WAV(K-1*HOKK,J  ,I-1)) * (1.0-DZF)))
               ENDIF
               
               IF(YHOMOG) THEN
                  DUDY = 0.0
                  DVDY = 0.0
                  DWDY = 0.0
               ELSE
               DUDY = RDDYPL *
     $       (((UAV(K  ,J+1,I  )-UAV(K  ,J-1,I  )) *      DXF
     $   +  (UAV(K  ,J+1,I-1*HOII)-UAV(K  ,J-1,I-1*HOII)) * (1.0-DXF)))
               DVDY = RDDY   * (VAV(K  ,J  ,I  )-VAV(K  ,J-1,I  ))
               DWDY = RDDYPL *
     $               (((WAV(K  ,J+1,I  )-WAV(K  ,J-1,I  )) *      DZF
     $  +  (WAV(K-1*HOKK,J+1,I  )-WAV(K-1*HOKK,J-1,I  )) * (1.0-DZF)))
               ENDIF

               IF(ZHOMOG) THEN
                  DUDZ = 0.0
                  DVDZ = 0.0
                  DWDZ = 0.0
               ELSE
               DUDZ = RDDZPL *
     $     (((UAV(K+1,J  ,I  )-UAV(K-1,J  ,I  )) *      DXF
     $  +  (UAV(K+1,J  ,I-1*HOII)-UAV(K-1,J  ,I-1*HOII)) * (1.0-DXF)))
               DVDZ = RDDZPL *
     $     (((VAV(K+1,J  ,I  )-VAV(K-1,J  ,I  )) *      DYF
     $  +  (VAV(K+1,J-1*HOJJ,I  )-VAV(K-1,J-1*HOJJ,I  )) * (1.0-DYF)))

               DWDZ = RDDZ   * (WAV(K,J,I) - WAV(K-1,J,I))
               ENDIF
               DFG(K,J,I) = GMOL*
     $              ((2.0*(DUDX*DUDX    + DVDY*DVDY    + DWDZ*DWDZ  )
     $              +     (DUDY+DVDX)*(DUDY+DVDX)
     $              +     (DUDZ+DWDX)*(DUDZ+DWDX)
     $              +     (DVDZ+DWDY)*(DVDZ+DWDY)))
            ENDDO
         ENDDO
      ENDDO
C

      RETURN
      END
