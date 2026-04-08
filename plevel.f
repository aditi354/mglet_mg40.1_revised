










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
      SUBROUTINE PLEVEL  (KK,JJ,II,KMX,JMX,IMX,P,B,
     $                    DDX, DDY, DDZ ,NFRO,NRGT,NBOT)
C*STARLET***************************************************************
C        P L E V E L      DAS DRUCKNIVEAU WIRD IN DER WEISE REDUZIERT,
C                         DASS DER VOLUMENMITTELWERT DES DRUCKES (UEBER
C                         DAS GESAMTE BERECHNUNGSGEBIET) ZU NULL WIRD.
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        P(KK,JJ,II)    + DRUCKFELD
C        DDX (II)       - KANTENLAENGE DER KONTROLLVOLUMEN IN X-RI.
C        DDY (JJ)       - KANTENLAENGE DER KONTROLLVOLUMEN IN Y-RI.
C        DDZ (KK)       - KANTENLAENGE DER KONTROLLVOLUMEN IN Z-RI.
C
C UPROG                 : KEINE
C
C DEFINE-DIREKTIVEN     : FRPER, RIPER
C
C        26.05.86 (HW)  : ORIGINAL
C        17.10.89 (HW)  : PLEVEL1 AUS PLEVEL ABGELEITET.
C        30.10.90 (HW)  : VOLUMENGEWICHTETE DRUCKWERTE
C        12. 4.93 (MM)  : B-FELD EINGEFUEHRT, 
C                         UEBERGABE DER RANDBEDINGUNGEN IM KOPF
C
C*STARLET***************************************************************
C
C
      REAL     P(KK,JJ,II),B(KK,JJ,II),
     $           DDX (II),  DDY (JJ),  DDZ (KK)

C
C                                 START- UND STOPINDIZES
C
      IM2    = IMX - 2
      JM2    = JMX - 2
      KM2    = KMX - 2

      ISTART = 3
      ISTOP  = IMX - 2
      IF(NFRO.EQ.1) ISTOP  = IMX - 1
C
      JSTART = 3
      JSTOP  = JMX - 2
      IF(NRGT.EQ.1) JSTOP  = JMX - 1
C
C                                 VOLUMENMITTELWERT DES DRUCKES
C                                 VGL. DISSERTATION (ABSCHNITT 4.5.3)
C
      PSUM   = 0.0
      VOLUME = 0.0
      DO  I = 3,IM2
      DO  J = 3,JM2
      DO  K = 3,KM2
               VOLIJK = DDX (I) * DDY (J) * DDZ (K)
     $                * B(K,J,I)
               VOLUME = VOLUME + VOLIJK
               PSUM   = PSUM + P(K,J,I) * VOLIJK
      ENDDO
      ENDDO
      ENDDO
      PMEAN  = PSUM / VOLUME

C
C
      DO  I = 3,IM2
      DO  J = 3,JM2  
      DO  K = 3,KM2
                  P(K,J,I) = P(K,J,I) - PMEAN*B(K,J,I)
      ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
