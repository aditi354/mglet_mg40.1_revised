










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
      SUBROUTINE HELICI  (KK,JJ,II,KMX,JMX,IMX,U,V,W,OX,OY,OZ,HE)
C*STARLET***************************************************************
C        H E L I C I      BERECHNUNG DER HELICITY:
C                         HELICITY = 1/2 * < U_I * OMEGA_I >;
C                         SUMMATIONSKONVENTION GILT !
C                         KUBUSFOERMIGE EINBAUTEN WERDEN BERUECKSICHTIGT
C*STARLET***************************************************************
C
C PARAM: KK, JJ, II     - ARRAYDIMENSIONEN
C        KMX, JMX, IMX  - GRENZEN DES BERECHNUNGSGEBIETES (MIT BOUND)
C        U(KK,JJ,II)    - X-KOMPONENTE DES GESCHWINDIGKEITSFELDES
C        V(KK,JJ,II)    - Y-KOMPONENTE DES GESCHWINDIGKEITSFELDES
C        W(KK,JJ,II)    - Z-KOMPONENTE DES GESCHWINDIGKEITSFELDES
C        OX(KK,JJ,II)   - X-KOMPONENTE DER VORTICITY
C        OY(KK,JJ,II)   - Y-KOMPONENTE DER VORTICITY
C        OZ(KK,JJ,II)   - Z-KOMPONENTE DER VORTICITY
C        HE(KK,JJ,II)   + HELICITY
C
C UPROG                 : KEINE
C
C DEFINE-DIREKTIVEN     : KEINE
C
C VERS:  29.06.88 (HW)  : ORIGINAL
C        03.01.89 (HW)  : DIE HELIZITAET IST JETZT DREIFACH IM MASCHEN-
C                         GITTER VERSCHOBEN (IN X-, Y- UND Z-RICHTUNG)
C                         VOELLIG NEUE REVISION !
C
C
C*STARLET***************************************************************
C
      REAL        U (KK,JJ,II),   V (KK,JJ,II),   W (KK,JJ,II),
     $           OX (KK,JJ,II),  OY (KK,JJ,II),  OZ (KK,JJ,II),
     $           HE (KK,JJ,II)
C
C
      KM2    = KMX-2
      JM2    = JMX-2
      IM2    = IMX-2
C
      DO 100 I = 2,IM2
         IKST   =      MAX0(I      ,3)
         IKSTPS = MIN0(MAX0(I+1    ,3),IM2)
         DO 100 J = 2,JM2
            JKST   =      MAX0(J      ,3)
            JKSTPS = MIN0(MAX0(J+1    ,3),JM2)
C
            KSTART = 3
C
            DO 100 K = KSTART,KM2
               HE (K,J,I) = 0.5 * (
     $                      (0.5 *(OX(K  ,J  ,I  ) + OX(K  ,J  ,I+1)))
     $                    * (0.25*(U (K  ,J  ,I  ) + U (K+1,J  ,I  )
     $                    +        U (K  ,J+1,I  ) + U (K+1,J+1,I  )))
     $                    + (0.5 *(OY(K  ,J  ,I  ) + OY(K  ,J+1,I  )))
     $                    * (0.25*(V (K  ,J  ,I  ) + V (K+1,J  ,I  )
     $                    +        V (K  ,J  ,I+1) + V (K+1,J  ,I+1)))
     $                    + (0.5 *(OZ(K  ,J  ,I  ) + OZ(K+1,J  ,I  )))
     $                    * (0.25*(W (K  ,J  ,I  ) + W (K  ,J+1,I  )
     $                    +        W (K  ,J  ,I+1) + W (K  ,J+1,I+1)))
     $                    )
  100 CONTINUE
C
      RETURN
      END
