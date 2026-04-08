










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
      SUBROUTINE MGBFTC (KKF,JJF,IIF,DXF,DYF,DZF,DDXF,DDYF,DDZF,FF,
     &                   KKC,JJC,IIC,DXC,DYC,DZC,DDXC,DDYC,DDZC,FC,
     &                   KPOS,JPOS,IPOS,CID,IPROCF,IPROCC)
C
C*MGLET*****************************************************************
C  M G B F T C         RESTRIKTION DER RANDBED.-BUFFER
C                      FEINES GITTER AUF GROBES GITTER
C
C   ORIGINAL:   25. 6. 93 (MM) AUS MGFTOC ABGELEITET
C                              NUR UFR IMPLEMENTIERT !!!!!
C               14. 7. 97 (MM) BOTTOM-BUFFER geht auch mit veraendertem
C                              Aufruf!
C               30.10. 97 (MM) MPI EINGEFUEHRT
C
C*MGLET*****************************************************************
C
C  PARAMETER
C             IIF,JJF,KKF    - ANZAHL DER FEINGITTERPUNKTE
C             DXF,DYF,DZF    - ABSTAND DER FEINGITTERPUNKTE
C            DDXF,DDYF,DDZF  - KANTENLAENGE DES FEINGITTERS
C             FF             - ZU RESTRINGIERENDE VARIABLE
C
C             IIC,JJC,KKC    - ANZAHL DER GROBGITTERPUNKTE
C             DXC,DYC,DZC    - ABSTAND DER GROBGITTERPUNKTE
C            DDXC,DDYC,DDZC  - KANTENLAENGE DES GROBGITTERS
C             FC             + RESTRINGIERTE VARIABLE
C
C             IPOS,JPOS,KPOS - "AUFHAENGEPUNKT" DES FEINEN GITTERS
C                              IM GROBEN GITTER
C
C             CID            - CHARACTER (LEN=1), GIBT AN WIE GROESSE
C                              IM MASCHENGITTER VERSETZT IST
C                              ("U","V","W" ODER "P")
C                               "U" BEDEUTET UFR
C
C***********************************************************************
C
      REAL  DXF(IIF), DYF(JJF), DZF(KKF),
     &     DDXF(IIF),DDYF(JJF),DDZF(KKF),
     &      DXC(IIC), DYC(JJC), DZC(KKC),
     &     DDXC(IIC),DDYC(JJC),DDZC(KKC)

      REAL   FF(KKF,JJF,IIF),
     &       FC(KKC,JJC,IIC)

      CHARACTER (LEN=1) CID
C
C
C
C-------------------------------------- GROESSE IN X-RICHTUNG  VERSETZT
C                                       MITTELWERT UEBER OBERFLAECHE
C
C                                       HIER: FRONT-BUFFER
C
C
      IF ( CID .EQ. 'U' ) THEN

C
C-------------------------------------------------- CHECK AUF KONSISTENZ
C
      IF ( KKC-4 .LT. (KKF-4)/2 ) CALL ERRR (502,' MGBFTC')
      IF ( JJC-4 .LT. (JJF-4)/2 ) CALL ERRR (503,' MGBFTC')
      IF ( IIC   .NE.     2     ) CALL ERRR (504,' MGBFTC')
      IF ( IIF   .NE.     2     ) CALL ERRR (505,' MGBFTC')

C301097         CALL SETS   (KKC,JJC, 2 ,KKC,JJC, 2 ,FC,0.0)

            IF = 2
            IC = 2
         DO JF = 3,JJF-2,2
            JC = JPOS + (JF-3)/2
         DO KF = 3,KKF-2,2
            KC = KPOS + (KF-3)/2


            FC (KC  ,JC  ,IC  ) = 1./(         DDYC(JC)*DDZC(KC))

     &      *(   FF(KF  ,JF  ,IF  )           *DDYF(JF  )*DDZF(KF  )
     &         + FF(KF  ,JF+1,IF  )           *DDYF(JF+1)*DDZF(KF  )
     &         + FF(KF+1,JF  ,IF  )           *DDYF(JF  )*DDZF(KF+1)
     &         + FF(KF+1,JF+1,IF  )           *DDYF(JF+1)*DDZF(KF+1))

         ENDDO
         ENDDO
C
C-------------------------------------- IM FEHLERFALL
C
      ELSE

         CALL ERRR (508,' MGFTOC')

      ENDIF
C
C
C
C----------------------------------------------------------------
C
      RETURN
      END


