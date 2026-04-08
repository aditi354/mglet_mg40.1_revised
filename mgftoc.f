










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
      SUBROUTINE MGFTOC (KKF,JJF,IIF,DXF,DYF,DZF,DDXF,DDYF,DDZF,FF,
     &                   KKC,JJC,IIC,DXC,DYC,DZC,DDXC,DDYC,DDZC,FC,
     &                   IPOS,JPOS,KPOS,CID,IPROCF,IPROCC,
     &                      NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,BP)
C
C*MGLET*****************************************************************
C  M G F T O C         RESTRIKTION
C                      VON WAHLWEISE U,V,W ODER P(DP,G USW.)
C                      FEINES GITTER AUF GROBES GITTER
C
C   ORIGINAL:   25. 6. 93 (MM)
C               18. 7. 95 (MM) MPI EINGEFUEHRT
C               12. 6. 96 (MM) _BTOPAR_FIXED_ EINGEFUEHRT
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
C
C***********************************************************************
C
      REAL  DXF(IIF), DYF(JJF), DZF(KKF),
     &     DDXF(IIF),DDYF(JJF),DDZF(KKF),
     &      DXC(IIC), DYC(JJC), DZC(KKC),
     &     DDXC(IIC),DDYC(JJC),DDZC(KKC)

      REAL   FF(KKF,JJF,IIF),
     &       FC(KKC,JJC,IIC),BP(KKF,JJF,IIF)

      CHARACTER (LEN=1) CID
C
C
C-------------------------------------------------- CHECK AUF KONSISTENZ
C
CTEST      IF ( KKC-4 .NE. (KKF-4)/2 ) CALL ERRR (501,' MGFTOC')
CTEST      IF ( JJC-4 .NE. (JJF-4)/2 ) CALL ERRR (502,' MGFTOC')
CTEST      IF ( IIC-4 .NE. (IIF-4)/2 ) CALL ERRR (503,' MGFTOC')
C
C-------------------------------------------------- ANPASSUNG DER RANDBED.
C
      NTW = 0
      IF ((NTOP.EQ.1).OR.(NTOP.EQ.7).OR.(NTOP.EQ.3))  NTW = 1
      IF (NTOP.EQ.8)  NTW = 1
C
      NBU = 0
      IF ((NBAC.EQ.1).OR.(NBAC.EQ.7).OR.(NBAC.EQ.3))  NBU = 1
      IF (NBAC.EQ.8)  NBU = 1
C
C
C
C----------------------------------------------- GROESSE NICHT VERSETZT
C                                                VOLUMENMITTELWERT
C
      IF ( CID .EQ. 'P' ) THEN

*poption parallel
         DO IF = 3,IIF-2,2
            IC = IPOS + (IF-3)/2
         DO JF = 3,JJF-2,2
            JC = JPOS + (JF-3)/2
         DO KF = 3,KKF-2,2
            KC = KPOS + (KF-3)/2


            FC (KC  ,JC  ,IC  ) = 1./(DDXC(IC)*DDYC(JC)*DDZC(KC))

     &      *(   FF(KF  ,JF  ,IF  )*DDXF(IF  )*DDYF(JF  )*DDZF(KF  )
     &         + FF(KF  ,JF  ,IF+1)*DDXF(IF+1)*DDYF(JF  )*DDZF(KF  )
     &         + FF(KF  ,JF+1,IF  )*DDXF(IF  )*DDYF(JF+1)*DDZF(KF  )
     &         + FF(KF  ,JF+1,IF+1)*DDXF(IF+1)*DDYF(JF+1)*DDZF(KF  )
     &         + FF(KF+1,JF  ,IF  )*DDXF(IF  )*DDYF(JF  )*DDZF(KF+1)
     &         + FF(KF+1,JF  ,IF+1)*DDXF(IF+1)*DDYF(JF  )*DDZF(KF+1)
     &         + FF(KF+1,JF+1,IF  )*DDXF(IF  )*DDYF(JF+1)*DDZF(KF+1)
     &         + FF(KF+1,JF+1,IF+1)*DDXF(IF+1)*DDYF(JF+1)*DDZF(KF+1))
         ENDDO
         ENDDO
         ENDDO
C
C-------------------------------------- GROESSE IN X-RICHTUNG  VERSETZT
C                                       MITTELWERT UEBER OBERFLAECHE
C
      ELSEIF ( CID .EQ. 'U' ) THEN

*poption parallel
         DO IF = 2,IIF-3+NBU,2
            IC = IPOS + (IF-2)/2 - 1
         DO JF = 3,JJF-2,2
            JC = JPOS + (JF-3)/2
         DO KF = 3,KKF-2,2
            KC = KPOS + (KF-3)/2


            FC (KC  ,JC  ,IC  ) = 1./(         DDYC(JC)*DDZC(KC))

     &      *(   FF(KF  ,JF  ,IF  )*BP(KF  ,JF  ,IF  )
     &                                  *DDYF(JF  )*DDZF(KF  )
     &         + FF(KF  ,JF+1,IF  )*BP(KF  ,JF+1,IF  )
     &                                  *DDYF(JF+1)*DDZF(KF  )
     &         + FF(KF+1,JF  ,IF  )*BP(KF+1,JF  ,IF  )
     &                                  *DDYF(JF  )*DDZF(KF+1)
     &         + FF(KF+1,JF+1,IF  )*BP(KF+1,JF+1,IF  )
     &                                  *DDYF(JF+1)*DDZF(KF+1))

         ENDDO
         ENDDO
         ENDDO
C
C-------------------------------------- GROESSE IN Y-RICHTUNG  VERSETZT
C                                       MITTELWERT UEBER OBERFLAECHE
C
      ELSEIF ( CID .EQ. 'V' ) THEN

*poption parallel
         DO IF = 3,IIF-2,2
            IC = IPOS + (IF-3)/2
         DO JF = 2,JJF-2,2
            JC = JPOS + (JF-2)/2 - 1
         DO KF = 3,KKF-2,2
            KC = KPOS + (KF-3)/2


            FC (KC  ,JC  ,IC  ) = 1./(DDXC(IC)         *DDZC(KC))

     &      *(   FF(KF  ,JF  ,IF  )*BP(KF  ,JF  ,IF  )
     &                             *DDXF(IF  )           *DDZF(KF  )
     &         + FF(KF  ,JF  ,IF+1)*BP(KF  ,JF  ,IF+1)
     &                             *DDXF(IF+1)           *DDZF(KF  )
     &         + FF(KF+1,JF  ,IF  )*BP(KF+1,JF  ,IF  )
     &                             *DDXF(IF  )           *DDZF(KF+1)
     &         + FF(KF+1,JF  ,IF+1)*BP(KF+1,JF  ,IF+1)
     &                             *DDXF(IF+1)           *DDZF(KF+1))

         ENDDO
         ENDDO
         ENDDO
C
C-------------------------------------- GROESSE IN Z-RICHTUNG  VERSETZT
C                                       MITTELWERT UEBER OBERFLAECHE
C
      ELSEIF ( CID .EQ. 'W' ) THEN

*poption parallel
         DO IF = 3,IIF-2,2
            IC = IPOS + (IF-3)/2
         DO JF = 3,JJF-2,2
            JC = JPOS + (JF-3)/2
         DO KF = 2,KKF-3+NTW,2
            KC = KPOS + (KF-2)/2 - 1


            FC (KC  ,JC  ,IC  ) = 1./(DDXC(IC)*DDYC(JC)         )

     &      *(   FF(KF  ,JF  ,IF  )*BP(KF  ,JF  ,IF  )
     &                             *DDXF(IF  )*DDYF(JF  )
     &         + FF(KF  ,JF  ,IF+1)*BP(KF  ,JF  ,IF+1)
     &                             *DDXF(IF+1)*DDYF(JF  )
     &         + FF(KF  ,JF+1,IF  )*BP(KF  ,JF+1,IF  )
     &                             *DDXF(IF  )*DDYF(JF+1)
     &         + FF(KF  ,JF+1,IF+1)*BP(KF  ,JF+1,IF+1)
     &                             *DDXF(IF+1)*DDYF(JF+1)           )

         ENDDO
         ENDDO
         ENDDO
      ELSEIF ( CID .EQ. 'A' ) THEN

*poption parallel
         DO IF = 3,IIF-2,2
            IC = IPOS + (IF-3)/2
         DO JF = 3,JJF-2,2
            JC = JPOS + (JF-3)/2
         DO KF = 3,KKF-2,2
            KC = KPOS + (KF-3)/2


            FC (KC  ,JC  ,IC  ) =
     &      MIN(1.0,
     &      (   FF(KF  ,JF  ,IF  )+
     &          FF(KF  ,JF  ,IF+1)+
     &          FF(KF  ,JF+1,IF  )+
     &          FF(KF  ,JF+1,IF+1)+
     &          FF(KF+1,JF  ,IF  )+
     &          FF(KF+1,JF  ,IF+1)+
     &          FF(KF+1,JF+1,IF  )+
     &          FF(KF+1,JF+1,IF+1)  ))
        END DO
        END DO
        END DO
      ELSEIF ( CID .EQ. 'B' ) THEN


*poption parallel
         DO IF = 2,IIF-3+NBU,2
            IC = IPOS + (IF-2)/2 - 1
         DO JF = 3,JJF-2,2
            JC = JPOS + (JF-3)/2
         DO KF = 3,KKF-2,2
            KC = KPOS + (KF-3)/2

            FC (KC  ,JC  ,IC  ) =
     &      MIN(1.0,
     &      (   FF(KF  ,JF  ,IF  )+
     &          FF(KF  ,JF+1,IF  )+
     &          FF(KF+1,JF  ,IF  )+
     &          FF(KF+1,JF+1,IF  )  ))
        END DO
        END DO
        END DO
      ELSEIF ( CID .EQ. 'C' ) THEN
*poption parallel
         DO IF = 3,IIF-2,2
            IC = IPOS + (IF-3)/2
         DO JF = 2,JJF-2,2
            JC = JPOS + (JF-2)/2 - 1
         DO KF = 3,KKF-2,2
            KC = KPOS + (KF-3)/2


            FC (KC  ,JC  ,IC  ) =
     &      MIN(1.0,
     &      (   FF(KF  ,JF  ,IF  )+
     &          FF(KF  ,JF  ,IF+1)+
     &          FF(KF+1,JF  ,IF  )+
     &          FF(KF+1,JF  ,IF+1)))
        END DO
        END DO
        END DO
      ELSEIF ( CID .EQ. 'D' ) THEN

*poption parallel
         DO IF = 3,IIF-2,2
            IC = IPOS + (IF-3)/2
         DO JF = 3,JJF-2,2
            JC = JPOS + (JF-3)/2
         DO KF = 2,KKF-3+NTW,2
            KC = KPOS + (KF-2)/2 - 1

            FC (KC  ,JC  ,IC  ) =
     &      MIN(1.0,
     &      (   FF(KF  ,JF  ,IF  )+
     &          FF(KF  ,JF  ,IF+1)+
     &          FF(KF  ,JF+1,IF  )+
     &          FF(KF  ,JF+1,IF+1)))


         ENDDO
         ENDDO
         ENDDO
      ELSEIF ( CID .EQ. 'E' ) THEN

*poption parallel
         DO IF = 3,IIF-2,2
            IC = IPOS + (IF-3)/2
         DO JF = 3,JJF-2,2
            JC = JPOS + (JF-3)/2
         DO KF = 3,KKF-2,2
            KC = KPOS + (KF-3)/2


            FC (KC  ,JC  ,IC  ) =
     &      (   FF(KF  ,JF  ,IF  )+
     &          FF(KF  ,JF  ,IF+1)+
     &          FF(KF  ,JF+1,IF  )+
     &          FF(KF  ,JF+1,IF+1)+
     &          FF(KF+1,JF  ,IF  )+
     &          FF(KF+1,JF  ,IF+1)+
     &          FF(KF+1,JF+1,IF  )+
     &          FF(KF+1,JF+1,IF+1)  )
        END DO
        END DO
        END DO

C
C-------------------------------------- IM FEHLERFALL
C
      ELSE

         CALL ERRR (504,' MGFTOC')

      ENDIF
C
C----------------------------------------------------------------
C

      RETURN
      END


