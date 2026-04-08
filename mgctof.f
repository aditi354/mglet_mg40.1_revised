










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
      SUBROUTINE MGCTOF (KKF,JJF,IIF,DXF,DYF,DZF,DDXF,DDYF,DDZF,FF,
     &                   KKC,JJC,IIC,DXC,DYC,DZC,DDXC,DDYC,DDZC,FC,
     &                   IPOS,JPOS,KPOS,CID,IPROCF,IPROCC)
C
C*MGLET*****************************************************************
C  M G C T O F         PROLONGATION
C                      VON WAHLWEISE U,V,W ODER P(DP,G USW.)
C                      GROBES GITTER AUF FEINES GITTER
C                      IST NUR FUER KORREKTURDRUCK ERLAUBT
C                      RANDZELLEN DES FEINGITTERS WERDEN MITBELEGT
C                      DESHALB MUESSEN RANDBEDINGUNGEN AUF GROBGITTER
C                      GESETZT SEIN
C
C   ORIGINAL:   28. 6. 93 (MM) AUS MGFTOC ABGELEITET
C               18. 7. 95 (MM) MPI EINGEFUEHRT
C               ....2. 96 (MM) INTERPOLATION 2. ORDNUNG FUER P EINGEFUEHRT
C                              NEUER CID: "B", FUER GROESSEN,
C                              DIR NICHT UEBER RAND INTERPOLIERT WERDEN
C                              DUERFEN
C               06.02.03  (TB) SCALAR FIELD PROCESSING ADDED
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
C                              ("U","V","W" ODER "P", "B")
C
C***********************************************************************
C
      REAL  DXF(IIF), DYF(JJF), DZF(KKF),
     &     DDXF(IIF),DDYF(JJF),DDZF(KKF),
     &      DXC(IIC), DYC(JJC), DZC(KKC),
     &     DDXC(IIC),DDYC(JJC),DDZC(KKC)

      REAL   FF(KKF,JJF,IIF),
     &       FC(KKC,JJC,IIC)

      REAL I3DO2,I2DO2 

      REAL STEIGUNGY1,STEIGUNGY2,STEUGUNGZ1,STEIGUNGZ2

      CHARACTER (LEN=1) CID
C
C
C-------------------------------------------------- CHECK AUF KONSISTENZ
C
CTEST      IF ( KKC-4 .NE. (KKF-4)/2 ) CALL ERRR (501,' MGCTOF')
CTEST      IF ( JJC-4 .NE. (JJF-4)/2 ) CALL ERRR (502,' MGCTOF')
CTEST      IF ( IIC-4 .NE. (IIF-4)/2 ) CALL ERRR (503,' MGCTOF')
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                 TEST OF DEFINE-DIRECTIVES

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C----------------------------------------------- GROESSE NICHT VERSETZT
C                                                VOLUMENMITTELWERT
C                                                KEINE INTERPOLATION ERLAUBT
C
      IF ( CID .EQ. 'B' ) THEN

         DO IF = 1,IIF-1,2
            IC = IPOS + (IF-3)/2
         DO JF = 1,JJF-1,2
            JC = JPOS + (JF-3)/2
         DO KF = 1,KKF-1,2
            KC = KPOS + (KF-3)/2


            FF (KF  ,JF  ,IF  ) = FC (KC  ,JC  ,IC  )
            FF (KF  ,JF  ,IF+1) = FC (KC  ,JC  ,IC  )
            FF (KF  ,JF+1,IF  ) = FC (KC  ,JC  ,IC  )
            FF (KF  ,JF+1,IF+1) = FC (KC  ,JC  ,IC  )
            FF (KF+1,JF  ,IF  ) = FC (KC  ,JC  ,IC  )
            FF (KF+1,JF  ,IF+1) = FC (KC  ,JC  ,IC  )
            FF (KF+1,JF+1,IF  ) = FC (KC  ,JC  ,IC  )
            FF (KF+1,JF+1,IF+1) = FC (KC  ,JC  ,IC  )

         ENDDO
         ENDDO
         ENDDO
C
C----------------------------------------------- GROESSE NICHT VERSETZT
C                                                VOLUMENMITTELWERT
C
      ELSEIF ( CID .EQ. 'P' ) THEN

         DO IF = 1,IIF-1,2
            IC = IPOS + (IF-3)/2
         DO JF = 1,JJF-1,2
            JC = JPOS + (JF-3)/2
         DO KF = 1,KKF-1,2
            KC = KPOS + (KF-3)/2


            FF (KF  ,JF  ,IF  ) = FC (KC  ,JC  ,IC  )
            FF (KF  ,JF  ,IF+1) = FC (KC  ,JC  ,IC  )
            FF (KF  ,JF+1,IF  ) = FC (KC  ,JC  ,IC  )
            FF (KF  ,JF+1,IF+1) = FC (KC  ,JC  ,IC  )
            FF (KF+1,JF  ,IF  ) = FC (KC  ,JC  ,IC  )
            FF (KF+1,JF  ,IF+1) = FC (KC  ,JC  ,IC  )
            FF (KF+1,JF+1,IF  ) = FC (KC  ,JC  ,IC  )
            FF (KF+1,JF+1,IF+1) = FC (KC  ,JC  ,IC  )




         ENDDO
         ENDDO
         ENDDO
C
C-------------------------------------- GROESSE IN X-RICHTUNG  VERSETZT
C                                       MITTELWERT UEBER OBERFLAECHE
C
      ELSEIF ( CID .EQ. 'U' ) THEN

         DO IF = 1,IIF-1,2
            IC = IPOS + (IF-3)/2
         DO JF = 1,JJF-1,2
            JC = JPOS + (JF-3)/2
         DO KF = 1,KKF-1,2
            KC = KPOS + (KF-3)/2

C
C                                       FEIN- UND GROBGITTERWERTE AUF
C                                       DERSELBEN I-POSITION

            FF (KF  ,JF  ,IF+1) = FC (KC  ,JC  ,IC  )
            FF (KF  ,JF+1,IF+1) = FC (KC  ,JC  ,IC  )
            FF (KF+1,JF  ,IF+1) = FC (KC  ,JC  ,IC  )
            FF (KF+1,JF+1,IF+1) = FC (KC  ,JC  ,IC  )
C
C                                       I-POSITION DER
C                                       FEINGITTERWERTE ZWISCHEN ZWEI
C                                       GROBGITTERWERTEN

            FF (KF  ,JF  ,IF  ) = 0.5*( FC(KC,JC,IC) + FC(KC,JC,IC-1))
            FF (KF  ,JF+1,IF  ) = 0.5*( FC(KC,JC,IC) + FC(KC,JC,IC-1))
            FF (KF+1,JF  ,IF  ) = 0.5*( FC(KC,JC,IC) + FC(KC,JC,IC-1))
            FF (KF+1,JF+1,IF  ) = 0.5*( FC(KC,JC,IC) + FC(KC,JC,IC-1))


         ENDDO
         ENDDO
         ENDDO


C
C-------------------------------------- GROESSE IN Y-RICHTUNG  VERSETZT
C                                       MITTELWERT UEBER OBERFLAECHE
C
      ELSEIF ( CID .EQ. 'V' ) THEN


         DO IF = 1,IIF-1,2
            IC = IPOS + (IF-3)/2
         DO JF = 1,JJF-1,2
            JC = JPOS + (JF-3)/2
         DO KF = 1,KKF-1,2
            KC = KPOS + (KF-3)/2

C
C                                       FEIN- UND GROBGITTERWERTE AUF
C                                       DERSELBEN J-POSITION

            FF (KF  ,JF+1,IF  ) = FC (KC  ,JC  ,IC  )
            FF (KF  ,JF+1,IF+1) = FC (KC  ,JC  ,IC  )
            FF (KF+1,JF+1,IF  ) = FC (KC  ,JC  ,IC  )
            FF (KF+1,JF+1,IF+1) = FC (KC  ,JC  ,IC  )
C
C                                       J-POSITION DER
C                                       FEINGITTERWERTE ZWISCHEN ZWEI
C                                       GROBGITTERWERTEN

            FF (KF  ,JF  ,IF  ) = 0.5*( FC(KC,JC,IC) + FC(KC,JC-1,IC))
            FF (KF  ,JF  ,IF+1) = 0.5*( FC(KC,JC,IC) + FC(KC,JC-1,IC))
            FF (KF+1,JF  ,IF  ) = 0.5*( FC(KC,JC,IC) + FC(KC,JC-1,IC))
            FF (KF+1,JF  ,IF+1) = 0.5*( FC(KC,JC,IC) + FC(KC,JC-1,IC))


         ENDDO
         ENDDO
         ENDDO

C
C-------------------------------------- GROESSE IN Z-RICHTUNG  VERSETZT
C                                       MITTELWERT UEBER OBERFLAECHE
C
      ELSEIF ( CID .EQ. 'W' ) THEN

         DO IF = 1,IIF-1,2
            IC = IPOS + (IF-3)/2
         DO JF = 1,JJF-1,2
            JC = JPOS + (JF-3)/2
         DO KF = 1,KKF-1,2
            KC = KPOS + (KF-3)/2

C
C                                       FEIN- UND GROBGITTERWERTE AUF
C                                       DERSELBEN K-POSITION

            FF (KF+1,JF  ,IF  ) = FC (KC  ,JC  ,IC  )
            FF (KF+1,JF  ,IF+1) = FC (KC  ,JC  ,IC  )
            FF (KF+1,JF+1,IF  ) = FC (KC  ,JC  ,IC  )
            FF (KF+1,JF+1,IF+1) = FC (KC  ,JC  ,IC  )
C
C                                       K-POSITION DER
C                                       FEINGITTERWERTE ZWISCHEN ZWEI
C                                       GROBGITTERWERTEN

            FF (KF  ,JF  ,IF  ) = 0.5*( FC(KC,JC,IC) + FC(KC-1,JC,IC))
            FF (KF  ,JF  ,IF+1) = 0.5*( FC(KC,JC,IC) + FC(KC-1,JC,IC))
            FF (KF  ,JF+1,IF  ) = 0.5*( FC(KC,JC,IC) + FC(KC-1,JC,IC))
            FF (KF  ,JF+1,IF+1) = 0.5*( FC(KC,JC,IC) + FC(KC-1,JC,IC))


         ENDDO
         ENDDO
         ENDDO

C-------------------------------------- GROESSE IN X-RICHTUNG VERSETZT
C                                       2-D INTERPOLATION (F.S.)
      ELSEIF ( CID .EQ. 'X' ) THEN

         DO IF = 1,IIF-1,2
            IC = IPOS + (IF-3)/2
         DO JF = 1,JJF-1,2
            JC = JPOS + (JF-3)/2
         DO KF = 3,KKF-5,2
            KC = KPOS + (KF-3)/2
C------------- INTERPOLATION IN YZ-EBENE FUER FINE AUF GLEICHER X-KORD
C              WIE GROB IF+1 = IC
        FF(KF+1,JF+1,IF+1) = FC(KC,JC,IC) 
     $   +DYF(JF)/2.0*(FC(KC,JC+1,IC)-FC(KC,JC,IC))/DYC(JC) 
        FF(KF+1,JF+2,IF+1) = FC(KC,JC,IC) 
     $   +(DYF(JF)/2.0+DYF(JF+1))*(FC(KC,JC+1,IC)-FC(KC,JC,IC))/DYC(JC) 
        FF(KF+2,JF+1,IF+1) = FC(KC+1,JC,IC) 
     $   +DYF(JF)/2.0*(FC(KC+1,JC+1,IC)-FC(KC+1,JC,IC))/DYC(JC) 
        FF(KF+2,JF+2,IF+1) = FC(KC+1,JC,IC) 
     $   +(DYF(JF)/2.0+DYF(JF+1))
     $   *(FC(KC+1,JC+1,IC)-FC(KC+1,JC,IC))/DYC(JC) 

         STEIGUNGZ1 = (FF(KF+2,JF+1,IF+1)-FF(KF+1,JF+1,IF+1))/DZC(KC)
         STEIGUNGZ2 = (FF(KF+2,JF+2,IF+1)-FF(KF+1,JF+2,IF+1))/DZC(KC)

        FF(KF+2,JF+1,IF+1) = FF(KF+1,JF+1,IF+1) + 
     $    (DZF(KF)/2.0+DZF(KF+1))*STEIGUNGZ1
        FF(KF+1,JF+1,IF+1) = FF(KF+1,JF+1,IF+1) + 
     $    (DZF(KF)/2.0)*STEIGUNGZ1
        FF(KF+2,JF+2,IF+1) = FF(KF+1,JF+2,IF+1) + 
     $    (DZF(KF)/2.0+DZF(KF+1))*STEIGUNGZ2
        FF(KF+1,JF+2,IF+1) = FF(KF+1,JF+2,IF+1) + 
     $    (DZF(KF)/2.0)*STEIGUNGZ2

         ENDDO
         ENDDO
         ENDDO
C-------------- Randbeingung NO-SLIP WALL TOP AND BOTTOM
C-------------- Bottom Index KF TOP Index KF2
         KF = 1
         KC = KPOS + (KF-3)/2
         KF2 = KKF-3
         KC2 = KPOS + (KF2-3)/2
         DO IF = 1,IIF-1,2
            IC = IPOS + (IF-3)/2
         DO JF = 1,JJF-1,2
            JC = JPOS + (JF-3)/2
C-------- bottom
        FF(KF+2,JF+1,IF+1) = FC(KC+1,JC,IC)
     $   +DYF(JF)/2.0*(FC(KC+1,JC+1,IC)-FC(KC+1,JC,IC))/DYC(JC)
        FF(KF+2,JF+2,IF+1) = FC(KC+1,JC,IC)
     $   +(DYF(JF)/2.0+DYF(JF+1))
     $   *(FC(KC+1,JC+1,IC)-FC(KC+1,JC,IC))/DYC(JC)
         STEIGUNGZ1 = (FF(KF+2,JF+1,IF+1)+FF(KF+2,JF+1,IF+1))/DZC(KC)
         STEIGUNGZ2 = (FF(KF+2,JF+2,IF+1)+FF(KF+2,JF+2,IF+1))/DZC(KC)
         FF(KF+2,JF+1,IF+1) = DZF(KF+1)/2.0*STEIGUNGZ1
         FF(KF+2,JF+2,IF+1) = DZF(KF+1)/2.0*STEIGUNGZ2
C-------- top
        FF(KF2+1,JF+1,IF+1) = FC(KC2,JC,IC)
     $   +DYF(JF)/2.0*(FC(KC2,JC+1,IC)-FC(KC2,JC,IC))/DYC(JC)
        FF(KF2+1,JF+2,IF+1) = FC(KC2,JC,IC)
     $   +(DYF(JF)/2.0+DYF(JF+1))
     $   *(FC(KC2,JC+1,IC)-FC(KC2,JC,IC))/DYC(JC)
         STEIGUNGZ1 = (FF(KF2+1,JF+1,IF+1)+FF(KF2+1,JF+1,IF+1))/DZC(KC2)
         STEIGUNGZ2 = (FF(KF2+1,JF+2,IF+1)+FF(KF2+1,JF+2,IF+1))/DZC(KC2)
         FF(KF2+1,JF+1,IF+1) = DZF(KF2+1)/2.0*STEIGUNGZ1
         FF(KF2+1,JF+2,IF+1) = DZF(KF2+1)/2.0*STEIGUNGZ2
         ENDDO
         ENDDO
C-------------- INTERPOLATION IN X-RICHTUNG FUER DIE UEBRIGEN FINE-KOORD
         DO IF = 3,IIF-1,2
          DO JF = 3,JJF-2
           DO KF = 3,KKF-2
            FF(KF,JF,IF) = 0.5*(FF(KF,JF,IF+1)+FF(KF,JF,IF-1))
           ENDDO
          ENDDO
         ENDDO


C-------------------------------------- GROESSE IN Y-RICHTUNG VERSETZT
C                                       2-D INTERPOLATION (F.S.)

      ELSEIF ( CID .EQ. 'Y' ) THEN

         DO IF = 1,IIF-1,2
            IC = IPOS + (IF-3)/2
         DO JF = 1,JJF-1,2
            JC = JPOS + (JF-3)/2
         DO KF = 3,KKF-5,2
            KC = KPOS + (KF-3)/2
C------------- INTERPOLATION IN XZ-EBENE FUER FINE AUF GLEICHER Y-KORD
C              WIE GROB JF+1 = JC
        FF(KF+1,JF+1,IF+1) = FC(KC,JC,IC)
     $   +DXF(IF)/2.0*(FC(KC,JC,IC+1)-FC(KC,JC,IC))/DXC(IC)
        FF(KF+1,JF+1,IF+2) = FC(KC,JC,IC)
     $   +(DXF(IF)/2.0+DXF(IF+1))*(FC(KC,JC,IC+1)-FC(KC,JC,IC))/DXC(IC)
        FF(KF+2,JF+1,IF+1) = FC(KC+1,JC,IC)
     $   +DXF(IF)/2.0*(FC(KC+1,JC,IC+1)-FC(KC+1,JC,IC))/DXC(IC)
        FF(KF+2,JF+1,IF+2) = FC(KC+1,JC,IC)
     $   +(DXF(IF)/2.0+DXF(IF+1))
     $   *(FC(KC+1,JC,IC+1)-FC(KC+1,JC,IC))/DXC(IC)

         STEIGUNGZ1 = (FF(KF+2,JF+1,IF+1)-FF(KF+1,JF+1,IF+1))/DZC(KC)
         STEIGUNGZ2 = (FF(KF+2,JF+1,IF+2)-FF(KF+1,JF+1,IF+2))/DZC(KC)

        FF(KF+2,JF+1,IF+1) = FF(KF+1,JF+1,IF+1) +
     $    (DZF(KF)/2.0+DZF(KF+1))*STEIGUNGZ1
        FF(KF+1,JF+1,IF+1) = FF(KF+1,JF+1,IF+1) +
     $    (DZF(KF)/2.0)*STEIGUNGZ1
        FF(KF+2,JF+1,IF+2) = FF(KF+1,JF+1,IF+2) +
     $    (DZF(KF)/2.0+DZF(KF+1))*STEIGUNGZ2
        FF(KF+1,JF+1,IF+2) = FF(KF+1,JF+1,IF+2) +
     $    (DZF(KF)/2.0)*STEIGUNGZ2

         ENDDO
         ENDDO
         ENDDO
C-------------- Randbeingung NO-SLIP WALL TOP AND BOTTOM
C-------------- Bottom Index KF TOP Index KF2
         KF = 1
         KC = KPOS + (KF-3)/2
         KF2 = KKF-3
         KC2 = KPOS + (KF2-3)/2
         DO IF = 1,IIF-1,2
            IC = IPOS + (IF-3)/2
         DO JF = 1,JJF-1,2
            JC = JPOS + (JF-3)/2
C-------- bottom
        FF(KF+2,JF+1,IF+1) = FC(KC+1,JC,IC)
     $   +DXF(IF)/2.0*(FC(KC+1,JC,IC+1)-FC(KC+1,JC,IC))/DXC(IC)
        FF(KF+2,JF+1,IF+2) = FC(KC+1,JC,IC)
     $   +(DXF(IF)/2.0+DXF(IF+1))
     $   *(FC(KC+1,JC,IC+1)-FC(KC+1,JC,IC))/DXC(IC)
         STEIGUNGZ1 = (FF(KF+2,JF+1,IF+1)+FF(KF+2,JF+1,IF+1))/DZC(KC)
         STEIGUNGZ2 = (FF(KF+2,JF+1,IF+2)+FF(KF+2,JF+1,IF+2))/DZC(KC)
         FF(KF+2,JF+1,IF+1) = DZF(KF+1)/2.0*STEIGUNGZ1
         FF(KF+2,JF+1,IF+2) = DZF(KF+1)/2.0*STEIGUNGZ2
C-------- top
        FF(KF2+1,JF+1,IF+1) = FC(KC2,JC,IC)
     $   +DXF(IF)/2.0*(FC(KC2,JC,IC+1)-FC(KC2,JC,IC))/DXC(IC)
        FF(KF2+1,JF+1,IF+2) = FC(KC2,JC,IC)
     $   +(DXF(IF)/2.0+DXF(IF+1))
     $   *(FC(KC2,JC,IC+1)-FC(KC2,JC,IC))/DXC(IC)
         STEIGUNGZ1 = (FF(KF2+1,JF+1,IF+1)+FF(KF2+1,JF+1,IF+1))/DZC(KC2)
         STEIGUNGZ2 = (FF(KF2+1,JF+1,IF+2)+FF(KF2+1,JF+1,IF+2))/DZC(KC2)
         FF(KF2+1,JF+1,IF+1) = DZF(KF2+1)/2.0*STEIGUNGZ1
         FF(KF2+1,JF+2,IF+1) = DZF(KF2+1)/2.0*STEIGUNGZ2
         ENDDO
         ENDDO
C-------------- INTERPOLATION IN Y-RICHTUNG FUER DIE UEBRIGEN FINE-KOORD
         DO IF = 3,IIF-2
          DO JF = 3,JJF-1,2
           DO KF = 3,KKF-2
            FF(KF,JF,IF) = 0.5*(FF(KF,JF+1,IF)+FF(KF,JF-1,IF))
           ENDDO
          ENDDO
         ENDDO
C-------------------------------------- GROESSE IN Z-RICHTUNG VERSETZT
C                                       2-D INTERPOLATION (F.S.)

      ELSEIF ( CID .EQ. 'Z' ) THEN

         DO IF = 1,IIF-1,2
            IC = IPOS + (IF-3)/2
         DO JF = 1,JJF-1,2
            JC = JPOS + (JF-3)/2
         DO KF = 3,KKF-5,2
            KC = KPOS + (KF-3)/2
C------------- INTERPOLATION IN XY-EBENE FUER FINE AUF GLEICHER Z-KORD
C              WIE GROB KF+1 = KC
        FF(KF+1,JF+1,IF+1) = FC(KC,JC,IC)
     $   +DXF(IF)/2.0*(FC(KC,JC,IC+1)-FC(KC,JC,IC))/DXC(IC)
        FF(KF+1,JF+1,IF+2) = FC(KC,JC,IC)
     $   +(DXF(IF)/2.0+DXF(IF+1))*(FC(KC,JC,IC+1)-FC(KC,JC,IC))/DXC(IC)
        FF(KF+1,JF+2,IF+1) = FC(KC,JC+1,IC)
     $   +DXF(IF)/2.0*(FC(KC,JC+1,IC+1)-FC(KC,JC+1,IC))/DXC(IC)
        FF(KF+1,JF+2,IF+2) = FC(KC+1,JC,IC)
     $   +(DXF(IF)/2.0+DXF(IF+1))
     $   *(FC(KC,JC+1,IC+1)-FC(KC,JC+1,IC))/DXC(IC)

         STEIGUNGY1 = (FF(KF+1,JF+2,IF+1)-FF(KF+1,JF+1,IF+1))/DYC(JC)
         STEIGUNGY2 = (FF(KF+1,JF+2,IF+2)-FF(KF+1,JF+1,IF+2))/DYC(JC)

        FF(KF+1,JF+2,IF+1) = FF(KF+1,JF+1,IF+1) +
     $    (DYF(JF)/2.0+DYF(JF+1))*STEIGUNGY1
        FF(KF+1,JF+1,IF+1) = FF(KF+1,JF+1,IF+1) +
     $    (DYF(JF)/2.0)*STEIGUNGY1
        FF(KF+1,JF+2,IF+2) = FF(KF+1,JF+1,IF+2) +
     $    (DYF(JF)/2.0+DYF(JF+1))*STEIGUNGY2
        FF(KF+1,JF+1,IF+2) = FF(KF+1,JF+1,IF+2) +
     $    (DZF(YF)/2.0)*STEIGUNGY2

         ENDDO
         ENDDO
         ENDDO
C-------------- INTERPOLATION IN Y-RICHTUNG FUER DIE UEBRIGEN FINE-KOORD
         DO IF = 3,IIF-1
          DO JF = 3,JJF-1
           DO KF = 3,KKF-3
            FF(KF,JF,IF) = 0.5*(FF(KF+1,JF,IF)+FF(KF-1,JF,IF))
           ENDDO
          ENDDO
         ENDDO

C
C----------------------------------------------- GROESSE NICHT VERSETZT
C                                                VOLUMENMITTELWERT
C                                                KEINE INTERPOLATION ERLAUBT
C
C
C-------------------------------------- IM FEHLERFALL
C
      ELSE

         CALL ERRR (506,' MGCTOF')

      ENDIF
C
C----------------------------------------------------------------
C
      RETURN
      END


