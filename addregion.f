










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

      SUBROUTINE ADDREGION (KK1,JJ1,II1,FELD1,IPROC1,
     $                       KA,KE,JA,JE,IA,IE,
     $                       KK2,JJ2,II2,FELD2,IPROC2,
     $                       KSHIFT,JSHIFT,ISHIFT,HILF,MYID,FAVERAGE)
C-MGLET---------------------------------------------------------------72
C
C                  COPYING OF A SUB-AREA FROM FIELD 1 INTO FIELD 2
C
C
C        25.10.95 (MM)  : ORIGINAL
C
C-MGLET---------------------------------------------------------------72

      REAL FELD1(KK1,JJ1,II1)
      REAL FELD2(KK2,JJ2,II2)
      REAL HILF (KK1*JJ1*II1)


C
C      	write (6,*)'addregion,1:',myid,iproc1,iproc2,ia,ie,ja,je,ka,ke,
C     $             ishift,jshift,kshift
C

C-MGLET--------------- FIELD ON SAME PROCESSORS ----------------------72

            
            DO I1=IA,IE
                     I2=I1+ISHIFT
               DO J1=JA,JE
                     J2=J1+JSHIFT
                  DO K1=KA,KE

                     K2=K1+KSHIFT

                     FELD2(K2,J2,I2) = FELD2(K2,J2,I2) + 
     $                                 FELD1(K1,J1,I1)*FAVERAGE
                     
                  ENDDO
               ENDDO
            ENDDO
            

C-MGLET---------------------------------------------------------------72
      RETURN
      END

