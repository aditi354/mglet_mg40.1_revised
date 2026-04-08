










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
      SUBROUTINE FDERFOWZ(KK,JJ,II,KSTART,JSTART,ISTART,
     $     KSTOP,JSTOP,ISTOP,DZ,DDZ,COEFDZ,RSGS,FAKTOR,
     $     F,FIN,FD,NBOT,NTOP,LCOL,DIAG,RCOL
     $     ,BP,GI,GJ,GK
     $     )
C*MGLET***************************************************************
C        F D E R F O W Z  (First DERive Fourth Order)
C        BERECHNUNG DER ERSTEN ABLEITUNG (DW/DZ)                    
C*MGLET***************************************************************
C
C PARAM: F(K,J,I)       - BEKANNTE GROESSE AN DEN GITTERKANTEN
C      : FIN(K,J,I)     - INTERPOLIERTE GROESSE
C                         ZWISCHEN DEN BEIDEN KANTEN       
C      : FD (K,J,I)       - ERSTE ABLEITUNG, FIRST DERINATIVE  DW/DZ
C
C      : LCOL           - ENTHAELT KOEFFIZIENTEN DER "LINKEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C      : DIAG           - ENTHAELT KOEFFIZIENTEN DER HAUPTDIAGONALEN
C      : RCOL           - ENTHAELT KOEFFIZIENTEN DER "RECHTEN SPALTE"
C                         PARALLEL ZUR HAUPTDIAGONALEN
C
C DEFINE DIREKTIVEN     : KEINE 
C
C UPROG                 : THOMASK     THOMAS-ALGORITHMUS IN Z-RICHTUNG
C
C VERS:  08.01.97(AM)  : ORIGINAL     (KOMPAKT 4. ORDNUNG)
C
C*MGLET***************************************************************

      IMPLICIT NONE

      INTEGER KK,JJ,II,K,J,I,KSTART,JSTART,ISTART,KSTOP,JSTOP,ISTOP,
     $        KA,KE,NBOT,NTOP

      REAL    DZ(KK),DDZ(KK),LCOL(KK),DIAG(KK),RCOL(KK),
     $        F(KK,JJ,II),FIN(KK,JJ,II),FD(KK,JJ,II),RSGS(KK,JJ,II),
     $        COEFDZ(KK,12*3),FAKTOR(KK,JJ,II),A,B,C,D,
     $        GI(KK,JJ,II),GJ(KK,JJ,II),GK(KK,JJ,II),BP(KK,JJ,II)
  

      IF (NBOT .EQ. 5) THEN
         KA = 3
      ELSE
         WRITE(6,*) 'RANDBEDINGUNG IN FDERFOWZ NOCH NICHT '
         STOP
      ENDIF  
      
      IF (NTOP .EQ. 5 .OR. NTOP .EQ. 3) THEN
         KE = KK-2
      ELSE
         WRITE(6,*) 'RANDBEDINGUNG IN FDERFOWZ NOCH NICHT '
         STOP
      ENDIF

C                                    ERSTER PHYSIKALISCHER PUNKT
      K = KA

      DO I = ISTART, ISTOP
         DO J = JSTART, JSTOP
            GI(K,J,I) = COEFDZ(K,1)
            GJ(K,J,I) = COEFDZ(K,2)
            GK(K,J,I) = COEFDZ(K,3)

            RSGS(K,J,I) =  COEFDZ(K,4) * F(K-1,J,I)
     $                   + COEFDZ(K,5) * F(K  ,J,I)
     $                   + COEFDZ(K,6) * F(K+1,J,I)
   
         ENDDO
      ENDDO

C                                    LETZTER PHYSIKALISCHER PUNKT

      K = KE

      DO I = ISTART, ISTOP
         DO J = JSTART, JSTOP
            GI(K,J,I) = COEFDZ(K,1)
            GJ(K,J,I) = COEFDZ(K,2)
            GK(K,J,I) = COEFDZ(K,3)
            RSGS(K,J,I) =  COEFDZ(K,4) * F(K  ,J,I)
     $                   + COEFDZ(K,5) * F(K-1,J,I)
     $                   + COEFDZ(K,6) * F(K-2,J,I)
C            RSGS(K,J,I) =
C     $                   (- F(K-1,J,I))*DZ(K)
         ENDDO
      ENDDO



C                                                IM GEBIET 

      DO I = ISTART, ISTOP
         DO J = JSTART, JSTOP
            DO K = KA+1 , KE-1

C          A =    BP(K,J,I)
C          B = (1-BP(K,J,I))*    BP(K+1,J,I)
C          C = (1-BP(K,J,I))*    BP(K-1,J,I)
C          D = (1-BP(K,J,I))* (1-BP(K-1,J,I))* (1-BP(K+1,J,I))

          A =    BP(K-1,J,I) *   BP(K,J,I)
          B = (1-BP(K-1,J,I))*   BP(K,J,I)
          C = (1-BP(K,J,I  ))*   BP(K-1,J,I)
          D = (1-BP(K,J,I  ))*(1-BP(K-1,J,I))

          GI(K,J,I) =  COEFDZ(K,1)    * A
     $               + COEFDZ(K,12+1) * B
     $               + COEFDZ(K,24+1) * C
          GJ(K,J,I) =  COEFDZ(K,2)    * A
     $               + COEFDZ(K,12+2) * B
     $               + COEFDZ(K,24+2) * C
     $               + 1.0            * D
          GK(K,J,I) =  COEFDZ(K,3)    * A
     $               + COEFDZ(K,12+3) * B
     $               + COEFDZ(K,24+3) * C

                RSGS(K,J,I) =
     $   (COEFDZ(K,4)    * F(K,J,I)+
     $    COEFDZ(K,5)    * FIN(K,J,I)+
     $    COEFDZ(K,6)    * F(K-1,J,I))*DDZ(K) * A
     $ + (COEFDZ(K,12+4) * F(K-1,J,I)+
     $    COEFDZ(K,12+5) * F(K,J,I)+
     $    COEFDZ(K,12+6) * F(K+1,J,I)) *        B
     $ + (COEFDZ(K,24+4) * F(K,J,I)+
     $    COEFDZ(K,24+5) * F(K-1,J,I)+
     $    COEFDZ(K,24+6) * F(K-2,J,I)) *        C


            ENDDO
         ENDDO
      ENDDO

      CALL THOMASK(KK,JJ,II,KA,KE,JSTART,JSTOP,ISTART,ISTOP,
     $     LCOL,DIAG,RCOL,RSGS,FAKTOR,FD
     $     ,BP,GI,GJ,GK
     $     )
C*********************  EXTRAPOLATION   *******************************

      IF (NTOP.EQ.3) THEN
         K = KE+1
         DO I = ISTART,ISTOP
            DO J = JSTART,JSTOP
               FD(K,J,I) = - COEFDZ(K,1) * FD(K-1,J,I)
     $                    +  COEFDZ(K,4) * F (K-1,J,I)
     $                    +  COEFDZ(K,5) * F (K-2,J,I)
     $                    +  COEFDZ(K,6) * F (K-3,J,I)
            ENDDO
         ENDDO
      ENDIF

      RETURN
      END

