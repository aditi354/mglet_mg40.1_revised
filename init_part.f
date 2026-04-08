










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
      SUBROUTINE INIT_PARTICLES (NXNPART,NPART_MAX,
     $     XPART,START_VEC,RUN,NPART,PFLOWCASE)

C--------------------------------------------------
C
C     INITIALISES PARTICLE POSITIONS:
C
C         PFLOWCASE == 1 > T-MISCHER
C         PFLOWCASE == 2 > CHANNEL
C
C     26.05.03 (FS): ORIGINAL
C
C                   
C--------------------------------------------------

      IMPLICIT NONE
      INTEGER NXNPART,NPART_MAX,NPART,K,J,I,F,PFLOWCASE
      INTEGER RUN(NPART_MAX)
      REAL    XPART( NPART_MAX,3),START_VEC(NPART_MAX,3)
      REAL    xo,yo,a,ro,Y,X,H,VOL(NPART_MAX),SUM,zo,Z

      IF (PFLOWCASE .EQ. 1) THEN
      WRITE(6,*)' NPART, NPART_MAX: ',NPART, NPART_MAX
      WRITE(6,*)'IMPLEMENTIERT FUER T-MIXER'
      IF (NXNPART**2/1.3 .ge. NPART_MAX) THEN
         CALL ERRR (501,' INIT_PART ')
      ENDIF
      xo = 0.25
      yo = 0.5
      K = 0
      SUM = 0.0
      DO I = 1,NXNPART
         X = 0.5/(NXNPART)*I
         DO J = 1,NXNPART
            Y = 0.25 + 0.5/(NXNPART)*J
            H = ((X-xo)**2+(Y-yo)**2)**0.5
            IF (H .lt. (0.25-(0.5*2**0.5/(2*NXNPART)))) THEN
               K = K + 1

               START_VEC(K,1) = X
               START_VEC(K,2) = Y
               START_VEC(K,3) = 0.0
C------------------------ CALCULATION OF VOLUMEFLUX OF PARTIKEL
               VOL(K) = 1.41*(1-(H/0.25)**2)*(0.5/NXNPART)**2
               SUM = SUM + 1.41*(1-(H/0.25)**2)*(0.5/NXNPART)**2
            ENDIF
         ENDDO
      ENDDO

      WRITE(6,*)'PUNKTE IM KREIS: ',K
      WRITE(6,*)'SUMME DES VOLUMENSTROMS: ',SUM
      NPART = K
      DO I=1,NPART
         XPART(I,1) =  START_VEC(I,1)
         XPART(I,2) =  START_VEC(I,2)
         XPART(I,3) =  START_VEC(I,3)
      ENDDO
      
      do i=1,NPART
         write (6,1000) 'PARTICLE: ',i,xpart(i,1),xpart(i,2),
     $        xpart(i,3),RUN(I),'VOL: ',VOL(I)
      enddo
C-------------------------- INITIALISATION FOR T-MIXER + FRED BODY
      ELSEIF (PFLOWCASE .EQ. 3) THEN
      WRITE(6,*)' NPART, NPART_MAX: ',NPART, NPART_MAX
      WRITE(6,*)'IMPLEMENTIERT FUER T-MIXER mit ZULAEUFEN'
      IF (NXNPART**2/1.3 .ge. NPART_MAX) THEN
         CALL ERRR (501,' INIT_PART ')
      ENDIF

      xo = 0.25
      yo = 0.8
      zo = 0.5
      K = 0
      SUM = 0.0
      DO I = 1,NXNPART
         X = 0.5/(NXNPART)*I
         DO J = 1,NXNPART
            Z = 0.25 + 0.5/(NXNPART)*J
            H = ((X-xo)**2+(Z-zo)**2)**0.5
            IF (H .lt. (0.25-(0.5*2**0.5/(2*NXNPART)))) THEN
               K = K + 1

               START_VEC(K,1) = X
               START_VEC(K,2) = yo
               START_VEC(K,3) = Z
C------------------------ CALCULATION OF VOLUMEFLUX OF PARTIKEL
               VOL(K) = 1.41*(1-(H/0.25)**2)*(0.5/NXNPART)**2
               SUM = SUM + 1.41*(1-(H/0.25)**2)*(0.5/NXNPART)**2
            ENDIF
         ENDDO
      ENDDO
      WRITE(6,*)'PUNKTE IM KREIS: ',K
      WRITE(6,*)'SUMME DES VOLUMENSTROMS: ',SUM
      NPART = K
      DO I=1,NPART
         XPART(I,1) =  START_VEC(I,1)
         XPART(I,2) =  START_VEC(I,2)
         XPART(I,3) =  START_VEC(I,3)
      ENDDO

      do i=1,NPART
         write (6,1000) 'PARTICLE: ',i,xpart(i,1),xpart(i,2),
     $        xpart(i,3),RUN(I),'VOL: ',VOL(I)
      enddo



C-------------------------- INITIALISATION FOR CHANNEL FLOW

      ELSEIF (PFLOWCASE .EQ. 2) THEN
       WRITE(6,*)'IMPLEMENTIERT FUER CHANNEL FLOW'
       F = 0
       DO I = 1,NXNPART
          X = 9.6/(NXNPART+1)*I 
        DO J = 1,NXNPART
           Y = 6.0/(NXNPART+1)*J
         DO K = 1,NXNPART
          F = F + 1
          START_VEC(F,1) = X
          START_VEC(F,2) = Y
          START_VEC(F,3) = 2.0/(NXNPART+1)*K
          XPART(F,1)=START_VEC(F,1)
          XPART(F,2)=START_VEC(F,2)
          XPART(F,3)=START_VEC(F,3)
         ENDDO
        ENDDO
       ENDDO
       NPART = F
      DO I=1,NPART
       WRITE(6,*)START_VEC(I,1),START_VEC(I,2),START_VEC(I,3)
      ENDDO
       WRITE(6,*)'GESAMTZAHL AN PARTIKEL: ',NPART   
      ELSE
         CALL ERRR (501,' INIT_PART ')
      ENDIF
       
      RETURN
 1000 FORMAT (A9,I3,3(2X,E12.4E2),I3,1X,A5,(E12.4E2))
      END
