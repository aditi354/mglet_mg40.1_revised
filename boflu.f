










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
      SUBROUTINE BOFLU(KK,JJ,II,UBO,U,HILF,AU,JJA,NBND,X,XBANF,XBEND)
C--MGLET----------------------------------------------------------
C
C     ADDS FLUCTUATIONS TO BOTTOM-BUFFERS FOR 
C     "STRUKTUR-PERIODISCHE" BOUNDARY CONDITION
C
C        UBO(X0)=UBO(X0) + (U(X) - <U>(X))
C
C    03.07.97 (MM.):    ORIGINAL AUS SVFLU ABGELEITET
C
C-------10--------20--------30--------40--------50--------60--------7072

CCCCCCADDED FOR CORRECTION BY TB, 280503 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      REAL
     $     UBO(JJ,II, 2),U(JJ,II),HILF(II,JJ),X(II)
      REAL AU(JJA,II)
	data icount /0/
	save icount

	icount = icount + 1

C---------------------------------------------------------------------72
      IF (NBND .LT. 1) CALL ERRR(501,"BOFLU")
C-------------------------------- AVERAGING IN Y ---------------------72

      FAC = 1./FLOAT(JJ - 2*NBND)

C                                        FIRST VALUE
      DO I=NBND+1,II-NBND+1
         HILF(I,NBND+1) = FAC*U(NBND+1,I)
      ENDDO

C                                        SUMMING UP


      DO I= NBND+1,II-NBND+1
         DO J=NBND+2,JJ-NBND
            HILF(I,J) = HILF(I,J-1) + FAC*U(J,I)
         ENDDO
      ENDDO
C                                      AVERAGE IS IN HILF(I,JJ-NBND)



C------------------------------------- NOW ADDING FLUCTUATIONS -------72

      DO I=NBND+1,II-NBND+1
         DISTANCE = MIN(0.0,(XBEND-X(I))/(XBEND-XBANF))
         FAK = EXP(DISTANCE)**4
CCC         WRITE (6,*) 'BOFLU:',I,X(I),DISTANCE,FAK
         DO J=   NBND+1,JJ-NBND

            UBO(J,I,2)=UBO(J,I,2) + (U(J,I)-HILF(I,JJ-NBND))*FAK

         ENDDO
      ENDDO

C---------------------------------------------------------------------72
       iwrite = -20000
       if (iwrite .eq. 1) then
       if (icount .lt.20) then
       if (jj .eq. 16) then
       do i=1,ii
        write (6,'(A10,i3,4(1X,E12.5E3))') 'boflu_grd1',icount,
     $              x(i),UBO(3,i, 2),U(3,i),HILF(I,JJ-NBND)
       enddo

       elseif (jj .eq. 28) then
       do i=1,ii
        write (6,'(A10,i3,4(1X,E12.5E3))') 'boflu_grd2',icount,
     $              x(i),UBO(3,i, 2),U(3,i),HILF(I,JJ-NBND)
       enddo
       elseif (jj .eq. 52) then
       do i=1,ii
        write (6,'(A10,i3,4(1X,E12.5E3))') 'boflu_grd3',icount,
     $              x(i),UBO(3,i, 2),U(3,i),HILF(I,JJ-NBND)
       enddo


       endif
       endif
       endif
C---------------------------------------------------------------------72

      RETURN
      END
