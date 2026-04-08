










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
      SUBROUTINE SEARCHINDEX(II,X,DX,DDX,
     $     NXGRAE,NPART,XPART,INDEX,HILF)
C----------------------------------------------------------
C
C         SUCHT DEN INDEX, INDEX SICH EIN PARTIKEL BEFINDET
C
C   1.9.99 (MM) ORIGINAL
C
C-----------------------------------------------------------
      
      INTEGER INDEX(NPART)

      REAL XPART(NPART)
      REAL X(II), DX(II), DDX(II), HILF(II)
      
      IF (NXGRAE .EQ. 1) THEN
C--------------------------------------- HIER: AEQUIDISTANTES GITTER

         XSTART = X(   2) + 0.5*DX(   2)
         XSTOP  = X(II-2) + 0.5*DX(II-2)

         XLENGTH= XSTOP - XSTART

         DO IPART = 1,NPART
            INDEX(IPART) = (XPART(IPART)-XSTART)/XLENGTH * (II-4) + 3
C            write (0,*) 'xpart:',xpart(ipart),index(ipart)
         ENDDO

      ELSE
C--------------------------------------- HIER: NICHT-AEQUIDISTANTES GITTER

         DO I=1,II
            HILF(I) = X(I) + 0.5*DX(I)
         ENDDO

         IALGO = 1
         IF (IALGO .EQ. 1) THEN
         
            DO I=2,II-1
               DO IPART = 1,NPART
                  IF (XPART(IPART) .GT. HILF(I-1) .AND.
     $                 XPART(IPART) .LE. HILF(I)        ) THEN
                     INDEX(IPART) = I
                  ENDIF
               ENDDO
            ENDDO
         
         ELSEIF (IALGO .EQ. 2) THEN

            DO IPART = 1,NPART
               DO I=2,II-1
                  IF ( XPART(IPART) .LE. HILF(I)) THEN
                     INDEX(IPART) = I
                     GOTO 1000
                  ENDIF
               ENDDO
 1000          CONTINUE
            ENDDO

         ELSE


            r1 = log(float(ii-4))/log(2.0)
            i1 = nint(r1)
            icontrol = 2**i1
            if (icontrol .ne. ii-4) then
               write (0,*) itest,' is not a power of 2!'
               stop
            endif
            
            DO IPART = 1,NPART

            IMARK = 0
            NX = II-4
            
            DO I=1,100

               NX = NX/2
               IMARK2 = IMARK + NX
               
C               write (0,*) i,nx,imark,imark2,imark2+2
               IF (XPART(IPART) .GE. HILF(IMARK2+2)) THEN
                  IMARK = IMARK + NX
               ELSE
                  IMARK = IMARK
               ENDIF
               
               IF (NX .EQ. 1) GOTO 2000
            ENDDO

 2000       CONTINUE

            INDEX(IPART) = IMARK + 3

            ENDDO
            
         ENDIF

      ENDIF

C         DO IPART = 1,NPART
C            WRITE (0,*) 'XPART, NON-EQU:',XPART(IPART),INDEX(IPART)
C         ENDDO


      RETURN
      END
