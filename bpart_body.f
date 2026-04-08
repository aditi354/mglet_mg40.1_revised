










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
      SUBROUTINE BPART_BODY (KK,JJ,II,X,Y,Z,DX,DY,DZ,
     $                   NPART,XP,YP,ZP,
     $                   NFRO,NBAC,NRGT,NLFT,NBOT,NTOP
     $    ,IINDEX,JINDEX,KINDEX,IINDEX_N,JINDEX_N,KINDEX_N,BP)
C------------------------------------------------------------
C
C     TREATS BOUNDARY CONDITIONS FOR PARTICLES
C     IN COMBINATION WITH FRED_BODY
C     12.10.04 (FS)  ORIGINAL
C
C------------------------------------------------------------

      REAL    X(II),    Y(II),    Z(II),
     $       DX(II),   DY(JJ),   DZ(KK)

      REAL XP(NPART), YP(NPART), ZP(NPART)
      INTEGER NI,NJ,NK,DN_X,DN_Y,DN_Z
      REAL BP(KK,JJ,II)
      INTEGER IINDEX(NPART),JINDEX(NPART),KINDEX(NPART)
      INTEGER IINDEX_N(NPART),JINDEX_N(NPART),KINDEX_N(NPART)
      LOGICAL DIR
C---------------------------------- BODY IN FLOWFIELD
      DO IPART = 1,NPART
C--------- IF PARTICLE POSITION IS INSIDE BODY
       NI = IINDEX(IPART)
       NJ = JINDEX(IPART)
       NK = KINDEX(IPART)
C--------- Debugging
C       IF(IPART.EQ.71) THEN
C          WRITE(6,*)'Number of PArtikel: ',IPART
C          WRITE(6,*)'Indizes, old: ',NI,NJ,NK
C          WRITE(6,*)'Indizes, new: ',IINDEX_N(IPART),
C     $    JINDEX_N(IPART),KINDEX_N(IPART)
C          WRITE(6,*)'Coordinates:',XP(IPART),YP(IPART),ZP(IPART)
C       ENDIF
C********************
       IF(BP(KINDEX_N(IPART),JINDEX_N(IPART),IINDEX_N(IPART)) .EQ. 0)
     $    THEN
C          WRITE(6,*)'Number of PArtikel: ',IPART
C          WRITE(6,*)'Indizes, old: ',NI,NJ,NK
C          WRITE(6,*)'Indizes, new: ',IINDEX_N(IPART),
C     $    JINDEX_N(IPART),KINDEX_N(IPART)
          DN_X = IINDEX_N(IPART) - NI
          DN_Y = JINDEX_N(IPART) - NJ
          DN_Z = KINDEX_N(IPART) - NK
C         WRITE(6,*)'DN_X: ',DN_X
C          WRITE(6,*)'DN_Y: ',DN_Y
C          WRITE(6,*)'DN_Z: ',DN_Z
C          WRITE(6,*)'BP(K,J,I+DNX): ',BP(NK,NJ,NI+DN_X)
C          WRITE(6,*)'BP(K,J+DNY,I): ',BP(NK,NJ+DN_Y,NI)
CC          WRITE(6,*)'BP(K+DNZ,J,I): ',BP(NK+DN_Z,NJ,NI)
C          WRITE(6,*)'BP(K,J+DNY,I+DNX): ',BP(NK,NJ+DN_Y,NI+DN_X)
C          WRITE(6,*)'BP(K+DNZ,J,I+DNX): ',BP(NK+DN_Z,NJ,NI+DN_X)
C          WRITE(6,*)'BP(K+DNZ,J+DNY,I): ',BP(NK+DN_Z,NJ+DN_Y,NI)
C          WRITE(6,*)'BP(K+DNZ,J+DNY,I+DNX): ',
C     $                               BP(NK+DN_Z,NJ+DN_Y,NI+DN_X)
          DIR = .FALSE.
c-------------
          IF(BP(NK,NJ,NI+DN_X).EQ.0)
     $    THEN
          DIR = .TRUE.
c------------- NNX is wether -1 or 0 depending on DN_X
          NNX = int(-0.5 + SIGN(0.5,float(DN_X)))
          XP(IPART) = (X(NI+NNX)+0.5*DX(NI+NNX))
     $                 -(XP(IPART)-(X(NI+NNX)+0.5*DX(NI+NNX)))
C          WRITE(6,*)'CASE 1'
          ENDIF
          IF(BP(NK,NJ+DN_Y,NI).EQ.0)
     $    THEN
          DIR = .TRUE.
          NNY = int(-0.5 + SIGN(0.5,float(DN_Y)))
          YP(IPART) = (Y(NJ+NNY)+0.5*DY(NJ+NNY))
     $                 -(YP(IPART)-(Y(NJ+NNY)+0.5*DY(NJ+NNY)))
C          WRITE(6,*)'CASE 2'
          ENDIF
          IF(BP(NK+DN_Z,NJ,NI).EQ.0)
     $    THEN
          DIR = .TRUE.
          NNZ = int(-0.5 + SIGN(0.5,float(DN_Z)))
          ZP(IPART) = (Z(NK+NNZ)+0.5*DZ(NK+NNZ))
     $                 -(ZP(IPART)-(Z(NK+NNZ)+0.5*DZ(NK+NNZ)))
C          WRITE(6,*)'CASE 3'
          ENDIF
c-------------
          IF(DIR.EQ.(.FALSE.)) THEN
           IF(BP(NK,NJ+DN_Y,NI+DN_X).EQ.0)
     $      THEN
C            WRITE(6,*)'CASE 4'
            NNX = int(-0.5 + SIGN(0.5,float(DN_X)))
            NNY = int(-0.5 + SIGN(0.5,float(DN_Y)))
          XP(IPART) = (X(NI+NNX)+0.5*DX(NI+NNX))
     $                 -(XP(IPART)-(X(NI+NNX)+0.5*DX(NI+NNX)))
          YP(IPART) = (Y(NJ+NNY)+0.5*DY(NJ+NNY))
     $                 -(YP(IPART)-(Y(NJ+NNY)+0.5*DY(NJ+NNY)))
           ELSEIF(BP(NK+DN_Z,NJ,NI+DN_X).EQ.0)
     $      THEN
C            WRITE(6,*)'CASE 5'
            NNX = int(-0.5 + SIGN(0.5,float(DN_X)))
            NNZ = int(-0.5 + SIGN(0.5,float(DN_Z)))
          XP(IPART) = (X(NI+NNX)+0.5*DX(NI+NNX))
     $                 -(XP(IPART)-(X(NI+NNX)+0.5*DX(NI+NNX)))
          ZP(IPART) = (Z(NK+NNZ)+0.5*DZ(NK+NNZ))
     $                 -(ZP(IPART)-(Z(NK+NNZ)+0.5*DZ(NK+NNZ)))
           ELSEIF(BP(NK+DN_Z,NJ+DN_Y,NI).EQ.0)
     $      THEN
C            WRITE(6,*)'CASE 6'
            NNY = int(-0.5 + SIGN(0.5,float(DN_Y)))
            NNZ = int(-0.5 + SIGN(0.5,float(DN_Z)))
          YP(IPART) = (Y(NJ+NNY)+0.5*DY(NJ+NNY))
     $                 -(YP(IPART)-(Y(NJ+NNY)+0.5*DY(NJ+NNY)))
          ZP(IPART) = (Z(NK+NNZ)+0.5*DZ(NK+NNZ))
     $                 -(ZP(IPART)-(Z(NK+NNZ)+0.5*DZ(NK+NNZ)))
c-------------
           ELSE
C            WRITE(6,*)'CASE 7'
            NNX = int(-0.5 + SIGN(0.5,float(DN_X)))
            NNY = int(-0.5 + SIGN(0.5,float(DN_Y)))
            NNZ = int(-0.5 + SIGN(0.5,float(DN_Z)))
C            WRITE(6,*)'CASE 7',NNX,NNY,NNZ
C            WRITE(6,*)'Xcor: ',XP(IPART)
          XP(IPART) = (X(NI+NNX)+0.5*DX(NI+NNX))
     $                 -(XP(IPART)-(X(NI+NNX)+0.5*DX(NI+NNX)))
C            WRITE(6,*)'Xcor: ',XP(IPART)
C            WRITE(6,*)'Ycor: ',YP(IPART)
          YP(IPART) = (Y(NJ+NNY)+0.5*DY(NJ+NNY))
     $                 -(YP(IPART)-(Y(NJ+NNY)+0.5*DY(NJ+NNY)))
C            WRITE(6,*)'Ycor: ',YP(IPART)
C            WRITE(6,*)'Zcor: ',ZP(IPART)
          ZP(IPART) = (Z(NK+NNZ)+0.5*DZ(NK+NNZ))
     $                 -(ZP(IPART)-(Z(NK+NNZ)+0.5*DZ(NK+NNZ)))
C            WRITE(6,*)'Zcor: ',ZP(IPART)
          ENDIF
          ENDIF

       ENDIF
       ENDDO



      RETURN
      END
